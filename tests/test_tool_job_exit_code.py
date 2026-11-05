"""A tool job that fails must write the `error` sentinel, not report success.

CLEAN, CRISPR-Copies and Mutagenesis pipe the tool through `tee` into the job log.
Without `set -o pipefail` the pipeline's status is tee's, so a crash (CLEAN job f289fd27
on 2026-11-04, `torch.cat(): expected a non-empty list of Tensors`) skipped the `error`
sentinel, the job template touched `finished`, and the job was marked completed with no
results.

These render the real job template through create_job with each tool's real command,
then run the job container's script with a stand-in `python` on PATH.

A failed run still exits 0: the job container is an init container, so a non-zero exit
would stop post-job from uploading the log. postjob.py fails the Job on `error` after
the upload (tests/test_postjob_exit_code.py).
"""
import os
import subprocess

import pytest

from services import kubejob_service
from services.clean_service import CleanService
from services.crispr_copies_service import CRISPRCopiesService
from services.mutagenesis_service import MutagenesisService

JOB_ID = "abc123"

# job type -> (command builder, a valid job_info). The CRISPR-Copies and Mutagenesis
# inputs are the known-good payloads from tests/test_tool_registry.py.
TOOLS = {
    "clean": (
        CleanService.build_clean_job_command,
        {"input_fasta": [{"header": "seq1", "sequence": "MKTAYIAKQR"}]},
    ),
    "crispr-copies": (
        CRISPRCopiesService.build_crispr_copies_job_command,
        {
            "organism": {"organismIdentifier": {"label": "E. coli", "value": "ecoli"}},
            "parameters": {
                "crisprSystem": "CRISPR/Cas9",
                "guideRNA": {
                    "pam": "NGG", "pamOrientation": "3prime",
                    "onTargetScore": "Doench et al. 2016", "guideLength": 20,
                    "seedLength": 12, "editDistance": 3,
                    "gcContentMin": 30, "gcContentMax": 70,
                    "distanceType": "hamming",
                },
                "homologyArm": {"lengthBp": 500},
                "intergenicRegion": {
                    "distance": 100, "distanceToMeasureGeneDensity": 5000,
                    "distalEndOfChromosomeToAvoid": 10000,
                },
            },
        },
    ),
    "mutagenesis": (
        MutagenesisService.build_mutagenesis_job_command,
        {"orfFile": {"filename": "orf.txt"}, "mutationList": {"filename": "mutations.csv"}},
    ),
}


@pytest.fixture(params=sorted(TOOLS))
def job_script(request, monkeypatch):
    job_type = request.param
    build, job_info = TOOLS[job_type]
    submitted = []

    class _FakeBatchV1Api:
        def create_namespaced_job(self, namespace, body):
            submitted.append(body)

    monkeypatch.setattr(kubejob_service, "api_batch_v1", _FakeBatchV1Api())
    monkeypatch.setattr(kubejob_service, "get_namespace", lambda: "test")

    command = build(job_id=JOB_ID, job_info=job_info)
    kubejob_service.create_job(job_type=job_type, job_id=JOB_ID, image_name="example/image:tag", command=command)
    pod_spec = submitted[0]["spec"]["template"]["spec"]
    # The workload runs as the init container named "job"
    job = next(c for c in pod_spec["initContainers"] if c["name"] == "job")
    assert job["command"][:2] == ["bash", "-c"]
    return job["command"][2]


def _run(script, tmp_path, python_body):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    python = bin_dir / "python"
    python.write_text("#!/bin/sh\n" + python_body)
    python.chmod(0o755)
    for d in ("data/inputs", "in", "out"):
        (tmp_path / d).mkdir(parents=True)
    env = {
        **os.environ,
        "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
        "JOB_INPUT_DIR": str(tmp_path / "in"),
        "JOB_OUTPUT_DIR": str(tmp_path / "out"),
    }
    result = subprocess.run(["bash", "-c", script], cwd=tmp_path, env=env, capture_output=True, text=True)
    return result.returncode, sorted(p.name for p in (tmp_path / "out").iterdir())


@pytest.mark.parametrize("failure", ["exit 1", "kill -9 $$"], ids=["exits-nonzero", "killed"])
def test_failed_run_writes_error_and_lets_post_job_run(job_script, tmp_path, failure):
    code, outputs = _run(job_script, tmp_path, f"echo Traceback\n{failure}\n")

    assert code == 0
    assert outputs == ["error", "finished", "log"]


def test_successful_run_finishes(job_script, tmp_path):
    code, outputs = _run(job_script, tmp_path, 'echo done\ntouch "$JOB_OUTPUT_DIR/result.csv"\n')

    assert code == 0
    assert outputs == ["finished", "log", "result.csv"]
