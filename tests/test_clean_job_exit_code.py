"""A CLEAN job whose inference fails must write the `error` sentinel, not report success.

The CLEAN command pipes the inference through `tee` into the job log. Without
`set -o pipefail` the pipeline's status is tee's, so a crash (job f289fd27 on 2026-11-04,
`torch.cat(): expected a non-empty list of Tensors`) skipped the `error` sentinel, the job
template touched `finished`, and the job was marked completed with no results.

These render the real job template through create_job with the real CLEAN command, then
run the job container's script with a stand-in `python` on PATH.

A failed inference still exits 0: the job container is an init container, so a non-zero
exit would stop post-job from uploading the log. postjob.py fails the Job on `error`
after the upload (tests/test_postjob_exit_code.py).
"""
import os
import subprocess

import pytest

from services import kubejob_service
from services.clean_service import CleanService

JOB_ID = "abc123"
JOB_INFO = {"input_fasta": [{"header": "seq1", "sequence": "MKTAYIAKQR"}]}


def _job_script(monkeypatch):
    submitted = []

    class _FakeBatchV1Api:
        def create_namespaced_job(self, namespace, body):
            submitted.append(body)

    monkeypatch.setattr(kubejob_service, "api_batch_v1", _FakeBatchV1Api())
    monkeypatch.setattr(kubejob_service, "get_namespace", lambda: "test")

    command = CleanService.build_clean_job_command(job_id=JOB_ID, job_info=JOB_INFO)
    kubejob_service.create_job(job_type="clean", job_id=JOB_ID, image_name="example/image:tag", command=command)
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
def test_failed_inference_writes_error_and_lets_post_job_run(monkeypatch, tmp_path, failure):
    script = _job_script(monkeypatch)

    code, outputs = _run(script, tmp_path, f"echo Traceback\n{failure}\n")

    assert code == 0
    assert outputs == ["error", "finished", "log"]


def test_successful_inference_finishes(monkeypatch, tmp_path):
    script = _job_script(monkeypatch)

    code, outputs = _run(script, tmp_path, f'echo done\ntouch "$JOB_OUTPUT_DIR/{JOB_ID}_maxsep.csv"\n')

    assert code == 0
    assert outputs == [f"{JOB_ID}_maxsep.csv", "finished", "log"]
