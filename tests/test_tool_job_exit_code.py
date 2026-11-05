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
import pathlib
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


def _init_container_scripts(monkeypatch, job_type):
    """Render the job and return {init container name: its bash -c script}."""
    build, job_info = TOOLS[job_type]
    submitted = []

    class _FakeBatchV1Api:
        def create_namespaced_job(self, namespace, body):
            submitted.append(body)

    monkeypatch.setattr(kubejob_service, "api_batch_v1", _FakeBatchV1Api())
    monkeypatch.setattr(kubejob_service, "get_namespace", lambda: "test")

    command = build(job_id=JOB_ID, job_info=job_info)
    kubejob_service.create_job(job_type=job_type, job_id=JOB_ID, image_name="example/image:tag", command=command)
    scripts = {}
    for container in submitted[0]["spec"]["template"]["spec"]["initContainers"]:
        assert container["command"][0].endswith("bash") and container["command"][1] == "-c"
        scripts[container["name"]] = container["command"][2]
    return scripts


@pytest.fixture(params=sorted(TOOLS))
def job_script(request, monkeypatch):
    # The workload runs as the init container named "job"
    return _init_container_scripts(monkeypatch, request.param)["job"]


def _run(script, tmp_path, python_body, dirs=("data/inputs", "in", "out")):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    python = bin_dir / "python"
    python.write_text("#!/bin/sh\n" + python_body)
    python.chmod(0o755)
    for d in dirs:
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


def _run_crispr_pod(monkeypatch, tmp_path, fetch="echo fetched", convert="echo converted", tool="echo done"):
    """Run CRISPR-Copies' `init` container (prejob.py + the setup step), then its `job`
    container, in order, sharing one output directory, as the pod does.

    Returns (init exit code, job exit code, output file names, contents of `log`).
    """
    scripts = _init_container_scripts(monkeypatch, "crispr-copies")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    python = bin_dir / "python"
    python.write_text(
        "#!/bin/sh\n"
        'case "$1" in\n'
        f"  fetch_organism.py) {fetch} ;;\n"
        f"  convert_annotations.py) {convert} ;;\n"
        f"  *) {tool} ;;\n"
        "esac\n"
    )
    python.chmod(0o755)
    prejob = tmp_path / "prejob.py"  # the backend's MinIO input sync
    prejob.write_text("#!/bin/sh\nexit 0\n")
    prejob.chmod(0o755)
    out = tmp_path / "out"
    env = {
        **os.environ,
        "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
        "JOB_INPUT_DIR": str(tmp_path / "in"),
        "JOB_OUTPUT_DIR": str(out),
    }
    codes = [
        subprocess.run(["bash", "-c", scripts[name]], cwd=tmp_path, env=env, capture_output=True, text=True).returncode
        for name in ("init", "job")
    ]
    log = (out / "log").read_text() if (out / "log").exists() else ""
    return codes[0], codes[1], sorted(p.name for p in out.iterdir()), log


@pytest.mark.parametrize("step", ["fetch", "convert"])
def test_a_failed_crispr_setup_step_keeps_its_log_and_skips_the_tool(monkeypatch, tmp_path, step):
    # fetch_organism.py and convert_annotations.py run in the backend's init container.
    # A failure there used to fail that container, so neither the tool nor post-job ran
    # and nothing was uploaded.
    failing = {step: "echo 'HTTPError: 404 from NCBI'; exit 1"}

    init_code, job_code, outputs, log = _run_crispr_pod(
        monkeypatch, tmp_path, tool='touch "$JOB_OUTPUT_DIR/tool-ran"', **failing)

    assert (init_code, job_code) == (0, 0)  # so post-job runs, uploads, then fails the Job
    assert outputs == ["error", "finished", "log"]
    assert "HTTPError: 404 from NCBI" in log


def test_a_crispr_run_keeps_the_setup_output_in_its_log(monkeypatch, tmp_path):
    init_code, job_code, outputs, log = _run_crispr_pod(
        monkeypatch, tmp_path, tool='echo done; touch "$JOB_OUTPUT_DIR/output.csv"')

    assert (init_code, job_code) == (0, 0)
    assert outputs == ["finished", "log", "output.csv"]
    assert "fetched" in log and "converted" in log and "done" in log


def test_a_failed_crispr_tool_after_setup_keeps_both_outputs(monkeypatch, tmp_path):
    init_code, job_code, outputs, log = _run_crispr_pod(
        monkeypatch, tmp_path, tool="echo Traceback; exit 1")

    assert (init_code, job_code) == (0, 0)
    assert outputs == ["error", "finished", "log"]
    assert "fetched" in log and "Traceback" in log


def test_the_chart_runs_the_same_crispr_setup_command_as_the_app_config():
    # The tests above render the job from app/cfg/config.yaml; the cluster uses the chart.
    import yaml
    from config import app_config

    chart = yaml.safe_load((pathlib.Path(__file__).resolve().parent.parent / "chart" / "values.yaml").read_text())
    assert (chart["config"]["kubernetes_jobs"]["crispr-copies"]["prejob_command"]
            == app_config["kubernetes_jobs"]["crispr-copies"]["prejob_command"])


def test_the_tool_does_not_run_when_an_earlier_step_fails(job_script, tmp_path):
    # The template runs `ls -al ${JOB_INPUT_DIR} && <command> && ...`. With no input
    # directory that first step fails, and the chain must stop there. A command that
    # starts with `set -o pipefail;` would run the tool anyway.
    code, outputs = _run(job_script, tmp_path, 'touch "$JOB_OUTPUT_DIR/tool-ran"\n', dirs=("data/inputs", "out"))

    assert code != 0
    assert "tool-ran" not in outputs
