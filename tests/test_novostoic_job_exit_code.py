"""A failed novostoic run must write the `error` sentinel and still let post-job run.

The novostoic commands in chart/values.yaml used to write a `fail` marker that nothing
read, so a failed run was marked 'completed' with no results (119 completed jobs on
prod, 2025-03-27 to 2026-10-26). They now write `error`, which postjob.py fails the Job
on after uploading the logs (tests/test_postjob_exit_code.py).

These render each deployed command through the real job template with create_job, then
run the job container's script with a stand-in `python` on PATH.
"""
import os
import pathlib
import subprocess

import pytest
import yaml

from services import kubejob_service

JOB_ID = "abc123"
VALUES = yaml.safe_load(
    (pathlib.Path(__file__).resolve().parent.parent / "chart" / "values.yaml").read_text())
JOBS = VALUES["config"]["kubernetes_jobs"]
# Built from the chart so that a novostoic job added there is covered too
NOVOSTOIC_TYPES = sorted(t for t in JOBS if t.startswith("novostoic-"))


def test_every_novostoic_type_is_covered():
    assert NOVOSTOIC_TYPES == [
        "novostoic-dgpredictor", "novostoic-enzrank", "novostoic-optstoic", "novostoic-pathways"]


@pytest.fixture(params=NOVOSTOIC_TYPES)
def job_script(request, monkeypatch):
    submitted = []

    class _FakeBatchV1Api:
        def create_namespaced_job(self, namespace, body):
            submitted.append(body)

    monkeypatch.setattr(kubejob_service, "api_batch_v1", _FakeBatchV1Api())
    monkeypatch.setattr(kubejob_service, "get_namespace", lambda: "test")

    kubejob_service.create_job(job_type=request.param, job_id=JOB_ID, image_name="example/image:tag",
                               command=JOBS[request.param]["command"])
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
    for d in ("in", "out"):
        (tmp_path / d).mkdir()
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
    code, outputs = _run(job_script, tmp_path, f"echo Traceback >&2\n{failure}\n")

    assert code == 0
    assert outputs == ["error", "error.log", "finished", "output.log"]


def test_successful_run_finishes(job_script, tmp_path):
    code, outputs = _run(job_script, tmp_path, 'echo done\ntouch "$JOB_OUTPUT_DIR/output.json"\n')

    assert code == 0
    assert outputs == ["error.log", "finished", "output.json", "output.log", "success"]
