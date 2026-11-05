"""postjob.py when uploading a job's outputs to MinIO fails.

postjob.py is the Job pod's main container; the tool itself runs as an init container.
It used to log the upload error and exit 0, so Kubernetes saw a successful Job, the
watcher marked it 'completed' and emailed the user, and the results were missing
(25 CLEAN and 3 ML-SimpleFold jobs on 2026-11-04, while MinIO was unreachable).
A non-zero exit makes the Job fail, so the watcher marks it 'error' instead.

postjob.py runs at import time, so it is executed with runpy, as coordinator.py is in
test_coordinator_step_error.py.
"""
import os
import pathlib
import runpy
import subprocess

import pytest
import yaml

from services import kubejob_service

POSTJOB = pathlib.Path(__file__).resolve().parent.parent / "app" / "postjob.py"


@pytest.fixture
def job_env(monkeypatch, tmp_path):
    out = tmp_path / "uws" / "jobs" / "clean" / "job1" / "out"
    out.mkdir(parents=True)
    (out / "finished").touch()
    monkeypatch.setenv("JOB_TYPE", "clean")
    monkeypatch.setenv("JOB_ID", "job1")
    monkeypatch.setenv("JOB_OUTPUT_DIR", str(out))
    return out


def run_postjob():
    runpy.run_path(str(POSTJOB), run_name="__main__")


def test_successful_upload_exits_zero(monkeypatch, job_env):
    calls = []
    monkeypatch.setattr(kubejob_service, "upload_local_directory_to_minio", lambda **kwargs: calls.append(kwargs))

    run_postjob()  # no SystemExit

    assert calls == [{"local_path": str(job_env), "bucket_name": "clean", "minio_prefix": "job1/out"}]


def test_failed_upload_exits_nonzero(monkeypatch, job_env):
    def unreachable(**kwargs):
        raise ConnectionError("Failed to establish a new connection: [Errno 111] Connection refused")

    monkeypatch.setattr(kubejob_service, "upload_local_directory_to_minio", unreachable)

    with pytest.raises(SystemExit) as exited:
        run_postjob()

    assert exited.value.code == 1


def test_missing_output_directory_exits_nonzero(monkeypatch, job_env, tmp_path):
    """The init container creates the output directory, so a missing one means nothing
    was uploaded. upload_local_directory_to_minio reports that by returning False."""
    monkeypatch.setenv("JOB_OUTPUT_DIR", str(tmp_path / "does-not-exist"))
    monkeypatch.setattr(kubejob_service.Minio, "fput_object", lambda *a, **k: pytest.fail("nothing to upload"))

    with pytest.raises(SystemExit) as exited:
        run_postjob()

    assert exited.value.code == 1


@pytest.mark.parametrize("sentinel", ["error", "fail"])
def test_error_sentinel_uploads_then_exits_nonzero(monkeypatch, job_env, sentinel):
    """A tool that fails writes `error` and exits 0 so its log is uploaded; postjob.py
    then fails the Job so the watcher marks it 'error'.

    The novostoic commands write `fail` instead. Without this check, 119 completed
    novostoic jobs on prod (2025-03-27 to 2026-10-26) were marked 'completed' with a
    `fail` marker and no output.json."""
    (job_env / sentinel).touch()
    calls = []
    monkeypatch.setattr(kubejob_service, "upload_local_directory_to_minio", lambda **kwargs: calls.append(kwargs))

    with pytest.raises(SystemExit) as exited:
        run_postjob()

    assert exited.value.code == 1
    assert calls == [{"local_path": str(job_env), "bucket_name": "clean", "minio_prefix": "job1/out"}]


NOVOSTOIC_TYPES = ["novostoic-dgpredictor", "novostoic-enzrank", "novostoic-optstoic", "novostoic-pathways"]


@pytest.mark.parametrize("job_type", NOVOSTOIC_TYPES)
def test_failed_novostoic_command_writes_the_fail_sentinel(tmp_path, job_type):
    """The `fail` check above only works while the deployed commands write `fail`."""
    values = yaml.safe_load((POSTJOB.parent.parent / "chart" / "values.yaml").read_text())
    command = values["config"]["kubernetes_jobs"][job_type]["command"]
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python").write_text("#!/bin/sh\necho Traceback >&2\nexit 1\n")
    (bin_dir / "python").chmod(0o755)
    out = tmp_path / "out"
    out.mkdir()
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}", "JOB_OUTPUT_DIR": str(out)}

    subprocess.run(["bash", "-c", command], cwd=tmp_path, env=env, check=True)

    assert (out / "fail").exists()
    assert not (out / "success").exists()
