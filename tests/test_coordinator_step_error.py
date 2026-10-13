"""coordinator.py when a step (subjob) ends in phase 'error'.

The step's error is logged, then the coordinator must exit with its own code. It used to
look up the step's Kubernetes Job first, with `namespace=None` (nothing sets NAMESPACE)
and the parent's job type in the Job name. That raised ApiValueError, which the outer
handler turned into the generic -500 with a traceback that hid the step's failure.

coordinator.py runs at import time, so it is executed with runpy against stubbed HTTP.
"""
import runpy
import pathlib

import pytest
import requests

from services import kubejob_service

COORDINATOR = pathlib.Path(__file__).resolve().parent.parent / "app" / "coordinator.py"


class _Response:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload


def test_failed_step_exits_with_step_error_code_without_kubernetes_lookup(monkeypatch, tmp_path):
    job_dir = tmp_path / "uws" / "jobs" / "ez-specificity" / "parent1"
    (job_dir / "in").mkdir(parents=True)
    (job_dir / "out").mkdir()
    monkeypatch.setenv("JOB_TYPE", "ez-specificity")
    monkeypatch.setenv("JOB_ID", "parent1")
    monkeypatch.setenv("JOB_INPUT_DIR", str(job_dir / "in"))
    monkeypatch.setenv("JOB_OUTPUT_DIR", str(job_dir / "out"))
    monkeypatch.setenv("MMLI_BACKEND_HOST", "http://backend")
    monkeypatch.setenv("EZSPEC_UNIDOCK_JOB_ID", "unidock1")
    monkeypatch.setenv("EZSPEC_INFERENCE_JOB_ID", "inference1")
    monkeypatch.delenv("NAMESPACE", raising=False)

    monkeypatch.setattr(kubejob_service, "upload_local_directory_to_minio", lambda **kwargs: None)
    monkeypatch.setattr(requests, "post", lambda url, **kwargs: _Response({"job_id": "unidock1"}))
    monkeypatch.setattr(requests, "get", lambda url, **kwargs: _Response([{"phase": "error"}]))

    def _read_namespaced_job(**kwargs):
        raise AssertionError("coordinator must not look up the failed step's Kubernetes Job")

    monkeypatch.setattr(kubejob_service.api_batch_v1, "read_namespaced_job", _read_namespaced_job, raising=False)

    with pytest.raises(SystemExit) as exited:
        runpy.run_path(str(COORDINATOR), run_name="__main__")

    assert exited.value.code == -800
