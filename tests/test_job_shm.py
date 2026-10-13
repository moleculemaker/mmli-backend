"""Tests for the optional memory-backed /dev/shm on job pods (`shmSizeLimit`).

Without it a job container gets the runtime's default 64Mi /dev/shm. EZSpecificity
inference overflows that once docked complexes are scored: the PyTorch DataLoader worker
dies with "Bus error ... out of shared memory".

These render the real job template through create_job and inspect the Job body that
would be submitted to the cluster.
"""
import pathlib

import pytest
import yaml

from config import app_config
from services import kubejob_service

_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_BASE_VALUES = "chart/values.yaml"
DEPLOYED_VALUES_FILES = [
    _BASE_VALUES,
    "chart/values.prod.yaml",
    "chart/values.staging.yaml",
    "chart/values.mmli2.prod.yaml",
    "chart/values.mmli2.staging.yaml",
]


def _coalesce(base, overlay):
    """Merge the way Helm does: maps merge key by key, anything else the overlay wins."""
    merged = dict(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _coalesce(merged[key], value)
        else:
            merged[key] = value
    return merged


def _deployed_config(values_file):
    """The `config:` block a pod deployed with `values_file` would actually be handed."""
    base = yaml.safe_load((_REPO_ROOT / _BASE_VALUES).read_text())
    if values_file == _BASE_VALUES:
        return base["config"]
    overlay = yaml.safe_load((_REPO_ROOT / values_file).read_text())
    return _coalesce(base, overlay)["config"]


def _submitted_pod_spec(monkeypatch, job_type, config):
    submitted = []

    class _FakeBatchV1Api:
        def create_namespaced_job(self, namespace, body):
            submitted.append(body)

    monkeypatch.setattr(kubejob_service, "api_batch_v1", _FakeBatchV1Api())
    monkeypatch.setattr(kubejob_service, "get_namespace", lambda: "test")
    monkeypatch.setattr(kubejob_service, "app_config", config)

    kubejob_service.create_job(job_type=job_type, job_id="abc123", image_name="example/image:tag")
    assert len(submitted) == 1
    return submitted[0]["spec"]["template"]["spec"]


def _job_container(pod_spec):
    # The workload runs as the init container named "job"; "post-job" is the only container
    return next(c for c in pod_spec["initContainers"] if c["name"] == "job")


@pytest.mark.parametrize("values_file", DEPLOYED_VALUES_FILES)
def test_deployed_ezspec_inference_gets_memory_backed_dev_shm(monkeypatch, values_file):
    pod_spec = _submitted_pod_spec(monkeypatch, "ezspec-inference", _deployed_config(values_file))

    dshm = [v for v in pod_spec["volumes"] if v["name"] == "dshm"]
    assert dshm == [{"name": "dshm", "emptyDir": {"medium": "Memory", "sizeLimit": "4Gi"}}]
    assert {"name": "dshm", "mountPath": "/dev/shm"} in _job_container(pod_spec)["volumeMounts"]


def test_job_type_without_shm_size_limit_has_no_dshm(monkeypatch):
    # app/cfg/config.yaml does not set shmSizeLimit for any job type
    assert "shmSizeLimit" not in app_config["kubernetes_jobs"]["ezspec-inference"]
    pod_spec = _submitted_pod_spec(monkeypatch, "ezspec-inference", app_config)

    assert all(v["name"] != "dshm" for v in pod_spec["volumes"])
    assert all(m["mountPath"] != "/dev/shm" for m in _job_container(pod_spec)["volumeMounts"])
