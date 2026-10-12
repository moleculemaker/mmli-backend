"""Tests for the optional memory-backed /dev/shm on job pods (`shmSizeLimit`).

Without it a job container gets the runtime's default 64Mi /dev/shm. EZSpecificity
inference overflows that once docked complexes are scored: the PyTorch DataLoader worker
dies with "Bus error ... out of shared memory".

These render the real job template through create_job and inspect the Job body that
would be submitted to the cluster.
"""
from config import app_config
from services import kubejob_service


def _submitted_job(monkeypatch, job_type, job_config):
    submitted = []

    class _FakeBatchV1Api:
        def create_namespaced_job(self, namespace, body):
            submitted.append(body)

    monkeypatch.setattr(kubejob_service, "api_batch_v1", _FakeBatchV1Api())
    monkeypatch.setattr(kubejob_service, "get_namespace", lambda: "test")
    monkeypatch.setitem(app_config["kubernetes_jobs"], job_type, job_config)

    kubejob_service.create_job(job_type=job_type, job_id="abc123", image_name="example/image:tag")
    assert len(submitted) == 1
    return submitted[0]


def _job_container(pod_spec):
    # The workload runs as the init container named "job"; "post-job" is the only container
    return next(c for c in pod_spec["initContainers"] if c["name"] == "job")


def test_shm_size_limit_mounts_memory_backed_dev_shm_on_job_container(monkeypatch):
    job_config = {**app_config["kubernetes_jobs"]["ezspec-inference"], "shmSizeLimit": "4Gi"}
    pod_spec = _submitted_job(monkeypatch, "ezspec-inference", job_config)["spec"]["template"]["spec"]

    dshm = [v for v in pod_spec["volumes"] if v["name"] == "dshm"]
    assert dshm == [{"name": "dshm", "emptyDir": {"medium": "Memory", "sizeLimit": "4Gi"}}]
    assert {"name": "dshm", "mountPath": "/dev/shm"} in _job_container(pod_spec)["volumeMounts"]


def test_no_shm_size_limit_leaves_pod_without_dshm(monkeypatch):
    job_config = {k: v for k, v in app_config["kubernetes_jobs"]["ezspec-inference"].items() if k != "shmSizeLimit"}
    pod_spec = _submitted_job(monkeypatch, "ezspec-inference", job_config)["spec"]["template"]["spec"]

    assert all(v["name"] != "dshm" for v in pod_spec["volumes"])
    assert all(m["mountPath"] != "/dev/shm" for m in _job_container(pod_spec)["volumeMounts"])
