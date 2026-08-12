"""Tests for job provenance: image digests, the parent/subjob link, and cancellation.

The watcher is exercised directly rather than through HTTP, because none of this is
reachable from an endpoint yet -- it is the data layer the versioned API will read.
"""
import json

import pytest
from sqlmodel import Session

from models.enums import JobStatus, JobType
from models.sqlmodel.models import Job
from services import kubejob_service


class _FakeContainerStatus:
    def __init__(self, image_id):
        self.image_id = image_id


class _FakePodStatus:
    def __init__(self, container_statuses):
        self.container_statuses = container_statuses


class _FakePod:
    def __init__(self, image_id):
        self.status = _FakePodStatus([_FakeContainerStatus(image_id)])


class _FakePodList:
    def __init__(self, pods):
        self.items = pods


@pytest.fixture
def fake_pods(monkeypatch):
    """Control what list_namespaced_pod returns, and count how often it is asked.

    The call count matters because the reconcile pass now runs every 600s over every
    listed Job, so "does the watcher stop asking" is a behavior worth asserting and not
    just an efficiency note.
    """
    state = {"pods": [], "calls": 0}

    class _FakeCoreApi:
        def list_namespaced_pod(self, namespace, label_selector=None):
            state["calls"] += 1
            return _FakePodList(state["pods"])

    monkeypatch.setattr(kubejob_service, "api_v1", _FakeCoreApi())
    return state


class TestImageDigestExtraction:
    def test_reads_a_containerd_style_image_id(self, fake_pods):
        fake_pods["pods"] = [_FakePod("moleculemaker/novostoic@sha256:9f2c4a")]

        assert kubejob_service.get_image_digest("novostoic-optstoic", "j1") == \
            "moleculemaker/novostoic@sha256:9f2c4a"

    def test_strips_the_docker_pullable_prefix(self, fake_pods):
        fake_pods["pods"] = [_FakePod("docker-pullable://moleculemaker/somn@sha256:abc123")]

        assert kubejob_service.get_image_digest("somn", "j1") == \
            "moleculemaker/somn@sha256:abc123"

    def test_ignores_an_image_id_with_no_digest(self, fake_pods):
        """Some runtimes report a bare tag before the image is pulled."""
        fake_pods["pods"] = [_FakePod("moleculemaker/somn:latest")]

        assert kubejob_service.get_image_digest("somn", "j1") is None

    def test_returns_none_when_no_pod_exists(self, fake_pods):
        fake_pods["pods"] = []

        assert kubejob_service.get_image_digest("somn", "j1") is None

    def test_returns_none_when_the_cluster_call_fails(self, monkeypatch):
        """Provenance is worth recording but never worth failing a job over."""
        class _Broken:
            def list_namespaced_pod(self, *a, **kw):
                raise RuntimeError("cluster unreachable")

        monkeypatch.setattr(kubejob_service, "api_v1", _Broken())

        assert kubejob_service.get_image_digest("somn", "j1") is None


class TestParentJobLink:
    def test_a_subjob_records_its_parent(self, client):
        parent_id = "parent-1"
        job_info = json.dumps({
            "parent_job_id": parent_id,
            "ezspec_unidock_job_id": "sub-dock",
            "ezspec_inference_job_id": "sub-infer",
        })

        resp = client.post("/ezspec-unidock/jobs", json={"job_info": job_info})

        assert resp.status_code == 201
        row = client.get("/ezspec-unidock/jobs/sub-dock").json()[0]
        assert row["parent_job_id"] == parent_id

    def test_the_subjob_id_comes_from_job_info_not_the_request(self, client):
        """coordinator.py sends the subjob id inside job_info; the row must use it."""
        job_info = json.dumps({
            "parent_job_id": "parent-1",
            "ezspec_unidock_job_id": "sub-dock",
            "ezspec_inference_job_id": "sub-infer",
        })

        client.post("/ezspec-inference/jobs", json={"job_info": job_info})

        assert client.get("/ezspec-inference/jobs/sub-infer").json()[0]["parent_job_id"] == "parent-1"

    def test_an_ordinary_job_has_no_parent(self, client):
        client.post("/somn/jobs", json={"job_id": "plain", "job_info": "[]"})

        assert client.get("/somn/jobs/plain").json()[0]["parent_job_id"] is None


class _FakeJobObject:
    """Minimal stand-in for a V1Job as the watcher sees it.

    conditions defaults to None, which _derive_phase reads as "still processing" -- as it
    does an empty list, since neither carries a terminal condition. Phase derivation is
    covered in detail in test_kube_phase.py; these tests only need a job object that
    lands on a given phase.
    """

    def __init__(self, job_id, job_type, conditions=None, succeeded=None, failed=None):
        self.metadata = type("meta", (), {
            "namespace": "mmli",
            "labels": {"type": "mmli-job", "jobId": job_id, "jobType": job_type},
        })()
        self.spec = type("spec", (), {"completions": 1})()
        self.status = type("status", (), {
            "conditions": conditions,
            "succeeded": succeeded,
            "failed": failed,
        })()


class _Condition:
    def __init__(self, type_, status="True"):
        self.type = type_
        self.status = status


@pytest.fixture
def watcher(monkeypatch, fake_pods, sync_db_url):
    """A KubeEventWatcher wired to the test database, with email suppressed.

    conftest stubs __init__ to a no-op, so the instance starts no thread and opens no
    connection; everything it needs is supplied here.
    """
    from sqlmodel import create_engine

    instance = kubejob_service.KubeEventWatcher()
    instance.logger = kubejob_service.log
    instance.engine = create_engine(sync_db_url)
    monkeypatch.setattr(instance, "send_notification_email", lambda *a, **kw: None)
    return instance


class TestCancellationIsTerminal:
    def _seed(self, watcher_instance, phase):
        with Session(watcher_instance.engine) as session:
            session.add(Job(
                job_id="j1", type=JobType.SOMN, phase=phase,
                time_created=0, user_agent="", deleted=0,
            ))
            session.commit()

    def _phase(self, watcher_instance):
        with Session(watcher_instance.engine) as session:
            return session.get(Job, "j1").phase

    def test_a_canceled_job_is_not_resurrected_by_a_stale_event(self, watcher):
        """Deleting a Job produces events that derive to 'processing'.

        Kubernetes never reports a terminal condition for a Job that was deleted, so
        without an explicit guard the watcher would move a job we deliberately stopped
        back to 'processing' and it would poll forever.
        """
        self._seed(watcher, JobStatus.CANCELED)

        watcher._reconcile_job_phase(
            _FakeJobObject("j1", "somn"), ignored_namespaces=[], required_labels={"type": "mmli-job"},
        )

        assert self._phase(watcher) == JobStatus.CANCELED

    def test_an_ordinary_job_still_advances(self, watcher):
        self._seed(watcher, JobStatus.QUEUED)

        watcher._reconcile_job_phase(
            _FakeJobObject("j1", "somn", conditions=[_Condition("SuccessCriteriaMet")]),
            ignored_namespaces=[], required_labels={"type": "mmli-job"},
        )

        assert self._phase(watcher) == JobStatus.COMPLETED


class TestDigestIsNotRetriedForever:
    """The digest read costs a list_namespaced_pod call, and the reconcile pass now
    revisits every listed Job every 600s rather than only on reconnect.

    A job that reached a terminal phase without yielding a digest never will: its pod is
    gone, and nothing about it will change. Re-reading it would cost one cluster call per
    job per pass, forever, for a value that cannot arrive.
    """

    def _seed(self, watcher_instance, phase):
        with Session(watcher_instance.engine) as session:
            session.add(Job(
                job_id="j1", type=JobType.SOMN, phase=phase,
                time_created=0, user_agent="", deleted=0,
            ))
            session.commit()

    def _reconcile(self, watcher_instance, conditions):
        watcher_instance._reconcile_job_phase(
            _FakeJobObject("j1", "somn", conditions=conditions),
            ignored_namespaces=[], required_labels={"type": "mmli-job"},
        )

    def test_a_settled_terminal_job_is_not_asked_again(self, watcher, fake_pods):
        """Already 'completed' in the DB and still completed in the cluster: no phase
        change, nothing left to learn, so the pod must not be listed."""
        self._seed(watcher, JobStatus.COMPLETED)
        fake_pods["calls"] = 0

        self._reconcile(watcher, [_Condition("Complete")])

        assert fake_pods["calls"] == 0

    def test_the_transition_into_a_terminal_phase_still_captures_it(self, watcher, fake_pods):
        """The completion event is usually the last moment the pod exists. Skipping it
        would lose the digest for exactly the short jobs whose whole lifecycle fits
        between two reconcile passes."""
        fake_pods["pods"] = [_FakePod("ianrinehart/somn@sha256:deadbeef")]
        self._seed(watcher, JobStatus.PROCESSING)

        self._reconcile(watcher, [_Condition("Complete")])

        with Session(watcher.engine) as session:
            assert session.get(Job, "j1").image_digest == "ianrinehart/somn@sha256:deadbeef"

    def test_a_running_job_is_still_retried(self, watcher, fake_pods):
        """While the job runs the digest may simply not be readable yet -- the image has
        to be pulled first -- so a null result must not stop later attempts."""
        self._seed(watcher, JobStatus.PROCESSING)
        fake_pods["calls"] = 0

        self._reconcile(watcher, None)
        self._reconcile(watcher, None)

        assert fake_pods["calls"] == 2


class TestDigestCapture:
    def test_the_watcher_records_the_digest(self, watcher, fake_pods):
        fake_pods["pods"] = [_FakePod("ianrinehart/somn@sha256:deadbeef")]
        with Session(watcher.engine) as session:
            session.add(Job(job_id="j1", type=JobType.SOMN, phase=JobStatus.QUEUED,
                            time_created=0, user_agent="", deleted=0))
            session.commit()

        watcher._reconcile_job_phase(
            _FakeJobObject("j1", "somn"), ignored_namespaces=[], required_labels={"type": "mmli-job"},
        )

        with Session(watcher.engine) as session:
            assert session.get(Job, "j1").image_digest == "ianrinehart/somn@sha256:deadbeef"

    def test_a_missing_digest_is_left_null_rather_than_guessed(self, watcher, fake_pods):
        fake_pods["pods"] = []
        with Session(watcher.engine) as session:
            session.add(Job(job_id="j1", type=JobType.SOMN, phase=JobStatus.QUEUED,
                            time_created=0, user_agent="", deleted=0))
            session.commit()

        watcher._reconcile_job_phase(
            _FakeJobObject("j1", "somn"), ignored_namespaces=[], required_labels={"type": "mmli-job"},
        )

        with Session(watcher.engine) as session:
            assert session.get(Job, "j1").image_digest is None
