"""Tests for reporting 'queued' until a job's container starts (#143).

The tool runs in the pod's `job` init container (services/templates/job.tpl.yaml), so the
Job object cannot tell a pod waiting for the GPU node from one computing: both report
active=1, ready=0, with no condition, from pod creation until the pod finishes. Only the
pod shows the difference. So the Job watch now reports every unfinished job as 'queued',
and a second watch, over pods, moves it to 'processing' when `job` starts.
"""
import logging
import time
from types import SimpleNamespace

import pytest
from sqlmodel import Session

from models.enums import JobStatus, JobType
from models.sqlmodel.models import Job
from services import kubejob_service
from services.kubejob_service import KubeEventWatcher, _job_container_started


def _state(running=False, terminated=False):
    return SimpleNamespace(
        running=SimpleNamespace(started_at="t") if running else None,
        terminated=SimpleNamespace(exit_code=0) if terminated else None,
        waiting=None if (running or terminated) else SimpleNamespace(reason="PodInitializing"),
    )


def _cs(name, **state):
    return SimpleNamespace(name=name, state=_state(**state))


def _pod(job_id="j1", init=None, containers=None, has_status=True, labels=None):
    return SimpleNamespace(
        metadata=SimpleNamespace(
            labels={"type": "mmli-job", "jobId": job_id, "jobType": "somn"} if labels is None else labels,
            resource_version="7",
        ),
        status=SimpleNamespace(init_container_statuses=init, container_statuses=containers)
        if has_status else None,
    )


# The shapes below are the ones a real pod passes through, in order.
WAITING_FOR_NODE = _pod()
INPUTS_DOWNLOADING = _pod(init=[_cs("init", running=True), _cs("job")], containers=[_cs("post-job")])
COMPUTING = _pod(init=[_cs("init", terminated=True), _cs("job", running=True)], containers=[_cs("post-job")])
UPLOADING = _pod(init=[_cs("init", terminated=True), _cs("job", terminated=True)],
                 containers=[_cs("post-job", running=True)])


class TestJobContainerStarted:
    @pytest.mark.parametrize("pod", [WAITING_FOR_NODE, INPUTS_DOWNLOADING, _pod(has_status=False)],
                             ids=["unscheduled", "init-running", "no-status"])
    def test_not_started_until_the_job_container_runs(self, pod):
        assert _job_container_started(pod) is False

    @pytest.mark.parametrize("pod", [COMPUTING, UPLOADING], ids=["computing", "uploading"])
    def test_started_once_the_job_container_runs_or_has_run(self, pod):
        assert _job_container_started(pod) is True

    def test_a_job_container_still_pulling_its_image_has_not_started(self):
        pod = _pod(init=[_cs("init", terminated=True), _cs("job")], containers=[_cs("post-job")])
        assert _job_container_started(pod) is False


@pytest.fixture
def watcher(sync_db_url):
    from sqlmodel import create_engine

    instance = KubeEventWatcher()  # __init__ is a no-op under conftest
    instance.logger = logging.getLogger("test-watcher")
    instance.engine = create_engine(sync_db_url)
    return instance


def _seed(w, phase):
    with Session(w.engine) as session:
        session.add(Job(job_id="j1", type=JobType.SOMN, phase=phase,
                        time_created=0, user_agent="", deleted=0))
        session.commit()


def _phase(w):
    with Session(w.engine) as session:
        return session.get(Job, "j1").phase


class TestPodReconcile:
    def test_a_queued_job_moves_to_processing_when_its_container_starts(self, watcher):
        _seed(watcher, JobStatus.QUEUED)
        watcher._reconcile_pod(COMPUTING)
        assert _phase(watcher) == JobStatus.PROCESSING

    def test_a_queued_job_stays_queued_while_its_pod_waits(self, watcher):
        _seed(watcher, JobStatus.QUEUED)
        watcher._reconcile_pod(WAITING_FOR_NODE)
        watcher._reconcile_pod(INPUTS_DOWNLOADING)
        assert _phase(watcher) == JobStatus.QUEUED

    @pytest.mark.parametrize("phase", [JobStatus.COMPLETED, JobStatus.ERROR, JobStatus.CANCELED])
    def test_a_late_pod_event_never_reopens_a_finished_job(self, watcher, phase):
        """The pod watch and the Job watch are separate streams, so a pod event can
        arrive after the Job watch has already written the terminal phase."""
        _seed(watcher, phase)
        watcher._reconcile_pod(UPLOADING)
        assert _phase(watcher) == phase

    def test_a_pod_without_a_job_id_is_ignored(self, watcher):
        _seed(watcher, JobStatus.QUEUED)
        watcher._reconcile_pod(_pod(labels={"type": "mmli-job"}, init=COMPUTING.status.init_container_statuses))
        assert _phase(watcher) == JobStatus.QUEUED

    def test_a_pod_for_an_unknown_job_is_ignored(self, watcher):
        _seed(watcher, JobStatus.QUEUED)
        watcher._reconcile_pod(_pod(job_id="nope", init=COMPUTING.status.init_container_statuses))
        with Session(watcher.engine) as session:
            assert session.get(Job, "nope") is None
        assert _phase(watcher) == JobStatus.QUEUED


class _JobObject:
    def __init__(self, conditions=None):
        self.metadata = SimpleNamespace(namespace="mmli",
                                        labels={"type": "mmli-job", "jobId": "j1", "jobType": "somn"})
        self.spec = SimpleNamespace(completions=1)
        self.status = SimpleNamespace(conditions=conditions, succeeded=None, failed=None)


class TestJobReconcileKeepsProcessing:
    @pytest.fixture(autouse=True)
    def _no_pods(self, monkeypatch):
        monkeypatch.setattr(kubejob_service, "get_image_digest", lambda *a: None)

    def _reconcile(self, w, job):
        w._reconcile_job_phase(job, ignored_namespaces=[], required_labels={"type": "mmli-job"})

    def test_an_unfinished_job_event_does_not_move_processing_back_to_queued(self, watcher):
        """The Job reads the same before and after its container starts, so its events
        and the 600s sweep must not undo what the pod watch wrote."""
        _seed(watcher, JobStatus.PROCESSING)
        self._reconcile(watcher, _JobObject())
        assert _phase(watcher) == JobStatus.PROCESSING

    def test_a_new_job_is_queued(self, watcher):
        _seed(watcher, JobStatus.QUEUED)
        self._reconcile(watcher, _JobObject())
        assert _phase(watcher) == JobStatus.QUEUED

    def test_a_finished_job_still_completes_from_either_phase(self, watcher):
        _seed(watcher, JobStatus.QUEUED)
        self._reconcile(watcher, _JobObject(conditions=[SimpleNamespace(type="Complete", status="True")]))
        assert _phase(watcher) == JobStatus.COMPLETED


class _Stop(BaseException):
    """Ends run_pods()'s endless loop. BaseException, so `except Exception` lets it out."""


class _FakeCoreApi:
    def __init__(self, pods):
        self.pods = pods
        self.list_calls = []

    def list_namespaced_pod(self, **kwargs):
        self.list_calls.append(kwargs)
        return SimpleNamespace(metadata=SimpleNamespace(resource_version="5"), items=self.pods)


class TestRunPods:
    @pytest.fixture
    def harness(self, monkeypatch, watcher):
        monkeypatch.setattr(kubejob_service.time, "sleep", lambda s: None)
        seen = []
        monkeypatch.setattr(watcher, "_reconcile_pod", seen.append)
        return watcher, seen

    def _watch(self, monkeypatch, streams):
        calls = []

        class _Watch:
            def stream(self, func, **kwargs):
                calls.append(kwargs)
                if not streams:
                    raise _Stop()
                return streams.pop(0)

        monkeypatch.setattr(kubejob_service.watch, "Watch", _Watch)
        return calls

    def test_lists_and_reconciles_then_watches_from_the_list(self, harness, monkeypatch):
        w, seen = harness
        api = _FakeCoreApi([WAITING_FOR_NODE])
        monkeypatch.setattr(kubejob_service, "api_v1", api)
        stream_calls = self._watch(monkeypatch, [iter([{"type": "MODIFIED", "object": COMPUTING}])])

        with pytest.raises(_Stop):
            w.run_pods()

        assert seen[:2] == [WAITING_FOR_NODE, COMPUTING]
        assert api.list_calls[0]["label_selector"] == "type=mmli-job"
        assert api.list_calls[0]["_request_timeout"] == kubejob_service.KUBE_REQUEST_TIMEOUT
        assert stream_calls[0]["label_selector"] == "type=mmli-job"
        assert stream_calls[0]["resource_version"] == "5"
        assert stream_calls[0]["timeout_seconds"] == kubejob_service.WATCH_TIMEOUT_SECONDS
        assert stream_calls[0]["_request_timeout"] == kubejob_service.WATCH_REQUEST_TIMEOUT

    def test_a_failing_watch_reconnects_and_relists(self, harness, monkeypatch):
        w, seen = harness
        api = _FakeCoreApi([])
        monkeypatch.setattr(kubejob_service, "api_v1", api)

        def _broken():
            raise RuntimeError("connection reset")
            yield  # pragma: no cover

        self._watch(monkeypatch, [_broken()])

        with pytest.raises(_Stop):
            w.run_pods()

        assert len(api.list_calls) == 2

    def test_one_bad_pod_does_not_stop_the_others(self, harness, monkeypatch):
        w, seen = harness

        def _reconcile(pod):
            seen.append(pod)
            if pod is WAITING_FOR_NODE:
                raise ValueError("bad pod")

        monkeypatch.setattr(w, "_reconcile_pod", _reconcile)
        monkeypatch.setattr(kubejob_service, "api_v1", _FakeCoreApi([WAITING_FOR_NODE, COMPUTING]))
        self._watch(monkeypatch, [])

        with pytest.raises(_Stop):
            w.run_pods()

        assert seen == [WAITING_FOR_NODE, COMPUTING]

    def test_records_progress_for_the_health_check(self, harness, monkeypatch):
        w, _ = harness
        w.last_pod_progress = 0
        monkeypatch.setattr(kubejob_service, "api_v1", _FakeCoreApi([]))
        self._watch(monkeypatch, [])

        with pytest.raises(_Stop):
            w.run_pods()

        assert w.last_pod_progress > 0


class TestHealth:
    def _watcher(self, job_progress, pod_progress):
        w = KubeEventWatcher()
        w.last_progress = job_progress
        w.last_pod_progress = pod_progress
        w.is_alive = lambda: True
        return w

    def test_a_stalled_pod_watch_is_unhealthy(self):
        """Without it, a job would stay 'queued' for its whole run and nothing would say so."""
        stale = time.monotonic() - kubejob_service.WATCHER_STALL_SECONDS - 1
        assert self._watcher(time.monotonic(), stale).is_healthy() is False

    def test_both_watches_progressing_is_healthy(self):
        assert self._watcher(time.monotonic(), time.monotonic()).is_healthy() is True


class _Thread:
    def __init__(self, alive):
        self.alive = alive

    def is_alive(self):
        return self.alive


class TestIsAlive:
    """conftest replaces is_alive for every test, so this calls the real one."""

    @pytest.mark.parametrize("job_alive,pod_alive,expected", [
        (True, True, True), (True, False, False), (False, True, False),
    ])
    def test_both_threads_must_be_alive(self, job_alive, pod_alive, expected, real_watcher_is_alive):
        w = KubeEventWatcher()
        w.thread, w.pod_thread = _Thread(job_alive), _Thread(pod_alive)
        assert real_watcher_is_alive(w) is expected
