"""Tests for the defenses against the 2026-10-15 watcher stall.

The watcher thread stopped making progress for two days while staying alive, so every
job's phase stayed at 'queued'. Nothing noticed: no call it made had a client-side
timeout, nothing checked whether it was still working, and the pod had no liveness
probe. These tests pin each piece of the fix:

  * the Kubernetes calls the watcher makes carry a client-side timeout, and a read
    timeout on the watch reconnects instead of ending the thread;
  * is_healthy() and /healthz report a dead OR stalled watcher.
"""
import logging
import time
from types import SimpleNamespace

import pytest
import urllib3

import main
from models.enums import JobType
from services import kubejob_service
from services.kubejob_service import KubeEventWatcher


class _Stop(BaseException):
    """Ends run()'s endless loop. BaseException, so run()'s `except Exception` lets it out."""


def _watcher(alive=True, last_progress=None):
    w = KubeEventWatcher()  # __init__ is a no-op under conftest
    w.logger = logging.getLogger("test-watcher")
    w.last_progress = time.monotonic() if last_progress is None else last_progress
    w.is_alive = lambda: alive
    return w


class TestIsHealthy:
    def test_a_running_watcher_that_made_progress_recently_is_healthy(self):
        assert _watcher().is_healthy() is True

    def test_a_running_but_stalled_watcher_is_not_healthy(self):
        """The 2026-10-15 failure: thread alive, blocked on a read, no progress."""
        stale = time.monotonic() - kubejob_service.WATCHER_STALL_SECONDS - 1
        assert _watcher(last_progress=stale).is_healthy() is False

    def test_a_dead_watcher_is_not_healthy(self):
        assert _watcher(alive=False).is_healthy() is False

    def test_the_stall_threshold_outlasts_a_quiet_watch_cycle(self):
        """A healthy watcher with no Job activity only records progress on reconnect,
        once per watch cycle; the threshold must not fire inside one of those."""
        assert kubejob_service.WATCHER_STALL_SECONDS > kubejob_service.WATCH_REQUEST_TIMEOUT[1]


class TestHealthz:
    def test_ok_while_the_watcher_is_healthy(self, client, monkeypatch):
        monkeypatch.setattr(main.watcher, "is_healthy", lambda: True, raising=False)
        assert client.get("/healthz").status_code == 200

    def test_503_once_the_watcher_stalls(self, client, monkeypatch):
        monkeypatch.setattr(main.watcher, "is_healthy", lambda: False, raising=False)
        assert client.get("/healthz").status_code == 503


class _FakeBatchApi:
    def __init__(self):
        self.list_calls = []

    def list_namespaced_job(self, **kwargs):
        self.list_calls.append(kwargs)
        return SimpleNamespace(metadata=SimpleNamespace(resource_version="1"), items=[])


def _fake_watch(streams):
    """A watch.Watch stand-in whose successive stream() calls yield from `streams`.

    Each entry is an iterator; once they run out, stream() raises _Stop to end run().
    """
    calls = []

    class _Watch:
        def stream(self, func, **kwargs):
            calls.append(kwargs)
            if not streams:
                raise _Stop()
            return streams.pop(0)

    return _Watch, calls


def _raising(exc):
    raise exc
    yield  # pragma: no cover - makes this a generator


class TestRunUsesClientTimeouts:
    @pytest.fixture
    def harness(self, monkeypatch):
        batch = _FakeBatchApi()
        monkeypatch.setattr(kubejob_service, "api_batch_v1", batch)
        monkeypatch.setattr(kubejob_service.time, "sleep", lambda s: None)
        return batch

    def test_the_list_and_the_watch_both_carry_a_client_timeout(self, harness, monkeypatch, real_watcher_run):
        watch_cls, stream_calls = _fake_watch([])
        monkeypatch.setattr(kubejob_service.watch, "Watch", watch_cls)

        with pytest.raises(_Stop):
            real_watcher_run(_watcher())

        assert harness.list_calls[0]["_request_timeout"] == kubejob_service.KUBE_REQUEST_TIMEOUT
        assert stream_calls[0]["_request_timeout"] == kubejob_service.WATCH_REQUEST_TIMEOUT
        assert stream_calls[0]["timeout_seconds"] == kubejob_service.WATCH_TIMEOUT_SECONDS

    def test_the_watch_read_timeout_outlasts_the_server_side_timeout(self):
        """Shorter, and a healthy stream with no Job changes would be torn down."""
        connect, read = kubejob_service.WATCH_REQUEST_TIMEOUT
        assert read > kubejob_service.WATCH_TIMEOUT_SECONDS

    def test_a_read_timeout_on_the_watch_reconnects_instead_of_killing_the_thread(
            self, harness, monkeypatch, real_watcher_run):
        timeout = urllib3.exceptions.ReadTimeoutError(None, "/apis/batch/v1/jobs", "Read timed out.")
        watch_cls, stream_calls = _fake_watch([_raising(timeout)])
        monkeypatch.setattr(kubejob_service.watch, "Watch", watch_cls)

        with pytest.raises(_Stop):
            real_watcher_run(_watcher())

        assert len(harness.list_calls) == 2  # re-listed (and reconciled) after the timeout

    def test_each_event_records_progress(self, harness, monkeypatch, real_watcher_run):
        event = {"object": SimpleNamespace(metadata=SimpleNamespace(resource_version="2"))}
        watch_cls, _ = _fake_watch([iter([event])])
        monkeypatch.setattr(kubejob_service.watch, "Watch", watch_cls)
        w = _watcher(last_progress=0)
        seen = []
        w._reconcile_job_phase = lambda obj, *a: seen.append(w.last_progress)

        with pytest.raises(_Stop):
            real_watcher_run(w)

        assert seen and w.last_progress > seen[0] > 0


def test_a_list_that_always_fails_is_reported_as_stalled(monkeypatch, real_watcher_run):
    """The loop retries every ~3s; progress must not be recorded for an attempt that
    failed, or /healthz stays green while no phase is ever updated."""
    class _Batch:
        def list_namespaced_job(self, **kwargs):
            raise kubejob_service.ApiException(status=401, reason="Unauthorized")

    clock = {"t": 1000.0, "sleeps": 0}

    def _sleep(seconds):  # advance a simulated clock by the real sleep durations
        clock["t"] += seconds
        clock["sleeps"] += 1
        if clock["t"] > 1000.0 + kubejob_service.WATCHER_STALL_SECONDS + 60:
            raise _Stop()

    monkeypatch.setattr(kubejob_service, "api_batch_v1", _Batch())
    monkeypatch.setattr(kubejob_service.time, "monotonic", lambda: clock["t"])
    monkeypatch.setattr(kubejob_service.time, "sleep", _sleep)
    w = _watcher(last_progress=1000.0)

    with pytest.raises(_Stop):
        real_watcher_run(w)

    assert clock["sleeps"] > 100  # hundreds of failed attempts, not one long sleep
    assert w.is_healthy() is False


def test_get_image_digest_carries_a_client_timeout(monkeypatch):
    calls = []

    class _CoreApi:
        def list_namespaced_pod(self, **kwargs):
            calls.append(kwargs)
            return SimpleNamespace(items=[])

    monkeypatch.setattr(kubejob_service, "api_v1", _CoreApi())
    kubejob_service.get_image_digest(str(JobType.CLEAN), "job-abc123")
    assert calls[0]["_request_timeout"] == kubejob_service.KUBE_REQUEST_TIMEOUT
