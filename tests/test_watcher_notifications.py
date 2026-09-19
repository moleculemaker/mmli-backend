"""Tests for when the watcher sends a completion email.

_reconcile_job_phase is run from two places: the live watch-event loop, and the
reconcile pass that lists every Job in the namespace on a 600s cadence. The second one
is what makes this worth pinning. A completed Job stays listed for
ttlSecondsAfterFinished (12h, app/cfg/config.yaml), so the sweep sees the same terminal
job roughly 72 times, and every sighting used to reach send_notification_email.

The MinIO marker was the only thing between that and 72 copies of "your result is
ready", and it fails open in both directions:

  * should_send_email treats every S3Error as "not sent yet"
    (MinIOService.check_file_exists returns False on all of them).
  * A stat_object against a bucket that does not exist is a bodyless 404, which minio
    7.1.17 synthesises as NoSuchKey -- indistinguishable from an absent marker
    (minio/api.py:361). mark_email_as_sent cannot repair that either, because put_object
    does not create buckets.

So the guard has to be the phase transition itself, not the marker.
"""
import pytest
from sqlmodel import Session, create_engine

from models.enums import JobStatus, JobType
from models.sqlmodel.models import Job
from services import kubejob_service


class _Condition:
    def __init__(self, type_, status="True"):
        self.type = type_
        self.status = status


class _FakeJobObject:
    """Minimal stand-in for a V1Job as the watcher sees it."""

    def __init__(self, job_id, job_type, conditions=None):
        self.metadata = type("meta", (), {
            "namespace": "mmli",
            "labels": {"type": "mmli-job", "jobId": job_id, "jobType": job_type},
        })()
        self.spec = type("spec", (), {"completions": 1})()
        self.status = type("status", (), {
            "conditions": conditions,
            "succeeded": None,
            "failed": None,
        })()


@pytest.fixture
def watcher(monkeypatch, sync_db_url):
    """A KubeEventWatcher wired to the test database, recording notification calls.

    conftest stubs __init__ to a no-op, so the instance starts no thread and opens no
    connection. get_image_digest is stubbed out because the digest path would otherwise
    reach for the cluster; it is covered in test_provenance.py.
    """
    instance = kubejob_service.KubeEventWatcher()
    instance.logger = kubejob_service.log
    instance.engine = create_engine(sync_db_url)

    monkeypatch.setattr(kubejob_service, "get_image_digest", lambda job_type, job_id: None)

    sent = []
    monkeypatch.setattr(
        instance,
        "send_notification_email",
        lambda job_id, job_type, updated_job, new_phase: sent.append((job_id, new_phase)),
    )
    instance.sent = sent
    return instance


def _seed(watcher_instance, phase):
    with Session(watcher_instance.engine) as session:
        session.add(Job(
            job_id="j1", type=JobType.SOMN, phase=phase, email="user@example.org",
            time_created=0, user_agent="", deleted=0,
        ))
        session.commit()


def _reconcile(watcher_instance, conditions):
    watcher_instance._reconcile_job_phase(
        _FakeJobObject("j1", "somn", conditions=conditions),
        ignored_namespaces=[], required_labels={"type": "mmli-job"},
    )


class TestNotificationIsBoundToTheTransition:
    def test_completing_notifies_once(self, watcher):
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Complete")])

        assert watcher.sent == [("j1", JobStatus.COMPLETED)]

    def test_a_job_already_completed_in_the_database_is_not_notified_again(self, watcher):
        """The sweep's steady state. This is the one that used to amplify: the row and
        the cluster already agree, so there is nothing to tell the user."""
        _seed(watcher, JobStatus.COMPLETED)

        _reconcile(watcher, [_Condition("Complete")])

        assert watcher.sent == []

    def test_repeated_sweeps_over_a_terminal_job_notify_once_in_total(self, watcher):
        """72 passes is what a 12h TTL against a 600s cadence actually produces."""
        _seed(watcher, JobStatus.QUEUED)

        for _ in range(72):
            _reconcile(watcher, [_Condition("Complete")])

        assert watcher.sent == [("j1", JobStatus.COMPLETED)]

    def test_failing_notifies_once(self, watcher):
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Failed")])

        assert watcher.sent == [("j1", JobStatus.ERROR)]

    def test_a_canceled_job_is_never_notified(self, watcher):
        """Cancellation returns before the phase logic; asserted here because the
        notification call moved, and that return is now what suppresses the email."""
        _seed(watcher, JobStatus.CANCELED)

        _reconcile(watcher, [_Condition("Complete")])

        assert watcher.sent == []
