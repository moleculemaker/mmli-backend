"""Tests for when the watcher sends a completion email.

_reconcile_job_phase runs from two places: the live watch-event loop, and the reconcile
pass that lists every Job in the namespace on a 600s cadence. The second one is what
makes this worth pinning. A completed Job stays listed for ttlSecondsAfterFinished (12h,
app/cfg/config.yaml), so the sweep sees the same terminal job roughly 72 times.

Three gates have stood here, and the first two each failed in one direction:

  * The MinIO `{job_id}/email-sent` marker. Retryable, but unbounded, because it fails
    OPEN: check_file_exists returned False for every S3Error, and a stat_object against
    a missing bucket is a bodyless 404 that minio 7.1.17 synthesizes as NoSuchKey
    (minio/api.py:361) -- indistinguishable from an absent marker. mark_email_as_sent
    could not repair it either, since put_object does not create buckets. 72 copies.
  * The phase transition. Bounded, but once-only: the phase is committed before the
    send, so a send that fails is never reattempted and the mail is simply lost.
  * job.notified_at. Bounded AND retryable, written only when an email was actually
    handed over.

The third is at-least-once, not exactly-once, and these tests do not claim otherwise:
the send precedes the record of it, in a separate transaction from the phase, so a
database failure in between leaves the row due and the next pass sends again. What is
pinned below is the part that is actually guaranteed -- sent once when nothing fails,
retried until it succeeds when something does, and never sent for a job that should not
get one.
"""
import logging

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


class _Notifier:
    """Stands in for send_notification_email, recording calls and its verdict.

    The return value is the contract under test: True means an email was handed to the
    email service and the caller must record that, anything else means it was not and
    the job stays due. `fails` makes every call raise instead, which is what a MinIO
    connection error or a 5xx does -- neither is an S3Error, so neither is caught inside
    the notifier.
    """

    def __init__(self, sends=True, fails=False):
        self.sends = sends
        self.fails = fails
        self.calls = []

    def __call__(self, job_id, job_type, updated_job, new_phase):
        self.calls.append((job_id, new_phase))
        if self.fails:
            raise ConnectionError("MinIO unreachable")
        return self.sends


@pytest.fixture
def watcher(monkeypatch, sync_db_url):
    """A KubeEventWatcher wired to the test database.

    conftest stubs __init__ to a no-op, so the instance starts no thread and opens no
    connection. get_image_digest is stubbed out because the digest path would otherwise
    reach for the cluster; it is covered in test_provenance.py.
    """
    instance = kubejob_service.KubeEventWatcher()
    instance.logger = kubejob_service.log
    instance.engine = create_engine(sync_db_url)

    monkeypatch.setattr(kubejob_service, "get_image_digest", lambda job_type, job_id: None)

    instance.notifier = _Notifier()
    instance.send_notification_email = instance.notifier
    return instance


def _seed(watcher_instance, phase, notified_at=None):
    with Session(watcher_instance.engine) as session:
        session.add(Job(
            job_id="j1", type=JobType.SOMN, phase=phase, email="user@example.org",
            notified_at=notified_at, time_created=0, user_agent="", deleted=0,
        ))
        session.commit()


def _seed_no_email(watcher_instance, phase):
    with Session(watcher_instance.engine) as session:
        session.add(Job(
            job_id="j1", type=JobType.SOMN, phase=phase, email=None,
            time_created=0, user_agent="", deleted=0,
        ))
        session.commit()


def _reconcile(watcher_instance, conditions):
    watcher_instance._reconcile_job_phase(
        _FakeJobObject("j1", "somn", conditions=conditions),
        ignored_namespaces=[], required_labels={"type": "mmli-job"},
    )


def _row(watcher_instance):
    with Session(watcher_instance.engine) as session:
        return session.get(Job, "j1")


class TestItIsSentOnce:
    def test_completing_notifies(self, watcher):
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Complete")])

        assert watcher.notifier.calls == [("j1", JobStatus.COMPLETED)]

    def test_the_send_is_recorded_on_the_row(self, watcher):
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Complete")])

        assert _row(watcher).notified_at is not None

    def test_repeated_sweeps_over_a_terminal_job_notify_once_in_total(self, watcher):
        """72 passes is what a 12h TTL against a 600s cadence actually produces."""
        _seed(watcher, JobStatus.QUEUED)

        for _ in range(72):
            _reconcile(watcher, [_Condition("Complete")])

        assert watcher.notifier.calls == [("j1", JobStatus.COMPLETED)]

    def test_a_job_already_notified_is_not_notified_again(self, watcher):
        """The steady state after a restart: the row says done, so the notifier is not
        even consulted -- there is no MinIO round trip left to get this wrong."""
        _seed(watcher, JobStatus.COMPLETED, notified_at=1_700_000_000)

        _reconcile(watcher, [_Condition("Complete")])

        assert watcher.notifier.calls == []

    def test_failing_notifies(self, watcher):
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Failed")])

        assert watcher.notifier.calls == [("j1", JobStatus.ERROR)]


class TestItIsRetriedUntilItSucceeds:
    """What the phase-transition guard could not do. The email is no longer tied to the
    single pass that carried the job into its terminal phase."""

    def test_a_raising_notifier_leaves_the_job_due(self, watcher):
        watcher.notifier.fails = True
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Complete")])

        assert _row(watcher).notified_at is None

    def test_a_notifier_that_sent_nothing_leaves_the_job_due(self, watcher):
        """False, not an exception: a send_email that failed inside the notifier, which
        catches and logs its own. Recording that as notified would lose the mail."""
        watcher.notifier.sends = False
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Complete")])

        assert _row(watcher).notified_at is None

    def test_the_next_pass_sends_it(self, watcher):
        """The phase does not change between the two passes -- that is the point. Under
        the transition guard this second call never happened."""
        watcher.notifier.fails = True
        _seed(watcher, JobStatus.PROCESSING)
        _reconcile(watcher, [_Condition("Complete")])

        watcher.notifier.fails = False
        _reconcile(watcher, [_Condition("Complete")])

        assert len(watcher.notifier.calls) == 2
        assert _row(watcher).notified_at is not None

    def test_recovery_does_not_then_re_send(self, watcher):
        watcher.notifier.fails = True
        _seed(watcher, JobStatus.PROCESSING)
        _reconcile(watcher, [_Condition("Complete")])
        watcher.notifier.fails = False
        _reconcile(watcher, [_Condition("Complete")])

        for _ in range(10):
            _reconcile(watcher, [_Condition("Complete")])

        assert len(watcher.notifier.calls) == 2

    def test_the_phase_is_still_recorded_when_the_notifier_raises(self, watcher):
        """The stranded-job bug must not come back through the email path: a job the
        cluster says is finished reads as finished whatever the mail does."""
        watcher.notifier.fails = True
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Complete")])

        assert _row(watcher).phase == JobStatus.COMPLETED

    def test_the_failure_is_logged_as_retryable(self, watcher, caplog):
        """run()'s per-job handler would report 'Reconcile skipped a job', which is
        untrue -- the phase was written."""
        watcher.notifier.fails = True
        _seed(watcher, JobStatus.PROCESSING)

        with caplog.at_level(logging.ERROR):
            _reconcile(watcher, [_Condition("Complete")])

        assert any("will retry on the next pass" in r.getMessage() for r in caplog.records)

    def test_the_exception_does_not_escape_into_the_watch_loop(self, watcher):
        watcher.notifier.fails = True
        _seed(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Complete")])  # must not raise


class TestSomeJobsAreNeverNotified:
    def test_a_canceled_job_is_never_notified(self, watcher):
        """Cancellation returns before any of this. Worth asserting from out here: the
        gate moved, and that early return is still what suppresses the email."""
        _seed(watcher, JobStatus.CANCELED)

        _reconcile(watcher, [_Condition("Complete")])

        assert watcher.notifier.calls == []

    def test_a_running_job_is_not_even_consulted(self, watcher):
        """The notifier returns False for a non-terminal phase anyway, so this changes
        no outcome -- but it is reached on every event and every sweep while notified_at
        is null, and running the dispatch for a job that cannot produce an email is work
        with no result."""
        _seed(watcher, JobStatus.QUEUED)

        _reconcile(watcher, None)

        assert watcher.notifier.calls == []
        assert _row(watcher).notified_at is None

    def test_a_running_job_stays_due_for_its_real_completion(self, watcher):
        """The guard must not settle the row early: the completion still has to send."""
        _seed(watcher, JobStatus.QUEUED)
        _reconcile(watcher, None)

        _reconcile(watcher, [_Condition("Complete")])

        assert watcher.notifier.calls == [("j1", JobStatus.COMPLETED)]
        assert _row(watcher).notified_at is not None

    def test_a_terminal_job_with_no_address_is_not_consulted(self, watcher):
        """Same reasoning, and this one is terminal forever: without the guard the
        dispatch would run on every sweep for the whole 12h the job stays listed."""
        _seed_no_email(watcher, JobStatus.PROCESSING)

        _reconcile(watcher, [_Condition("Complete")])

        assert watcher.notifier.calls == []
        assert _row(watcher).notified_at is None
