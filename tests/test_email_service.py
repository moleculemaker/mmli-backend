"""EmailService must time out, and must report a failed send.

Both matter because the caller is usually the kube watcher thread. Without a timeout an
unresponsive mail server blocks it indefinitely. And while send_email swallowed errors,
send_notification_email reported every send as a success, so the watcher stamped
job.notified_at on mail that was never sent and never retried it.
"""
import smtplib
from types import SimpleNamespace

import pytest

from models.enums import JobStatus, JobType
from services import email_service
from services.kubejob_service import KubeEventWatcher


class _FakeSMTP:
    instances = []
    fail = None

    def __init__(self, host, port=0, timeout=None):
        self.timeout = timeout
        _FakeSMTP.instances.append(self)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def sendmail(self, *args):
        if _FakeSMTP.fail:
            raise _FakeSMTP.fail


@pytest.fixture(autouse=True)
def fake_smtp(monkeypatch):
    _FakeSMTP.instances = []
    _FakeSMTP.fail = None
    monkeypatch.setattr(email_service.smtplib, "SMTP", _FakeSMTP)


def test_the_smtp_connection_has_a_timeout():
    email_service.EmailService().send_email("a@illinois.edu", "s", "b")
    assert _FakeSMTP.instances[0].timeout == email_service.SMTP_TIMEOUT_SECONDS


def test_a_failed_send_raises_instead_of_being_swallowed():
    _FakeSMTP.fail = smtplib.SMTPServerDisconnected("gone")
    with pytest.raises(smtplib.SMTPServerDisconnected):
        email_service.EmailService().send_email("a@illinois.edu", "s", "b")


def test_through_the_real_email_service_a_failed_send_stays_retryable():
    """The existing notifier tests stub EmailService, so they could not see that the
    real one never raised."""
    _FakeSMTP.fail = TimeoutError("timed out")
    w = KubeEventWatcher()  # __init__ is a no-op under conftest
    w.email_service = email_service.EmailService()
    job = SimpleNamespace(job_id="job-abc123", email="a@illinois.edu")

    assert w.send_notification_email("job-abc123", str(JobType.CLEANDB_MEPESM),
                                     job, JobStatus.COMPLETED) is False
