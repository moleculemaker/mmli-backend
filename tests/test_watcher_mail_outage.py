"""A mail server outage must not stall the watcher.

Notification sends run one at a time on the watcher thread, and now that a failed send
stays due, every pass retries every due job. With the mail server unreachable, each of
those sends waits out the SMTP timeout: 20 due jobs cost 600s of every pass before the
watch even starts, and phase updates fall behind by as much. So once a send finds the
server unreachable, the rest of that pass skips sending; the next pass tries once more.

A refused recipient is different: it fails only that job, so it must not stop the
others. Every smtplib exception subclasses OSError, which is why the distinction is
drawn explicitly in _mail_server_unavailable.
"""
import logging
import smtplib
import socket
from types import SimpleNamespace

import pytest
from sqlmodel import Session, create_engine

from models.enums import JobStatus, JobType
from models.sqlmodel.models import Job
from services import email_service, kubejob_service
from services.kubejob_service import KubeEventWatcher, _mail_server_unavailable


class _Stop(BaseException):
    pass


def _job_object(job_id):
    return SimpleNamespace(
        metadata=SimpleNamespace(namespace="mmli", resource_version="1",
                                 labels={"type": "mmli-job", "jobId": job_id, "jobType": "clean"}),
        spec=SimpleNamespace(completions=1),
        status=SimpleNamespace(conditions=[SimpleNamespace(type="Complete", status="True")],
                               succeeded=1, failed=None),
    )


class _SMTP:
    """smtplib.SMTP stand-in. `behavior(recipient)` is called per send: it may raise."""
    behavior = None
    connects = 0
    delivered = []

    def __init__(self, host, port=0, timeout=None):
        _SMTP.connects += 1
        if _SMTP.behavior == "unreachable":
            raise TimeoutError("timed out")

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def sendmail(self, from_addr, to_addr, msg):
        if callable(_SMTP.behavior):
            _SMTP.behavior(to_addr)
        _SMTP.delivered.append(to_addr)

    def quit(self):
        pass

    def close(self):
        pass


@pytest.fixture
def harness(monkeypatch, sync_db_url, real_watcher_run):
    _SMTP.behavior, _SMTP.connects, _SMTP.delivered = None, 0, []
    monkeypatch.setattr(email_service.smtplib, "SMTP", _SMTP)
    monkeypatch.setattr(kubejob_service, "get_image_digest", lambda *a: None)
    monkeypatch.setattr(kubejob_service.time, "sleep", lambda s: None)

    w = KubeEventWatcher()  # __init__ is a no-op under conftest
    w.logger = logging.getLogger("test-watcher")
    w.last_progress = 0
    w.engine = create_engine(sync_db_url)
    w.email_service = email_service.EmailService()

    def run(job_ids, passes):
        with Session(w.engine) as s:
            for i, job_id in enumerate(job_ids):
                s.add(Job(job_id=job_id, type=JobType.CLEAN, phase=JobStatus.COMPLETED,
                          email=f"user{i}@example.org", time_created=0, user_agent="", deleted=0))
            s.commit()
        lists = []

        class _Batch:
            def list_namespaced_job(self, **kw):
                lists.append(1)
                if len(lists) > passes:
                    raise _Stop()
                return SimpleNamespace(metadata=SimpleNamespace(resource_version="1"),
                                       items=[_job_object(j) for j in job_ids])

        class _Watch:
            def stream(self, func, **kw):
                return iter([])

        monkeypatch.setattr(kubejob_service, "api_batch_v1", _Batch())
        monkeypatch.setattr(kubejob_service.watch, "Watch", _Watch)
        with pytest.raises(_Stop):
            real_watcher_run(w)
        with Session(w.engine) as s:
            return {j: s.get(Job, j).notified_at for j in job_ids}

    return run


def test_an_unreachable_mail_server_costs_one_attempt_per_pass(harness):
    _SMTP.behavior = "unreachable"
    notified = harness([f"j{i}" for i in range(20)], passes=3)

    assert _SMTP.connects == 3
    assert all(v is None for v in notified.values())  # still due, retried next pass


def test_a_refused_recipient_does_not_stop_the_other_jobs(harness):
    def refuse_first(to_addr):
        if to_addr == "user0@example.org":
            raise smtplib.SMTPRecipientsRefused({to_addr: (550, b"no such user")})

    _SMTP.behavior = refuse_first
    notified = harness(["j0", "j1", "j2"], passes=1)

    assert notified["j0"] is None
    assert notified["j1"] is not None and notified["j2"] is not None


def test_the_mail_server_gets_a_fresh_attempt_each_pass(harness):
    """An outage that ends mid-run must not leave sending switched off."""
    attempts = {"n": 0}
    real_init = _SMTP.__init__

    def flaky_init(self, *a, **kw):
        attempts["n"] += 1
        if attempts["n"] == 1:
            raise TimeoutError("timed out")
        real_init(self, *a, **kw)

    _SMTP.__init__ = flaky_init
    try:
        notified = harness(["j0", "j1"], passes=2)
    finally:
        _SMTP.__init__ = real_init

    assert all(v is not None for v in notified.values())


@pytest.mark.parametrize("exc", [
    TimeoutError("timed out"),
    ConnectionRefusedError(111, "refused"),
    socket.gaierror(-2, "Name or service not known"),
    smtplib.SMTPServerDisconnected("gone"),
    smtplib.SMTPConnectError(421, b"busy"),
])
def test_transport_failures_mean_the_server_is_unavailable(exc):
    assert _mail_server_unavailable(exc) is True


@pytest.mark.parametrize("exc", [
    smtplib.SMTPRecipientsRefused({"a@b": (550, b"no")}),
    smtplib.SMTPSenderRefused(553, b"no", "from@b"),
    smtplib.SMTPDataError(554, b"rejected"),
    ValueError("not an smtp problem"),
])
def test_per_message_failures_do_not(exc):
    assert _mail_server_unavailable(exc) is False
