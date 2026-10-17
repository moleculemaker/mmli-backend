"""A message the server accepted must not be reported as a failed send.

`with smtplib.SMTP(...)` closes with SMTP.__exit__, which sends QUIT and raises
SMTPResponseException if the reply is not 221. That happens after DATA was answered 250,
i.e. after the mail was delivered. Now that send_email raises, send_notification_email
returns False, job.notified_at stays null, and the watcher sends the same mail again on
every pass.

A timeout on QUIT does not raise: getreply turns the OSError into SMTPServerDisconnected,
which __exit__ swallows (smtplib.py in Python 3.10). The timeout case is pinned below so
that stays true.

test_email_service.py replaces smtplib.SMTP wholesale, so it cannot see __exit__. These
tests use the real smtplib against a minimal SMTP server on a local socket.
"""
import socket
import threading
from types import SimpleNamespace

import pytest

from models.enums import JobStatus, JobType
from services import email_service
from services.kubejob_service import KubeEventWatcher


class _SMTPServer:
    """Accepts every message. `quit_reply` is sent in answer to QUIT; None sends nothing."""

    def __init__(self, quit_reply):
        self.quit_reply = quit_reply
        self.delivered = 0
        self._sock = socket.socket()
        self._sock.bind(("127.0.0.1", 0))
        self._sock.listen()
        self.port = self._sock.getsockname()[1]
        threading.Thread(target=self._serve, daemon=True).start()

    def _serve(self):
        while True:
            try:
                conn, _ = self._sock.accept()
            except OSError:
                return
            threading.Thread(target=self._session, args=(conn,), daemon=True).start()

    def _session(self, conn):
        f = conn.makefile("rb")
        conn.sendall(b"220 test ESMTP\r\n")
        try:
            for line in f:
                cmd = line.strip().upper()
                if cmd == b"DATA":
                    conn.sendall(b"354 go ahead\r\n")
                    for body in f:
                        if body == b".\r\n":
                            break
                    self.delivered += 1
                    conn.sendall(b"250 queued\r\n")
                elif cmd == b"QUIT":
                    if self.quit_reply is not None:
                        conn.sendall(self.quit_reply)
                    else:
                        threading.Event().wait(5)  # outlast the client's timeout
                    return
                else:
                    conn.sendall(b"250 ok\r\n")
        finally:
            conn.close()

    def close(self):
        self._sock.close()


def _serve(monkeypatch, quit_reply):
    server = _SMTPServer(quit_reply)
    real_smtp = email_service.smtplib.SMTP

    def smtp_on_test_port(host, port=0, **kwargs):
        return real_smtp("127.0.0.1", server.port, **kwargs)

    monkeypatch.setattr(email_service.smtplib, "SMTP", smtp_on_test_port)
    monkeypatch.setattr(email_service, "SMTP_TIMEOUT_SECONDS", 1)
    return server


def _notify():
    w = KubeEventWatcher()  # __init__ is a no-op under conftest
    w.email_service = email_service.EmailService()
    job = SimpleNamespace(job_id="job-abc123", email="a@illinois.edu")
    return w.send_notification_email("job-abc123", str(JobType.CLEANDB_MEPESM),
                                     job, JobStatus.COMPLETED)


QUIT_REPLIES = {
    # Control: the harness itself reports a clean session as delivered.
    "221": b"221 bye\r\n",
    "non_221": b"421 closing\r\n",
    "timeout": None,
}


@pytest.mark.parametrize("quit_reply", QUIT_REPLIES.values(), ids=QUIT_REPLIES.keys())
def test_a_delivered_message_is_reported_as_sent(monkeypatch, quit_reply):
    server = _serve(monkeypatch, quit_reply)
    try:
        sent = _notify()
    finally:
        server.close()
    assert server.delivered == 1
    assert sent is True


def test_a_delivered_message_is_not_sent_again_on_the_next_pass(monkeypatch):
    """What the watcher does with the False: the job stays due, and each pass resends."""
    server = _serve(monkeypatch, b"421 closing\r\n")
    try:
        passes = 0
        while passes < 3 and not _notify():
            passes += 1
    finally:
        server.close()
    assert server.delivered == 1
