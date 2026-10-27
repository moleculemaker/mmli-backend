"""request_with_retry lets coordinator.py ride out a backend restart.

On 2026-10-17 one ConnectionError during a restart ended 17 ez-specificity coordinators.
"""
from types import SimpleNamespace

import pytest
import requests

from services import backend_retry


@pytest.fixture
def responses(monkeypatch):
    """Queue of outcomes for requests.request: an int is a status code, an exception is
    raised."""
    queue = []
    calls = []

    def _request(method, url, **kwargs):
        calls.append(kwargs)
        outcome = queue.pop(0)
        if isinstance(outcome, BaseException):
            raise outcome
        return SimpleNamespace(status_code=outcome)

    monkeypatch.setattr(backend_retry.requests, "request", _request)
    monkeypatch.setattr(backend_retry.time, "sleep", lambda s: None)
    return SimpleNamespace(queue=queue, calls=calls)


def test_rides_out_a_backend_restart(responses):
    responses.queue += [requests.ConnectionError("refused"), 503, 200]
    resp = backend_retry.request_with_retry("GET", "http://backend/x/jobs/1")
    assert resp.status_code == 200
    assert len(responses.calls) == 3


def test_rides_out_a_reset_partway_through_the_response_body(responses):
    responses.queue += [requests.exceptions.ChunkedEncodingError("Connection broken"), 200]
    resp = backend_retry.request_with_retry("GET", "http://backend/x/jobs/1")
    assert resp.status_code == 200


def test_every_attempt_has_a_timeout(responses):
    responses.queue += [200]
    backend_retry.request_with_retry("GET", "http://backend/x/jobs/1")
    assert responses.calls[0]["timeout"] == backend_retry.REQUEST_TIMEOUT_SECONDS


def test_a_4xx_is_returned_without_retrying(responses):
    responses.queue += [404]
    assert backend_retry.request_with_retry("GET", "http://backend/x/jobs/1").status_code == 404
    assert len(responses.calls) == 1


def test_gives_up_once_the_budget_is_spent(responses):
    responses.queue += [requests.ConnectionError("refused")]
    with pytest.raises(requests.ConnectionError):
        backend_retry.request_with_retry("GET", "http://backend/x/jobs/1", budget_seconds=0)
