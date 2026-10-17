import time

import requests

from config import get_logger

log = get_logger(__name__)

# How long a caller keeps retrying an unreachable backend before giving up. A rolling
# restart of the backend takes well under a minute; this leaves room for a slow one.
RETRY_BUDGET_SECONDS = 600
RETRY_INTERVAL_SECONDS = 10
REQUEST_TIMEOUT_SECONDS = 30


def request_with_retry(method, url, budget_seconds=RETRY_BUDGET_SECONDS,
                       interval_seconds=RETRY_INTERVAL_SECONDS, **kwargs):
    """requests.request(), retried while the backend is unreachable or answering 5xx.

    For callers that outlive a backend restart, such as coordinator.py. On 2026-10-17 a
    restart made one GET fail with ConnectionError, and that single failure killed 17
    ez-specificity jobs that had been waiting on their sub-jobs for up to two days.

    Only for requests that are safe to repeat. Not job creation: a POST that fails after
    the backend committed the row is retried into "duplicate, 200", for a job that was
    never launched, and its caller then waits on a row stuck at 'queued' forever.

    Only transport failures and 5xx are retried. A 4xx is the caller's mistake and is
    returned as-is, as is the last 5xx once the budget runs out; a transport failure
    after the budget re-raises.
    """
    kwargs.setdefault('timeout', REQUEST_TIMEOUT_SECONDS)
    deadline = time.monotonic() + budget_seconds
    while True:
        try:
            resp = requests.request(method, url, **kwargs)
            if resp.status_code < 500 or time.monotonic() >= deadline:
                return resp
            log.warning(f'{method} {url} answered {resp.status_code}; retrying in {interval_seconds}s')
        except (requests.ConnectionError, requests.Timeout) as e:
            if time.monotonic() >= deadline:
                raise
            log.warning(f'{method} {url} failed ({e}); retrying in {interval_seconds}s')
        time.sleep(interval_seconds)
