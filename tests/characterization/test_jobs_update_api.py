"""The legacy PUT and PATCH job endpoints are gone.

Both wrote `phase`, `time_start` and `time_end` with no authentication, to a job looked
up by the `job_id` in the request body. The watcher is the only thing that should set
those fields. Neither had a working caller: the frontends' generated client builds the
PATCH URL as `/{job_type}/jobs//`, which never reached the handler. See #138.
"""
from models.enums import JobStatus

DEFAULTS = "defaults"


def _create(client, job_id="j1", run_id="1"):
    return client.post(f"/{DEFAULTS}/jobs", json={"job_id": job_id, "run_id": run_id})


def _phase(client):
    return client.get(f"/{DEFAULTS}/jobs/j1/1").json()["phase"]


class TestPatchJob:
    def test_patch_is_not_allowed(self, client):
        _create(client)
        before = _phase(client)

        resp = client.patch(
            f"/{DEFAULTS}/jobs/j1/1",
            json={"job_id": "j1", "run_id": 1, "phase": JobStatus.CANCELED},
        )

        assert resp.status_code == 405
        assert _phase(client) == before


class TestPutJob:
    def test_put_is_not_allowed(self, client):
        _create(client)
        before = _phase(client)

        resp = client.put(
            f"/{DEFAULTS}/jobs/j1/1",
            json={"job_id": "j1", "run_id": "1", "phase": JobStatus.CANCELED, "time_start": 1, "time_end": 2},
        )

        assert resp.status_code == 405
        assert _phase(client) == before
