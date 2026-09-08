"""Characterization tests for the PUT and PATCH job endpoints.

These were missing from the original characterization pass. They matter because
`JobUpdate.run_id` is declared `int` while `Job.run_id` is `Optional[str]`, and how
that mismatch is coerced is exactly the kind of thing a Pydantic major-version upgrade
changes.
"""

DEFAULTS = "defaults"


def _create(client, job_id="j1", run_id="1"):
    return client.post(f"/{DEFAULTS}/jobs", json={"job_id": job_id, "run_id": run_id})


class TestPatchJob:
    def test_patch_updates_the_phase(self, client):
        _create(client)

        resp = client.patch(
            f"/{DEFAULTS}/jobs/j1/1",
            json={"job_id": "j1", "run_id": 1, "phase": "completed"},
        )

        assert resp.status_code == 200
        assert resp.json()["phase"] == "completed"

    def test_patch_accepts_a_numeric_string_run_id(self, client):
        """run_id is declared int here but str on the Job model."""
        _create(client)

        resp = client.patch(
            f"/{DEFAULTS}/jobs/j1/1",
            json={"job_id": "j1", "run_id": "1", "phase": "completed"},
        )

        assert resp.status_code == 200

    def test_patch_updates_timestamps(self, client):
        _create(client)

        resp = client.patch(
            f"/{DEFAULTS}/jobs/j1/1",
            json={"job_id": "j1", "run_id": 1, "time_start": 111, "time_end": 222},
        )

        assert resp.status_code == 200
        assert resp.json()["time_start"] == 111
        assert resp.json()["time_end"] == 222

    def test_patch_leaves_omitted_fields_alone(self, client):
        _create(client)
        client.patch(f"/{DEFAULTS}/jobs/j1/1", json={"job_id": "j1", "run_id": 1, "phase": "processing"})

        resp = client.patch(f"/{DEFAULTS}/jobs/j1/1", json={"job_id": "j1", "run_id": 1, "time_start": 5})

        assert resp.json()["phase"] == "processing"

    def test_patch_unknown_job_returns_404(self, client):
        resp = client.patch(
            f"/{DEFAULTS}/jobs/nope/1",
            json={"job_id": "nope", "run_id": 1, "phase": "completed"},
        )

        assert resp.status_code == 404

    def test_patch_without_run_id_is_rejected(self, client):
        """run_id has no default on JobUpdate, so it is required."""
        _create(client)

        resp = client.patch(f"/{DEFAULTS}/jobs/j1/1", json={"job_id": "j1", "phase": "completed"})

        assert resp.status_code == 422


class TestPutJob:
    def test_put_updates_phase_and_timestamps(self, client):
        _create(client)

        resp = client.put(
            f"/{DEFAULTS}/jobs/j1/1",
            json={"job_id": "j1", "run_id": "1", "phase": "completed", "time_start": 1, "time_end": 2},
        )

        assert resp.status_code == 200
        body = resp.json()
        assert body["phase"] == "completed"
        assert body["time_start"] == 1
        assert body["time_end"] == 2

    def test_put_unknown_job_returns_404(self, client):
        resp = client.put(
            f"/{DEFAULTS}/jobs/nope/1",
            json={"job_id": "nope", "run_id": "1", "phase": "completed"},
        )

        assert resp.status_code == 404
