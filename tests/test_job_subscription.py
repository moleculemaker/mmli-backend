"""POST /{job_type}/jobs/{job_id}/subscription: add a notification address after submission.

Anyone holding a job id can call this, so it only fills an empty address and never
replaces one. See #141.
"""
import pytest
import sqlmodel

import main
from models.sqlmodel.models import Job

SOMN = "somn"


def _row(client, job_id):
    engine = sqlmodel.create_engine(client.sync_db_url)
    with sqlmodel.Session(engine) as session:
        return session.get(Job, job_id)


def _set_stored_email(client, job_id, value):
    engine = sqlmodel.create_engine(client.sync_db_url)
    with sqlmodel.Session(engine) as session:
        job = session.get(Job, job_id)
        job.email = value
        session.add(job)
        session.commit()


def _create(client, job_id="j1", **body):
    resp = client.post(f"/{SOMN}/jobs", json={"job_id": job_id, **body})
    assert resp.status_code == 201
    return resp


def _subscribe(client, email, job_type=SOMN, job_id="j1"):
    return client.post(f"/{job_type}/jobs/{job_id}/subscription", json={"email": email})


class TestSubscribe:
    def test_sets_the_address_when_the_job_has_none(self, client):
        _create(client)

        resp = _subscribe(client, "  Someone@EXAMPLE.org ")

        assert resp.status_code == 204
        assert _row(client, "j1").email == "Someone@example.org"

    def test_an_empty_stored_address_counts_as_none(self, client):
        """Rows from before #137 can hold '' rather than null."""
        _create(client)
        _set_stored_email(client, "j1", "")

        assert _subscribe(client, "a@example.org").status_code == 204
        assert _row(client, "j1").email == "a@example.org"

    def test_resubmitting_the_same_address_succeeds(self, client):
        _create(client)
        _subscribe(client, "a@example.org")

        assert _subscribe(client, "a@EXAMPLE.org").status_code == 204
        assert _row(client, "j1").email == "a@example.org"

    def test_a_different_address_never_replaces_the_stored_one(self, client):
        _create(client)
        _subscribe(client, "owner@example.org")

        resp = _subscribe(client, "someone-else@example.org")

        assert resp.status_code == 409
        assert _row(client, "j1").email == "owner@example.org"
        # Says that a subscriber exists, never who.
        assert "owner@example.org" not in resp.text

    def test_an_address_given_at_submission_is_not_replaced(self, client):
        _create(client, email="owner@example.org")

        assert _subscribe(client, "someone-else@example.org").status_code == 409
        assert _row(client, "j1").email == "owner@example.org"


class TestRejected:
    @pytest.mark.parametrize("address", [
        pytest.param("victim@example.org, other@example.com", id="list"),
        pytest.param('"other@example.com" <victim@example.org>', id="display-name"),
        pytest.param("victim@example.org\x00other", id="nul"),
        pytest.param("not-an-address", id="no-at-sign"),
        pytest.param("", id="empty"),
        pytest.param("   ", id="blank"),
    ])
    def test_an_unusable_address_is_422(self, client, address):
        _create(client)

        assert _subscribe(client, address).status_code == 422
        assert _row(client, "j1").email is None

    def test_a_missing_email_field_is_422(self, client):
        _create(client)

        assert client.post(f"/{SOMN}/jobs/j1/subscription", json={}).status_code == 422

    def test_an_unknown_job_is_404(self, client):
        assert _subscribe(client, "a@example.org", job_id="nope").status_code == 404

    def test_a_job_of_another_type_is_404(self, client):
        """The path's job_type is checked against the row, not ignored."""
        _create(client)

        assert _subscribe(client, "a@example.org", job_type="defaults").status_code == 404
        assert _row(client, "j1").email is None

    def test_an_invalid_job_type_is_400(self, client):
        assert _subscribe(client, "a@example.org", job_type="not-a-real-tool").status_code == 400


class TestGeneratedClients:
    def test_the_path_parameters_are_declared(self):
        """The old PATCH route did not declare its path parameters, so the generated
        Angular client sent `/{job_type}/jobs//` and never reached the handler."""
        operation = main.app.openapi()["paths"]["/{job_type}/jobs/{job_id}/subscription"]["post"]

        declared = {p["name"] for p in operation["parameters"] if p["in"] == "path"}
        assert declared == {"job_type", "job_id"}
