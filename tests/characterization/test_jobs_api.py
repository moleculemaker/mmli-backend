"""Characterization tests for the legacy Jobs API.

These pin the CURRENT behavior of `/{job_type}/jobs`, including behavior that is
wrong. Tests covering known defects are marked with a `DEFECT:` comment naming what is
wrong; when a later change fixes one, the test diff is the record of that decision.

Nothing here should be read as an endorsement of the behavior it asserts.
"""
import json

import pytest

from models.enums import JobStatus
from services import kubejob_service


DEFAULTS = "defaults"


def _post_job(client, job_type=DEFAULTS, **body):
    return client.post(f"/{job_type}/jobs", json=body)


class TestCreateJob:
    def test_new_job_returns_201_with_stringified_fields(self, client):
        resp = _post_job(client, email="a@b.edu", job_info='{"x": 1}')

        assert resp.status_code == 201
        body = resp.json()
        assert set(body) == {"job_id", "run_id", "email", "job_info"}
        assert body["email"] == "a@b.edu"
        assert body["job_info"] == '{"x": 1}'
        assert len(body["job_id"]) == 32  # uuid4().hex, dashes stripped

    def test_absent_optional_fields_serialize_as_the_string_None(self, client):
        # DEFECT: every field is passed through str(), so a null run_id/email is
        # returned as the four-character string "None" rather than JSON null.
        resp = _post_job(client)

        assert resp.status_code == 201
        assert resp.json()["run_id"] == "None"
        assert resp.json()["email"] == "None"

    def test_job_info_defaults_to_empty_object_string(self, client):
        assert _post_job(client).json()["job_info"] == "{}"

    def test_client_may_choose_its_own_job_id(self, client):
        resp = _post_job(client, job_id="client-chosen-id")

        assert resp.status_code == 201
        assert resp.json()["job_id"] == "client-chosen-id"

    def test_a_kubernetes_job_is_submitted(self, client):
        _post_job(client, job_id="k8s-1")

        assert len(client.k8s_jobs) == 1
        assert client.k8s_jobs[0]["job_id"] == "k8s-1"
        assert client.k8s_jobs[0]["job_type"] == DEFAULTS

    def test_invalid_job_type_returns_400(self, client):
        resp = _post_job(client, job_type="not-a-real-tool")

        assert resp.status_code == 400
        assert "Invalid job type" in resp.json()["detail"]


class TestCreateJobWithExistingId:
    """The collision branch: POSTing a job_id that already exists.

    This endpoint is unauthenticated, so naming an existing id proves only that the
    caller can name it. Each test here previously asserted a defect; see the diff of
    this file in the "legacy fixes" change for what the behavior used to be.
    """

    def test_does_not_disclose_the_stored_email(self, client):
        _post_job(client, job_id="victim", email="owner@illinois.edu", job_info='{"secret": 1}')

        resp = _post_job(client, job_id="victim")

        # Was: returned the owner's address, making the endpoint an oracle for
        # harvesting researchers' emails from guessed job_ids.
        assert resp.status_code == 200
        assert resp.json()["email"] is None

    def test_does_not_disclose_the_stored_job_info(self, client):
        _post_job(client, job_id="victim", email="owner@illinois.edu", job_info='{"secret": 1}')

        resp = _post_job(client, job_id="victim")

        # Was: returned the owner's submitted inputs.
        assert resp.json()["job_info"] is None

    def test_still_returns_the_job_id(self, client):
        """coordinator.py reads job_id from a 200 or a 201, so it must survive."""
        _post_job(client, job_id="victim")

        resp = _post_job(client, job_id="victim")

        assert resp.status_code == 200
        assert resp.json()["job_id"] == "victim"
        assert set(resp.json()) == {"job_id", "run_id", "email", "job_info"}

    def test_does_not_overwrite_the_stored_job_info(self, client):
        _post_job(client, job_id="victim", email="owner@illinois.edu", job_info='{"secret": 1}')

        _post_job(client, job_id="victim", job_info='{"tampered": true}')

        # Was: the caller's job_info replaced the owner's, destroying the record of
        # what was actually run.
        stored = client.get(f"/{DEFAULTS}/jobs/victim").json()[0]["job_info"]
        assert stored == '{"secret": 1}'

    def test_does_not_submit_a_second_kubernetes_job(self, client):
        _post_job(client, job_id="victim")
        _post_job(client, job_id="victim")

        # Was: create_job ran before the existence check, so a duplicate POST forked a
        # second run of the same job.
        assert len(client.k8s_jobs) == 1

    def test_does_not_disturb_the_existing_jobs_phase(self, client):
        _post_job(client, job_id="victim")

        _post_job(client, job_id="victim")

        assert client.get(f"/{DEFAULTS}/jobs/victim").json()[0]["phase"] == JobStatus.QUEUED


class TestKubernetesRejection:
    """create_job swallows ApiException and reports failure via its return value."""

    def test_a_rejected_submission_returns_400(self, client, monkeypatch):
        monkeypatch.setattr(
            kubejob_service,
            "create_job",
            lambda **kw: {"status": "error", "message": "jobs.batch already exists"},
        )

        resp = _post_job(client, job_id="rejected")

        # Was: 201, because the exception never propagated and the return value was
        # ignored -- the caller was told the job started when no pod existed.
        assert resp.status_code == 400
        assert "already exists" in resp.json()["detail"]

    def test_a_rejected_submission_leaves_the_job_in_error_not_queued(self, client, monkeypatch):
        monkeypatch.setattr(
            kubejob_service,
            "create_job",
            lambda **kw: {"status": "error", "message": "exceeded quota"},
        )

        _post_job(client, job_id="rejected")

        # Was: stranded at 'queued' forever, with nothing running to advance it.
        assert client.get(f"/{DEFAULTS}/jobs/rejected").json()[0]["phase"] == JobStatus.ERROR

    def test_a_successful_submission_is_unaffected(self, client):
        resp = _post_job(client, job_id="fine")

        assert resp.status_code == 201
        assert client.get(f"/{DEFAULTS}/jobs/fine").json()[0]["phase"] == JobStatus.QUEUED


class TestReadJobs:
    def test_list_by_type_returns_every_users_job_including_email(self, client):
        _post_job(client, job_id="j1", email="one@illinois.edu")
        _post_job(client, job_id="j2", email="two@illinois.edu")

        resp = client.get(f"/{DEFAULTS}/jobs")

        # DEFECT: unauthenticated, unpaginated, and returns other users' addresses.
        assert resp.status_code == 200
        emails = sorted(job["email"] for job in resp.json())
        assert emails == ["one@illinois.edu", "two@illinois.edu"]

    def test_list_by_job_id_returns_a_list_not_an_object(self, client):
        _post_job(client, job_id="j1")

        body = client.get(f"/{DEFAULTS}/jobs/j1").json()

        # A single job is still wrapped in a list by this endpoint.
        assert isinstance(body, list)
        assert len(body) == 1
        assert body[0]["job_id"] == "j1"

    def test_new_jobs_start_queued(self, client):
        _post_job(client, job_id="j1")

        assert client.get(f"/{DEFAULTS}/jobs/j1").json()[0]["phase"] == JobStatus.QUEUED

    def test_unknown_job_id_returns_an_empty_list_not_404(self, client):
        resp = client.get(f"/{DEFAULTS}/jobs/does-not-exist")

        assert resp.status_code == 200
        assert resp.json() == []

    def test_invalid_job_type_returns_400(self, client):
        assert client.get("/not-a-real-tool/jobs").status_code == 400


class TestDeleteJob:
    def test_delete_removes_the_row_and_the_kubernetes_job(self, client, deleted_k8s_jobs):
        _post_job(client, job_id="j1", run_id="r1")

        resp = client.delete(f"/{DEFAULTS}/jobs/j1/r1")

        assert resp.status_code == 200
        assert client.get(f"/{DEFAULTS}/jobs/j1").json() == []
        # Was: the row was dropped but the pod kept running, consuming cluster
        # resources with nothing left to observe or stop it.
        assert deleted_k8s_jobs == [{"job_type": DEFAULTS, "job_id": "j1"}]

    def test_delete_succeeds_even_if_the_cluster_call_fails(self, client, monkeypatch):
        """A cluster problem must not leave a user unable to delete their own record."""
        def _boom(**kwargs):
            raise RuntimeError("cluster unreachable")

        monkeypatch.setattr(kubejob_service, "delete_job", _boom)
        _post_job(client, job_id="j1", run_id="r1")

        resp = client.delete(f"/{DEFAULTS}/jobs/j1/r1")

        assert resp.status_code == 200
        assert client.get(f"/{DEFAULTS}/jobs/j1").json() == []

    def test_delete_unknown_job_returns_404(self, client, deleted_k8s_jobs):
        assert client.delete(f"/{DEFAULTS}/jobs/nope/nope").status_code == 404
        assert deleted_k8s_jobs == []

    def test_delete_with_an_invalid_job_type_returns_400(self, client):
        """Every other route on this router validates the type before querying.

        Whether filtering on an unrecognized type is harmless depends on the schema.
        The deployed schema (built by alembic) stores `type` as a varchar and simply
        matches nothing. A schema built from the SQLModel metadata - which is what
        SQLModel.metadata.create_all() produces, and what this suite uses - maps
        JobType to a native Postgres enum, where the driver raises
        InvalidTextRepresentationError and the request 500s.

        Validating first makes the route behave the same either way. See the README
        note on schema divergence.
        """
        resp = client.delete("/not-a-real-tool/jobs/j1/r1")

        assert resp.status_code == 400
        assert "Invalid job type" in resp.json()["detail"]


class TestPerToolInputValidation:
    def test_chemscraper_requires_input_file(self, client):
        resp = _post_job(client, job_type="chemscraper", job_info="{}")

        assert resp.status_code == 400
        assert "input_file" in resp.json()["detail"]

    def test_ml_simplefold_requires_fasta(self, client):
        resp = _post_job(client, job_type="ml-simplefold", job_info="{}")

        assert resp.status_code == 400
        assert "fasta" in resp.json()["detail"]

    def test_ez_specificity_requires_enzymes_and_substrates(self, client):
        resp = _post_job(client, job_type="ez-specificity", job_info="{}")

        assert resp.status_code == 400
        assert "enzymes" in resp.json()["detail"]

    @pytest.mark.parametrize(
        "job_info, expected",
        [
            ('{"enzymes": {}, "substrates": []}', "must be lists"),
            ('{"enzymes": [], "substrates": ["c"]}', "enzyme"),
            ('{"enzymes": [{"filename": "a.pdb"}], "substrates": []}', "substrate"),
        ],
    )
    def test_ez_specificity_rejects_malformed_input(self, client, job_info, expected):
        resp = _post_job(client, job_type="ez-specificity", job_info=job_info)

        assert resp.status_code == 400
        assert expected in resp.json()["detail"]

    def test_ez_specificity_rejects_an_enzyme_that_was_never_uploaded(self, client):
        resp = _post_job(
            client,
            job_type="ez-specificity",
            job_info='{"enzymes": [{"filename": "missing.pdb"}], "substrates": ["CCO"]}',
        )

        assert resp.status_code == 400
        assert "not found in uploads" in resp.json()["detail"]


class TestJobInfoEscaping:
    def test_job_info_is_stored_as_an_opaque_string(self, client):
        """job_info is a string field, not a modelled object.

        Nothing about its per-tool shape appears in the OpenAPI schema; this is the
        core FAIR problem the /v1 API exists to solve.
        """
        _post_job(client, job_id="j1", job_info='{"nested": {"a": [1, 2]}}')

        stored = client.get(f"/{DEFAULTS}/jobs/j1").json()[0]["job_info"]

        assert isinstance(stored, str)
        assert json.loads(stored) == {"nested": {"a": [1, 2]}}
