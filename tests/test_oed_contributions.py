"""Tests for POST /oed/contributions, the OED About page's data contribution form."""
import inspect
import json

import pytest

from routers import oed


FORM = {"email": "contributor@example.org", "license_accepted": "true", "rights_confirmed": "true"}


def _submit(client, filename="kinetics.xlsx", content=b"PK\x03\x04data", **form):
    return client.post("/oed/contributions",
                       data={**FORM, **form},
                       files={"file": (filename, content, "application/octet-stream")})


def _stored(client):
    return client.minio.objects.get(oed.CONTRIBUTIONS_BUCKET, {})


class TestSubmitContribution:
    def test_stores_file_and_submission_record(self, client):
        response = _submit(client)

        assert response.status_code == 201
        contribution_id = response.json()["contributionID"]
        stored = _stored(client)
        assert stored[f"{contribution_id}/kinetics.xlsx"] == b"PK\x03\x04data"
        submission = json.loads(stored[f"{contribution_id}/submission.json"])
        assert submission["email"] == "contributor@example.org"
        assert submission["filename"] == "kinetics.xlsx"
        assert submission["license"] == "CC-BY-4.0"
        assert submission["rights_confirmed"] is True
        assert submission["submitted_at"] == response.json()["submitted_at"]

    def test_email_is_not_in_the_object_key(self, client):
        _submit(client)

        assert not any("contributor@example.org" in key for key in _stored(client))

    def test_each_submission_gets_its_own_id(self, client):
        first = _submit(client).json()["contributionID"]
        second = _submit(client).json()["contributionID"]

        assert first != second
        assert len(_stored(client)) == 4

    def test_client_path_is_reduced_to_base_name(self, client):
        response = _submit(client, filename="C:\\Users\\me\\..\\data.csv")

        assert response.status_code == 201
        assert set(_stored(client)) == {f"{response.json()['contributionID']}/{name}"
                                        for name in ("data.csv", "submission.json")}

    @pytest.mark.parametrize("field", ["license_accepted", "rights_confirmed"])
    def test_requires_both_agreements(self, client, field):
        response = _submit(client, **{field: "false"})

        assert response.status_code == 400
        assert _stored(client) == {}

    @pytest.mark.parametrize("filename", ["notes.pdf", "kinetics.xlsx.exe", ""])
    def test_rejects_other_file_types(self, client, filename):
        response = _submit(client, filename=filename)

        assert response.status_code in (400, 422)
        assert _stored(client) == {}

    def test_rejects_invalid_email(self, client):
        response = _submit(client, email="not-an-email")

        assert response.status_code == 422
        assert _stored(client) == {}

    def test_rejects_empty_file(self, client):
        response = _submit(client, content=b"")

        assert response.status_code == 400
        assert _stored(client) == {}

    def test_rejects_file_over_size_limit(self, client, monkeypatch):
        monkeypatch.setattr(oed, "CONTRIBUTION_MAX_BYTES", 10)

        assert _submit(client, content=b"x" * 10).status_code == 201
        response = _submit(client, content=b"x" * 11)

        assert response.status_code == 413
        assert len(_stored(client)) == 2

    def test_storage_failure_is_a_500(self, client, monkeypatch):
        monkeypatch.setattr(client.minio, "upload_file", lambda *args: False)

        assert _submit(client).status_code == 500

    def test_failed_data_upload_leaves_the_submission_record(self, client, monkeypatch):
        upload_file = client.minio.upload_file
        monkeypatch.setattr(client.minio, "upload_file",
                            lambda bucket, name, data: name.endswith("submission.json") and upload_file(bucket, name, data))

        assert _submit(client).status_code == 500
        [key] = _stored(client)
        assert key.endswith("/submission.json")

    def test_upload_route_still_rejects_the_contributions_bucket(self, client):
        response = client.post(f"/{oed.CONTRIBUTIONS_BUCKET}/upload",
                               files={"file": ("kinetics.xlsx", b"data", "application/octet-stream")})

        assert response.status_code == 400

    def test_handler_is_sync_so_minio_calls_run_off_the_event_loop(self):
        # The MinIO client blocks; FastAPI only moves a handler to its threadpool if it is a plain def.
        assert not inspect.iscoroutinefunction(oed.submit_contribution)
