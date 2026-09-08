"""Tests for the versioned API.

These are not characterization tests: /v1 is new, so they describe intended behavior
rather than pinning existing behavior.
"""
import json

import pytest

from models.enums import JobStatus
from services import kubejob_service


SOMN_INPUT = [{
    "reactant_pair_name": "pair1",
    "el": "CC(C)c1ccc(Br)cc1", "el_name": "el1", "el_input_type": "smi", "el_idx": "0",
    "nuc": "NCc1ccccc1", "nuc_name": "nuc1", "nuc_input_type": "smi", "nuc_idx": "0",
}]
OPTSTOIC_INPUT = {"primary_precursor": "MNXM1137670", "target_molecule": "MNXM26"}


def _submit(client, tool="novostoic-optstoic", payload=None, **kwargs):
    return client.post(f"/v1/tools/{tool}/jobs",
                       json=OPTSTOIC_INPUT if payload is None else payload, **kwargs)


class TestServiceInfo:
    def test_describes_the_service(self, client):
        body = client.get("/v1/service-info").json()

        assert body["id"] == "org.moleculemaker.alphasynthesis"
        assert body["toolCount"] > 0

    def test_states_the_access_model_explicitly(self, client):
        """FAIR requires access conditions be stated, not that they be open."""
        auth = client.get("/v1/service-info").json()["auth"]

        assert auth["type"] == "none"
        assert auth["model"] == "capability-url"
        assert "secret" in auth["description"]

    def test_links_are_absolute_and_usable(self, client):
        links = client.get("/v1/service-info").json()["links"]

        assert client.get(links["tools"]).status_code == 200


class TestDiscovery:
    def test_lists_tools_as_json_ld(self, client):
        body = client.get("/v1/tools").json()

        assert body["@context"] == "https://schema.org"
        assert body["@type"] == "ItemList"
        assert body["numberOfItems"] == len(body["itemListElement"])

    def test_internal_subjobs_are_not_listed(self, client):
        ids = {t["identifier"] for t in client.get("/v1/tools").json()["itemListElement"]}

        assert "ez-specificity" in ids
        assert "ezspec-unidock" not in ids

    def test_a_tool_descriptor_carries_provenance_fields(self, client):
        body = client.get("/v1/tools/somn").json()

        assert body["@type"] == "SoftwareApplication"
        assert body["identifier"] == "somn"
        # Published as null rather than omitted: nobody has supplied them yet.
        assert "license" in body and "citation" in body
        assert set(body["edam"]) == {"operation", "topic"}

    def test_unknown_tool_returns_a_problem_listing_what_exists(self, client):
        resp = client.get("/v1/tools/not-a-real-tool")

        assert resp.status_code == 404
        assert resp.headers["content-type"].startswith("application/problem+json")
        body = resp.json()
        assert body["type"].endswith("/unknown-tool")
        assert "somn" in body["known_tools"]

    def test_input_schema_is_served_with_an_id(self, client):
        body = client.get("/v1/tools/somn/input-schema").json()

        assert body["$schema"] == "https://json-schema.org/draft/2020-12/schema"
        assert body["$id"].endswith("/v1/tools/somn/input-schema")

    def test_a_tool_without_a_schema_says_so(self, client):
        """ezspec subjobs take no job_info of their own."""
        resp = client.get("/v1/tools/ezspec-unidock/input-schema")

        assert resp.status_code == 404


class TestSubmission:
    def test_a_valid_submission_returns_201_with_links(self, client):
        resp = _submit(client)

        assert resp.status_code == 201
        body = resp.json()
        assert body["status"] == JobStatus.QUEUED
        assert body["tool"] == "novostoic-optstoic"
        assert body["links"]["self"].endswith(f"/v1/jobs/{body['job_id']}")

    def test_the_server_assigns_the_job_id(self, client):
        """A client-supplied id would be authority over a job in an unauthenticated API."""
        resp = client.post("/v1/tools/novostoic-optstoic/jobs",
                           json={**OPTSTOIC_INPUT, "job_id": "chosen-by-client"})

        # job_id is not part of the schema, so it is rejected outright rather than
        # silently ignored.
        assert resp.status_code == 422

    def test_a_kubernetes_job_is_submitted(self, client):
        resp = _submit(client)

        assert len(client.k8s_jobs) == 1
        assert client.k8s_jobs[0]["job_id"] == resp.json()["job_id"]

    def test_timestamps_are_iso_8601(self, client):
        body = _submit(client).json()

        assert body["submitted_at"].endswith("Z")
        # Unset timestamps are null, not 1970.
        assert body["started_at"] is None
        assert body["finished_at"] is None

    def test_unknown_tool_is_rejected(self, client):
        resp = client.post("/v1/tools/not-a-real-tool/jobs", json={})

        assert resp.status_code == 404

    def test_internal_subjob_types_cannot_be_submitted(self, client):
        resp = client.post("/v1/tools/ezspec-unidock/jobs", json={})

        assert resp.status_code == 404


class TestInputValidation:
    def test_invalid_input_reports_a_json_pointer(self, client):
        resp = client.post("/v1/tools/novostoic-optstoic/jobs",
                           json={"primary_precursor": "MNXM1"})

        assert resp.status_code == 422
        body = resp.json()
        assert body["type"].endswith("/invalid-input")
        assert body["schema_url"] == "/v1/tools/novostoic-optstoic/input-schema"
        assert body["errors"]

    def test_every_error_is_reported_not_just_the_first(self, client):
        resp = client.post("/v1/tools/somn/jobs", json=[{"reactant_pair_name": "x"}])

        assert resp.status_code == 422
        assert len(resp.json()["errors"]) > 1

    def test_pointer_locates_the_offending_element(self, client):
        bad = [dict(SOMN_INPUT[0]), dict(SOMN_INPUT[0])]
        bad[1]["el_input_type"] = "not-a-format"

        resp = client.post("/v1/tools/somn/jobs", json=bad)

        assert resp.status_code == 422
        assert any(e["pointer"] == "/1/el_input_type" for e in resp.json()["errors"])

    def test_no_kubernetes_job_is_started_for_invalid_input(self, client):
        client.post("/v1/tools/novostoic-optstoic/jobs", json={})

        assert client.k8s_jobs == []

    def test_malformed_json_body_is_a_problem_not_a_crash(self, client):
        resp = client.post("/v1/tools/novostoic-optstoic/jobs",
                           content=b"{not json", headers={"content-type": "application/json"})

        assert resp.status_code == 400
        assert resp.headers["content-type"].startswith("application/problem+json")


class TestIdempotency:
    def test_replaying_a_key_returns_the_original_job(self, client):
        first = _submit(client, headers={"Idempotency-Key": "abc"})
        second = _submit(client, headers={"Idempotency-Key": "abc"})

        assert first.json()["job_id"] == second.json()["job_id"]

    def test_replaying_a_key_does_not_start_a_second_job(self, client):
        _submit(client, headers={"Idempotency-Key": "abc"})
        _submit(client, headers={"Idempotency-Key": "abc"})

        assert len(client.k8s_jobs) == 1

    def test_different_keys_produce_different_jobs(self, client):
        first = _submit(client, headers={"Idempotency-Key": "abc"})
        second = _submit(client, headers={"Idempotency-Key": "xyz"})

        assert first.json()["job_id"] != second.json()["job_id"]
        assert len(client.k8s_jobs) == 2

    def test_the_same_key_on_a_different_tool_is_independent(self, client):
        a = _submit(client, headers={"Idempotency-Key": "shared"})
        b = _submit(client, tool="somn", payload=SOMN_INPUT,
                    headers={"Idempotency-Key": "shared"})

        assert a.json()["job_id"] != b.json()["job_id"]

    def test_no_key_means_every_submission_is_a_new_job(self, client):
        first = _submit(client)
        second = _submit(client)

        assert first.json()["job_id"] != second.json()["job_id"]


class TestMultipartSubmission:
    def test_files_are_stored_under_the_assigned_job_id(self, client):
        resp = client.post(
            "/v1/tools/molli/jobs",
            data={"inputs": json.dumps({"CORES_FILE_NAME": "cores.cdxml",
                                        "SUBS_FILE_NAME": "subs.cdxml"})},
            files=[("files", ("cores.cdxml", b"<cores/>", "application/xml")),
                   ("files", ("subs.cdxml", b"<subs/>", "application/xml"))],
        )

        assert resp.status_code == 201
        job_id = resp.json()["job_id"]
        # The whole point of multipart: the client never had to know the id in advance.
        assert client.minio.get_file("molli", f"{job_id}/in/cores.cdxml") == b"<cores/>"

    def test_inputs_part_is_still_schema_validated(self, client):
        resp = client.post(
            "/v1/tools/molli/jobs",
            data={"inputs": json.dumps({"CORES_FILE_NAME": "cores.cdxml"})},
            files=[("files", ("cores.cdxml", b"<cores/>", "application/xml"))],
        )

        assert resp.status_code == 422


class TestJobRetrieval:
    def test_a_job_is_retrievable_by_id(self, client):
        job_id = _submit(client).json()["job_id"]

        body = client.get(f"/v1/jobs/{job_id}").json()

        assert body["job_id"] == job_id

    def test_unknown_job_returns_a_problem(self, client):
        resp = client.get("/v1/jobs/never-existed")

        assert resp.status_code == 404
        assert resp.json()["type"].endswith("/unknown-job")

    def test_provenance_is_reported_including_what_is_unknown(self, client):
        job_id = _submit(client).json()["job_id"]

        provenance = client.get(f"/v1/jobs/{job_id}").json()["provenance"]

        assert provenance["image"]
        # Null until the watcher reads it from a running pod; never guessed.
        assert provenance["image_digest"] is None

    def test_there_is_no_endpoint_that_lists_jobs(self, client):
        """Deliberate: listing would hand out other people's job ids and emails."""
        _submit(client)

        assert client.get("/v1/jobs").status_code == 404
        assert client.get("/v1/tools/novostoic-optstoic/jobs").status_code == 405


class TestResults:
    def _finish(self, client, job_id, phase=JobStatus.COMPLETED):
        import sqlmodel
        from models.sqlmodel.models import Job as JobModel
        engine = sqlmodel.create_engine(client.sync_db_url)
        with sqlmodel.Session(engine) as session:
            job = session.get(JobModel, job_id)
            job.phase = phase
            session.add(job)
            session.commit()

    def test_running_job_returns_409_with_retry_after(self, client):
        job_id = _submit(client).json()["job_id"]

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 409
        assert resp.headers["Retry-After"] == "10"
        assert resp.json()["type"].endswith("/job-not-finished")

    def test_unknown_job_returns_404_not_409(self, client):
        """The legacy endpoint cannot tell these apart; this one must."""
        resp = client.get("/v1/jobs/never-existed/results")

        assert resp.status_code == 404

    def test_failed_job_says_so_rather_than_returning_nothing(self, client):
        job_id = _submit(client).json()["job_id"]
        self._finish(client, job_id, JobStatus.ERROR)

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 409
        assert resp.json()["type"].endswith("/job-failed")

    def test_completed_job_with_no_output_is_distinguishable(self, client):
        job_id = _submit(client).json()["job_id"]
        self._finish(client, job_id)

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 409
        assert "no result file" in resp.json()["detail"]

    def test_completed_job_returns_its_results(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put("novostoic-optstoic", f"{job_id}/out/output.json", b"[]")
        self._finish(client, job_id)

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 200
        assert resp.json()["job_id"] == job_id


class TestArtifacts:
    def test_lists_output_files_with_download_links(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put("novostoic-optstoic", f"{job_id}/out/output.json", b"{}")

        body = client.get(f"/v1/jobs/{job_id}/artifacts").json()

        assert body["count"] == 1
        assert body["items"][0]["name"] == "output.json"

    def test_logs_are_hidden_unless_requested(self, client):
        job_id = _submit(client).json()["job_id"]
        for name in ("output.json", "error.log", "output.log", "success"):
            client.minio.put("novostoic-optstoic", f"{job_id}/out/{name}", b"x")

        default = client.get(f"/v1/jobs/{job_id}/artifacts").json()
        with_logs = client.get(f"/v1/jobs/{job_id}/artifacts?include=logs").json()

        # error.log carries whatever the container printed, and anyone holding the job
        # id can read it.
        assert [i["name"] for i in default["items"]] == ["output.json"]
        assert with_logs["count"] == 4

    def test_an_artifact_can_be_downloaded(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put("novostoic-optstoic", f"{job_id}/out/output.json", b'{"a":1}')

        resp = client.get(f"/v1/jobs/{job_id}/artifacts/output.json")

        assert resp.status_code == 200
        assert resp.content == b'{"a":1}'

    def test_path_traversal_is_rejected(self, client):
        job_id = _submit(client).json()["job_id"]

        resp = client.get(f"/v1/jobs/{job_id}/artifacts/..%2F..%2Fother/out/secret.json")

        assert resp.status_code in (400, 404)

    def test_unknown_artifact_returns_404(self, client):
        job_id = _submit(client).json()["job_id"]

        assert client.get(f"/v1/jobs/{job_id}/artifacts/nope.json").status_code == 404


class TestCancel:
    def test_canceling_stops_the_job_and_records_it(self, client, deleted_k8s_jobs):
        job_id = _submit(client).json()["job_id"]

        resp = client.post(f"/v1/jobs/{job_id}/cancel")

        assert resp.status_code == 200
        assert resp.json()["status"] == JobStatus.CANCELED
        assert deleted_k8s_jobs == [{"job_type": "novostoic-optstoic", "job_id": job_id}]

    def test_cancel_is_idempotent(self, client):
        job_id = _submit(client).json()["job_id"]
        client.post(f"/v1/jobs/{job_id}/cancel")

        resp = client.post(f"/v1/jobs/{job_id}/cancel")

        assert resp.status_code == 200
        assert resp.json()["status"] == JobStatus.CANCELED

    def test_a_finished_job_cannot_be_canceled(self, client):
        job_id = _submit(client).json()["job_id"]
        TestResults()._finish(client, job_id)

        resp = client.post(f"/v1/jobs/{job_id}/cancel")

        assert resp.status_code == 409
        assert resp.json()["type"].endswith("/not-cancelable")

    def test_cancel_succeeds_even_if_the_cluster_call_fails(self, client, monkeypatch):
        def _boom(**kwargs):
            raise RuntimeError("cluster unreachable")

        monkeypatch.setattr(kubejob_service, "delete_job", _boom)
        job_id = _submit(client).json()["job_id"]

        resp = client.post(f"/v1/jobs/{job_id}/cancel")

        assert resp.status_code == 200
        assert resp.json()["status"] == JobStatus.CANCELED


class TestErrorFormat:
    @pytest.mark.parametrize("path", [
        "/v1/tools/nope", "/v1/jobs/nope", "/v1/jobs/nope/results",
    ])
    def test_errors_are_problem_json(self, client, path):
        resp = client.get(path)

        assert resp.headers["content-type"].startswith("application/problem+json")
        body = resp.json()
        assert {"type", "title", "status"} <= set(body)
        assert body["status"] == resp.status_code


class TestRouteOrdering:
    """Regression guard on a shadowing bug that is easy to reintroduce.

    The legacy router declares /{job_type}/jobs/{job_id}, which matches /v1/jobs/<id>
    with job_type="v1". Starlette matches in registration order, so if the /v1 mount is
    added after the legacy routers, most of the versioned API disappears behind
    "Invalid job type: v1" -- with a 400 and the legacy error shape, not a problem
    document. Every /v1 test would fail, but only after someone reordered main.py.
    """

    def test_v1_is_mounted_before_the_legacy_catch_all_routes(self, client):
        import main

        paths = [getattr(r, "path", None) for r in main.app.routes]
        mount_index = paths.index("/v1")
        legacy_index = min(
            i for i, p in enumerate(paths) if p == "/{job_type}/jobs/{job_id}"
        )

        assert mount_index < legacy_index

    def test_a_job_path_under_v1_is_not_answered_by_the_legacy_router(self, client):
        resp = client.get("/v1/jobs/does-not-exist")

        assert resp.status_code == 404
        assert "Invalid job type" not in resp.text

    def test_the_legacy_route_still_works(self, client):
        """The reorder must not cost the legacy surface anything."""
        client.post("/defaults/jobs", json={"job_id": "legacy-1"})

        assert client.get("/defaults/jobs/legacy-1").json()[0]["job_id"] == "legacy-1"


class TestOpenApi:
    def test_v1_publishes_openapi_31(self, client):
        """3.1 is what makes the published documents true JSON Schema."""
        assert client.get("/v1/openapi.json").json()["openapi"].startswith("3.1")

    def test_the_mount_path_is_advertised_as_a_server(self, client):
        """Swagger UI resolves request URLs against this.

        Without it, every "Try it out" in the docs is sent to the unmounted path and
        answers 404 -- which is invisible in a static reading of the page and only shows
        up when a request is actually executed.
        """
        assert client.get("/v1/openapi.json").json()["servers"] == [{"url": "/v1"}]

    def test_submission_documents_both_json_and_multipart(self, client):
        content = (client.get("/v1/openapi.json").json()
                   ["paths"]["/tools/{tool}/jobs"]["post"]["requestBody"]["content"])

        assert set(content) == {"application/json", "multipart/form-data"}

    def test_json_is_listed_first_so_the_docs_default_to_it(self, client):
        """Swagger UI selects whichever content type comes first.

        Only the multipart form is declared through function parameters, so that is all
        FastAPI infers; left alone the docs present the file-upload path as the default
        and the simple case looks like the complicated one.
        """
        content = (client.get("/v1/openapi.json").json()
                   ["paths"]["/tools/{tool}/jobs"]["post"]["requestBody"]["content"])

        assert list(content)[0] == "application/json"

    def test_the_json_body_points_at_the_per_tool_schema(self, client):
        """It is polymorphic per tool, so this is as specific as it can honestly be."""
        content = (client.get("/v1/openapi.json").json()
                   ["paths"]["/tools/{tool}/jobs"]["post"]["requestBody"]["content"])

        assert "input-schema" in content["application/json"]["schema"]["description"]

    def test_legacy_openapi_is_still_served_separately(self, client):
        assert client.get("/openapi.json").status_code == 200
