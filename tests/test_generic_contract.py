"""Tests for the generic job contract (services/generic_job.py).

No real tool uses the contract yet, so a stand-in tool is registered for each test. Its
schemas mirror moleculemaker/job-template's sample tool.
"""
import json

import pytest
import yaml

from config import app_config
from models.enums import JobStatus
from services import tool_registry

TOOL = "sample"

INPUT_SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "type": "object",
    "required": ["sequences"],
    "additionalProperties": False,
    "properties": {
        "sequences": {"type": "array", "minItems": 1, "items": {"type": "string"},
                      "description": "Protein sequences."},
    },
}

OUTPUT_SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "type": "object",
    "required": ["rows", "artifacts"],
    "additionalProperties": False,
    "properties": {
        "rows": {"type": "array", "items": {
            "type": "object", "required": ["length"],
            "properties": {"length": {"type": "integer"}},
        }},
        "artifacts": {"type": "array"},
    },
}

INPUT = {"sequences": ["ACD", "EFG"]}
RESULTS = {"rows": [{"length": 3}, {"length": 3}],
           "artifacts": [{"name": "results.csv", "media_type": "text/csv"}]}


@pytest.fixture(autouse=True)
def generic_tool(monkeypatch):
    registry = dict(tool_registry.get_registry())
    registry[TOOL] = {
        "identifier": TOOL, "name": "Sample", "summary": "s", "description": "d",
        "execution": "kubernetes-job", "schema_status": "verified", "internal": False,
        "input_schema": INPUT_SCHEMA, "output_schema": OUTPUT_SCHEMA,
    }
    monkeypatch.setattr(tool_registry, "_REGISTRY", registry)
    monkeypatch.setitem(app_config["kubernetes_jobs"], TOOL,
                        {"image": "ghcr.io/moleculemaker/job-template:sha-3a5c0da"})


def _submit(client, payload=INPUT):
    return client.post(f"/v1/tools/{TOOL}/jobs", json=payload)


def _finish(client, job_id):
    import sqlmodel
    from models.sqlmodel.models import Job as JobModel
    engine = sqlmodel.create_engine(client.sync_db_url)
    with sqlmodel.Session(engine) as session:
        job = session.get(JobModel, job_id)
        job.phase = JobStatus.COMPLETED
        session.add(job)
        session.commit()


class TestSubmission:
    def test_validated_input_is_written_as_job_json(self, client):
        job_id = _submit(client).json()["job_id"]

        assert json.loads(client.minio.get_file(TOOL, f"{job_id}/in/job.json")) == INPUT

    def test_the_image_entrypoint_runs_with_no_command_or_tool_env(self, client):
        _submit(client)

        launched = client.k8s_jobs[0]
        assert launched["command"] == ""
        assert launched["environment"] == []
        assert launched["image_name"] == "ghcr.io/moleculemaker/job-template:sha-3a5c0da"

    def test_invalid_input_is_rejected_before_anything_is_written(self, client):
        resp = _submit(client, {"sequences": []})

        assert resp.status_code == 422
        assert client.minio.objects.get(TOOL, {}) == {}
        assert client.k8s_jobs == []

    def test_uploaded_files_sit_next_to_job_json(self, client):
        resp = client.post(
            f"/v1/tools/{TOOL}/jobs",
            data={"inputs": json.dumps(INPUT)},
            files=[("files", ("structure.pdb", b"ATOM", "chemical/x-pdb"))],
        )

        job_id = resp.json()["job_id"]
        assert client.minio.get_file(TOOL, f"{job_id}/in/structure.pdb") == b"ATOM"
        assert client.minio.get_file(TOOL, f"{job_id}/in/job.json") is not None

    def test_an_upload_named_job_json_is_refused(self, client):
        """It would be silently overwritten by the input document."""
        resp = client.post(
            f"/v1/tools/{TOOL}/jobs",
            data={"inputs": json.dumps(INPUT)},
            files=[("files", ("job.json", b"{}", "application/json"))],
        )

        assert resp.status_code == 422
        assert resp.json()["type"].endswith("/invalid-input")
        assert client.k8s_jobs == []

    def test_the_legacy_route_cannot_overwrite_an_existing_jobs_input(self, client, monkeypatch):
        """The legacy route takes a client-chosen job_id and skips schema validation."""
        # The legacy route only accepts JobType values, so make a real one generic
        registry = dict(tool_registry.get_registry())
        registry["molli"] = {**registry["molli"], "input_schema": INPUT_SCHEMA,
                             "output_schema": OUTPUT_SCHEMA}
        monkeypatch.setattr(tool_registry, "_REGISTRY", registry)
        job_id = client.post("/v1/tools/molli/jobs", json=INPUT).json()["job_id"]

        resp = client.post("/molli/jobs", json={"job_id": job_id, "job_info": '{"evil": 1}'})

        assert resp.status_code == 400
        assert "/v1/tools/molli/jobs" in resp.json()["detail"]
        assert json.loads(client.minio.get_file("molli", f"{job_id}/in/job.json")) == INPUT
        assert len(client.k8s_jobs) == 1

    def test_the_legacy_upload_route_cannot_overwrite_an_existing_jobs_input(self, client, monkeypatch):
        """Same attack through the other legacy route that takes a client-chosen job_id."""
        registry = dict(tool_registry.get_registry())
        registry["molli"] = {**registry["molli"], "input_schema": INPUT_SCHEMA,
                             "output_schema": OUTPUT_SCHEMA}
        monkeypatch.setattr(tool_registry, "_REGISTRY", registry)
        job_id = client.post("/v1/tools/molli/jobs", json=INPUT).json()["job_id"]

        resp = client.post(f"/molli/upload?job_id={job_id}",
                           files={"file": ("job.json", b'{"evil": 1}', "application/json")})

        assert resp.status_code == 400
        assert json.loads(client.minio.get_file("molli", f"{job_id}/in/job.json")) == INPUT

    @pytest.mark.parametrize("name", ["x/../job.json", "./job.json", "../job.json",
                                      "sub/file.txt", "..\\job.json"])
    def test_upload_names_with_path_segments_are_refused(self, client, name):
        resp = client.post(
            f"/v1/tools/{TOOL}/jobs",
            data={"inputs": json.dumps(INPUT)},
            files=[("files", (name, b"{}", "application/json"))],
        )

        assert resp.status_code == 422
        assert client.k8s_jobs == []


class TestResults:
    def test_results_json_is_served_as_written(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put(TOOL, f"{job_id}/out/results.json", json.dumps(RESULTS).encode())
        _finish(client, job_id)

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 200
        assert resp.json()["results"] == RESULTS

    def test_results_that_break_the_output_schema_are_a_502(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put(TOOL, f"{job_id}/out/results.json", b'{"rows": [{"length": "3"}]}')
        _finish(client, job_id)

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 502
        body = resp.json()
        assert body["type"].endswith("/upstream-failure")
        assert body["schema_url"] == f"/v1/tools/{TOOL}/output-schema"
        # Same {pointer, detail} shape as every other /v1 problem
        assert {"pointer": "", "detail": "'artifacts' is a required property"} in body["errors"]
        assert {"pointer": "/rows/0/length", "detail": "'3' is not of type 'integer'"} in body["errors"]

    def test_results_that_are_not_json_are_a_502(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put(TOOL, f"{job_id}/out/results.json", b"NaN-ish {")
        _finish(client, job_id)

        assert client.get(f"/v1/jobs/{job_id}/results").status_code == 502

    @pytest.mark.parametrize("constant", [b"NaN", b"Infinity", b"-Infinity", b"1e999", b"-1e999"])
    def test_non_json_numbers_are_a_502_not_a_500(self, client, constant):
        """Python's json accepts these and jsonschema calls them numbers, but they are
        not JSON and cannot be rendered into a response."""
        job_id = _submit(client).json()["job_id"]
        # In a field the schema leaves unconstrained, so only the JSON parse can catch it
        results = b'{"rows": [{"length": 3, "score": ' + constant + b'}], "artifacts": []}'
        client.minio.put(TOOL, f"{job_id}/out/results.json", results)
        _finish(client, job_id)

        assert client.get(f"/v1/jobs/{job_id}/results").status_code == 502
        assert client.get(f"/{TOOL}/results/{job_id}").status_code == 502

    def test_deeply_nested_results_are_a_502_not_a_500(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put(TOOL, f"{job_id}/out/results.json", b"[" * 100_000 + b"]" * 100_000)
        _finish(client, job_id)

        assert client.get(f"/v1/jobs/{job_id}/results").status_code == 502
        assert client.get(f"/{TOOL}/results/{job_id}").status_code == 502

    def test_validation_too_deep_for_a_recursive_schema_is_a_502(self, client, monkeypatch):
        """Parses fine, but a $ref-recursive schema recurses once per level when validating."""
        registry = dict(tool_registry.get_registry())
        registry[TOOL] = {**registry[TOOL], "output_schema": {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "$ref": "#/$defs/node",
            "$defs": {"node": {"type": "array", "items": {"$ref": "#/$defs/node"}}},
        }}
        monkeypatch.setattr(tool_registry, "_REGISTRY", registry)
        job_id = _submit(client).json()["job_id"]
        client.minio.put(TOOL, f"{job_id}/out/results.json", b"[" * 500 + b"]" * 500)
        _finish(client, job_id)

        assert client.get(f"/v1/jobs/{job_id}/results").status_code == 502
        assert client.get(f"/{TOOL}/results/{job_id}").status_code == 502

    def test_one_huge_failing_value_is_truncated(self, client):
        job_id = _submit(client).json()["job_id"]
        rows = [{"length": "x" * 100_000}]
        client.minio.put(TOOL, f"{job_id}/out/results.json",
                         json.dumps({"rows": rows, "artifacts": []}).encode())
        _finish(client, job_id)

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 502
        assert len(resp.content) < 5_000

    def test_a_large_invalid_document_is_not_echoed_back(self, client):
        """jsonschema messages repr the failing value, which can be the whole document."""
        job_id = _submit(client).json()["job_id"]
        rows = [{"length": "x" * 100} for _ in range(2000)]
        client.minio.put(TOOL, f"{job_id}/out/results.json", json.dumps({"rows": rows}).encode())
        _finish(client, job_id)

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 502
        assert len(resp.content) < 20_000

    def test_a_completed_job_without_results_json_is_a_409(self, client):
        job_id = _submit(client).json()["job_id"]
        _finish(client, job_id)

        resp = client.get(f"/v1/jobs/{job_id}/results")

        assert resp.status_code == 409
        assert "no result file" in resp.json()["detail"]

    def test_the_legacy_results_route_serves_the_same_document(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put(TOOL, f"{job_id}/out/results.json", json.dumps(RESULTS).encode())

        resp = client.get(f"/{TOOL}/results/{job_id}")

        assert resp.status_code == 200
        assert resp.json() == RESULTS

    def test_listed_artifacts_download_by_name(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put(TOOL, f"{job_id}/out/results.csv", b"length\n3\n3\n")

        resp = client.get(f"/v1/jobs/{job_id}/artifacts/results.csv")

        assert resp.status_code == 200
        assert resp.content == b"length\n3\n3\n"


class TestDiscovery:
    def test_output_schema_is_served_with_an_id(self, client):
        body = client.get(f"/v1/tools/{TOOL}/output-schema").json()

        assert body["$id"].endswith(f"/v1/tools/{TOOL}/output-schema")
        assert body["required"] == ["rows", "artifacts"]

    def test_descriptor_links_the_output_schema(self, client):
        links = client.get(f"/v1/tools/{TOOL}").json()["links"]

        assert links["output_schema"].endswith(f"/v1/tools/{TOOL}/output-schema")

    def test_tools_without_the_contract_publish_no_output_schema(self, client):
        assert client.get("/v1/tools/somn").json()["links"]["output_schema"] is None
        assert client.get("/v1/tools/somn/output-schema").status_code == 404


class TestRegistryLoading:
    def test_an_output_schema_without_an_input_schema_is_rejected(self, tmp_path, monkeypatch):
        (tmp_path / "out.json").write_text(json.dumps(OUTPUT_SCHEMA))
        (tmp_path / "tools.yaml").write_text(yaml.safe_dump({"tools": {
            "molli": {"name": "m", "output_schema": "out.json"},
        }}))
        monkeypatch.setattr(tool_registry, "TOOLS_FILEPATH", tmp_path / "tools.yaml")
        monkeypatch.setattr(tool_registry, "SCHEMAS_DIR", tmp_path)

        with pytest.raises(tool_registry.ToolRegistryError, match="no input_schema"):
            tool_registry._build()
