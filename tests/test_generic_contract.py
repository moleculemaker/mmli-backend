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
        assert any("artifacts" in e for e in body["errors"])
        assert any("'3' is not of type 'integer'" in e for e in body["errors"])

    def test_results_that_are_not_json_are_a_502(self, client):
        job_id = _submit(client).json()["job_id"]
        client.minio.put(TOOL, f"{job_id}/out/results.json", b"NaN-ish {")
        _finish(client, job_id)

        assert client.get(f"/v1/jobs/{job_id}/results").status_code == 502

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
