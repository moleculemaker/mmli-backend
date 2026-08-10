"""Tests for the MCP server.

Driven over the wire as JSON-RPC against the mounted /mcp endpoint, rather than by
calling the handlers directly, so that what is verified is what an agent would actually
receive -- protocol framing included.
"""
import json

import pytest

from routers import mcp_app
from services import tool_registry

MCP_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json, text/event-stream",
}

OPTSTOIC = {"primary_precursor": "MNXM1137670", "target_molecule": "MNXM26"}


def _rpc(client, method, params=None, request_id=1):
    payload = {"jsonrpc": "2.0", "id": request_id, "method": method}
    if params is not None:
        payload["params"] = params
    response = client.post("/mcp", json=payload, headers=MCP_HEADERS)
    return response.json()


def _call(client, name, arguments):
    return _rpc(client, "tools/call", {"name": name, "arguments": arguments})


def _payload(result):
    """Unwrap the JSON an MCP tool returns inside its text content block."""
    return json.loads(result["result"]["content"][0]["text"])


@pytest.fixture(autouse=True)
def reset_rate_limit():
    mcp_app._submit_history.clear()
    yield
    mcp_app._submit_history.clear()


class TestProtocol:
    def test_initialize_succeeds(self, client):
        body = _rpc(client, "initialize", {
            "protocolVersion": "2025-06-18", "capabilities": {},
            "clientInfo": {"name": "test", "version": "1"},
        })

        assert body["result"]["serverInfo"]["name"] == "alphasynthesis-tools"

    def test_server_advertises_tool_support(self, client):
        body = _rpc(client, "initialize", {
            "protocolVersion": "2025-06-18", "capabilities": {},
            "clientInfo": {"name": "test", "version": "1"},
        })

        assert "tools" in body["result"]["capabilities"]


class TestToolListing:
    def test_lists_a_submit_tool_per_registered_tool(self, client):
        tools = _rpc(client, "tools/list")["result"]["tools"]
        names = {t["name"] for t in tools}

        expected = {
            mcp_app._submit_tool_name(identifier)
            for identifier, entry in tool_registry.list_tools().items()
            if entry.get("input_schema") is not None
        }
        assert expected <= names

    def test_internal_subjobs_are_not_exposed(self, client):
        names = {t["name"] for t in _rpc(client, "tools/list")["result"]["tools"]}

        assert "submit_ezspec_unidock" not in names
        assert "submit_ezspec_inference" not in names

    def test_lifecycle_tools_are_present(self, client):
        names = {t["name"] for t in _rpc(client, "tools/list")["result"]["tools"]}

        assert {"get_job_status", "get_job_results", "list_job_artifacts",
                "cancel_job", "describe_tool"} <= names

    def test_every_tool_carries_an_object_input_schema(self, client):
        """MCP requires an object at the root of inputSchema."""
        for tool in _rpc(client, "tools/list")["result"]["tools"]:
            assert tool["inputSchema"]["type"] == "object", tool["name"]

    def test_a_submit_tool_carries_the_published_schema(self, client):
        tools = {t["name"]: t for t in _rpc(client, "tools/list")["result"]["tools"]}

        schema = tools["submit_novostoic_optstoic"]["inputSchema"]
        assert set(schema["required"]) == {"primary_precursor", "target_molecule"}
        # Property descriptions are what an agent reads to fill the call in.
        assert schema["properties"]["primary_precursor"]["description"]

    @pytest.mark.parametrize("tool_name,expected", [
        # SOMN's published schema is a top-level array.
        ("submit_somn", "array"),
        # cleandb-mepesm's is an anyOf accepting an object or a bare string, so it has
        # no top-level type at all.
        ("submit_cleandb_mepesm", None),
    ])
    def test_a_non_object_schema_is_nested_rather_than_reshaped(self, client, tool_name, expected):
        """MCP needs an object at the root; nesting keeps the published schema authoritative."""
        tools = {t["name"]: t for t in _rpc(client, "tools/list")["result"]["tools"]}

        schema = tools[tool_name]["inputSchema"]
        assert schema["type"] == "object"
        assert schema["required"] == ["inputs"]
        assert schema["properties"]["inputs"].get("type") == expected

    def test_a_nested_schema_still_submits_correctly(self, client):
        body = _payload(_call(client, "submit_cleandb_mepesm",
                              {"inputs": {"sequence": "MEDIPDTSRPPLKYVK"}}))

        assert body["tool"] == "cleandb-mepesm"

    def test_tools_needing_file_uploads_say_so(self, client):
        """MCP has no upload mechanism, so an agent must be told before it tries."""
        tools = {t["name"]: t for t in _rpc(client, "tools/list")["result"]["tools"]}

        assert "use the HTTP API" in tools["submit_molli"]["description"]

    def test_descriptions_tell_an_agent_the_job_is_asynchronous(self, client):
        tools = {t["name"]: t for t in _rpc(client, "tools/list")["result"]["tools"]}

        assert "get_job_status" in tools["submit_aceretro"]["description"]


class TestSubmission:
    def test_submitting_returns_a_job_id(self, client):
        body = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))

        assert body["status"] == "queued"
        assert body["job_id"]

    def test_submitting_starts_a_kubernetes_job(self, client):
        _call(client, "submit_novostoic_optstoic", OPTSTOIC)

        assert len(client.k8s_jobs) == 1

    def test_array_shaped_input_is_unwrapped_before_submission(self, client):
        somn = [{
            "reactant_pair_name": "p1",
            "el": "CC(C)c1ccc(Br)cc1", "el_name": "e", "el_input_type": "smi", "el_idx": "0",
            "nuc": "NCc1ccccc1", "nuc_name": "n", "nuc_input_type": "smi", "nuc_idx": "0",
        }]

        body = _payload(_call(client, "submit_somn", {"inputs": somn}))

        assert body["tool"] == "somn"

    def test_invalid_input_surfaces_the_validation_problem(self, client):
        """The agent needs to see which field was wrong so it can correct itself."""
        result = _call(client, "submit_novostoic_optstoic", {"primary_precursor": "X"})

        assert result["result"]["isError"] is True
        text = result["result"]["content"][0]["text"]
        assert "target_molecule" in text

    def test_invalid_input_starts_no_job(self, client):
        _call(client, "submit_novostoic_optstoic", {})

        assert client.k8s_jobs == []


class TestJobLifecycle:
    def test_status_is_retrievable(self, client):
        job_id = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))["job_id"]

        body = _payload(_call(client, "get_job_status", {"job_id": job_id}))

        assert body["job_id"] == job_id
        assert body["status"] == "queued"

    def test_results_report_that_the_job_is_unfinished(self, client):
        job_id = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))["job_id"]

        result = _call(client, "get_job_results", {"job_id": job_id})

        assert result["result"]["isError"] is True
        assert "job-not-finished" in result["result"]["content"][0]["text"]

    def test_artifacts_can_be_listed(self, client):
        job_id = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))["job_id"]
        client.minio.put("novostoic-optstoic", f"{job_id}/out/output.json", b"{}")

        body = _payload(_call(client, "list_job_artifacts", {"job_id": job_id}))

        assert [i["name"] for i in body["items"]] == ["output.json"]

    def test_logs_are_excluded_unless_requested(self, client):
        job_id = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))["job_id"]
        client.minio.put("novostoic-optstoic", f"{job_id}/out/error.log", b"boom")

        without = _payload(_call(client, "list_job_artifacts", {"job_id": job_id}))
        with_logs = _payload(_call(client, "list_job_artifacts",
                                   {"job_id": job_id, "include_logs": True}))

        assert without["count"] == 0
        assert with_logs["count"] == 1

    def test_a_job_can_be_canceled(self, client, deleted_k8s_jobs):
        job_id = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))["job_id"]

        body = _payload(_call(client, "cancel_job", {"job_id": job_id}))

        assert body["status"] == "canceled"

    def test_unknown_job_is_an_error_not_a_silent_empty_result(self, client):
        result = _call(client, "get_job_status", {"job_id": "never-existed"})

        assert result["result"]["isError"] is True
        assert "unknown-job" in result["result"]["content"][0]["text"]


class TestDescribeTool:
    def test_returns_metadata_including_provenance_fields(self, client):
        body = _payload(_call(client, "describe_tool", {"tool": "somn"}))

        assert body["identifier"] == "somn"
        # Null rather than absent: nobody has supplied them yet.
        assert "license" in body and "citation" in body

    def test_unknown_tool_is_an_error(self, client):
        result = _call(client, "describe_tool", {"tool": "not-a-real-tool"})

        assert result["result"]["isError"] is True


class TestRateLimiting:
    def test_submissions_are_limited(self, client):
        """The ingress rule cannot see individual MCP calls, so this is enforced here."""
        for _ in range(mcp_app.SUBMIT_LIMIT):
            assert _call(client, "submit_novostoic_optstoic", OPTSTOIC)["result"].get("isError") is not True

        result = _call(client, "submit_novostoic_optstoic", OPTSTOIC)

        assert result["result"]["isError"] is True
        assert "Rate limit exceeded" in result["result"]["content"][0]["text"]

    def test_the_limit_does_not_apply_to_reads(self, client):
        job_id = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))["job_id"]
        mcp_app._submit_history.clear()

        # A polling agent must never be throttled.
        for _ in range(mcp_app.SUBMIT_LIMIT * 3):
            result = _call(client, "get_job_status", {"job_id": job_id})
            assert result["result"].get("isError") is not True

    def test_the_limit_message_tells_the_agent_what_to_do_instead(self, client):
        for _ in range(mcp_app.SUBMIT_LIMIT + 1):
            result = _call(client, "submit_novostoic_optstoic", OPTSTOIC)

        assert "Poll get_job_status" in result["result"]["content"][0]["text"]


class TestAdapterFidelity:
    def test_mcp_and_v1_produce_the_same_job_document(self, client):
        """MCP is an adapter over /v1, so the two must not drift apart."""
        via_mcp = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))
        via_http = client.post("/v1/tools/novostoic-optstoic/jobs", json=OPTSTOIC).json()

        assert set(via_mcp) == set(via_http)
        assert via_mcp["tool"] == via_http["tool"]

    def test_mcp_submission_is_visible_over_http(self, client):
        job_id = _payload(_call(client, "submit_novostoic_optstoic", OPTSTOIC))["job_id"]

        assert client.get(f"/v1/jobs/{job_id}").status_code == 200
