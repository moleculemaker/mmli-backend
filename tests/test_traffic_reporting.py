"""Tests for caller attribution and the usage report."""
import json
import time

import pytest
import sqlmodel

from models.sqlmodel.models import Job
from services import analytics
from routers import reports

MCP_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json, text/event-stream",
}
OPTSTOIC = {"primary_precursor": "MNXM1137670", "target_molecule": "MNXM26"}


def _row(client, job_id):
    engine = sqlmodel.create_engine(client.sync_db_url)
    with sqlmodel.Session(engine) as session:
        return session.get(Job, job_id)


def _submit_v1(client, **kwargs):
    return client.post("/v1/tools/novostoic-optstoic/jobs", json=OPTSTOIC, **kwargs)


def _submit_mcp(client):
    response = client.post("/mcp", json={
        "jsonrpc": "2.0", "id": 1, "method": "tools/call",
        "params": {"name": "submit_novostoic_optstoic", "arguments": OPTSTOIC},
    }, headers=MCP_HEADERS)
    return json.loads(response.json()["result"]["content"][0]["text"])["job_id"]


class TestSurfaceAttribution:
    def test_legacy_submissions_are_labeled_legacy(self, client):
        client.post("/defaults/jobs", json={"job_id": "j1"})

        assert _row(client, "j1").client_surface == "legacy"

    def test_v1_submissions_are_labeled_v1(self, client):
        job_id = _submit_v1(client).json()["job_id"]

        assert _row(client, job_id).client_surface == "v1"

    def test_mcp_submissions_are_labeled_mcp(self, client):
        """Without this the adapter's own httpx client is all that gets recorded."""
        job_id = _submit_mcp(client)

        assert _row(client, job_id).client_surface == "mcp"

    def test_an_external_caller_cannot_claim_to_be_an_agent(self, client):
        """The surface header is only honored with the process-local token.

        Otherwise anyone could label their traffic as agent traffic and quietly corrupt
        the adoption figures this exists to produce.
        """
        job_id = _submit_v1(client, headers={
            analytics.SURFACE_HEADER: "mcp",
            analytics.SURFACE_TOKEN_HEADER: "guessed",
        }).json()["job_id"]

        assert _row(client, job_id).client_surface == "v1"

    def test_a_missing_token_is_not_honored(self, client):
        job_id = _submit_v1(client, headers={analytics.SURFACE_HEADER: "mcp"}).json()["job_id"]

        assert _row(client, job_id).client_surface == "v1"


class TestCallerFields:
    def test_user_agent_is_recorded_on_legacy_submissions(self, client):
        """It was a TODO, so every historical legacy job has an empty string."""
        client.post("/defaults/jobs", json={"job_id": "j1"},
                    headers={"User-Agent": "curl/8.0"})

        assert _row(client, "j1").user_agent == "curl/8.0"

    def test_origin_is_recorded_when_sent(self, client):
        job_id = _submit_v1(client, headers={
            "Origin": "https://somn.frontend.mmli1.ncsa.illinois.edu"}).json()["job_id"]

        assert _row(client, job_id).client_origin.endswith("ncsa.illinois.edu")

    def test_a_script_sends_no_origin_and_that_is_the_signal(self, client):
        job_id = _submit_v1(client).json()["job_id"]

        assert _row(client, job_id).client_origin is None


class TestFingerprint:
    def test_disabled_by_default(self, client):
        """No salt configured means no fingerprint, rather than a weak one."""
        job_id = _submit_v1(client).json()["job_id"]

        assert _row(client, job_id).client_fingerprint is None

    def test_is_derived_when_a_salt_is_configured(self, monkeypatch):
        monkeypatch.setattr(analytics, "ANALYTICS_SALT", "a-real-secret")

        assert analytics.client_fingerprint("192.0.2.10") is not None

    def test_the_address_itself_is_never_the_stored_value(self, monkeypatch):
        monkeypatch.setattr(analytics, "ANALYTICS_SALT", "a-real-secret")

        fingerprint = analytics.client_fingerprint("192.0.2.10")

        assert "192.0.2.10" not in fingerprint

    def test_the_same_client_is_stable_within_a_period(self, monkeypatch):
        """Stability is what makes a distinct-client count possible at all."""
        monkeypatch.setattr(analytics, "ANALYTICS_SALT", "a-real-secret")
        january = time.mktime(time.strptime("2026-01-05", "%Y-%m-%d"))
        december = time.mktime(time.strptime("2026-12-20", "%Y-%m-%d"))

        assert (analytics.client_fingerprint("192.0.2.10", january)
                == analytics.client_fingerprint("192.0.2.10", december))

    def test_it_changes_between_years(self, monkeypatch):
        """Annual rotation bounds how long any pseudonym stays linkable."""
        monkeypatch.setattr(analytics, "ANALYTICS_SALT", "a-real-secret")
        this_year = time.mktime(time.strptime("2026-06-01", "%Y-%m-%d"))
        next_year = time.mktime(time.strptime("2027-06-01", "%Y-%m-%d"))

        assert (analytics.client_fingerprint("192.0.2.10", this_year)
                != analytics.client_fingerprint("192.0.2.10", next_year))

    def test_different_clients_differ(self, monkeypatch):
        monkeypatch.setattr(analytics, "ANALYTICS_SALT", "a-real-secret")

        assert (analytics.client_fingerprint("192.0.2.10")
                != analytics.client_fingerprint("192.0.2.11"))


class TestNotifyEmail:
    def test_v1_records_the_notify_email(self, client):
        """Without it a /v1 job runs for hours and tells nobody it finished."""
        job_id = _submit_v1(client, headers={
            "X-Notify-Email": "a@illinois.edu"}).json()["job_id"]

        assert _row(client, job_id).email == "a@illinois.edu"

    def test_it_is_optional(self, client):
        job_id = _submit_v1(client).json()["job_id"]

        assert _row(client, job_id).email is None

    def test_it_is_not_part_of_the_input_document(self, client):
        """Sent as a header so the body stays exactly the published schema."""
        response = client.post("/v1/tools/novostoic-optstoic/jobs",
                               json={**OPTSTOIC, "notify_email": "a@illinois.edu"})

        assert response.status_code == 422


class TestReportAccessControl:
    def _authorize(self, monkeypatch, groups):
        monkeypatch.setattr(reports, "REPORTING_GROUP", "mmli-reporting")
        monkeypatch.setattr(reports, "validate_auth_cookie",
                            lambda request: {"email": "a@b.edu", "groups": groups})

    def test_group_members_are_allowed(self, client, monkeypatch):
        self._authorize(monkeypatch, ["mmli-reporting"])

        assert client.get("/internal/reports/usage").status_code == 200

    def test_authenticated_non_members_are_refused(self, client, monkeypatch):
        """Resolving a user proves only that somebody is logged in."""
        self._authorize(monkeypatch, ["some-other-group"])

        assert client.get("/internal/reports/usage").status_code == 403

    def test_unauthenticated_callers_are_refused(self, client, monkeypatch):
        monkeypatch.setattr(reports, "REPORTING_GROUP", "mmli-reporting")
        monkeypatch.setattr(reports, "validate_auth_cookie", lambda request: None)

        assert client.get("/internal/reports/usage").status_code == 401

    def test_an_unconfigured_group_refuses_rather_than_falls_open(self, client, monkeypatch):
        """An unset group must not silently mean "any authenticated user"."""
        monkeypatch.setattr(reports, "REPORTING_GROUP", "")
        monkeypatch.setattr(reports, "validate_auth_cookie",
                            lambda request: {"email": "a@b.edu", "groups": ["anything"]})

        assert client.get("/internal/reports/usage").status_code == 503


class TestReportContent:
    @pytest.fixture
    def authorized(self, monkeypatch):
        monkeypatch.setattr(reports, "REPORTING_GROUP", "mmli-reporting")
        monkeypatch.setattr(reports, "validate_auth_cookie",
                            lambda request: {"email": "a@b.edu", "groups": ["mmli-reporting"]})

    def test_counts_jobs_by_tool(self, client, authorized):
        _submit_v1(client)
        _submit_v1(client)

        body = client.get("/internal/reports/usage").json()

        by_tool = {row["tool"]: row for row in body["by_tool"]}
        assert by_tool["novostoic-optstoic"]["jobs"] == 2
        assert body["total_jobs"] == 2

    def test_counts_by_surface(self, client, authorized):
        _submit_v1(client)
        _submit_mcp(client)
        client.post("/defaults/jobs", json={"job_id": "legacy-1"})

        body = client.get("/internal/reports/usage").json()

        assert body["by_surface"] == {"v1": 1, "mcp": 1, "legacy": 1}

    def test_counts_identified_users_separately_from_jobs(self, client, authorized):
        _submit_v1(client, headers={"X-Notify-Email": "a@illinois.edu"})
        _submit_v1(client, headers={"X-Notify-Email": "a@illinois.edu"})
        _submit_v1(client)

        row = {r["tool"]: r for r in client.get("/internal/reports/usage").json()["by_tool"]}

        assert row["novostoic-optstoic"]["jobs"] == 3
        assert row["novostoic-optstoic"]["identified_users"] == 1

    def test_returns_no_raw_rows_and_no_addresses(self, client, authorized):
        """The constraint that keeps this from becoming a user export."""
        _submit_v1(client, headers={"X-Notify-Email": "secret@illinois.edu"})

        body = client.get("/internal/reports/usage").text

        assert "secret@illinois.edu" not in body
        assert "job_id" not in body

    def test_rejects_a_malformed_date(self, client, authorized):
        assert client.get("/internal/reports/usage?from=last-tuesday").status_code == 400

    def test_rejects_an_inverted_range(self, client, authorized):
        response = client.get("/internal/reports/usage?from=2026-06-01&to=2026-01-01")

        assert response.status_code == 400

    def test_a_period_with_no_traffic_reports_zero(self, client, authorized):
        _submit_v1(client)

        body = client.get("/internal/reports/usage?from=2020-01-01&to=2020-12-31").json()

        assert body["total_jobs"] == 0
