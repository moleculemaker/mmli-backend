"""Characterization tests for the deprecated ChemScraper analyze endpoint.

`/chemscraper/analyze` predates ChemScraper becoming a normal Kubernetes job and is
marked deprecated, but it is still mounted and still reachable.
"""
import pytest

from services.chemscraper_service import ChemScraperService


@pytest.fixture
def stub_analysis(monkeypatch):
    """Neutralize the background task.

    TestClient drains background tasks synchronously before returning, so without this
    the endpoint's response depends on the whole extraction pipeline: it would look for
    an uploaded PDF and then POST it to the external ChemScraper service. These tests
    are about what the endpoint accepts, not about that pipeline.
    """
    async def _noop(self, *args, **kwargs):
        return True

    monkeypatch.setattr(ChemScraperService, "runChemscraperOnDocument", _noop)


def _analyze(client, **overrides):
    body = {"jobId": "j1", "user_email": "a@b.edu", "fileList": ["paper.pdf"]}
    body.update(overrides)
    return client.post("/chemscraper/analyze", json=body)


class TestAnalyzeValidation:
    def test_empty_file_list_returns_400(self, client):
        resp = _analyze(client, fileList=[])

        # Was: the function fell out of its `if` with no return statement at all, which
        # FastAPI renders as HTTP 200 with a null body -- reporting success for a
        # request that was never accepted, for a job that would never run.
        assert resp.status_code == 400
        assert "fileList" in resp.json()["detail"]

    def test_empty_job_id_returns_400(self, client):
        resp = _analyze(client, jobId="")

        assert resp.status_code == 400
        assert "jobId" in resp.json()["detail"]

    def test_a_valid_request_is_accepted(self, client, stub_analysis):
        resp = _analyze(client)

        assert resp.status_code == 202
        assert resp.json()["jobId"] == "j1"

    def test_an_accepted_request_creates_a_processing_job(self, client, stub_analysis):
        _analyze(client)

        jobs = client.get("/chemscraper/jobs/j1").json()

        assert len(jobs) == 1
        assert jobs[0]["phase"] == "processing"
