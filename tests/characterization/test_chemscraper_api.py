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


class TestAnalyzeLeavesExistingJobsAlone:
    """#150's threat model on this route: it has its own public ingress, takes a
    caller-chosen jobId, and its background task rewrites the row's phase and times."""

    def _seed(self, sync_db_url, job_type, phase):
        from sqlmodel import Session, create_engine
        from models.sqlmodel.models import Job
        engine = create_engine(sync_db_url)
        with Session(engine) as session:
            session.add(Job(job_id="victim", type=job_type, phase=phase,
                            time_created=0, time_start=1, time_end=2, user_agent="", deleted=0))
            session.commit()
        engine.dispose()

    def _row(self, sync_db_url):
        from sqlmodel import Session, create_engine
        from models.sqlmodel.models import Job
        engine = create_engine(sync_db_url)
        with Session(engine) as session:
            job = session.get(Job, "victim")
            row = (job.type, job.phase, job.time_start, job.time_end)
        engine.dispose()
        return row

    @pytest.mark.parametrize("job_type", ["aceretro", "chemscraper"])
    def test_an_existing_job_is_refused_and_untouched(self, client, sync_db_url, job_type):
        """Unstubbed: the real background task is what flipped a completed job to error."""
        from models.enums import JobStatus
        self._seed(sync_db_url, job_type, JobStatus.COMPLETED)

        try:
            resp = _analyze(client, jobId="victim", fileList=["nope.pdf"], user_email="")
        except Exception:
            # TestClient re-raises the background task's failure (no such PDF); what
            # matters is what it did to the row on the way.
            resp = None

        assert self._row(sync_db_url) == (job_type, JobStatus.COMPLETED, 1, 2)
        assert resp is not None and resp.status_code == 409
