"""Notification addresses are validated at submission, on every route that stores one.

The address becomes the SMTP envelope recipient and the `To:` header. These inputs are
the ones issue #137 sent through the real EmailService to see what reached the server.
"""
import pytest
import sqlmodel

from models.sqlmodel.models import Job
from services.chemscraper_service import ChemScraperService
from services.notify_email import InvalidNotifyEmail, normalize_notify_email

OPTSTOIC = {"primary_precursor": "MNXM1137670", "target_molecule": "MNXM26"}

BAD_ADDRESSES = [
    pytest.param("victim@example.org\r\nBcc: x@example.com\r\nSubject: injected", id="header-injection"),
    pytest.param("victim@example.org\r\n other@example.com", id="folded"),
    pytest.param("victim@example.org, other@example.com", id="list"),
    pytest.param('"other@example.com" <victim@example.org>', id="display-name"),
    pytest.param("victim@example.org\x00other", id="nul"),
    pytest.param("vic\x7ftim@example.org", id="del"),
    pytest.param("a@" + "b" * 250 + ".org", id="over-254"),
    pytest.param("not-an-address", id="no-at-sign"),
]


def _row(client, job_id):
    engine = sqlmodel.create_engine(client.sync_db_url)
    with sqlmodel.Session(engine) as session:
        return session.get(Job, job_id)


def _rows(client):
    engine = sqlmodel.create_engine(client.sync_db_url)
    with sqlmodel.Session(engine) as session:
        return session.exec(sqlmodel.select(Job)).all()


class TestNormalize:
    @pytest.mark.parametrize("address", BAD_ADDRESSES)
    def test_rejects(self, address):
        with pytest.raises(InvalidNotifyEmail):
            normalize_notify_email(address)

    @pytest.mark.parametrize("blank", [None, "", "   "])
    def test_blank_means_no_address(self, blank):
        assert normalize_notify_email(blank) is None

    def test_trims_and_lowercases_the_domain_only(self):
        assert normalize_notify_email("  Victim@EXAMPLE.Org ") == "Victim@example.org"

    def test_internationalized_domain_is_stored_as_ascii(self):
        assert normalize_notify_email("a@münchen.de") == "a@xn--mnchen-3ya.de"

    def test_non_ascii_local_part_is_rejected(self):
        with pytest.raises(InvalidNotifyEmail):
            normalize_notify_email("mü@example.org")


class TestLegacyJobs:
    @pytest.mark.parametrize("address", BAD_ADDRESSES)
    def test_rejects_before_anything_is_created(self, client, address):
        resp = client.post("/defaults/jobs", json={"job_id": "j1", "email": address})

        assert resp.status_code == 422
        assert _row(client, "j1") is None
        assert client.k8s_jobs == []

    def test_stores_the_normalized_address(self, client):
        client.post("/defaults/jobs", json={"job_id": "j1", "email": " a@B.EDU "})

        assert _row(client, "j1").email == "a@b.edu"

    def test_blank_is_stored_as_no_address(self, client):
        client.post("/defaults/jobs", json={"job_id": "j1", "email": ""})

        assert _row(client, "j1").email is None


class TestV1Jobs:
    @pytest.mark.parametrize("address", BAD_ADDRESSES)
    def test_rejects_before_anything_is_created(self, client, address):
        resp = client.post("/v1/tools/novostoic-optstoic/jobs", json=OPTSTOIC,
                           headers={"X-Notify-Email": address})

        assert resp.status_code == 422
        assert resp.json()["type"].endswith("/invalid-input")
        assert _rows(client) == []
        assert client.k8s_jobs == []

    def test_stores_the_normalized_address(self, client):
        job_id = client.post("/v1/tools/novostoic-optstoic/jobs", json=OPTSTOIC,
                             headers={"X-Notify-Email": "a@B.EDU"}).json()["job_id"]

        assert _row(client, job_id).email == "a@b.edu"


class TestChemscraperAnalyze:
    @pytest.fixture(autouse=True)
    def stub_analysis(self, monkeypatch):
        async def _noop(self, *args, **kwargs):
            return True

        monkeypatch.setattr(ChemScraperService, "runChemscraperOnDocument", _noop)

    def _analyze(self, client, address):
        return client.post("/chemscraper/analyze", json={
            "jobId": "j1", "user_email": address, "fileList": ["paper.pdf"]})

    @pytest.mark.parametrize("address", BAD_ADDRESSES)
    def test_rejects_before_anything_is_created(self, client, address):
        assert self._analyze(client, address).status_code == 422
        assert _row(client, "j1") is None

    def test_stores_the_normalized_address(self, client):
        self._analyze(client, "a@B.EDU")

        assert _row(client, "j1").email == "a@b.edu"
