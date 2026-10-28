"""Characterization tests for the legacy Files API.

As with the Jobs API tests, these pin CURRENT behavior. Assertions covering known
defects carry a `DEFECT:` comment.

Note the route parameter is named `bucket_name`: the public URL exposes the storage
bucket, which happens to equal the job type. That coupling is itself part of what /v1
replaces.
"""

import pytest

NOVOSTOIC_OPTSTOIC = "novostoic-optstoic"


class TestUpload:
    def test_upload_mints_a_job_id_when_none_is_given(self, client):
        resp = client.post(
            f"/{NOVOSTOIC_OPTSTOIC}/upload",
            files={"file": ("input.txt", b"hello", "text/plain")},
        )

        assert resp.status_code == 200
        body = resp.json()
        assert len(body["jobID"]) == 32
        assert "uploaded_at" in body

    def test_upload_honours_a_caller_supplied_job_id(self, client):
        resp = client.post(
            f"/{NOVOSTOIC_OPTSTOIC}/upload?job_id=my-job",
            files={"file": ("input.txt", b"hello", "text/plain")},
        )

        assert resp.json()["jobID"] == "my-job"
        assert client.minio.get_file(NOVOSTOIC_OPTSTOIC, "my-job/in/input.txt") == b"hello"

    def test_upload_for_an_existing_job_is_refused(self, client):
        """#150: any caller who knew a job id could replace that job's input files."""
        client.post("/mutagenesis/upload?job_id=victim",
                    files={"file": ("mutations.csv", b"A1G", "text/csv")})
        client.post("/defaults/jobs", json={"job_id": "victim"})

        resp = client.post("/mutagenesis/upload?job_id=victim",
                           files={"file": ("mutations.csv", b"EVIL", "text/csv")})

        assert resp.status_code == 409
        assert client.minio.get_file("mutagenesis", "victim/in/mutations.csv") == b"A1G"

    def test_several_uploads_before_the_job_exists_still_work(self, client):
        """The legitimate flow: every frontend uploads all its files, then creates."""
        first = client.post("/mutagenesis/upload",
                            files={"file": ("a.csv", b"A1G", "text/csv")}).json()["jobID"]
        second = client.post(f"/mutagenesis/upload?job_id={first}",
                             files={"file": ("b.csv", b"B2C", "text/csv")})

        assert second.status_code == 200
        assert client.minio.get_file("mutagenesis", f"{first}/in/b.csv") == b"B2C"

    @pytest.mark.parametrize("job_id", ["victim/in/x", "victim/in/.."])
    def test_a_job_id_that_is_a_path_is_refused(self, client, job_id):
        """The existence check matches job_id exactly; MinIO sees {job_id}/in/{name}."""
        client.post("/mutagenesis/upload?job_id=victim", files={"file": ("f", b"A1G", "text/csv")})
        client.post("/defaults/jobs", json={"job_id": "victim"})

        resp = client.post("/mutagenesis/upload", params={"job_id": job_id},
                           files={"file": ("f", b"EVIL", "text/csv")})

        assert resp.status_code == 400
        assert client.minio.get_file("mutagenesis", "victim/in/f") == b"A1G"
        assert not [k for k in client.minio.objects["mutagenesis"] if k != "victim/in/f"]

    @pytest.mark.parametrize("name", ["../in/f", "x/../f", "a\\b", ".."])
    def test_a_file_name_that_is_a_path_is_refused(self, client, name):
        resp = client.post("/mutagenesis/upload?job_id=fresh",
                           files={"file": (name, b"EVIL", "text/csv")})

        assert resp.status_code == 400
        assert not client.minio.objects.get("mutagenesis")

    def test_chemscraper_uploads_must_be_pdfs(self, client):
        resp = client.post(
            "/chemscraper/upload",
            files={"file": ("not-a.pdf", b"just text", "application/pdf")},
        )

        # DEFECT: a rejected upload is a client error, but this reports 500.
        assert resp.status_code == 500
        assert resp.json() == {"error": "Unable to upload file"}

    def test_chemscraper_accepts_a_real_pdf(self, client):
        resp = client.post(
            "/chemscraper/upload",
            files={"file": ("real.pdf", b"%PDF-1.4 ...", "application/pdf")},
        )

        assert resp.status_code == 200

    def test_any_other_bucket_accepts_any_content(self, client):
        """Only chemscraper validates content; every other bucket takes anything."""
        resp = client.post(
            f"/{NOVOSTOIC_OPTSTOIC}/upload",
            files={"file": ("anything.bin", b"\x00\x01\x02", "application/octet-stream")},
        )

        assert resp.status_code == 200


class TestResults:
    def test_missing_output_returns_200_with_a_null_body(self, client):
        client.post(f"/{NOVOSTOIC_OPTSTOIC}/jobs", json={"job_id": "j1"})

        resp = client.get(f"/{NOVOSTOIC_OPTSTOIC}/results/j1")

        # DEFECT: a polling client cannot distinguish "still running" from "finished
        # with no results" from "wrong job id" -- all three are 200 null.
        assert resp.status_code == 200
        assert resp.json() is None

    def test_unknown_job_id_returns_404(self, client):
        resp = client.get(f"/{NOVOSTOIC_OPTSTOIC}/results/never-existed")

        # Was: the service `return`ed (rather than `raise`d) an HTTPException, so
        # FastAPI serialized the exception object as an ordinary response body -- HTTP
        # 200, with the 404 buried inside the payload where no standard client looks.
        assert resp.status_code == 404
        assert resp.json()["detail"] == "Job not found"

    def test_a_known_job_with_no_output_still_returns_200_null(self, client):
        """Deliberately unchanged.

        Frontends poll this endpoint while a job runs and treat a null body as "not
        ready yet". Turning that into a 404 or a 409 is the right answer, but it is a
        contract change for every existing client, so it belongs with the versioned
        API rather than in a security fix.
        """
        client.post(f"/{NOVOSTOIC_OPTSTOIC}/jobs", json={"job_id": "j1"})

        resp = client.get(f"/{NOVOSTOIC_OPTSTOIC}/results/j1")

        assert resp.status_code == 200
        assert resp.json() is None

    def test_invalid_bucket_returns_400(self, client):
        resp = client.get("/not-a-real-tool/results/j1")

        assert resp.status_code == 400
        assert "Invalid job type" in resp.json()["detail"]


class TestInputsAndErrors:
    def test_inputs_returns_404_when_nothing_was_uploaded(self, client):
        resp = client.get(f"/{NOVOSTOIC_OPTSTOIC}/inputs/j1")

        assert resp.status_code == 404

    def test_inputs_lists_uploaded_files(self, client):
        client.minio.put(NOVOSTOIC_OPTSTOIC, "j1/in/a.txt", b"a")

        resp = client.get(f"/{NOVOSTOIC_OPTSTOIC}/inputs/j1")

        assert resp.status_code == 200
        assert len(resp.json()) == 1

    def test_errors_returns_404_when_absent(self, client):
        assert client.get(f"/{NOVOSTOIC_OPTSTOIC}/errors/j1").status_code == 404

    def test_errors_returns_the_error_file(self, client):
        client.minio.put(NOVOSTOIC_OPTSTOIC, "j1/errors.txt", b"boom")

        resp = client.get(f"/{NOVOSTOIC_OPTSTOIC}/errors/j1")

        assert resp.status_code == 200
