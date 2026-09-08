"""Regression tests for the ChemScraper results CSV round-trip.

A Molecule is never read back from the database. It is written to a results CSV in
MinIO and rebuilt from that CSV on every read, so its field types survive only as far
as pandas' dtype inference: doc_no, page_no, PubChemCID and molecularWeight are all
declared `str` and all come back as int64 or float64.

Pydantic 1 coerced those to strings. Pydantic 2 rejects them by default, which made
GET /chemscraper/results/{job_id} return 500 for every job whose results CSV had
already been written.
"""
import pytest
from pydantic import ValidationError

from models.molecule import Molecule

CHEMSCRAPER = "chemscraper"

# Column order as ChemScraperService.fetchExternalDataAndStoreResults writes it: a
# DataFrame built from Molecule.model_dump(), so the leading unnamed column is the
# index. doc_no, page_no, PubChemCID and OtherInstances hold digits and molecularWeight
# holds a decimal, which is what makes pandas infer numbers for them on the way back in.
RESULTS_CSV = (
    ",id,flagged,atom_count,doc_no,file_path,page_no,name,SMILE,structure,minX,minY,"
    "width,height,PubChemCID,molecularFormula,molecularWeight,chemicalSafety,"
    "Description,Location,OtherInstances,fingerprint\n"
    "0,0,False,3,1,/inputs/paper.pdf,1,Ethanol,CCO,<svg/>,1609,1949,201,108,702,"
    "C2H6O,46.07,,Unavailable, | page: 1,1,abc\n"
    "1,1,False,1,1,/inputs/paper.pdf,2,Methanol,CO,<svg/>,10,20,30,40,887,"
    "CH4O,32.04,,Unavailable, | page: 2,2,def\n"
)


class TestResultsEndpoint:
    def test_a_results_csv_with_numeric_columns_is_served(self, client):
        client.minio.put(CHEMSCRAPER, "j1/out/j1-results.csv", RESULTS_CSV.encode())

        resp = client.get(f"/{CHEMSCRAPER}/results/j1")

        # Was: 500, "2 validation errors for Molecule -- doc_no, page_no: Input should
        # be a valid string [input_type=int]".
        assert resp.status_code == 200
        assert len(resp.json()) == 2

    def test_numeric_columns_are_served_as_the_strings_they_were_written_as(self, client):
        client.minio.put(CHEMSCRAPER, "j1/out/j1-results.csv", RESULTS_CSV.encode())

        first, second = client.get(f"/{CHEMSCRAPER}/results/j1").json()

        assert (first["doc_no"], first["page_no"]) == ("1", "1")
        assert (second["doc_no"], second["page_no"]) == ("1", "2")
        assert (first["PubChemCID"], first["molecularWeight"]) == ("702", "46.07")

    def test_a_missing_results_csv_is_still_a_404(self, client):
        resp = client.get(f"/{CHEMSCRAPER}/results/j1")

        assert resp.status_code == 404


class TestMoleculeCoercion:
    def test_numbers_are_accepted_for_string_fields(self):
        molecule = _molecule(doc_no=1, page_no=2, PubChemCID=702, molecularWeight=46.07)

        assert (molecule.doc_no, molecule.page_no) == ("1", "2")
        assert (molecule.PubChemCID, molecule.molecularWeight) == ("702", "46.07")

    def test_coercion_does_not_extend_to_other_types(self):
        """The model is lenient about numbers, not about types in general."""
        with pytest.raises(ValidationError):
            _molecule(doc_no=["1"])

        with pytest.raises(ValidationError):
            _molecule(SMILE={"smile": "CCO"})


def _molecule(**overrides):
    fields = {
        "id": 0,
        "flagged": False,
        "atom_count": 3,
        "doc_no": "1",
        "file_path": "/inputs/paper.pdf",
        "page_no": "1",
        "SMILE": "CCO",
        "structure": "<svg/>",
        "minX": 1609,
        "minY": 1949,
        "width": 201,
        "height": 108,
        "Location": " | page: 1",
        "OtherInstances": [],
        "fingerprint": "abc",
    }
    fields.update(overrides)
    return Molecule(**fields)
