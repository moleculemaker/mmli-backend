from pydantic import BaseModel, ConfigDict, Field
from typing import List


class Molecule(BaseModel):
    # Molecules are not read back from the database. They are written to a results CSV
    # in MinIO and rebuilt from it on every read (see
    # ChemScraperService.fetchExternalDataAndStoreResults and .resultPostProcess), so
    # each field's type survives only as far as pandas' dtype inference. A doc_no of
    # "1" is written as 1 and comes back as int64; PubChemCID and molecularWeight do
    # the same once PubChem supplies numbers for every row.
    #
    # Pydantic 1 coerced numbers to str, so this went unnoticed for years. Pydantic 2
    # rejects them, which turns GET /chemscraper/results/{job_id} into a 500. Restoring
    # the v1 coercion here keeps the endpoint's output byte-for-byte what it was; the
    # alternative -- declaring dtypes at the pd.read_csv call -- would change what an
    # empty cell deserializes to, and the frontend has always been served the "nan"
    # that the old coercion produced.
    model_config = ConfigDict(coerce_numbers_to_str=True)

    id: int
    flagged: bool
    atom_count: int
    doc_no: str
    file_path: str
    page_no: str
    name: str = Field(default='Unavailable')
    SMILE: str
    structure: str
    minX: float
    minY: float
    width: float
    height: float
    PubChemCID: str = Field(default='Unavailable')
    molecularFormula: str = Field(default='Unavailable')
    molecularWeight: str = Field(default='Unavailable')
    # Not supporting after pubchem batching
    chemicalSafety: List[str] = Field(default=[])
    # Not supporting after pubchem batching
    Description: str = Field(default='Unavailable')
    Location: str
    OtherInstances: List[str]
    fingerprint: str
