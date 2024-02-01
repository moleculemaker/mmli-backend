from pydantic import BaseModel, Field
from typing import List

class Molecule(BaseModel):
    id: int
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

