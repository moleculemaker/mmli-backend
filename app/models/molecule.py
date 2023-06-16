from pydantic import BaseModel
from typing import List

class Molecule(BaseModel):
    doc_no: str
    file_path: str
    page_no: str
    name: str
    SMILE: str
    structure: str
    minX: float
    minY: float
    width: float
    height: float
    PubChemCID: str
    molecularFormula: str
    molecularWeight: str
    chemicalSafety: List[str]
    Description: str

