from typing import List
from pydantic import BaseModel

class ExportRequestBody(BaseModel):
    jobId: str
    pdf_filename: str
    cdxml: bool
    cdxml_filter: str
    cdxml_selected_pages: List[int]
    csv: bool
    csv_filter: str
    csv_molecules: List[int]