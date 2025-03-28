from typing import List, Optional
from pydantic import BaseModel

class ExportRequestBody(BaseModel):
    jobId: str
    input_filename: str
    cdxml: bool
    cdxml_filter: str
    cdxml_selected_pages: List[int]
    csv: bool
    csv_filter: str
    csv_molecules: List[int]