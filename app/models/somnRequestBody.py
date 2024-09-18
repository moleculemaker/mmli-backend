from typing import Optional
from pydantic import BaseModel

class SomnRequestBody(BaseModel):
    jobId: str
    user_email: str
    reactant_pair_name: str
    el: str
    el_name: str
    nuc: str
    nuc_name: str
    nuc_idx: str
    el_idx: str
