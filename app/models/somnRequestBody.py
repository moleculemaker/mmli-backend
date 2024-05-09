from typing import List
from pydantic import BaseModel

class SomnRequestBody(BaseModel):
    jobId: str
    user_email: str
    reactant_pair_name: str
    el: str
    el_name: str
    nuc: str
    nuc_name: str