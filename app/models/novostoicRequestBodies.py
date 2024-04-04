from typing import List
from pydantic import BaseModel

class OptstoicRequestBody(BaseModel):
    jobId: str
    user_email: str
    primary_precursor: str
    target_molecule: str