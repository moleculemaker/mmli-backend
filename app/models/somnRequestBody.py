from typing import List
from pydantic import BaseModel

class SomnRequestBody(BaseModel):
    jobId: str
    user_email: str
    reactant_pair_name: str
    aryl_halide_name: str
    aryl_halide_smiles: str
    amine_name: str
    amine_smiles: str