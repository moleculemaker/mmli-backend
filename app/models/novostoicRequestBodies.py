from typing import Dict, List, Union
from pydantic import BaseModel

class OptstoicRequestBody(BaseModel):
    jobId: str
    user_email: str
    primary_precursor: str
    target_molecule: str


class NovostoicRequestBody(BaseModel):
    jobId: str
    user_email: str
    
    stoic: str
    substrate: str
    product: str
    max_steps: int
    iterations: int

    # BG: comment out as those are not used
    # is_thermodynamic_feasible: bool
    # thermodynamical_feasible_reaction_only: bool
    # use_enzyme_selection: bool
    # num_enzyme_candidates: int


class DgPredictorRequestBody(BaseModel):
    class ReactionSmiles(BaseModel):
        type: str
        smiles: str

    class ReactionKeggId(BaseModel):
        type: str
        molecule_number: str
        molecule_inchi_or_smiles: str
        reaction_keggid: str
    
    jobId: str
    user_email: str

    ph: float
    ionic_strength: float
    reactions: List[Union[ReactionSmiles, ReactionKeggId]]


class EnzRankRequestBody(BaseModel):
    jobId: str
    user_email: str

    primary_precursor: str
    enzyme_sequence: str