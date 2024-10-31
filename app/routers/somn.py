import re
from typing import List, Literal
from traceback import format_exc

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from openbabel import pybel as pb

from config import get_logger
from services.somn_service import SomnService, SomnException

from models.sqlmodel.models import JobType

from services.shared import draw_chemical_svg

router = APIRouter()  

log = get_logger(__name__)

class CheckReactionSiteRequest(BaseModel):
    input: str
    role: Literal['el', 'nuc']
    input_type: Literal['smi', 'cml', 'cdxml']

class CheckReactionSiteResponse(BaseModel):
    reaction_site_idxes: List[int]
    smiles: str
    has_chiral: bool
    num_heavy_atoms: int
    svg: str

@router.post(f"/{JobType.SOMN}/all-reaction-sites", tags=['Somn'], response_model=CheckReactionSiteResponse)
async def check_reaction_sites(
    request: CheckReactionSiteRequest
):
    input = request.input.replace('\"', '"')
    role = request.role
    input_type = request.input_type
    
    input = input.strip()
    if input_type == 'cml':
        input = re.sub(r'> +<', '><', input)
    elif input_type == 'cdxml':
        input = re.sub(r' +', ' ', input)
    
    try:
        reactionSiteIdxes, _, has_chiral, num_heavy_atoms = SomnService.check_user_input_substrates(input, input_type, role)
    except SomnException as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        stack_trace = format_exc()
        log.error(f"Error when checking reaction sites for {input}: {stack_trace}")
        raise HTTPException(status_code=400, detail=str('Invalid user input'))
    
    def beforeDraw(d2d):
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.continuousHighlight=False
        dopts.highlightBondWidthMultiplier = 24
        dopts.highlightRadius = .4
        dopts.prepareMolsBeforeDrawing=True
        dopts.fillHighlights=False
        
    if input_type == 'cdxml' or input_type == 'cml':
        mol = pb.readstring(input_type, input)
        input = mol.write('smi').strip()

    return {
        "reaction_site_idxes": reactionSiteIdxes,
        "has_chiral": has_chiral,
        "num_heavy_atoms": num_heavy_atoms,
        "smiles": input,
        "svg": draw_chemical_svg(input, 
                             width=450,
                             height=300,
                             beforeDraw=beforeDraw, 
                             highlightAtoms=reactionSiteIdxes)
    }