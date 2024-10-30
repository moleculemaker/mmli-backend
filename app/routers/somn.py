from typing import List, Literal

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from openbabel import pybel as pb

from services.somn_service import SomnService, SomnException

from models.sqlmodel.models import JobType

from services.shared import draw_chemical_svg

router = APIRouter()  

class CheckReactionSiteResponse(BaseModel):
    reaction_site_idxes: List[int]
    svg: str

@router.get(f"/{JobType.SOMN}/all-reaction-sites", tags=['Somn'], response_model=CheckReactionSiteResponse)
async def check_reaction_sites(
    input: str, 
    role: Literal['el', 'nuc'],
    input_type: Literal['smi', 'cml', 'cdxml']
):
    try:
        reaction_sites_idxes = SomnService.check_user_input_substrates(input, input_type, role)
    except SomnException as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
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
        input = mol.write('smi')

    return {
        "reaction_site_idxes": reaction_sites_idxes,
        "smiles": input,
        "svg": draw_chemical_svg(input, 
                                 width=450,
                                 height=300,
                                 beforeDraw=beforeDraw, 
                                 highlightAtoms=reaction_sites_idxes)
    }