from typing import List, Literal
from traceback import format_exc

from fastapi import APIRouter, HTTPException, logger
from pydantic import BaseModel

from services.somn_service import SomnService, SomnException

from models.sqlmodel.models import JobType

from services.shared import draw_chemical_svg

router = APIRouter()

class CheckReactionSiteResponse(BaseModel):
    reaction_site_idxes: List[int]
    has_chiral: bool
    svg: str

@router.get(f"/{JobType.SOMN}/all-reaction-sites", tags=['Somn'], response_model=CheckReactionSiteResponse)
async def check_reaction_sites(smiles: str, role: Literal['el', 'nuc']):
    try:
        reactionSiteIdxes, _, has_chiral = SomnService.check_user_input_substrates(smiles, role)
    except SomnException as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        stack_trace = format_exc()
        logger.logger.error(f"Error when checking reaction sites for {smiles}: {stack_trace}")
        raise HTTPException(status_code=400, detail=str('Invalid user input'))
    
    def beforeDraw(d2d):
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.continuousHighlight=False
        dopts.highlightBondWidthMultiplier = 24
        dopts.highlightRadius = .4
        dopts.prepareMolsBeforeDrawing=True
        dopts.fillHighlights=False

    return {
        "reaction_site_idxes": reactionSiteIdxes,
        "has_chiral": has_chiral,
        "svg": draw_chemical_svg(smiles, 
                             width=450,
                             height=300,
                             beforeDraw=beforeDraw, 
                             highlightAtoms=reactionSiteIdxes)
    }