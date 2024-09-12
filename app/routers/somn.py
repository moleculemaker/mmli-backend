from datetime import datetime
from typing import List, Literal

from fastapi import APIRouter, Depends, BackgroundTasks, HTTPException, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from services.minio_service import MinIOService
from services.email_service import EmailService
from services.somn_service import SomnService, SomnException

from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.somnRequestBody import SomnRequestBody
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobType

from services.shared import draw_chemical_svg

router = APIRouter()

class CheckReactionSiteResponse(BaseModel):
    reaction_site_idxes: List[int]
    svg: str

@router.get(f"/{JobType.SOMN}/all-reaction-sites", tags=['Somn'], response_model=CheckReactionSiteResponse)
async def check_reaction_sites(smiles: str, role: Literal['el', 'nuc']):
    try:
        reactionSiteIdxes = SomnService.check_user_input_substrates(smiles, role)
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

    return {
        "reaction_site_idxes": reactionSiteIdxes,
        "svg": draw_chemical_svg(smiles, 
                             width=450,
                             height=300,
                             beforeDraw=beforeDraw, 
                             highlightAtoms=reactionSiteIdxes)
    }