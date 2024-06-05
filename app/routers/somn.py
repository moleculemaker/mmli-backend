from datetime import datetime
from typing import List

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

from rdkit import Chem

router = APIRouter()

class CheckReactionSiteResponse(BaseModel):
    reaction_site_idxes: List[int]
    svg: str


@router.get(f"/{JobType.SOMN}/all-reaction-sites", tags=['Somn'], response_model=CheckReactionSiteResponse)
async def check_reaction_sites(smiles: str, role: str, hightlight_idxes: List[int] = None):
    try:
        reactionSiteIdxes = SomnService.check_user_input_substrates(smiles, role)
    except SomnException as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=400, detail=str('Invalid user input'))

    mol = Chem.MolFromSmiles(smiles)
    d2d = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(450, 300)

    # Settings for SVG Image coloring via draw options
    dopts = d2d.drawOptions()
    dopts.useBWAtomPalette()
    dopts.continuousHighlight=False
    dopts.highlightBondWidthMultiplier = 24
    dopts.highlightRadius = .4
    dopts.prepareMolsBeforeDrawing=True
    dopts.fillHighlights=False

    Chem.rdDepictor.Compute2DCoords(mol)
    
    reactionSiteIdxes = [int(i) for i in reactionSiteIdxes] if hightlight_idxes is None else hightlight_idxes
    d2d.DrawMolecule(mol, highlightAtoms=reactionSiteIdxes)
    d2d.FinishDrawing()

    return {
        "reaction_site_idxes": reactionSiteIdxes,
        "svg": d2d.GetDrawingText()
    }