import json
import time
from datetime import datetime

from fastapi import APIRouter, Depends, BackgroundTasks, status
from fastapi.responses import JSONResponse

from services.minio_service import MinIOService
from services.email_service import EmailService
from services.somn_service import SomnService

from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.somnRequestBody import SomnRequestBody
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobType

from rdkit import Chem

router = APIRouter()

@router.post(f"/{JobType.SOMN}/run", tags=['Somn'])
async def start_somn(
    requestBody: SomnRequestBody, 
    background_tasks: BackgroundTasks, 
    service: MinIOService = Depends(), 
    db: AsyncSession = Depends(get_session), 
    email_service: EmailService = Depends()
):
    #ASSUMPTION: a job of novostoic/optstoic has been created in the db already
    if requestBody.jobId == "":
        content = {"jobId": requestBody.jobId, "error_message": "Job ID not provided."}
        return JSONResponse(content=content, status_code=400)
    
    existing_job = await db.execute(select(Job).where(Job.job_id == requestBody.jobId))
    if not existing_job.first():
        content = {"jobId": requestBody.jobId, "error_message": "Job not found."}
        return JSONResponse(content=content, status_code=404)

    #TODO: add validation for the input
    payload = {
        "job_id": requestBody.jobId,
    }

    somnService = SomnService(db=db)
    background_tasks.add_task(
        somnService.runSomn, 
        JobType.SOMN,
        payload,
        service, 
        email_service
    )
    content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)


@router.get(f"/{JobType.SOMN}/all-reaction-sites", tags=['Somn'])
async def check_reaction_sites(smiles: str, role: str):
    reactionSiteIdxes = SomnService.check_user_input_substrates(smiles, role)

    mol = Chem.MolFromSmiles(smiles)
    d2d = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(300, 150)

    # Settings for SVG Image coloring via draw options
    dopts = d2d.drawOptions()
    dopts.continuousHighlight=False
    dopts.highlightBondWidthMultiplier = 24
    dopts.highlightRadius = .4
    dopts.setHighlightColour((.9,.9,0))
    dopts.prepareMolsBeforeDrawing=True
    dopts.fillHighlights=False

    def get_highlight_atom_color(idxes):
        return { idx: (0, 0, 0, 0) for idx in idxes }

    Chem.rdDepictor.Compute2DCoords(mol)

    d2d.DrawMolecule(mol, 
                     highlightAtoms=reactionSiteIdxes, 
                     highlightAtomColors=get_highlight_atom_color(reactionSiteIdxes))
    d2d.FinishDrawing()

    return {
        "reaction_site_idxes": reactionSiteIdxes,
        "svg": d2d.GetDrawingText()
    }