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

from services.shared import draw_chemical_svg

router = APIRouter()

class CheckReactionSiteResponse(BaseModel):
    reaction_site_idxes: List[int]
    svg: str

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


@router.get(f"/{JobType.SOMN}/all-reaction-sites", tags=['Somn'], response_model=CheckReactionSiteResponse)
async def check_reaction_sites(smiles: str, role: str):
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