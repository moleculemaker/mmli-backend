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