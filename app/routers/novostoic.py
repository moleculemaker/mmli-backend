import time
from datetime import datetime

from fastapi import APIRouter, Depends, BackgroundTasks, status
from fastapi.responses import JSONResponse

from services.minio_service import MinIOService
from services.email_service import EmailService
from services.novostoic_service import NovostoicService

from sqlmodel.ext.asyncio.session import AsyncSession
from models.novostoicRequestBodies import OptstoicRequestBody
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobType


router = APIRouter()

@router.post("/novostoic/optstoic/run", tags=['Novostoic'])
async def start_optstoic(
    requestBody: OptstoicRequestBody, 
    background_tasks: BackgroundTasks, 
    service: MinIOService = Depends(), 
    db: AsyncSession = Depends(get_session), 
    email_service: EmailService = Depends()
):
    if len(requestBody.primary_precursor) > 0 \
        and len(requestBody.target_molecule) > 0 \
        and requestBody.jobId != "":
        # Create a new job and add to the DB
        curr_time = time.time()
        db_job = Job(
            email=requestBody.user_email,
            job_info=requestBody.primary_precursor + "|" + requestBody.target_molecule,
            job_id=requestBody.jobId,
            # Run ID takes the default value since the app doesn't support multiple runs right now
            type=JobType.NOVOSTOIC_OPTSTOIC,
            user_agent='',
            time_created=int(curr_time)
        )

        try:
            db.add(db_job)
            await db.commit()
        except Exception as e:
            content = {"jobId": requestBody.jobId, "error_message": "Database Error Occured.", "error_details": str(e)}
            return JSONResponse(content=content, status_code=400)
    
        novostoicService = NovostoicService(db=db)
        background_tasks.add_task(
            novostoicService.runOptstoic, 
            'novostoic',
            requestBody.primary_precursor, 
            requestBody.target_molecule, 
            requestBody.jobId, 
            service, 
            email_service
        )
        content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
        return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)