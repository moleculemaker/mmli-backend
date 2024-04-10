import json
import time
from datetime import datetime

from fastapi import APIRouter, Depends, BackgroundTasks, status
from fastapi.responses import JSONResponse

from services.minio_service import MinIOService
from services.email_service import EmailService
from services.novostoic_service import NovostoicService

from sqlmodel.ext.asyncio.session import AsyncSession
from models.somnRequestBody import SomnRequestBody
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobType


router = APIRouter()

@router.post("/somn/run", tags=['Somn'])
async def start_somn(
    requestBody: SomnRequestBody, 
    background_tasks: BackgroundTasks, 
    service: MinIOService = Depends(), 
    db: AsyncSession = Depends(get_session), 
    email_service: EmailService = Depends()
):
    #TODO: add validation for the input
    if requestBody.jobId != "":
        # Create a new job and add to the DB
        curr_time = time.time()
        db_job = Job(
            email=requestBody.user_email,
            job_info=json.dumps(requestBody), #TODO: use appropriate data structure
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
    
        somnService = somnService(db=db)
        background_tasks.add_task(
            somnService.runSomn, 
            'somn',
            {}, #TODO: use appropriate data structure
            service, 
            email_service
        )
        content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
        return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)