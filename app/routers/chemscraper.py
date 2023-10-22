
from datetime import datetime
from fastapi import APIRouter, Depends, status, BackgroundTasks
from fastapi.responses import JSONResponse
import time

from services.minio_service import MinIOService

from services.chemscraper_service import ChemScraperService

from models.analyzeRequestBody import AnalyzeRequestBody

from sqlmodel.ext.asyncio.session import AsyncSession
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job
from models.enums import JobType

router = APIRouter()


@router.post("/chemscraper/analyze", tags=['ChemScraper'])
async def analyze_documents(requestBody: AnalyzeRequestBody, background_tasks: BackgroundTasks, service: MinIOService = Depends(), db: AsyncSession = Depends(get_session)):
    # Analyze only one document for NSF demo
    if len(requestBody.fileList) > 0 and requestBody.jobId != "":
        # Create a new job and add to the DB
        separator = "|"
        db_job = Job(
            email=requestBody.user_email,
            job_info=separator.join(requestBody.fileList),
            job_id=requestBody.jobId,
            # Run ID takes the default value since the app doesn't support multiple runs right now
            type=JobType.CHEMSCRAPER,
            user_agent='',
            time_created=int(time.time())
        )

        try: 
            db.add(db_job)
            await db.commit()
        except:
            content = {"jobId": requestBody.jobId, "error_message": "Job ID already exists with the same ID. Please retry using another job ID."}
            return JSONResponse(content=content, status_code=400) 

        filename = requestBody.fileList[0]
        chemscraperService = ChemScraperService(db=db)
        objectPath = f"inputs/{requestBody.jobId}/{filename}"
        background_tasks.add_task(chemscraperService.runChemscraperOnDocument, 'chemscraper', filename, objectPath, requestBody.jobId, service)
        content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
        return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED) 


