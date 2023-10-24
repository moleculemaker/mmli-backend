
from datetime import datetime
from fastapi import APIRouter, Depends, status, BackgroundTasks
from fastapi.responses import JSONResponse
import time

from services.minio_service import MinIOService

from services.chemscraper_service import ChemScraperService

from models.analyzeRequestBody import AnalyzeRequestBody
from models.sqlmodel.models import FlaggedMolecule, FlaggedMoleculeDelete

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
        curr_time = time.time()
        db_job = Job(
            email=requestBody.user_email,
            job_info=separator.join(requestBody.fileList),
            job_id=requestBody.jobId,
            # Run ID takes the default value since the app doesn't support multiple runs right now
            type=JobType.CHEMSCRAPER,
            user_agent='',
            time_created=int(curr_time)
        )

        try: 
            db.add(db_job)
            await db.commit()
        except Exception as e:
            content = {"jobId": requestBody.jobId, "error_message": "Database Error Occured.", "error_details": str(e)}
            return JSONResponse(content=content, status_code=400) 

        filename = requestBody.fileList[0]
        chemscraperService = ChemScraperService(db=db)
        objectPath = f"inputs/{requestBody.jobId}/{filename}"
        background_tasks.add_task(chemscraperService.runChemscraperOnDocument, 'chemscraper', filename, objectPath, requestBody.jobId, service)
        content = {"jobId": requestBody.jobId, "submitted_at": datetime.fromtimestamp(curr_time).isoformat()}
        return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED) 

@router.post("/chemscraper/flag", tags=['ChemScraper'])
async def flag_molecule(requestBody: FlaggedMolecule, db: AsyncSession = Depends(get_session)):
    flagged_molecule = FlaggedMolecule(
        smile=requestBody.smile,
        job_id=requestBody.job_id,
        doc_id=requestBody.doc_id,
        time_created=int(time.time())
    )
    try:
        db.add(flagged_molecule)
        await db.commit()
    except Exception as e:
            content = {"jobId": requestBody.job_id, "error_message": "Unable to flag molecule. Database Error Occured.", "error_details": str(e)}
            return JSONResponse(content=content, status_code=400) 
    
    content = {"success_message": "Molecule Flag Successful"}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED) 

@router.delete("/chemscraper/flag", tags=['ChemScraper'])
async def delete_flagged_molecule(requestBody: FlaggedMoleculeDelete, db: AsyncSession = Depends(get_session)):
    flagged_molecule = db.get(FlaggedMolecule, (requestBody.smile, requestBody.job_id))
    try: 
        if flagged_molecule:
            db.delete(flagged_molecule)
            await db.commit()
        else:
            raise Exception("Flagged Molecule with Job ID and SMILE representation Not Found.")
    except Exception as e:
                content = {"jobId": requestBody.job_id, "error_message": "Unable to delete flagged molecule. Database Error Occured.", "error_details": str(e)}
                return JSONResponse(content=content, status_code=400) 
         
    content = {"success_message": "Flagged Molecule Delete Successful"}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED) 
