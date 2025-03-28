
from datetime import datetime
import os
from fastapi import APIRouter, Depends, status, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse
import time

from sqlmodel import select

from services.minio_service import MinIOService
from services.rdkit_service import RDKitService
from services.chemscraper_service import ChemScraperService
from services.email_service import EmailService

from models.analyzeRequestBody import AnalyzeRequestBody
from models.sqlmodel.models import FlaggedMolecule, FlaggedMoleculeDelete

from sqlmodel.ext.asyncio.session import AsyncSession
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job
from models.enums import JobType, JobStatus
import pandas as pd
import io

router = APIRouter()

@router.post("/chemscraper/analyze", tags=['ChemScraper'])
async def analyze_documents(requestBody: AnalyzeRequestBody, background_tasks: BackgroundTasks, service: MinIOService = Depends(), db: AsyncSession = Depends(get_session), email_service: EmailService = Depends()):
    # Analyze only one document for NSF demo
    if len(requestBody.fileList) > 0 and requestBody.jobId != "":
        # Create a new job and add to the DB
        separator = "|"
        curr_time = time.time()
        db_job: Job = await db.get(Job, requestBody.jobId)
        if not db_job:
            db_job = Job(
                email=requestBody.user_email,
                job_info=separator.join(requestBody.fileList),
                job_id=requestBody.jobId,
                phase=JobStatus.PROCESSING,
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
        objectPath = f"{requestBody.jobId}/in/{filename}"
        background_tasks.add_task(chemscraperService.runChemscraperOnDocument, 'chemscraper', filename, objectPath, requestBody.jobId, service, email_service)
        content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
        return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)

@router.get("/chemscraper/similarity-sorted-order/{job_id}")
def get_similarity_sorted_order(job_id: str, smile_string: str, service: MinIOService = Depends()):
    bucket_name = "chemscraper"
    csv_content = service.get_file(bucket_name, f"{job_id}/out/{job_id}-results.csv")
    if csv_content is None:
        filename = f"{job_id}/out/{job_id}-results.csv"
        raise HTTPException(status_code=404, detail=f"File {filename} not found")
    df = pd.read_csv(io.BytesIO(csv_content))
    # Check if sort_column exists in DataFrame
    sort_column = "fingerprint"
    if sort_column not in df.columns:
        raise HTTPException(status_code=400, detail=f"Column {sort_column} not found in CSV file")
    rdkitService = RDKitService()
    try:
        input_fingerprint =  rdkitService.getFingerprint(smile_string)
    except Exception as e:
        content = {"smile_string": smile_string, "error_message": "Could not generate fingerprint", "error_details": str(e)}
        return JSONResponse(content=content, status_code=400)
    # Calculate Tanimoto similarity for each molecule in the DataFrame
    similarity_scores = []
    for index, row in df.iterrows():
        fingerprint = row['fingerprint']
        if fingerprint == "0" or input_fingerprint == "0":
            similarity = 0.0
        else:
             similarity = rdkitService.getTanimotoSimilarity(input_fingerprint, fingerprint)
        similarity_scores.append(similarity)
    # Add similarity scores as a new column in the DataFrame
    df['similarity'] = similarity_scores
    # Sort DataFrame by sort_column and similarity
    df_sorted = df.sort_values(by='similarity', ascending=False)
    # Return the IDs from each row as a list
    return df_sorted['id'].tolist()

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
    results = await db.execute(select(FlaggedMolecule).where(
        FlaggedMolecule.smile == requestBody.smile
        and FlaggedMolecule.job_id == requestBody.job_id))

    flagged_molecule = results.scalars().first()

    try:
        if flagged_molecule:
            await db.delete(flagged_molecule)
            await db.commit()
        else:
            raise Exception("Flagged Molecule with Job ID and SMILE representation Not Found.")
    except Exception as e:
                content = {"jobId": requestBody.job_id, "error_message": "Unable to delete flagged molecule. Database Error Occured.", "error_details": str(e)}
                return JSONResponse(content=content, status_code=400) 
         
    content = {"success_message": "Flagged Molecule Delete Successful"}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)