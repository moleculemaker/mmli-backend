import time
from datetime import datetime

from fastapi import APIRouter, Depends, BackgroundTasks, status
from fastapi.responses import JSONResponse
from sqlalchemy import or_

from services.minio_service import MinIOService
from services.email_service import EmailService
from services.novostoic_service import NovostoicService

from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.novostoicRequestBodies import DgPredictorRequestBody, NovostoicRequestBody, EnzRankRequestBody, OptstoicRequestBody
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobType, ChemicalIdentifier


router = APIRouter()

@router.post(f"/{JobType.NOVOSTOIC_OPTSTOIC}/run", tags=['Novostoic'])
async def start_optstoic(
    requestBody: OptstoicRequestBody, 
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

    novostoicService = NovostoicService(db=db)
    background_tasks.add_task(
        novostoicService.runOptstoic, 
        JobType.NOVOSTOIC_OPTSTOIC,
        payload,
        service, 
        email_service
    )
    content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)

@router.get(f"/chemical/auto-complete", tags=['Novostoic'])
async def get_chemical_auto_complete(search: str, db: AsyncSession = Depends(get_session)):
    existing_chemicals = await db.execute(
        select(ChemicalIdentifier)
        .filter(or_(
            ChemicalIdentifier.name.like(f"%{search}%"),
            ChemicalIdentifier.smiles.like(f"%{search}%"),
            ChemicalIdentifier.inchi.like(f"%{search}%"),
            ChemicalIdentifier.inchi_key.like(f"%{search}%"),
            ChemicalIdentifier.metanetx_id.like(f"%{search}%"),
            ChemicalIdentifier.kegg_id.like(f"%{search}%")
        )
        ).limit(20)
    )
    return existing_chemicals.all()

@router.post(f"/{JobType.NOVOSTOIC_NOVOSTOIC}/run", tags=['Novostoic'])
async def start_novostoic(
    requestBody: NovostoicRequestBody, 
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

    novostoicService = NovostoicService(db=db)
    background_tasks.add_task(
        novostoicService.runNovostoic, 
        JobType.NOVOSTOIC_NOVOSTOIC,
        payload,
        service, 
        email_service
    )
    content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)


@router.post(f"/{JobType.NOVOSTOIC_ENZRANK}/run", tags=['Novostoic'])
async def start_enzrank(
    requestBody: EnzRankRequestBody, 
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

    novostoicService = NovostoicService(db=db)
    background_tasks.add_task(
        novostoicService.runEnzRank, 
        JobType.NOVOSTOIC_ENZRANK,
        payload,
        service, 
        email_service
    )
    content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)


@router.post(f"/{JobType.NOVOSTOIC_DGPREDICTOR}/run", tags=['Novostoic'])
async def start_dgpredictor(
    requestBody: DgPredictorRequestBody, 
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

    novostoicService = NovostoicService(db=db)
    background_tasks.add_task(
        novostoicService.runDgPredictor, 
        JobType.NOVOSTOIC_DGPREDICTOR,
        payload,
        service, 
        email_service
    )
    content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)