import time
from datetime import datetime
from typing import Optional, Union

from fastapi import APIRouter, Depends, BackgroundTasks, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from sqlalchemy import or_

from services.minio_service import MinIOService
from services.email_service import EmailService
from services.novostoic_service import NovostoicService

from sqlmodel import select, func
from sqlmodel.ext.asyncio.session import AsyncSession

from models.novostoicRequestBodies import DgPredictorRequestBody, NovostoicRequestBody, EnzRankRequestBody, OptstoicRequestBody
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobType, ChemicalIdentifier

from rdkit.Chem import CanonSmiles

class ChemicalAutoCompleteResponse(BaseModel):
    name: str
    smiles: str
    inchi: str
    inchi_key: str
    metanetx_id: str
    kegg_id: str
    structure: Optional[str]

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

@router.get(f"/chemical/auto-complete", 
            tags=['Novostoic'], 
            response_model=list[ChemicalAutoCompleteResponse], 
            description="Returns a list of chemicals that match the search string limited to 20 results.")
async def get_chemical_auto_complete(search: str, db: AsyncSession = Depends(get_session)):
    existing_chemicals = await db.execute(
        select(ChemicalIdentifier.name, 
               ChemicalIdentifier.smiles, 
               ChemicalIdentifier.inchi, 
               ChemicalIdentifier.inchi_key, 
               ChemicalIdentifier.metanetx_id, 
               ChemicalIdentifier.kegg_id
        ).filter(or_(
            ChemicalIdentifier.name.like(f"{search}%"),
            ChemicalIdentifier.smiles.like(f"{search}%"),
            ChemicalIdentifier.inchi.like(f"{search}%"),
            ChemicalIdentifier.inchi_key.like(f"{search}%"),
            ChemicalIdentifier.metanetx_id.like(f"{search}%"),
            ChemicalIdentifier.kegg_id.like(f"{search}%")
        )).limit(20)
    )

    return [
        {
            "name": chemical[0].lower(),
            "smiles": chemical[1],
            "inchi": chemical[2],
            "inchi_key": chemical[3],
            "metanetx_id": chemical[4],
            "kegg_id": chemical[5]
        } for chemical in existing_chemicals.all()
    ]

@router.get(f"/chemical/validate", tags=['Novostoic'], response_model=Union[ChemicalAutoCompleteResponse, None])
async def validate_chemical(search: str, db: AsyncSession = Depends(get_session)):
    try:
        smiles = CanonSmiles(search)
    except:
        smiles = search

    existing_chemicals = (await db.execute(
        select(ChemicalIdentifier.name, 
               ChemicalIdentifier.smiles, 
               ChemicalIdentifier.inchi, 
               ChemicalIdentifier.inchi_key, 
               ChemicalIdentifier.metanetx_id, 
               ChemicalIdentifier.kegg_id
        ).filter(or_(
            ChemicalIdentifier.smiles == smiles,
            ChemicalIdentifier.name == search,
            ChemicalIdentifier.inchi == search,
            ChemicalIdentifier.inchi_key == search,
            ChemicalIdentifier.metanetx_id == search,
            ChemicalIdentifier.kegg_id == search)
        )
    )).all()

    if not len(existing_chemicals):
        return 
    
    chemical = existing_chemicals[0]
    return {
        "name": chemical[0],
        "smiles": chemical[1],
        "inchi": chemical[2],
        "inchi_key": chemical[3],
        "metanetx_id": chemical[4],
        "kegg_id": chemical[5]
    }

@router.post(f"/{JobType.NOVOSTOIC_PATHWAYS}/run", tags=['Novostoic'])
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
        JobType.NOVOSTOIC_PATHWAYS,
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