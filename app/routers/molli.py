from datetime import datetime
import os
from fastapi import APIRouter, Depends, status
from fastapi.responses import JSONResponse
import time

from models.sqlmodel.models import SavedMolecule, SavedMoleculeDelete
from sqlmodel.ext.asyncio.session import AsyncSession
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job
from models.enums import JobType, JobStatus

from sqlmodel import select

router = APIRouter()


@router.get("/molli/saved_molecules", tags=['Molli'])
async def save_molecule(email: str, job_id: str, db: AsyncSession = Depends(get_session)):
    try:
        result = await db.execute(select(SavedMolecule).where(SavedMolecule.email == email).where(SavedMolecule.job_id == job_id))
    except Exception as e:
        content = {"jobId": job_id, "message": "Unable to get saved molecules. Database Error Occured.", "error_details": str(e)}
        return JSONResponse(content=content, status_code=400) 
    save_molecules = result.scalars().all()
    return save_molecules

@router.post("/molli/saved_molecules", tags=['Molli'])
async def save_molecule(requestBody: SavedMolecule, db: AsyncSession = Depends(get_session)):
    saved_molecule = SavedMolecule(
        email=requestBody.email,
        job_id=requestBody.job_id,
        molecule_id=requestBody.molecule_id,
        time_created=int(time.time())
    )
    try:
        db.add(saved_molecule)
        await db.commit()
    except Exception as e:
        content = {"jobId": requestBody.job_id, "molecule_id": requestBody.molecule_id, "message": "Unable to save molecule. Database Error Occured.", "error_details": str(e)}
        return JSONResponse(content=content, status_code=400) 
    
    content = {"message": "Molecule Save Successful."}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED) 

@router.delete("/molli/saved_molecules", tags=['Molli'])
async def delete_saved_molecule(requestBody: SavedMoleculeDelete, db: AsyncSession = Depends(get_session)):
    results = await db.execute(select(SavedMolecule).where(
        SavedMolecule.email == requestBody.email
        and SavedMolecule.job_id == requestBody.job_id
        and SavedMolecule.molecule_id == requestBody.molecule_id))

    saved_molecule = results.scalars().first()

    try:
        if saved_molecule:
            await db.delete(saved_molecule)
            await db.commit()
        else:
            raise Exception("Saved Molecule with Job ID and Molecule ID representation Not Found.")
    except Exception as e:
            content = {"jobId": requestBody.job_id, "molecule_id": requestBody.molecule_id, "message": "Unable to unsave molecule. Database Error Occured.", "error_details": str(e)}
            return JSONResponse(content=content, status_code=400) 
         
    content = {"message": "Molecule UnSave Successful."}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)