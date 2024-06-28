from datetime import datetime
import os
from services.userinfo_service import validate_auth_cookie
from fastapi import APIRouter, Depends, status, Request
from fastapi.responses import JSONResponse
import time

from models.sqlmodel.models import SavedMolecule, SavedMoleculeDelete
from sqlmodel.ext.asyncio.session import AsyncSession
from models.sqlmodel.db import get_session

from sqlmodel import select

router = APIRouter()

def is_example_job(job_id):
    return job_id.startswith("example")

@router.get("/molli/saved_molecule", tags=['Molli'])
async def save_molecule(request: Request, job_id: str, db: AsyncSession = Depends(get_session)):
    user = validate_auth_cookie(request)
    email = user["email"] if user is not None else None

    if email == None and is_example_job(job_id):
        return []
    
    try:
        result = await db.execute(select(SavedMolecule).where(
            (SavedMolecule.email == email) & 
            (SavedMolecule.job_id == job_id))
        )
    except Exception as e:
        content = {"jobId": job_id, "message": "Unable to get saved molecules.", "error_details": str(e)}
        return JSONResponse(content=content, status_code=400) 
    save_molecules = result.scalars().all()
    return save_molecules

@router.post("/molli/saved_molecule", tags=['Molli'])
async def save_molecule(request: Request, requestBody: SavedMolecule, db: AsyncSession = Depends(get_session)):
    user = validate_auth_cookie(request)   
    email = user["email"] if user is not None else None
    
    if email == None and is_example_job(requestBody.job_id):
        return JSONResponse(content={"message": "Unauthorized. Please log in to save molecules."}, status_code=401)

    saved_molecule = SavedMolecule(
        email=email,
        job_id=requestBody.job_id,
        molecule_id=requestBody.molecule_id,
        time_created=int(time.time())
    )

    try:
        db.add(saved_molecule)
        await db.commit()
    except Exception as e:
        content = {"jobId": requestBody.job_id, "molecule_id": requestBody.molecule_id, "message": "Unable to save molecule.", "error_details": str(e)}
        return JSONResponse(content=content, status_code=400) 
    
    content = {"message": "Molecule Save Successful."}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED) 

@router.delete("/molli/saved_molecule", tags=['Molli'])
async def delete_saved_molecule(request: Request, requestBody: SavedMoleculeDelete, db: AsyncSession = Depends(get_session)):
    user = validate_auth_cookie(request)
    email = user["email"] if user is not None else None

    results = await db.execute(select(SavedMolecule).where((SavedMolecule.id == requestBody.id)))
    saved_molecule = results.scalars().first()
    
    if saved_molecule:
        if email == saved_molecule.email:
            await db.delete(saved_molecule)
            await db.commit()
        else:
            return JSONResponse(content={"message": "Unauthorized User"}, status_code=401)
            
    content = {"message": "Molecule UnSave Successful."}
    return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)