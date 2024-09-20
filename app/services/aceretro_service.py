from typing import List
from fastapi import HTTPException
import requests
import io
import os
import csv
import zipfile
import time
import pandas as pd
from io import StringIO

from config import get_logger, app_config
from models.sqlmodel.models import Job
from models.enums import JobStatus

from services.minio_service import MinIOService
from services.rdkit_service import RDKitService
from services.pubchem_service import PubChemService
from services.email_service import EmailService
from models.molecule import Molecule
from fastapi import Depends
from fastapi.responses import JSONResponse

from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.models import FlaggedMolecule

log = get_logger(__name__)

class ACERetroService:
    chemscraper_api_baseURL = app_config['external']['chemscraper']['apiBaseUrl']
    chemscraper_frontend_baseURL = app_config['external']['chemscraper']['frontendBaseUrl']

    def __init__(self, db) -> None:
        self.db = db
    
    async def update_job_phase(self, jobObject, phase: JobStatus):
        jobObject.phase = phase
        if phase == JobStatus.PROCESSING:
            jobObject.time_start = int(time.time())
        else:
            jobObject.time_end = int(time.time())
        self.db.add(jobObject)
        await self.db.commit()

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        # rdkitService = RDKitService()
        filepath = f"/{job_id}/out/output.json"
        content = service.get_file(bucket_name, filepath)
        if content is None:
            raise HTTPException(status_code=404, detail=f"File {filepath} not found")
        return content