import asyncio
import json
import os
import time

from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.models import Job
from models.enums import JobStatus

from services.minio_service import MinIOService
from services.email_service import EmailService

class NovostoicService:
    novostoic_frontend_baseURL = os.environ.get("NOVOSTOIC_FRONTEND_URL")

    def __init__(self, db) -> None:
        self.db = db
    
    @staticmethod
    async def optstoicResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        return service.get_file(bucket_name, f"results/{job_id}/{job_id}.csv")
    
    @staticmethod
    async def novostoicResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        return await NovostoicService.optstoicResultPostProcess(bucket_name, job_id, service, db)
    
    @staticmethod
    async def enzRankResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        return await NovostoicService.optstoicResultPostProcess(bucket_name, job_id, service, db)
    
    @staticmethod
    async def dgPredictorResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        file = service.get_file(bucket_name, f"{job_id}/out/output.json")
        if not file:
            return None
        return json.loads(file)
