import json
import os

from fastapi import HTTPException
from models.enums import JobType
from models.sqlmodel.models import Job
from services.minio_service import MinIOService
from sqlmodel.ext.asyncio.session import AsyncSession
from services.shared import get_iupac_name
from config import log

class OEDService:
    def __init__(self, db) -> None:
        self.db = db
    
    @staticmethod
    async def chemInfoResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        file = service.get_file(bucket_name, f"{job_id}/out/job.json")
        if not file:
            return None
        data = json.loads(file)

        job = await db.get(Job, job_id)
        if not job:
            return HTTPException(status_code=404, detail="Job not found")
        
        job_info = json.loads(job.job_info)
        data['query_smiles'] = {
            "iupac_name": get_iupac_name(job_info['query_smiles']),
            "smiles": job_info['query_smiles'],
            "matches_return": {
                "mcs": len(data['mcs']),
                "fragment": len(data['fragment']),
                "tanimoto": len(data['tanimoto']),
            }
        }

        return data

    @staticmethod
    async def propertyPredictionResultPostProcess(bucket_name: str, job_id: str, service: MinIOService):
        """
        Outputs stored in Minio: /{job_id}/out/*  Bucket name: oed-*
        """
        
        filename_map = {
            JobType.OED_DLKCAT: "dlkcat-output.json",
            JobType.OED_UNIKP: "unikp-output.json",
            JobType.OED_CATPRED: "catpred-output.json",
        }
        
        content = service.get_file(bucket_name, f"{job_id}/out/{filename_map[bucket_name]}")
        if not content:
            return HTTPException(status_code=404, detail=f"File {filename_map[bucket_name]} not found")
        
        return json.loads(content)
        
