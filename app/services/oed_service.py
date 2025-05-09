import json

from fastapi import HTTPException
from models.sqlmodel.models import Job
from services.minio_service import MinIOService
from sqlmodel.ext.asyncio.session import AsyncSession
from services.shared import get_iupac_name

class OEDService:
    
    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
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