import json
import os

from fastapi import HTTPException
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
    async def propertyPredictionResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """
        Outputs stored in Minio: /{job_id}/out/*  Bucket name: oed-*
        """
        folder_path = f"{job_id}/out/"
        objects = service.list_files(bucket_name, folder_path, recursive=True)

        # Iterate over folder and add all contents to a dictionary
        content = {}
        for obj in objects:
            file_name = os.path.basename(obj.object_name).split('/')[-1]
            if file_name.endswith('.json'):
                content[file_name] = json.loads(service.get_file(bucket_name=bucket_name, object_name=obj.object_name))
            elif file_name.endswith('.csv'):
                content[file_name] = service.get_file(bucket_name=bucket_name, object_name=obj.object_name)
            else:
                log.warning(f'Skipping unrecognized file extension: ' + str(file_name))

        # Return the dictionary if it has contents
        if not content:
            raise HTTPException(status_code=404, detail=f"No output files were found")

        return content
