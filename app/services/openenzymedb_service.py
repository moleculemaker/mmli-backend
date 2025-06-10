import json
import os
import time

from fastapi import HTTPException
from sqlmodel.ext.asyncio.session import AsyncSession

from config import log
from models.enums import JobStatus
from services.minio_service import MinIOService
from services.reactionminer_search.core import txt_eval, smi_eval, mm_eval


class OpenEnzymeDBService:
    def __init__(self, db) -> None:
        self.db = db

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """
        Outputs stored in Minio: /{job_id}/out/*  Bucket name: oed-*
        """
        folder_path = f"/{job_id}/out/"
        objects = service.list_files(bucket_name, folder_path)

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

