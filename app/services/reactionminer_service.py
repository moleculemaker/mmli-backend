import json
import os
import time

from fastapi import HTTPException
from sqlmodel.ext.asyncio.session import AsyncSession

from config import log
from models.enums import JobStatus
from services.minio_service import MinIOService


class ReactionMinerService:
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
        """
        Inputs stored in Minio:  /{job_id}/in/[name].pdf    Bucket name: reactionminer
        Outputs stored in Minio: /{job_id}/out/[name].json  Bucket name: reactionminer
        """
        folder_path = f"/{job_id}/out/"
        objects = service.list_files(bucket_name, folder_path)

        # Iterate over folder and add all contents to a dictionary
        content = {}
        for obj in objects:
            file_name = os.path.basename(obj.object_name).split('/')[-1]
            content[file_name] = json.loads(service.get_file(bucket_name=bucket_name, object_name=obj.object_name))

        # Return the dictionary if it has contents
        if not content:
            raise HTTPException(status_code=404, detail=f"No output files were found")

        return content
