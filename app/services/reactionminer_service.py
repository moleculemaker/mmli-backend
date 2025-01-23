import time

from fastapi import HTTPException
from sqlmodel.ext.asyncio.session import AsyncSession

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
        Inputs stored in Minio:  /{job_id}/in/input.json    Bucket name: reactionminer
        Outputs stored in Minio: /{job_id}/out/output.json  Bucket name: reactionminer
        """
        filepath = f"/{job_id}/out/output.json"
        content = service.get_file(bucket_name, filepath)
        if content is None:
            raise HTTPException(status_code=404, detail=f"File {filepath} not found")
        return content