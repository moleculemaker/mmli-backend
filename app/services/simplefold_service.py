from fastapi import HTTPException

from config import get_logger
from services.minio_service import MinIOService
from sqlmodel.ext.asyncio.session import AsyncSession

log = get_logger(__name__)


class SimpleFoldService:

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """
        Inputs stored in Minio:  /{job_id}/in/<fasta_file>    Bucket name: mlsimplefold
        Outputs stored in Minio: /{job_id}/out/               Bucket name: mlsimplefold
        """
        output_urls = service.get_file_urls(bucket_name, f"{job_id}/out/")
        if output_urls is None:
            raise HTTPException(status_code=404, detail=f"Output files not found for job {job_id}")
        return output_urls
