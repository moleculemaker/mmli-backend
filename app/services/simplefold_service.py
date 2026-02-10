from fastapi import HTTPException

from config import get_logger
from services.minio_service import MinIOService
from sqlmodel.ext.asyncio.session import AsyncSession

log = get_logger(__name__)


class SimpleFoldService:

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """
        Inputs stored in Minio:  /{job_id}/in/<fasta_file>    Bucket name: ml-simplefold
        Outputs stored in Minio: /{job_id}/out/               Bucket name: ml-simplefold
        """
        cif_path = f"{job_id}/out/predictions_simplefold_100M/input_sampled_0.cif"
        output_urls = service.get_file_urls(bucket_name, cif_path)
        if output_urls is None or len(output_urls) == 0:
            raise HTTPException(status_code=404, detail=f"Output CIF file not found for job {job_id}")
        return output_urls
