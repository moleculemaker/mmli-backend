import json

from fastapi import HTTPException
from fastapi.responses import JSONResponse

from config import get_logger
from services.minio_service import MinIOService
from sqlmodel.ext.asyncio.session import AsyncSession

log = get_logger(__name__)


class EzSpecificityService:

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """
        Post-process ez-specificity docking results.

        Inputs stored in MinIO:  /{job_id}/in/   (job_config.json + PDB files)
        Outputs stored in MinIO: /{job_id}/out/  (results.json + complexes/*.pdb)
        """
        results_path = f"{job_id}/out/results.json"
        content = service.get_file(bucket_name, results_path)
        if content is None:
            raise HTTPException(status_code=404, detail=f"Results file not found for job {job_id}")

        try:
            results = json.loads(content.decode('utf-8'))
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            raise HTTPException(status_code=500, detail=f"Failed to parse results for job {job_id}: {str(e)}")

        # Convert complex filenames to presigned MinIO URLs
        complex_urls = service.get_file_urls(bucket_name, f"{job_id}/out/complexes/")
        url_map = {}
        if complex_urls:
            for url in complex_urls:
                # Extract filename from URL path
                filename = url.split('/')[-1].split('?')[0]
                url_map[filename] = url

        for result in results:
            complex_filename = result.get('complex', '')
            if complex_filename and complex_filename in url_map:
                result['complex'] = url_map[complex_filename]

        return results
