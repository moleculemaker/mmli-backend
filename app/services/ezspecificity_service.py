import csv
import io
import json


from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.enums import JobType
from models.sqlmodel.models import Job

from config import get_logger
from services.minio_service import MinIOService

log = get_logger(__name__)

class EzSpecificityService:

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        # TODO: what format is expected for the results in the frontend?
        # JSON array of objects w/ fields: rank, enzyme, substrate, ez_score
        #   Alt (if inputMode != pairs): rank, complexLabel, ez_score

        # TODO: fetch MinIO files and present in the correct format
        job_status = await db.get(Job, job_id)
        job_info = job_status.job_info
        job_config = json.loads(job_info.replace('\\"', '"'))

        unidock_id = job_config['ezspec_unidock_job_id']
        inference_id = job_config['ezspec_inference_job_id']

        unidock_results = await EzSpecificityService.unidockResultPostProcess(
            bucket_name=JobType.EZSPEC_UNIDOCK,
            job_id=unidock_id,
            service=service,
            db=db
        )
        inference_results = await EzSpecificityService.inferenceResultPostProcess(
            bucket_name=JobType.EZSPEC_INFERENCE,
            job_id=inference_id,
            service=service,
            db=db
        )

        return unidock_results | inference_results | { 'job_config': job_config }

    @staticmethod
    async def unidockResultPostProcess(bucket_name, job_id, service, db):
        # TODO: what format is expected for the results in the frontend?

        # TODO: fetch MinIO files and present in the correct format
        job_status = await db.get(Job, job_id)
        job_info = job_status.job_info
        job_config = json.loads(job_info.replace('\\"', '"'))

        unidock_id = job_config['ezspec_unidock_job_id']

        # Fetch ezspec-unidock results
        enzyme_bytes = service.get_file(bucket_name, f"{job_id}/out/{unidock_id}/out/Enzymes.csv")
        enzymes_stream = io.StringIO(enzyme_bytes.decode('utf-8'))
        substrate_bytes = service.get_file(bucket_name, f"{job_id}/out/{unidock_id}/out/Substrates.csv")
        substrates_stream = io.StringIO(substrate_bytes.decode('utf-8'))
        data_bytes = service.get_file(bucket_name, f"{job_id}/out/{unidock_id}/out/data.csv")
        data_stream = io.StringIO(data_bytes.decode('utf-8'))

        return {
            'job_config': job_config,
            'enzymes': list(csv.reader(enzymes_stream)),
            'substrates': list(csv.reader(substrates_stream)),
            'data': list(csv.reader(data_stream)),
        }

    @staticmethod
    async def inferenceResultPostProcess(bucket_name, job_id, service, db):
        # TODO: what format is expected for the results in the frontend?

        # TODO: fetch MinIO files and present in the correct format
        job_status = await db.get(Job, job_id)
        job_info = job_status.job_info
        job_config = json.loads(job_info.replace('\\"', '"'))

        inference_id = job_config['ezspec_inference_job_id']

        # Fetch ezspec-inference results
        results_bytes = service.get_file(bucket_name, f"{job_id}/out/{inference_id}/out/results.csv")
        results_stream = io.StringIO(results_bytes.decode('utf-8'))

        return {
            'job_config': job_config,
            'results': list(csv.reader(results_stream)),
        }
