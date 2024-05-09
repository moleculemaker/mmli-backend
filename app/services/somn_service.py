import os
import time
import asyncio
import json

from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.models import Job
from models.enums import JobStatus

from services.minio_service import MinIOService
from services.email_service import EmailService

class SomnService:
    somn_frontend_baseURL = os.environ.get("SOMN_FRONTEND_URL")

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

    async def runSomn(self, bucket_name: str, payload, service: MinIOService, email_service: EmailService):
        job_id = payload["job_id"]
        db_job = await self.db.get(Job, job_id)

        await self.update_job_phase(db_job, JobStatus.PROCESSING)

        #TODO: use somn image
        await asyncio.sleep(3)

        sample_response = "somn:pong"
        upload_result = service.upload_file(bucket_name, f"results/{job_id}/{job_id}.csv", bytes(sample_response, 'utf-8'))

        if upload_result:
            await self.update_job_phase(db_job, JobStatus.COMPLETED)
            if(db_job.email):
                try:
                    email_service.send_email(
                        db_job.email, 
                        f'''Result for your Somn Job ({db_job.job_id}) is ready''', 
                        f'''The result for your Somn Job is available at {self.somn_frontend_baseURL}/results/{db_job.job_id}''')
                except Exception as e:
                    print(e)
                return True
        else:
            error_content = 'error'
            upload_result = service.upload_file(bucket_name, f"results/{job_id}/error.txt", bytes(error_content, 'utf-8'))
            await self.update_job_phase(db_job, JobStatus.ERROR)
            if db_job.email:
                try:
                    email_service.send_email(
                        db_job.email, 
                        f'''Somn Job ({db_job.job_id}) failed''', 
                        f'''An error occurred in computing the result for your Somn job.''')
                except Exception as e:
                    print(e)
            
        return False
    
    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        job = await db.get(Job, job_id)
        job_info = json.loads(job.job_info)
        file_name = f'{job_info['nuc_name']}_{job_info['el_name']}_processed.csv'
        
        return service.get_file(bucket_name, file_name)
