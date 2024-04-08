import asyncio
import os
import time
from models.sqlmodel.models import Job
from models.enums import JobStatus

from services.minio_service import MinIOService
from services.email_service import EmailService

class NovostoicService:
    novostoic_frontend_baseURL = os.environ.get("NOVOSTOIC_FRONTEND_URL")

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

    async def runOptstoic(self, bucket_name: str, payload, service: MinIOService, email_service: EmailService):
        job_id = payload["job_id"]
        db_job = await self.db.get(Job, job_id)

        await self.update_job_phase(db_job, JobStatus.PROCESSING)

        #TODO: use optstoic image
        await asyncio.sleep(3)

        sample_response = "message\npong"
        upload_result = service.upload_file(bucket_name, f"results/{job_id}/{job_id}.csv", bytes(sample_response, 'utf-8'))

        if upload_result:
            await self.update_job_phase(db_job, JobStatus.COMPLETED)
            if(db_job.email):
                try:
                    email_service.send_email(
                        db_job.email, 
                        f'''Result for your Optstoic Job ({db_job.job_id}) is ready''', 
                        f'''The result for your ChemScraper Job is available at {self.novostoic_frontend_baseURL}/results/{db_job.job_id}''')
                except Exception as e:
                    print(e)
                return True
        else:
            error_content = 'error'
            upload_result = service.upload_file(bucket_name, f"results/{job_id}/error.txt", error_content)
            await self.update_job_phase(db_job, JobStatus.ERROR)
            if db_job.email:
                try:
                    email_service.send_email(
                        db_job.email, 
                        f'''Novostoic/Optstoic Job ({db_job.job_id}) failed''', 
                        f'''An error occurred in computing the result for your Novostoic/Optstoic job.''')
                except Exception as e:
                    print(e)
            
        return False

    async def runNovostoic(self, bucket_name: str, payload, service: MinIOService, email_service: EmailService):
        #TODO: use novostoic image
        return self.runOptstoic(bucket_name, payload, service, email_service) 
    
    async def runEnzRank(self, bucket_name: str, payload, service: MinIOService, email_service: EmailService):
        #TODO: use enzRank image
        return self.runOptstoic(bucket_name, payload, service, email_service) 
    
    async def runDgPredictor(self, bucket_name: str, payload, service: MinIOService, email_service: EmailService):
        #TODO: use dGPredictor image
        return self.runOptstoic(bucket_name, payload, service, email_service) 
