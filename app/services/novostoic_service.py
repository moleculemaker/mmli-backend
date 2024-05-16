import asyncio
import os
import time
import json
import re

from fastapi import HTTPException

from sqlmodel.ext.asyncio.session import AsyncSession
from sqlmodel import select
from sqlalchemy import or_

from models.sqlmodel.models import Job, ChemicalIdentifier
from models.enums import JobStatus

from services.minio_service import MinIOService
from services.email_service import EmailService

from rdkit.Chem import CanonSmiles

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

        sample_response = f"{bucket_name}:pong"
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
        return await self.runOptstoic(bucket_name, payload, service, email_service) 
    
    async def runEnzRank(self, bucket_name: str, payload, service: MinIOService, email_service: EmailService):
        #TODO: use enzRank image
        return await self.runOptstoic(bucket_name, payload, service, email_service) 
    
    async def runDgPredictor(self, bucket_name: str, payload, service: MinIOService, email_service: EmailService):
        #TODO: use dGPredictor image
        return await self.runOptstoic(bucket_name, payload, service, email_service) 
    
    @staticmethod
    async def getChemicalByMetanetXId(metanetx_id: str, db: AsyncSession):
        chemical = (await db.execute(select(ChemicalIdentifier).where(ChemicalIdentifier.metanetx_id == metanetx_id))).all()
        return chemical[0][0] if chemical else None
    
    @staticmethod
    async def optstoicResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        # get job info
        job = await db.get(Job, job_id)
        if not job:
            raise HTTPException(status_code=404, detail="Job not found")
        
        job_info = json.loads(job.job_info)
        primary_precursor = await NovostoicService.getChemicalByMetanetXId(job_info['primary_precursor'], db)
        target_molecule = await NovostoicService.getChemicalByMetanetXId(job_info['target_molecule'], db)
        
        # get all files of the job
        start_idx = 1
        file_name = f"results/{job_id}/metanetx_id/soln_dict_metanetx_id{start_idx}.json"
        files = []
        while service.check_file_exists(bucket_name, file_name):
            file = service.get_file(bucket_name, file_name)
            files.append(json.loads(file))
            start_idx += 1
            file_name = f"results/{job_id}/metanetx_id/soln_dict_metanetx_id{start_idx}.json"
            
        reactions = []
        for file in files:
            reactants = []
            products = []
            for key in file:
                if key == 'dG_Range' or \
                    key == target_molecule.metanetx_id or \
                    key == primary_precursor.metanetx_id:
                    continue
                
                y = int(file[key])
                target = reactants if y < 0 else products
                chemical = await NovostoicService.getChemicalByMetanetXId(key, db)
                target.append({
                    "molecule": chemical if chemical else key,
                    "amount": abs(y)
                })
            deltaG = re.search(r"dG range = (-?\d+)", file['dG_Range'])
            reactions.append({
                "stoichiometry": {
                    "reactants": reactants,
                    "products": products,
                },
                "yield": file[target_molecule.metanetx_id],
                "deltaG": deltaG.group(1) if deltaG else None,
            })
            
        return {
            "primaryPrecursor": primary_precursor,
            "targetMolecule": target_molecule,
            "results": reactions
        }
    
    @staticmethod
    async def novostoicResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        return await NovostoicService.optstoicResultPostProcess(bucket_name, job_id, service, db)
    
    @staticmethod
    async def enzRankResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        return await NovostoicService.optstoicResultPostProcess(bucket_name, job_id, service, db)
    
    @staticmethod
    async def dgPredictorResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        return await NovostoicService.optstoicResultPostProcess(bucket_name, job_id, service, db)
