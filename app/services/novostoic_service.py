import asyncio
import json
import os
import time

from fastapi import HTTPException
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.models import Job
from models.enums import JobStatus

from services.minio_service import MinIOService
from services.email_service import EmailService

from routers.novostoic import validate_chemical

class NovostoicService:
    novostoic_frontend_baseURL = os.environ.get("NOVOSTOIC_FRONTEND_URL")

    def __init__(self, db) -> None:
        self.db = db
    
    @staticmethod
    async def optstoicResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        file = service.get_file(bucket_name, f"{job_id}/out/output.json")
        if not file:
            return None
        data = json.loads(file)
        cache = {}
        for stoic in data:
            for reactant in stoic['stoichiometry']['reactants']:
                if reactant['molecule'] not in cache:
                    cache[reactant['molecule']] = await validate_chemical(reactant['molecule'], db)
                reactant['molecule'] = cache[reactant['molecule']]
            for product in stoic['stoichiometry']['products']:
                if product['molecule'] not in cache:
                    cache[product['molecule']] = await validate_chemical(product['molecule'], db)
                product['molecule'] = cache[product['molecule']]
        return data
    
    @staticmethod
    async def novostoicResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        job = await db.get(Job, job_id)
        if not job:
            return HTTPException(status_code=404, detail="Job not found")
        
        file = service.get_file(bucket_name, f"{job_id}/out/output.json")
        if not file:
            return None
        pathways = json.loads(file)
        cache = {}
        result = {}
        config = json.loads(job.job_info)
        if config['substrate'] not in cache:
            cache[config['substrate']] = await validate_chemical(config['substrate'], db)
        if config['product'] not in cache:
            cache[config['product']] = await validate_chemical(config['product'], db)
        result['primaryPrecursor'] = cache[config['substrate']]
        result['targetMolecule'] = cache[config['product']]
        
        for reactant in config['reactants']:
            if reactant['molecule'] not in cache:
                cache[reactant['molecule']] = await validate_chemical(reactant['molecule'], db)
            reactant['molecule'] = cache[reactant['molecule']]
        
        for product in config['products']:
            if product['molecule'] not in cache:
                cache[product['molecule']] = await validate_chemical(product['molecule'], db)
            product['molecule'] = cache[product['molecule']]
            
        result['stoichiometry'] = {
            'reactants': config['reactants'],
            'products': config['products']
        }
        
        for pathway in pathways:
            if pathway['primaryPrecursor'] not in cache:
                cache[pathway['primaryPrecursor']] = await validate_chemical(pathway['primaryPrecursor'], db)
                
            if pathway['targetMolecule'] not in cache:
                cache[pathway['targetMolecule']] = await validate_chemical(pathway['targetMolecule'], db)
            
            pathway['primaryPrecursor'] = cache[pathway['primaryPrecursor']]
            pathway['targetMolecule'] = cache[pathway['targetMolecule']]
            
            for reactant in pathway['reactants']:
                if reactant['molecule'] not in cache:
                    cache[reactant['molecule']] = await validate_chemical(reactant['molecule'], db)
                reactant['molecule'] = cache[reactant['molecule']]
            
            for product in pathway['products']:
                if product['molecule'] not in cache:
                    cache[product['molecule']] = await validate_chemical(product['molecule'], db)
                product['molecule'] = cache[product['molecule']]
                
        result['pathways'] = pathways
        return result
    
    @staticmethod
    async def enzRankResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        file = service.get_file(bucket_name, f"{job_id}/out/output.json")
        if not file:
            return None
        return json.loads(file)
    
    @staticmethod
    async def dgPredictorResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        file = service.get_file(bucket_name, f"{job_id}/out/output.json")
        if not file:
            return None
        return json.loads(file)
