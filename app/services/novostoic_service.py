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
from services.shared import draw_chemical_svg

from routers.novostoic import validate_chemical

class NovostoicService:
    novostoic_frontend_baseURL = os.environ.get("NOVOSTOIC_FRONTEND_URL")

    def __init__(self, db) -> None:
        self.db = db
    
    @staticmethod
    async def optstoicResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        job = await db.get(Job, job_id)
        if not job:
            return HTTPException(status_code=404, detail="Job not found")
        
        file = service.get_file(bucket_name, f"{job_id}/out/output.json")
        if not file:
            return None
        result = {}
        data = json.loads(file)
        cache = {}
        config = json.loads(job.job_info)
        if config['primary_precursor'] not in cache:
            cache[config['primary_precursor']] = await validate_chemical(config['primary_precursor'], db)
        if config['target_molecule'] not in cache:
            cache[config['target_molecule']] = await validate_chemical(config['target_molecule'], db)
        result['primaryPrecursor'] = cache[config['primary_precursor']]
        result['targetMolecule'] = cache[config['target_molecule']]
        
        for stoic in data:
            for reactant in stoic['stoichiometry']['reactants']:
                if reactant['molecule'] not in cache:
                    cache[reactant['molecule']] = await validate_chemical(reactant['molecule'], db)
                reactant['molecule'] = cache[reactant['molecule']]
            for product in stoic['stoichiometry']['products']:
                if product['molecule'] not in cache:
                    cache[product['molecule']] = await validate_chemical(product['molecule'], db)
                product['molecule'] = cache[product['molecule']]
        result['results'] = data        
        
        return result
    
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
            for step in pathway:
                if step['primaryPrecursor'] not in cache:
                    cache[step['primaryPrecursor']] = await validate_chemical(step['primaryPrecursor'], db)
                    
                if step['targetMolecule'] not in cache:
                    cache[step['targetMolecule']] = await validate_chemical(step['targetMolecule'], db)
                
                step['primaryPrecursor'] = cache[step['primaryPrecursor']]
                step['targetMolecule'] = cache[step['targetMolecule']]
                
                for reactant in step['reactants']:
                    if reactant['molecule'] not in cache:
                        cache[reactant['molecule']] = await validate_chemical(reactant['molecule'], db)
                    reactant['molecule'] = cache[reactant['molecule']]
                
                for product in step['products']:
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
        data = json.loads(file)
        
        # enzyme rank runs on smiles, which passed through add_info field with format <component_id>:<smiles>
        _, smiles = data['addInfo'].split(':')
        data['primaryPrecursor'] = await validate_chemical(smiles, db)

        return data
    
    @staticmethod
    async def dgPredictorResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        results = json.loads(service.get_file(bucket_name, f"{job_id}/out/output.json"))
        
        molecule_cache = {}
        for result in results:
            # get info for each molecule
            for molecule in result['molecules']:
                y = result['molecules'][molecule]
                result['molecules'][molecule] = { 'amount': y }
                if molecule not in result['novelMolecules']:
                    if not molecule in molecule_cache:
                        molecule_cache[molecule] = await validate_chemical(molecule, db)
                    result['molecules'][molecule].update(molecule_cache[molecule])
                else:
                    value = result['novelMolecules'][molecule]
                    if not value in molecule_cache:
                        molecule_cache[value] = draw_chemical_svg(value)
                        
                    result['molecules'][molecule].update({
                        'is_cofactor': False,
                        'type': 'inchi' if value.lower().startswith('inchi') else 'smiles',
                        'value': value,
                        'structure': molecule_cache[value],
                    })

        return results
