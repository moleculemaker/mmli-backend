import base64
import os
import json
import re
import time
import uuid
import csv
import io
import traceback

from typing import List

from fastapi import Depends, HTTPException, APIRouter, UploadFile
from fastapi.openapi.models import Response
from fastapi.params import Path, Body, File
from pydantic.fields import Annotated, Optional
from sqlalchemy import delete
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession
from starlette import status
from starlette.responses import JSONResponse

from config import get_logger, app_config
from models.enums import JobType, JobStatus, JobTypes
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobCreate, JobUpdate

from services import kubejob_service
from services.clean_service import CleanService
from services.molli_service import MolliService
from services.minio_service import MinIOService
from services.somn_service import SomnService

router = APIRouter()

log = get_logger(__name__)


@router.post("/{job_type}/jobs", response_model=Job, tags=['Jobs'], description="Create a new run for a new or existing Job")
async def create_job(
        job_id: Optional[str] = Body(default=None),
        run_id: Optional[str] = Body(default=None),
        email: Optional[str] = Body(default=None),
        job_info: Optional[str] = Body(default="{}"),
        job_type: str = Path(),
        service: MinIOService = Depends(),
        db: AsyncSession = Depends(get_session)
):
    job_id = job_id if job_id else str(kubejob_service.generate_uuid())
    #run_id = run_id if run_id else str(kubejob_service.generate_uuid())

    # Check if this job_id already exists
    statement = select(Job).where(Job.job_id == job_id)
    existing_jobs = await db.exec(statement)
    db_job: Job = existing_jobs.first()
    #if db_job:
    #    raise HTTPException(status_code=409, detail=f"Job already exists with job_id={job_id}")

    # Validate Job type
    # TODO: Set command+image based on job_type
    if job_type == JobType.CHEMSCRAPER:
        raise HTTPException(status_code=501, detail=f"job_type={job_type} is not yet supported by the Jobs API")

        # TODO: Create a simple job container for executing ChemScraper HTTP requests

        #log.debug("Creating CHEMSCRAPER job")
        # runs as a background_task
        #command = 'N/A'
        #image_name = 'N/A'
    elif job_type in JobTypes:
        log.debug(f"Creating Kubernetes job: {job_type}")
        # runs in Kubernetes, read Docker image name from config
        image_name = app_config['kubernetes_jobs'][job_type]['image']

        # Command + environment are set differently for each job (see below)
        command = ''
        environment = []

        # Mount in secrets/volumes at runtime
        volumes = []
        secrets = []

        if job_type == JobType.DEFAULT:
            command = app_config['kubernetes_jobs'][job_type]['command']
            #command = f'ls -al /uws/jobs/{job_type}/{job_id}'

        elif job_type == JobType.ACERETRO:
            # ACERetro jobs
            # Example usage:
            # curl -X POST https://mmli.kastan.ai/aceretro/jobs \
            #   -H "Content-Type: application/json" \
            #   -d '{
            #     "job_id": "123",
            #     "email": "user@gmail.com",
            #     "job_info": "{\"smiles\": \"O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO\"}"
            #   }'
            log.info(f"------------------ STARTING ACERETRO JOB ------------------  job[{job_type}]: " + job_id)
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            
            from config import MINIO_SERVER, MINIO_ACCESS_KEY, MINIO_SECRET_KEY
            environment = [{
                'name': 'MINIO_URL',
                # 'value': app_config['minio']['apiBaseUrl']
                'value': MINIO_SERVER
            },
            {
                'name': 'MINIO_ACCESS_KEY',
                # 'value': app_config['minio']['accessKey']
                'value': MINIO_ACCESS_KEY
            },
            {
                'name': 'MINIO_SECRET_ACCESS_KEY',
                # 'value': app_config['minio']['accessKey']
                'value': MINIO_SECRET_KEY
            },
            {
                'name': 'MINIO_SECURE',
                # 'value': app_config['minio']['secretKey']
                'value': False # HARD CODED, like everywhere else in this repo
            }]

            command = f"python entrypoint.py --job_id {job_id}"
            # Job is created at end of function
        elif job_type == JobType.REACTIONMINER:
            log.debug(f'Running ReactionMiner job: {job_id}')
            environment = app_config['kubernetes_jobs']['reactionminer']['env']

        elif job_type == JobType.OED_DLKCAT or job_type == JobType.OED_UNIKP or job_type == JobType.OED_CATPRED:
            # Example: "job_info": "{\"input_pairs\":[{\"name\":\"example\",\"sequence\":\"MEDIPDTSRPPLKYVK...\",\"type\":\"FASTA\",\"smiles\":\"OC1=CC=C(C[C@@H](C(O)=O)N)C=C1\"}]}"
            log.debug(f'Running OpenEnzymeDB job: {job_type} - {job_id}')
            environment = [{'OED_JOB_INFO': json.dumps(job_info)}]

        elif job_type == JobType.SOMN:
            #  Build up example_request.csv from user input, upload to MinIO?
            job_config = json.loads(job_info.replace('\"', '"'))
            
            # Canonicalize SMILES and update names from reference files
            for config in job_config:
                if config['el_input_type'] == 'smi':
                    config['el'] = SomnService.canonicalize_smiles(config['el'])
                if config['nuc_input_type'] == 'smi':
                    config['nuc'] = SomnService.canonicalize_smiles(config['nuc'])
            
            # Generate unique name mappings
            job_config, el_name_map = SomnService.generate_name_mapping(job_config, 'el')
            job_config, nuc_name_map = SomnService.generate_name_mapping(job_config, 'nuc')
            
            # job_config, el_name_map = SomnService.update_names_from_reference(job_config, el_name_map, 'el')
            # job_config, nuc_name_map = SomnService.update_names_from_reference(job_config, nuc_name_map, 'nuc')
            
            print('updated config: ', job_config, el_name_map, nuc_name_map)
            
            # job_info is automatically stored in Postgres to retain user input
            job_info = json.dumps({
                'info': job_config,
                'el_name_map': el_name_map,
                'nuc_name_map': nuc_name_map
            })
            
            # Some SMILESs fail to pass 3d generation test, so we process them as mol2 input
            for config in job_config:
                config['el'] = SomnService.process_molecule_input(config, 'el')
                config['nuc'] = SomnService.process_molecule_input(config, 'nuc')
            
            
            file = io.StringIO()
            writer = csv.writer(file)
            writer.writerow([
                "user", 
                "nuc", 
                "el", 
                "nuc_name", 
                "el_name",
                "nuc_idx",
                "el_idx"
            ])
            for config in job_config:
                writer.writerow([
                    config["reactant_pair_name"], 
                    config["nuc"],
                    config["el"], 
                    config["nuc_name"], 
                    config["el_name"],
                    config["nuc_idx"],
                    config["el_idx"]
                ])
            
            upload_result = service.upload_file(job_type, f"/{job_id}/in/example_request.csv", file.getvalue().encode('utf-8'))
            if not upload_result:
                raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")

            # We assume that file has already been uploaded to MinIO
            somn_project_dir = app_config['kubernetes_jobs']['somn']['projectDirectory']

            command = app_config['kubernetes_jobs'][job_type]['command']
            #command = f"ls -al && whoami && somn predict {project_id} {model_set} {new_predictions_name}"

            environment = [{
                'name': 'SOMN_PROJECT_DIR',
                'value': somn_project_dir
            }]
        
        elif job_type == JobType.NOVOSTOIC_OPTSTOIC:
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']

            # environment = [{
            #     # TBD... 
            #     # 'name': 'SOMN_PROJECT_DIR',
            #     # 'value': somn_project_dir
            # }]

            # Run a Kubernetes Job with the given image + command + environment
            try:
                log.debug(f"Creating Kubernetes job[{job_type}]: " + job_id)
            except Exception as ex:
                log.error("Failed to create Job: " + str(ex))
                raise HTTPException(status_code=400, detail="Failed to create Job: " + str(ex))
            
        elif job_type == JobType.NOVOSTOIC_PATHWAYS:
            if service.ensure_bucket_exists(job_type):
                job_info = json.loads(job_info.replace('\"', '"'))
                stoic = f'{job_info["substrate"]["amount"]} {job_info["substrate"]["molecule"]}'
                for coReactant in job_info['reactants']:
                    stoic += f' + {coReactant["amount"]} {coReactant["molecule"]}'
                stoic += " <=> "
                for coProduct in job_info['products']:
                    stoic += f'{coProduct["amount"]} {coProduct["molecule"]} + '
                stoic += f'{job_info["product"]["amount"]} {job_info["product"]["molecule"]}'
                
                job_info['stoic'] = stoic
                job_info['substrate'] = job_info['substrate']['molecule']
                job_info['product'] = job_info['product']['molecule']
                job_info['num_enzymes'] = job_info['num_enzymes'] if 'num_enzymes' in job_info else 0
                
                job_info = json.dumps(job_info)
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']
            
        elif job_type == JobType.NOVOSTOIC_ENZRANK:
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']
            
        elif job_type == JobType.NOVOSTOIC_DGPREDICTOR:
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']

        elif job_type == JobType.CLEAN:
            # Build up input.FASTA from user input
            job_config = json.loads(job_info.replace('\"', '"'))
            command = CleanService.build_clean_job_command(job_id=job_id, job_info=job_config)
        elif job_type == JobType.MOLLI:
            # Pass path to CORES/SUBS files into the container
            command = app_config['kubernetes_jobs'][job_type]['command']
            job_config = json.loads(job_info.replace('\"', '"'))
            environment = MolliService.build_molli_job_environment(job_id=job_id, job_info=job_config)

        # Run a Kubernetes Job with the given image + command + environment
        try:
            log.debug(f"Creating Kubernetes job[{job_type}]: " + job_id)
            kubejob_service.create_job(job_type=job_type, job_id=job_id, run_id=run_id, image_name=image_name, command=command, environment=environment)
        except Exception as ex:
            log.error("Failed to create Job: " + str(ex))
            log.error(traceback.format_exc())
            raise HTTPException(status_code=400, detail="Failed to create Job: " + str(ex))

    else:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    # TODO: Set user_agent based on requestor
    user_agent = ''

    # TODO: Validation
    # TODO: Set internal metadata / fields
    if not db_job:
        # Create a new DB job from user input
        db_job: Job = Job(
            # User input
            email=email,
            job_info=job_info,
            job_id=job_id,
            run_id=run_id,

            type=job_type,
            command=command,
            image=image_name,

            # Job metadata
            deleted=0,
            time_created=int(time.time()),

            # Set ser metadata
            user_agent=user_agent,
        )

        db.add(db_job)
        await db.commit()
        await db.refresh(db_job)

        return JSONResponse(status_code=status.HTTP_201_CREATED, content={
            'job_id': str(db_job.job_id),
            'run_id': str(db_job.run_id),
            'email': str(db_job.email),
            'job_info': str(db_job.job_info),
        })

    return JSONResponse(status_code=status.HTTP_200_OK, content={
        'job_id': str(db_job.job_id),
        'run_id': str(db_job.run_id),
        'email': str(db_job.email),
        'job_info': str(db_job.job_info),
    })


@router.get("/{job_type}/jobs", tags=['Jobs'], description="Get a list of all job runs by type")
async def list_jobs_by_type(job_type: str, db: AsyncSession = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    results = await db.execute(select(Job).where(Job.type == job_type))
    jobs = results.scalars().all()
    return jobs


@router.get("/{job_type}/jobs/{job_id}", tags=['Jobs'], description="Get a list of all job runs by type and job_id")
async def list_jobs_by_type_and_job_id(job_type: str, job_id: str, db: AsyncSession = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    results = await db.execute(select(Job).where(Job.type == job_type).where(Job.job_id == job_id))
    jobs = results.scalars().all()
    return jobs


@router.get("/{job_type}/jobs/{job_id}/{run_id}", tags=['Jobs'], description="Get a single job by type, job_id, and run_id")
async def get_job_by_type_and_job_id_and_run_id(job_type: str, job_id: str, run_id: str, db: AsyncSession = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    result = await db.execute(select(Job).where(Job.job_id == job_id).where(Job.run_id == run_id))
    job = result.first()
    if not job:
        raise HTTPException(status_code=404, detail=f"No job found with job_id={job_id} and run_id={run_id}")
    return job.Job


@router.put("/{job_type}/jobs/{job_id}/{run_id}", response_model=Job, tags=['Jobs'], description="Overwrite all writeable fields of an existing Job")
async def update_existing_job(job: Job, job_type: str, db: AsyncSession = Depends(get_session)):
    # Check if this job_id already exists
    db_job = await db.get(Job, job.job_id)
    if not db_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job.job_id} and run_id={job.run_id}")

    # TODO: Validation
    # TODO: Set internal metadata / fields
    db_job.phase = job.phase
    db_job.time_start = job.time_start
    db_job.time_end = job.time_end

    # Update existing DB job from user input
    db.add(db_job)
    await db.commit()
    await db.refresh(db_job)

    return db_job


@router.patch("/{job_type}/jobs/{job_id}/{run_id}", response_model=Job, tags=['Jobs'], description="Update one or more fields of an existing Job")
async def patch_existing_job(job: JobUpdate, job_type: str, db: AsyncSession = Depends(get_session)):
    # Check if this job_id already exists
    db_job = await db.get(Job, job.job_id)
    if not db_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job.job_id} and run_id={job.run_id}")

    # TODO: Validation
    # TODO: Set internal metadata / fields
    db_job.phase = job.phase if job.phase is not None else db_job.phase
    db_job.time_start = job.time_start if job.time_start is not None else db_job.time_start
    db_job.time_end = job.time_end if job.time_end is not None else db_job.time_end

    # Update existing DB job from user input
    db.add(db_job)
    await db.commit()
    await db.refresh(db_job)

    return db_job


@router.delete("/{job_type}/jobs/{job_id}/{run_id}", tags=['Jobs'], description="Delete a single Job by type, job_id, and run_id")
async def delete_job_by_type_and_job_id_and_run_id(job_type: str, job_id: str, run_id: str, db: AsyncSession = Depends(get_session)):
    # Check if this job_id already exists
    result = await db.execute(select(Job).where(Job.type == job_type).where(Job.job_id == job_id).where(Job.run_id == run_id))
    db_job = result.first()
    if not db_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job_id} and run_id={run_id}")

    await db.execute(delete(Job).where(Job.type == job_type).where(Job.job_id == job_id).where(Job.run_id == run_id))
    await db.commit()

    return db_job.Job
