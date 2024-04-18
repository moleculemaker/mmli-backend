import base64
import re
import time
import uuid
from typing import List

from fastapi import Depends, HTTPException, APIRouter, File, UploadFile
from sqlalchemy import delete
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from config import get_logger, app_config
from models.enums import JobType, JobStatus, JobTypes
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobCreate, JobUpdate

from services import kubejob_service, clean_service, molli_service

router = APIRouter()

log = get_logger(__name__)


@router.post("/{job_type}/jobs", response_model=Job, tags=['Jobs'], description="Create a new run for a new or existing Job")
async def create_job(job: JobCreate, job_type: str, files: List[UploadFile] = File(...), db: AsyncSession = Depends(get_session)):
    job_id = str(uuid.uuid4()) if job.job_id is None else job.job_id

    # Check if this job_id already exists
    existing_jobs = await db.execute(select(Job).where(Job.job_id == job_id))
    if existing_jobs.first():
        raise HTTPException(status_code=409, detail="Job already exists with job_id=" + job_id)

    # Validate Job type
    # TODO: Set command+image based on job_type
    if job_type == JobType.CHEMSCRAPER:
        log.debug("Creating CHEMSCRAPER job")
        # runs as a background_task
        command = 'N/A'
        image_name = 'N/A'
    elif job_type in JobTypes:
        log.debug(f"Creating Kubernetes job: {job_type}")
        # runs in Kubernetes, read Docker image name from config
        image_name = app_config['kubernetes_jobs'][job_type]['image']

        # Command + environment are set differently for each job (see below)
        command = ''
        environment = []

        if job_type == JobType.DEFAULT:
            command = app_config['kubernetes_jobs'][job_type]['command']
        elif job_type == JobType.SOMN:
            project_id = '44eb8d94effa11eea46f18c04d0a4970'
            model_set = 'apr-2024'
            new_predictions_name = 'asdf'

            # TODO: Build up example_request.csv from user input
            input_file_path = f'somn_scratch/{project_id}/scratch/example_request.csv'

            command = None
            #command = f"ls -al && whoami && somn predict {project_id} {model_set} {new_predictions_name}"

        # TODO: support NOVOSTOIC/CLEAN/MOLLI job types
        elif job_type == JobType.CLEAN or job_type == JobType.MOLLI or 'novostoic' in job_type:
            raise HTTPException(status_code=501, detail=f'{job_type} not yet implemented')

        elif job_type == JobType.CLEAN:
            # TODO: Build up input.FASTA from user input
            command = clean_service.build_clean_job_command(job_id=job_id, job_info=job.job_info)
        elif job_type == JobType.MOLLI:
            # TODO: Map/pass/mount CORES/SUBS files into the container
            command = app_config['kubernetes_jobs'][job_type]['command']
            environment = molli_service.build_molli_job_environment(job_id=job_id, files=files)

        # Run a Kubernetes Job with the given image + command + environment
        log.debug("Creating Kubernetes job: " + job_type)
        kubejob_service.create_job(job_type=job_type, job_id=job_id, image_name=image_name, command=command, environment=environment)

    else:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    # TODO: Set user_agent based on requestor
    user_agent = ''

    # TODO: Validation
    # TODO: Set internal metadata / fields
    # Create a new DB job from user input
    db_job = Job(
        # User input
        email=job.email,
        job_info=job.job_info,
        job_id=job_id,
        run_id=job.run_id,

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

    return db_job


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
        raise HTTPException(status_code=404, detail="No job found with job_id=" + job_id)
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
