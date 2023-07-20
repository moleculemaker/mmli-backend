import time
import uuid

from fastapi import Depends, HTTPException, APIRouter
from sqlmodel import Session, select

from models.enums import JobType, JobStatus
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobCreate, JobUpdate

router = APIRouter()


@router.post("/{job_type}/jobs", response_model=Job, tags=['Jobs'])
async def create_job(job: JobCreate, job_type: str, db: Session = Depends(get_session)):
    # Check that body was provided
    if not job:
        raise HTTPException(status_code=400, detail="'job' is required")

    # Check if this job_id already exists
    existing_job = db.execute(select(Job).where(Job.job_id == job.job_id).where(Job.run_id == job.run_id)).first()
    if existing_job:
        raise HTTPException(status_code=409, detail="Job already exists with job_id=" + job.job_id + " and run_id=" + job.run_id)

    # Validate Job type
    # TODO: Set command+image based on job_type
    if job_type == JobType.CHEMSCRAPER:
        print("Creating CHEMSCRAPER job")
        # runs as a background_task
        command = ''
        image = ''
    elif job_type == JobType.MOLLI:
        print("Creating MOLLI job")
        # runs in Kubernetes
        image = 'moleculemaker/molli:ncsa-workflow'
        command = ''  # insert reference to input file
    elif job_type == JobType.CLEAN:
        print("Creating CLEAN job")
        # runs in Kubernetes
        image = 'moleculemaker/clean-image-amd64'
        command = ''  # insert reference to input file
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
        job_id=str(uuid.uuid4()) if job.job_id is None else job.job_id,
        run_id=str(uuid.uuid4()) if job.run_id is None else job.run_id,

        type=job_type,
        command=command,
        image=image,

        # Job metadata
        phase=JobStatus.PENDING,
        deleted=0,
        time_created=int(time.time()),
        time_start=0,
        time_end=0,

        # Set ser metadata
        user_agent=user_agent,
    )

    db.add(db_job)
    db.commit()

    return db_job


@router.post("/{job_type}/jobs", response_model=Job, tags=['Jobs'])
async def create_job(job: JobCreate, job_type: str, db: Session = Depends(get_session)):
    # Check if this job_id already exists
    existing_job = db.execute(select(Job).where(Job.job_id == job.job_id).where(Job.run_id == job.run_id)).first()
    if existing_job:
        raise HTTPException(status_code=409, detail="Job already exists with job_id=" + job.job_id + " and run_id=" + job.run_id)

    # Validate Job type
    # TODO: Set command+image based on job_type
    if job_type == JobType.CHEMSCRAPER:
        print("Creating CHEMSCRAPER job")
        # runs as a background_task
        command = ''
        image = ''
    elif job_type == JobType.MOLLI:
        print("Creating MOLLI job")
        # runs in Kubernetes
        image = 'moleculemaker/molli:ncsa-workflow'
        command = ''  # insert reference to input file
    elif job_type == JobType.CLEAN:
        print("Creating CLEAN job")
        # runs in Kubernetes
        image = 'moleculemaker/clean-image-amd64'
        command = ''  # insert reference to input file
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
        job_id=str(uuid.uuid4()) if job.job_id is None else job.job_id,
        run_id=str(uuid.uuid4()) if job.run_id is None else job.run_id,

        type=job_type,
        command=command,
        image=image,

        # Job metadata
        phase=JobStatus.PENDING,
        deleted=0,
        time_created=int(time.time()),
        time_start=0,
        time_end=0,

        # Set ser metadata
        user_agent=user_agent,
    )

    db.add(db_job)
    db.commit()

    return db_job

@router.get("/{job_type}/jobs", tags=['Jobs'])
async def list_jobs_by_type(job_type: str, db: Session = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    results = db.execute(select(Job).where(Job.type == job_type))
    jobs = results.scalars().all()
    return jobs


@router.get("/{job_type}/jobs/{job_id}", tags=['Jobs'])
async def list_jobs_by_type_and_jobID(job_type: str, job_id: str, db: Session = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    #job = db.get(Job, job_id)

    results = db.execute(select(Job).where(Job.type == job_type).where(Job.job_id == job_id))
    jobs = results.scalars().all()
    return jobs


@router.get("/{job_type}/jobs/{job_id}/{run_id}", tags=['Jobs'])
async def get_job_by_type_and_job_id_and_run_id(job_type: str, job_id: str, run_id: str, db: Session = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    #job = db.get(Job, job_id)
    job = db.execute(select(Job).where(Job.job_id == job_id).where(Job.run_id == run_id)).first()
    if not job:
        raise HTTPException(status_code=404, detail="No job found with job_id=" + job_id)
    return job.Job


@router.put("/{job_type}/jobs/{job_id}/{run_id}", response_model=Job, tags=['Jobs'])
async def update_job(job: Job, job_type: str, db: Session = Depends(get_session)):
    # Check if this job_id already exists
    existing_job = db.execute(select(Job)
                              .where(Job.type == job_type)
                              .where(Job.job_id == job.job_id)
                              .where(Job.run_id == job.run_id)).first()
    if not existing_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job.job_id} and run_id={job.run_id}")

    # TODO: Validation
    # TODO: Set internal metadata / fields

    # Update existing DB job from user input
    db.add(job)
    db.commit()

    return job


@router.patch("/{job_type}/jobs/{job_id}/{run_id}", response_model=Job, tags=['Jobs'])
async def patch_job(job: JobUpdate, job_type: str, db: Session = Depends(get_session)):
    # Check if this job_id already exists
    db_job = db.get(Job, job.id)
    if not db_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job.job_id} and run_id={job.run_id}")

    # TODO: Validation
    job_data = job.dict(exclude_unset=True)

    # TODO: Set internal metadata / fields
    for key, value in job_data.items():
        setattr(db_job, key, value)

    # Update existing DB job from user input
    db.add(db_job)
    db.commit()
    db.refresh(db_job)

    return db_job
