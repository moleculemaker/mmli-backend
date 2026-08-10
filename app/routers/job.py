import time
import traceback

from typing import List, Optional

# Body/Path come from `fastapi`, not `fastapi.params`. The latter holds the underlying
# parameter classes; the public helpers are what the framework expects as defaults, and
# only they accept the keyword form used below.
from fastapi import Body, Depends, HTTPException, APIRouter, Path, UploadFile
from sqlalchemy import delete
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession
from starlette import status
from starlette.responses import JSONResponse

from config import get_logger, app_config, STATUS_ERROR
from models.enums import JobType, JobStatus, JobTypes
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobCreate, JobUpdate

from services import job_builder, kubejob_service
from services.minio_service import MinIOService

router = APIRouter()

log = get_logger(__name__)

async def _mark_job_failed(db: AsyncSession, db_job: Job) -> None:
    """Move a just-created job to 'error' when nothing will ever run to advance it.

    Only ever called for a row this request created. A job that already existed must
    keep whatever phase it legitimately reached.
    """
    db_job.phase = JobStatus.ERROR
    db.add(db_job)
    await db.commit()


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

    # Validate Job type
    # TODO: Set command+image based on job_type
    if job_type in JobTypes:
        # Look the job up only once the type is known to be valid, so this cannot
        # depend on how strictly the database enforces the column.
        #
        # The deployed schema (built by alembic) stores `type` as a varchar and would
        # simply match no rows. A schema built from the SQLModel metadata instead - which
        # is what SQLModel.metadata.create_all() produces, and what the test suite uses -
        # maps JobType to a native Postgres enum, where an unrecognized value makes the
        # driver raise InvalidTextRepresentationError and the request 500s.
        #
        # Validating first makes the behavior identical under both, which matters
        # because those two schemas are not currently guaranteed to agree.
        statement = select(Job).where(Job.type == job_type).where(Job.job_id == job_id)
        existing_jobs = await db.exec(statement)
        db_job: Job = existing_jobs.first()
        #if db_job:
        #    raise HTTPException(status_code=409, detail=f"Job already exists with job_id={job_id}")

        # Per-tool preparation (image, command, environment, and any input files a
        # tool needs written before it starts) lives in services/job_builder.py, because
        # the versioned API needs exactly the same work. Some tools rewrite job_id or
        # job_info, so both come back out.
        prepared = job_builder.prepare_job(
            job_type=job_type, job_id=job_id, job_info=job_info, service=service,
        )
        job_id = prepared.job_id
        job_info = prepared.job_info
        image_name = prepared.image_name
        command = prepared.command
        environment = prepared.environment
        parent_job_id = prepared.parent_job_id

        # TODO: Set user_agent based on requestor
        user_agent = ''

        # TODO: Validation
        # TODO: Set internal metadata / fields
        #
        # Write the DB row BEFORE creating the Kubernetes Job. KubeWatcher looks this row up
        # by jobId for every event it receives and drops events whose row does not exist yet;
        # since the watch never replays old events, a Job that got created (and, for the short
        # ones, finished) inside that window is stranded at 'queued' even after Kubernetes
        # reports it complete.
        created_new_job = not db_job
        if created_new_job:
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
                parent_job_id=parent_job_id,

                # Job metadata
                deleted=0,
                time_created=int(time.time()),

                # Set ser metadata
                user_agent=user_agent,
            )

            db.add(db_job)
            await db.commit()
            await db.refresh(db_job)

        # A POST naming a job_id that already exists is not a resubmission, and this
        # endpoint has no authentication: the caller has proved only that they can name
        # an id, which is not authority over the job that holds it. Stop here.
        #
        # Previously this path fell through and (a) launched a second Kubernetes Job
        # under the same name, (b) replaced the stored job_info with the caller's, and
        # (c) returned the owner's email address and inputs in the response body.
        # The response keys and status code are unchanged so existing clients - notably
        # coordinator.py, which reads job_id from a 200 or a 201 - keep working.
        if not created_new_job:
            log.warning(f'Ignoring duplicate job submission for existing job_id={job_id} (type={job_type})')
            return JSONResponse(status_code=status.HTTP_200_OK, content={
                'job_id': str(db_job.job_id),
                'run_id': str(db_job.run_id),
                # Withheld: naming an id is not authorisation to read its owner's
                # address or the inputs they submitted.
                'email': None,
                'job_info': None,
            })

        # Run a Kubernetes Job with the given image + command + environment
        try:
            log.debug(f"Creating Kubernetes job[{job_type}]: " + job_id)
            create_response = kubejob_service.create_job(
                job_type=job_type,
                job_id=job_id,
                run_id=run_id,
                image_name=image_name,
                command=command,
                environment=environment
            )
        except Exception as ex:
            log.error("Failed to create Job: " + str(ex))
            log.error(traceback.format_exc())
            # The row we just wrote is already visible to clients, and nothing will ever run
            # to move it along - fail it here rather than leaving it polling 'queued' forever.
            await _mark_job_failed(db, db_job)
            raise HTTPException(status_code=400, detail="Failed to create Job: " + str(ex))

        # create_job catches ApiException internally, logs it, and returns a response
        # carrying status=error rather than propagating it. Kubernetes rejections -
        # an exceeded quota, an unpullable image, a name that already exists - are
        # therefore invisible to the except clause above. Without this check the caller
        # is told 201 while the row sits at 'queued' forever with no pod behind it.
        if isinstance(create_response, dict) and create_response.get('status') == STATUS_ERROR:
            message = create_response.get('message', 'unknown error')
            log.error(f'Kubernetes rejected job[{job_type}] {job_id}: {message}')
            await _mark_job_failed(db, db_job)
            raise HTTPException(status_code=400, detail="Failed to create Job: " + str(message))

        return JSONResponse(status_code=status.HTTP_201_CREATED, content={
            'job_id': str(db_job.job_id),
            'run_id': str(db_job.run_id),
            'email': str(db_job.email),
            'job_info': str(db_job.job_info),
        })

    else:
        log.error("Failed to create job - invalid job type: " + str(job_type))
        raise HTTPException(status_code=400, detail="Invalid job type: " + str(job_type))


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
    # Validate before querying, as every other route on this router does. Filtering on
    # an unrecognized type is only harmless if the column is a plain varchar; against a
    # schema built from the SQLModel metadata it is a native enum and the driver raises.
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    result = await db.execute(select(Job).where(Job.type == job_type).where(Job.job_id == job_id).where(Job.run_id == run_id))
    db_job = result.first()
    if not db_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job_id} and run_id={run_id}")

    # Delete the Kubernetes Job before dropping the row. Removing only the row leaves a
    # pod running with nothing tracking it: KubeWatcher looks each event up by jobId and
    # discards events whose row is gone, so the work continues, consumes cluster
    # resources, and writes its output to MinIO with no way to observe or stop it.
    #
    # Deliberately best-effort. delete_job already swallows ApiException (a Job that
    # Kubernetes has since reaped is a 404, which is the expected case for anything
    # older than ttlSecondsAfterFinished), and a cluster problem must not leave the
    # caller unable to delete their own record.
    try:
        kubejob_service.delete_job(job_type=job_type, job_id=job_id)
    except Exception as ex:
        log.error(f'Failed to delete Kubernetes job[{job_type}] {job_id}, deleting DB row anyway: {ex}')

    await db.execute(delete(Job).where(Job.type == job_type).where(Job.job_id == job_id).where(Job.run_id == run_id))
    await db.commit()

    return db_job.Job
