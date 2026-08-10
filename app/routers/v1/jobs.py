"""Job submission and retrieval for the versioned API.

Differences from the legacy Jobs API that matter to a scripted client:

  job ids are server-assigned      naming an id is not authority over it, and this API
                                   is unauthenticated
  inputs are validated             against the tool's published JSON Schema, with JSON
                                   Pointers to whatever was wrong
  files ride along with submission a single multipart request, rather than a separate
                                   upload that has to guess the job id first
  retries are safe                 via Idempotency-Key
  status codes are meaningful      results are 409 while running, not 200 with a null
                                   body, so a polling loop can branch on the status line
  responses carry links            so clients never build URLs by string concatenation
"""
import json
import time
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, File, Form, Header, Request, UploadFile
from jsonschema import Draft202012Validator
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from config import STATUS_ERROR, app_config, get_logger
from models.enums import JobStatus, JobTypes
from models.sqlmodel.db import get_session
from models.sqlmodel.models import IdempotencyKey, Job
from routers.v1.problems import (
    IDEMPOTENCY_CONFLICT, INVALID_INPUT, JOB_FAILED, JOB_NOT_FINISHED, NOT_CANCELABLE,
    UNKNOWN_JOB, UNKNOWN_TOOL, UPSTREAM_FAILURE, ProblemException,
)
from services import job_builder, kubejob_service, tool_registry
from services.minio_service import MinIOService

router = APIRouter()
log = get_logger(__name__)

TERMINAL_PHASES = {JobStatus.COMPLETED, JobStatus.ERROR, JobStatus.CANCELED}

# Files a tool writes for operational reasons rather than as scientific output. Hidden
# from the artifact listing unless asked for: error logs carry whatever the container
# chose to print, which can include internal paths and configuration.
LOG_ARTIFACTS = {'output.log', 'error.log', 'errors.txt', 'success', 'fail'}

# How long a client should wait before polling again. Jobs here run for minutes to
# hours, so a tight loop helps nobody.
RETRY_AFTER_SECONDS = 10


def _base_url(request: Request) -> str:
    return str(request.base_url).rstrip('/') + '/v1'


def _iso(epoch: Optional[int]) -> Optional[str]:
    """Epoch seconds to ISO-8601 UTC. 0 is the schema's "unset", not 1970."""
    if not epoch:
        return None
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(epoch))


def _job_document(job: Job, base: str) -> Dict[str, Any]:
    """Render a job row as the /v1 job resource."""
    links = {
        'self': f'{base}/jobs/{job.job_id}',
        'results': f'{base}/jobs/{job.job_id}/results',
        'artifacts': f'{base}/jobs/{job.job_id}/artifacts',
        'tool': f'{base}/tools/{job.type}',
        'cancel': f'{base}/jobs/{job.job_id}/cancel',
    }
    if job.parent_job_id:
        links['parent'] = f'{base}/jobs/{job.parent_job_id}'

    return {
        'job_id': job.job_id,
        'tool': job.type,
        'status': job.phase,
        'parent_job_id': job.parent_job_id,
        'submitted_at': _iso(job.time_created),
        'started_at': _iso(job.time_start),
        'finished_at': _iso(job.time_end),
        'inputs': _decode_job_info(job.job_info),
        'provenance': {
            # What was configured, which may be an untagged or mutable reference...
            'image': job.image,
            # ...and what actually ran. Null when the pod was reaped before the watcher
            # could read it, or for jobs predating digest capture.
            'image_digest': job.image_digest,
            'command': job.command,
        },
        'links': links,
    }


def _decode_job_info(job_info: Optional[str]) -> Any:
    """Return job_info as JSON when it parses, otherwise as the raw string.

    Several tools store a derived form rather than what the user sent, and one stores a
    bare sequence string, so this cannot assume valid JSON.
    """
    if job_info is None:
        return None
    try:
        return json.loads(job_info)
    except (ValueError, TypeError):
        return job_info


def _validate_against_schema(tool: str, inputs: Any) -> None:
    """Validate submitted inputs, reporting every failure with a JSON Pointer."""
    schema = tool_registry.get_input_schema(tool)
    if schema is None:
        return

    errors = [
        {
            'pointer': '/' + '/'.join(str(p) for p in error.absolute_path),
            'detail': error.message,
        }
        for error in sorted(Draft202012Validator(schema).iter_errors(inputs),
                            key=lambda e: list(e.absolute_path))
    ]
    if errors:
        raise ProblemException(
            422, 'Input failed schema validation',
            f'The submitted inputs do not match the schema published for "{tool}".',
            type_uri=INVALID_INPUT, errors=errors,
            schema_url=f'/v1/tools/{tool}/input-schema',
        )


def _require_registered_tool(tool: str) -> Dict[str, Any]:
    entry = tool_registry.get_tool(tool)
    if entry is None or entry.get('internal'):
        raise ProblemException(
            404, 'Unknown tool',
            f'No tool is registered with identifier "{tool}".',
            type_uri=UNKNOWN_TOOL, known_tools=sorted(tool_registry.list_tools()),
        )
    return entry


async def _get_job_or_404(db: AsyncSession, job_id: str) -> Job:
    job = await db.get(Job, job_id)
    if job is None:
        raise ProblemException(
            404, 'Unknown job',
            f'No job exists with id "{job_id}".',
            type_uri=UNKNOWN_JOB,
        )
    return job


# The handler returns a plain dict rather than a response model, so FastAPI infers
# nothing about its shape and the docs would show the 201 body as the bare string
# "string". An example is used rather than a model: _job_document is the single place
# that decides this shape, and a parallel Pydantic model would be one more thing to keep
# in step with it.
_JOB_EXAMPLE = {
    'job_id': '4b3a85ee4efa42c3b5e5c1641e0588f4',
    'tool': 'novostoic-optstoic',
    'status': 'queued',
    'parent_job_id': None,
    'submitted_at': '2026-08-10T01:40:32Z',
    'started_at': None,
    'finished_at': None,
    'inputs': {'primary_precursor': 'MNXM1137670', 'target_molecule': 'MNXM26'},
    'provenance': {
        'image': 'moleculemaker/novostoic',
        'image_digest': None,
        'command': 'python ./novostoic-job.py optstoic',
    },
    'links': {
        'self': '/v1/jobs/4b3a85ee4efa42c3b5e5c1641e0588f4',
        'results': '/v1/jobs/4b3a85ee4efa42c3b5e5c1641e0588f4/results',
        'artifacts': '/v1/jobs/4b3a85ee4efa42c3b5e5c1641e0588f4/artifacts',
        'tool': '/v1/tools/novostoic-optstoic',
        'cancel': '/v1/jobs/4b3a85ee4efa42c3b5e5c1641e0588f4/cancel',
    },
}

_JOB_RESPONSE = {
    'content': {'application/json': {'example': _JOB_EXAMPLE}},
}


@router.post('/tools/{tool}/jobs', tags=['Jobs'], status_code=201,
             summary='Submit a job',
             responses={201: {'description': 'Job accepted and queued.', **_JOB_RESPONSE}})
async def submit_job(
    tool: str,
    request: Request,
    inputs: Optional[str] = Form(default=None),
    files: Optional[List[UploadFile]] = File(default=None),
    idempotency_key: Optional[str] = Header(default=None, alias='Idempotency-Key'),
    service: MinIOService = Depends(),
    db: AsyncSession = Depends(get_session),
):
    """Submit a job, as JSON or as multipart with files.

    JSON:       Content-Type: application/json, body is the tool's input document
    Multipart:  an `inputs` part holding that same document as JSON text, plus any
                number of `files` parts

    Multipart exists because several tools need their input files present before the
    job starts. The legacy flow required uploading first to learn a job id and then
    submitting with it, which is what made client-chosen ids necessary. Here the server
    assigns the id and writes the files itself, so neither is needed.
    """
    _require_registered_tool(tool)

    # Read the input document from whichever way it was sent.
    content_type = request.headers.get('content-type', '')
    if content_type.startswith('application/json'):
        raw = await request.body()
        try:
            submitted = json.loads(raw) if raw else {}
        except ValueError as ex:
            raise ProblemException(
                400, 'Body is not valid JSON', str(ex), type_uri=INVALID_INPUT,
            )
        uploads: List[UploadFile] = []
    else:
        try:
            submitted = json.loads(inputs) if inputs else {}
        except ValueError as ex:
            raise ProblemException(
                400, 'The "inputs" part is not valid JSON', str(ex),
                type_uri=INVALID_INPUT,
            )
        uploads = files or []

    _validate_against_schema(tool, submitted)

    # Replaying a key returns the original job rather than starting a second container.
    if idempotency_key:
        existing = (await db.exec(
            select(IdempotencyKey)
            .where(IdempotencyKey.key == idempotency_key)
            .where(IdempotencyKey.job_type == tool)
        )).first()
        if existing:
            job = await db.get(Job, existing.job_id)
            if job is not None:
                return _job_document(job, _base_url(request))
            # The key outlived the job it named, which means the job was deleted. Say so
            # rather than silently minting a new one under a key the client believes
            # already identifies something.
            raise ProblemException(
                409, 'Idempotency key refers to a job that no longer exists',
                f'Key "{idempotency_key}" was used for job {existing.job_id}, which has '
                f'since been deleted. Retry with a new key.',
                type_uri=IDEMPOTENCY_CONFLICT,
            )

    job_id = str(kubejob_service.generate_uuid())

    # Files first: preparation for several tools reads what was uploaded.
    if uploads:
        service.ensure_bucket_exists(tool)
        for upload in uploads:
            content = await upload.read()
            if not service.upload_file(tool, f'{job_id}/in/{upload.filename}', content):
                raise ProblemException(
                    502, 'Could not store uploaded file',
                    f'Storing "{upload.filename}" failed.', type_uri=UPSTREAM_FAILURE,
                )

    job_info = json.dumps(submitted)
    prepared = job_builder.prepare_job(
        job_type=tool, job_id=job_id, job_info=job_info, service=service,
    )

    db_job = Job(
        job_id=prepared.job_id,
        job_info=prepared.job_info,
        type=tool,
        command=prepared.command,
        image=prepared.image_name,
        parent_job_id=prepared.parent_job_id,
        email=None,
        run_id=None,
        deleted=0,
        time_created=int(time.time()),
        user_agent=request.headers.get('user-agent', ''),
    )
    db.add(db_job)
    await db.commit()
    await db.refresh(db_job)

    # Row first, then the cluster: KubeWatcher drops events whose row does not exist yet
    # and the watch never replays, so the reverse order strands short jobs at 'queued'.
    try:
        response = kubejob_service.create_job(
            job_type=tool, job_id=prepared.job_id, run_id=None,
            image_name=prepared.image_name, command=prepared.command,
            environment=prepared.environment,
        )
    except Exception as ex:
        log.error(f'Failed to create Kubernetes job for {prepared.job_id}: {ex}')
        await _fail(db, db_job)
        raise ProblemException(
            502, 'Could not start job', str(ex), type_uri=UPSTREAM_FAILURE,
        )

    # create_job catches ApiException internally and signals failure by return value, so
    # a rejection never reaches the except above.
    if isinstance(response, dict) and response.get('status') == STATUS_ERROR:
        message = response.get('message', 'unknown error')
        log.error(f'Kubernetes rejected job {prepared.job_id}: {message}')
        await _fail(db, db_job)
        raise ProblemException(
            502, 'Could not start job', str(message), type_uri=UPSTREAM_FAILURE,
        )

    if idempotency_key:
        db.add(IdempotencyKey(key=idempotency_key, job_type=tool,
                              job_id=prepared.job_id, time_created=int(time.time())))
        await db.commit()

    return _job_document(db_job, _base_url(request))


async def _fail(db: AsyncSession, db_job: Job) -> None:
    db_job.phase = JobStatus.ERROR
    db.add(db_job)
    await db.commit()


@router.get('/jobs/{job_id}', tags=['Jobs'], summary='Get a job',
            responses={200: {'description': 'The job.', **_JOB_RESPONSE}})
async def get_job(job_id: str, request: Request, db: AsyncSession = Depends(get_session)):
    job = await _get_job_or_404(db, job_id)
    return _job_document(job, _base_url(request))


@router.get('/jobs/{job_id}/results', tags=['Jobs'], summary='Get a job\'s results')
async def get_results(job_id: str, request: Request,
                      service: MinIOService = Depends(),
                      db: AsyncSession = Depends(get_session)):
    """Return post-processed results once the job has finished.

    Status codes carry the meaning, so a polling client need not parse the body:

      200  finished, results attached
      409  still queued or running (Retry-After), or finished without usable results
      404  no such job

    The legacy endpoint answers 200 with a null body for "running", "finished with
    nothing" and "no such job" alike, which a script cannot tell apart.
    """
    job = await _get_job_or_404(db, job_id)

    if job.phase not in TERMINAL_PHASES:
        raise ProblemException(
            409, 'Job has not finished',
            f'Job {job_id} is {job.phase}. Retry after {RETRY_AFTER_SECONDS} seconds.',
            type_uri=JOB_NOT_FINISHED,
            headers={'Retry-After': str(RETRY_AFTER_SECONDS)},
            status_of_job=str(job.phase),
        )

    if job.phase in (JobStatus.ERROR, JobStatus.CANCELED):
        raise ProblemException(
            409, 'Job did not complete successfully',
            f'Job {job_id} finished with status "{job.phase}", so it has no results.',
            type_uri=JOB_FAILED, status_of_job=str(job.phase),
        )

    from routers.files import get_results as legacy_get_results
    results = await legacy_get_results(job.type, job_id, service, db)
    if results is None:
        raise ProblemException(
            409, 'Job completed but produced no results',
            f'Job {job_id} reports "completed" but no result file is present. This '
            f'usually means the tool exited without writing its expected output.',
            type_uri=JOB_FAILED, status_of_job=str(job.phase),
        )

    return {'job_id': job_id, 'tool': job.type, 'results': results}


@router.get('/jobs/{job_id}/artifacts', tags=['Jobs'],
            summary='List the files a job produced')
async def list_artifacts(job_id: str, request: Request, include: Optional[str] = None,
                         service: MinIOService = Depends(),
                         db: AsyncSession = Depends(get_session)):
    """List raw output files.

    Logs and status-marker files are excluded unless `?include=logs`. Their content is
    whatever the container printed, which can carry internal paths and configuration,
    and anyone holding a job id can read them.
    """
    job = await _get_job_or_404(db, job_id)
    base = _base_url(request)
    want_logs = include == 'logs'

    prefix = f'{job_id}/out/'
    names = service.list_files(job.type, prefix, recursive=True) or []
    items = []
    for name in names:
        # list_files returns either object names or objects, depending on caller.
        name = getattr(name, 'object_name', name)
        short = name[len(prefix):] if name.startswith(prefix) else name
        if not short:
            continue
        if not want_logs and short.split('/')[-1] in LOG_ARTIFACTS:
            continue
        items.append({
            'name': short,
            'href': f'{base}/jobs/{job_id}/artifacts/{short}',
        })

    return {
        'job_id': job_id,
        'count': len(items),
        'items': sorted(items, key=lambda i: i['name']),
        'logs_included': want_logs,
    }


@router.get('/jobs/{job_id}/artifacts/{artifact_path:path}', tags=['Jobs'],
            summary='Download one output file')
async def get_artifact(job_id: str, artifact_path: str,
                       service: MinIOService = Depends(),
                       db: AsyncSession = Depends(get_session)):
    from fastapi.responses import Response

    job = await _get_job_or_404(db, job_id)

    # Confine reads to this job's output directory. Without this, a path such as
    # ../../other-job/out/x would read another job's results.
    if '..' in artifact_path or artifact_path.startswith('/'):
        raise ProblemException(
            400, 'Invalid artifact path',
            'Artifact paths are relative to the job output directory.',
            type_uri=INVALID_INPUT,
        )

    content = service.get_file(job.type, f'{job_id}/out/{artifact_path}')
    if content is None:
        raise ProblemException(
            404, 'Unknown artifact',
            f'Job {job_id} has no output file "{artifact_path}".',
            type_uri=UNKNOWN_JOB,
        )
    return Response(content=content, media_type='application/octet-stream')


@router.post('/jobs/{job_id}/cancel', tags=['Jobs'], summary='Cancel a running job')
async def cancel_job(job_id: str, request: Request, db: AsyncSession = Depends(get_session)):
    """Stop a job and mark it canceled.

    Idempotent: canceling an already-canceled job succeeds. A job that has already
    finished cannot be canceled, and says so rather than pretending.
    """
    job = await _get_job_or_404(db, job_id)

    if job.phase == JobStatus.CANCELED:
        return _job_document(job, _base_url(request))

    if job.phase in (JobStatus.COMPLETED, JobStatus.ERROR):
        raise ProblemException(
            409, 'Job has already finished',
            f'Job {job_id} is "{job.phase}" and cannot be canceled.',
            type_uri=NOT_CANCELABLE, status_of_job=str(job.phase),
        )

    try:
        kubejob_service.delete_job(job_type=job.type, job_id=job_id)
    except Exception as ex:
        # Best-effort, as in the legacy delete path: a cluster problem must not leave a
        # caller unable to stop their own job as far as this API is concerned.
        log.error(f'Failed to delete Kubernetes job for {job_id}: {ex}')

    job.phase = JobStatus.CANCELED
    job.time_end = int(time.time())
    db.add(job)
    await db.commit()
    await db.refresh(job)

    return _job_document(job, _base_url(request))
