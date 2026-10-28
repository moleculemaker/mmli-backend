import json
import uuid
from datetime import datetime, timezone
from pathlib import PureWindowsPath

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile, status
from fastapi.responses import JSONResponse
from pydantic import EmailStr

from config import get_logger
from services.minio_service import MinIOService

router = APIRouter()

log = get_logger(__name__)

# Contributions live in their own bucket, outside the job-type buckets that
# POST /{bucket_name}/upload is restricted to. Each one is stored as
# {CONTRIBUTIONS_BUCKET}/{contribution_id}/{filename} plus a submission.json beside it.
CONTRIBUTIONS_BUCKET = 'oed-contributions'
CONTRIBUTION_EXTENSIONS = ('.xlsx', '.xls', '.csv')
CONTRIBUTION_MAX_BYTES = 50 * 1024 * 1024
CONTRIBUTION_LICENSE = 'CC-BY-4.0'


@router.post(
    "/oed/contributions", tags=['OED'],
    summary="Submit a dataset to the Open Enzyme Database",
    description=(
        "Accept a data contribution from the OED About page. The file must be "
        f"{', '.join(CONTRIBUTION_EXTENSIONS)} and at most {CONTRIBUTION_MAX_BYTES // (1024 * 1024)} MB. "
        f"The contributor must agree to release it under {CONTRIBUTION_LICENSE} and confirm "
        "they have the rights to share it. Stored for review at MinIO "
        f"`{CONTRIBUTIONS_BUCKET}/{{contribution_id}}/`, with the contributor details in "
        "`submission.json`. Returns `{\"contributionID\": ..., \"submitted_at\": ...}`."
    ),
    status_code=status.HTTP_201_CREATED,
)
# Sync, so FastAPI runs it in its threadpool: the MinIO client blocks, and in an async
# handler a 50 MB put_object would stall every other request on the event loop.
def submit_contribution(
    file: UploadFile = File(...),
    email: EmailStr = Form(...),
    license_accepted: bool = Form(...),
    rights_confirmed: bool = Form(...),
    minio: MinIOService = Depends(),
):
    if not license_accepted:
        raise HTTPException(status_code=400, detail=f'Contributions must be released under {CONTRIBUTION_LICENSE}')
    if not rights_confirmed:
        raise HTTPException(status_code=400, detail='Contributor must confirm they have the rights to share this data')

    # Browsers may send a full client path (old IE/Edge send C:\...\name.xlsx); keep the
    # base name only, so the object key is always {contribution_id}/{name}.
    filename = PureWindowsPath(file.filename or '').name
    if not filename.lower().endswith(CONTRIBUTION_EXTENSIONS):
        raise HTTPException(status_code=400, detail=f'File must be one of: {", ".join(CONTRIBUTION_EXTENSIONS)}')

    # Read one byte past the limit so an oversized upload is detected without holding it all.
    content = file.file.read(CONTRIBUTION_MAX_BYTES + 1)
    if len(content) > CONTRIBUTION_MAX_BYTES:
        raise HTTPException(status_code=413, detail=f'File exceeds {CONTRIBUTION_MAX_BYTES // (1024 * 1024)} MB')
    if not content:
        raise HTTPException(status_code=400, detail='File is empty')

    contribution_id = uuid.uuid4().hex
    submitted_at = datetime.now(timezone.utc).isoformat()
    submission = {
        'contribution_id': contribution_id,
        'submitted_at': submitted_at,
        'email': email,
        'filename': filename,
        'size_bytes': len(content),
        'license': CONTRIBUTION_LICENSE,
        'rights_confirmed': True,
    }

    if (minio.ensure_bucket_exists(CONTRIBUTIONS_BUCKET)
            and minio.upload_file(CONTRIBUTIONS_BUCKET, f'{contribution_id}/{filename}', content)
            and minio.upload_file(CONTRIBUTIONS_BUCKET, f'{contribution_id}/submission.json', json.dumps(submission, indent=2))):
        return JSONResponse(content={'contributionID': contribution_id, 'submitted_at': submitted_at},
                            status_code=status.HTTP_201_CREATED)

    log.error(f'Failed to store contribution {contribution_id}: {filename}')
    return JSONResponse(content={'error': 'Unable to store contribution'}, status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)
