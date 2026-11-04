#!/bin/env python3
import os
import sys

from config import get_logger, MINIO_SERVER
from services.kubejob_service import upload_local_directory_to_minio

bucket_name = os.getenv('JOB_TYPE')
job_id = os.getenv('JOB_ID')
job_output_dir = os.getenv('JOB_OUTPUT_DIR')

log = get_logger(__name__)

# This is the Job pod's main container, so its exit code decides whether the Job
# succeeded. Exiting 0 after a failed upload marks the job 'completed' with no results.
try:
    log.info(f'Uploading to MinIO: {job_output_dir}')
    uploaded = upload_local_directory_to_minio(
        local_path=job_output_dir,
        bucket_name=bucket_name,
        minio_prefix=f"{job_id}/out")
except Exception as ex:
    log.error(f'Failed to upload output files from job[{bucket_name}]: {job_id} - {str(ex)}')
    sys.exit(1)

# upload_local_directory_to_minio returns False, without raising, when the output
# directory is missing
if uploaded is False:
    log.error(f'Failed to upload output files from job[{bucket_name}]: {job_id} - no output directory {job_output_dir}')
    sys.exit(1)

log.info(f'Uploaded successfully to MinIO: {job_output_dir}')
