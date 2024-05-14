#!/bin/env python3
import os

from config import get_logger
from services.kubejob_service import upload_local_directory_to_minio

bucket_name = os.getenv('JOB_TYPE')
job_id = os.getenv('JOB_ID')
job_output_dir = os.getenv('JOB_OUTPUT_DIR')

log = get_logger(__name__)

try:
    log.info(f'Uploading to MinIO: {job_output_dir}')
    upload_local_directory_to_minio(local_path=job_output_dir, bucket_name=bucket_name)
    log.info(f'Uploaded successfully to MinIO: {job_output_dir}')
except Exception as ex:
    log.error(f'Failed to upload output files from job[{bucket_name}]: {job_id} - {str(ex)}')
