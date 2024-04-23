#!/bin/env python3
import os

from config import get_logger
from services.kubejob_service import download_remote_directory_from_minio

bucket_name = os.getenv('JOB_TYPE')
job_id = os.getenv('JOB_ID')
job_output_dir = os.getenv('JOB_OUTPUT_DIR')

target_directory = os.sep.join(job_output_dir.split(os.sep)[0:-2])

remote_path = os.getenv('JOB_INPUT_DIR')

log = get_logger(__name__)

try:
    download_remote_directory_from_minio(
        remote_path=remote_path,
        bucket_name=bucket_name,
        target_directory=target_directory,
    )
except Exception as ex:
    log.error(f'Failed to download input files from job[{bucket_name}]: {job_id} - {str(ex)}')

