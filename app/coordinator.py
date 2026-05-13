#!/bin/env python3
import json
import os
import shutil
import sys
import traceback
from time import sleep

import requests
from fastapi import HTTPException
from requests import Response

from config import get_logger, app_config
from services.kubejob_service import download_remote_directory_from_minio, upload_local_directory_to_minio, \
    api_batch_v1, get_job_name_from_id

job_type = bucket_name = os.getenv('JOB_TYPE')
job_id = os.getenv('JOB_ID')

job_input_dir = os.getenv('JOB_INPUT_DIR')
job_output_dir = os.getenv('JOB_OUTPUT_DIR')
namespace = os.getenv('NAMESPACE')
target_directory = os.sep.join(job_input_dir.split(os.sep)[0:-2])

remote_path = os.path.join(job_id, 'in')

hostname = os.getenv('MMLI_BACKEND_HOST')

log = get_logger(__name__)


def create_job_step(step_type: str, subjob_id: str, job_config=None) -> str:
    req_body = {
        "job_id": subjob_id,
        # "run_id": step_type,
        # "email": "",
        "job_info": json.dumps(job_config)   # Any additional JSON inputs?
    } if job_config is not None else { "job_id": subjob_id }

    # Create a subjob that will run this step on the input files
    resp: Response = requests.post(
        f'{hostname}/{step_type}/jobs',
        data=json.dumps(req_body),
        headers={
            "Content-Type": "application/json",
        }
    )

    # Ensure that job was submitted successfully
    try:
        resp.raise_for_status()
        if resp.status_code == 200 or resp.status_code == 201:
            job_status = resp.json()
            step_id = job_status['job_id']

            log.debug(
                f'job_type={job_type}  job_id={job_id}  step_type={step_type}  step_id={subjob_id}  |  Job step created! status={str(job_status)}')
            log.info(f'job_type={job_type}  job_id={job_id}  step_type={step_type}  step_id={subjob_id}  |  Waiting for step to be completed...')

            # Return job_id / step_id
            return step_id
    except Exception as e:
        log.error(f'job_type={job_type}  job_id={job_id}  step={step_type}  step_id={subjob_id}  |  Failed to get status for step: {str(e)}')
        log.error(traceback.format_exc())
        sys.exit(resp.status_code)


def get_job_step_status(step_name: str, step_id: str) -> str:
    resp = requests.get(f'{hostname}/{step_name}/jobs/{step_id}')
    resp.raise_for_status()
    job_list = resp.json()
    log.debug(f'job_type={job_type}  job_id={job_id}  step_type={step_name}  step_id={step_id}  |  job_list={str(job_list)}')
    job_status = job_list[0] if len(job_list) > 0 else None
    log.debug(f'job_type={job_type}  job_id={job_id}  step={step_name}  step_id={step_id}  |  job_status={str(job_status)}')
    return job_status


try:
    scratch_dir = job_output_dir
    log.info(f'job_type={job_type}  job_id={job_id}  |  Running and coordinating job steps...')
    log.debug(f'job_type={job_type}  job_id={job_id}  |     job_input_dir={job_input_dir}')
    log.debug(f'job_type={job_type}  job_id={job_id}  |     job_output_dir={job_output_dir}')
    log.debug(f'job_type={job_type}  job_id={job_id}  |     scratch_dir={scratch_dir}')

    # See cfg/config.yaml (defaults) and chart/values.yaml (per-environment configs)
    job_defn = app_config['kubernetes_jobs'][job_type]
    keys = job_defn['steps'] if 'steps' in job_defn else None   # ['ezspec-docking', 'ezspec-inference']

    if keys is None:
        log.error(f'Failed to coordinate job_type={job_type}: unsupported job_type')
        sys.exit(-400)
    if len(keys) == 0:
        log.error(f'Failed to coordinate job_type={job_type}: job_type has no steps defined')
        sys.exit(-500)

    ezspec_unidock_job_id = os.getenv('EZSPEC_UNIDOCK_JOB_ID')
    ezspec_inference_job_id = os.getenv('EZSPEC_INFERENCE_JOB_ID')
    steps = dict(zip(keys, [
        ezspec_unidock_job_id,
        ezspec_inference_job_id
    ]))

    log.info(f'job_type={job_type}  job_id={job_id}  |  Running job steps: {str(steps)}')

    # To cut down on minio spam / size, we could do all of these as local "copy" from the NFS
    # But this would make the prejob/postjob unnecessary for these steps, and that's currently not optional
    # So we may want to weigh or options before investing too heavily in one solution or the other
    # This may seem wasteful, but should work for now and provide us with accounting/logging if stuff goes wrong

    # Our prejob should have already downloaded user inputs to job_input_dir
    # Copy these to our scratch folder, which is passed along to all subjobs
    shutil.copytree(job_input_dir, scratch_dir, dirs_exist_ok=True)

    # All subjobs can use the same job_info/job_config as the parent job
    job_config = {
        "parent_job_id": job_id,
        "ezspec_unidock_job_id": ezspec_unidock_job_id,
        "ezspec_inference_job_id": ezspec_inference_job_id,
    }

    # Our steps (subjobs) are now ready to run!

    for subjob_type, subjob_id in steps.items():
        log.info(f'job_type={job_type}  job_id={job_id}  step={subjob_type}  step_id={subjob_id}  |  Running job step: {subjob_id}')

        # Upload our local scratch directory to MinIO for processing
        upload_local_directory_to_minio(
            local_path=scratch_dir,
            bucket_name=subjob_type,
            minio_prefix=f'{subjob_id}/in'
        )

        # Create our Job step
        step_id = create_job_step(step_type=subjob_type, subjob_id=subjob_id, job_config=job_config)
        if step_id is None:
            log.error(f'job_type={job_type}  job_id={job_id}  step={subjob_type}  step_id={subjob_id}  |  FATAL: step_id=None. Aborting...')
            log.error(traceback.format_exc())
            sys.exit(-600)

        if subjob_id != step_id:
            log.error(
                f'job_type={job_type}  job_id={job_id}  step={subjob_type}  step_id={subjob_id}  |  ERROR: Step / SubJob ID mismatch: {str(e)}')
            log.error(traceback.format_exc())
            sys.exit(-700)

        # Wait for job status to be Completed
        try:
            log.debug(f'Polling for job status to be "completed"...')
            job_status = get_job_step_status(step_name=subjob_type, step_id=step_id)
            while job_status['phase'] != 'completed' and job_status['phase'] != 'error':
                job_status = get_job_step_status(step_name=subjob_type, step_id=step_id)
                log.debug(f'Waiting 10 seconds before polling again...')
                sleep(10)
        except Exception as e:
            log.error(f'job_type={job_type}  job_id={job_id}  step={subjob_type}  step_id={subjob_id}  |  Failed to get status for step: {str(e)}')
            log.error(traceback.format_exc())
            sys.exit(-400)

        # Ensure no error in job execution
        if job_status['phase'] == 'error':
            log.error(f'job_type={job_type}  job_id={job_id}  step={subjob_type}  step_id={subjob_id}  |  Job step failed with error: {str(job_status)}')
            api_response = api_batch_v1.read_namespaced_job(
                namespace=namespace,
                name=get_job_name_from_id(job_type=job_type, job_id=subjob_id),
            )
            sys.exit(api_response['status'])

        log.info(f'job_type={job_type}  job_id={job_id}  step={subjob_type}  step_id={subjob_id}  |  Job step completed successfully!')
        log.info(f'job_type={job_type}  job_id={job_id}  step={subjob_type}  step_id={subjob_id}  |  Collecting step outputs...')

        # Output files should now be present in MinIO: e.g. bucket=step  path={job_id}/out
        # Download them to our local scratch directory
        download_remote_directory_from_minio(
            remote_path=f'{subjob_id}/out',
            bucket_name=subjob_type,
            target_directory=scratch_dir
        )

        log.info(f'job_type={job_type}  job_id={job_id}  step={subjob_type}  step_id={subjob_id}  |  Job step finished!')
        continue

    # All per-step output files should now be present in our local scratch directory + MinIO
    # upload_local_directory_to_minio(
    #     local_path=target_directory,
    #     bucket_name=job_type,
    #     minio_prefix=f'{job_id}/out',
    # )
    # If this is too noisy, we can copy only those output files that the frontend needs
    # shutil.copy(f'{scratch_dir}/example.out', f'{job_output_dir}/')
    log.info(f'job_type={job_type}  job_id={job_id}  |  All job steps completed!')

    # After the job completes, our postjob should now upload all files
    # from job_output_dir to our parent job's output folder in MinIO

    sys.exit(0)
except Exception as ex:
    log.error(f'Failed to coordinate job_type={job_type}: job_id={job_id} - {str(ex)}')
    traceback.print_exc()
    sys.exit(-500)

