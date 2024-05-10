from config import get_logger

log = get_logger(__name__)


def build_molli_job_environment(job_info, job_id):
    environment = []

    input_dir = f'/uws/jobs/molli/{job_id}/in'
    CORES_FILE_NAME = job_info['CORES_FILE_NAME']
    SUBS_FILE_NAME = job_info['SUBS_FILE_NAME']

    # Tell the job to use these files when it runs
    environment.append({'name': 'CORES_INPUT_FILE', 'value': f'{input_dir}/{CORES_FILE_NAME}'})
    environment.append({'name': 'SUBS_INPUT_FILE', 'value': f'{input_dir}/{SUBS_FILE_NAME}'})

    return environment
