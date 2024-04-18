from config import get_logger

log = get_logger(__name__)


def build_molli_job_environment(files, job_id):
    environment = []

    # Write user input files to the shared folder
    input_dir_in_apiserver = '/uws/job/input'
    for name in ['cores', 'subs']:
        # Read uploaded file
        fileinfo = files[name][0]
        # log.debug(f"File info: {fileinfo}")

        # Write file to dest_folder
        filepath = f'{input_dir_in_apiserver}/{job_id}.{name}.cdxml'
        log.debug(f'[molli:{job_id}]  Writing {name} file to: {filepath}')
        fh = open(filepath, 'wb')
        fh.write(fileinfo['body'])
        fh.close()

    # Tell the job to use these files when it runs
    input_dir_in_container = '/app/data/inputs'
    environment.append({'name': 'CORES_INPUT_FILE', 'value': f'{input_dir_in_container}/{job_id}.cores.cdxml'})
    environment.append({'name': 'SUBS_INPUT_FILE', 'value': f'{input_dir_in_container}/{job_id}.subs.cdxml'})

    return environment
