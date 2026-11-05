"""Run a backend job helper with its own log messages copied to the job's log file.

The CRISPR-Copies setup step (fetch_organism.py, convert_annotations.py) runs in the
init container. Its log file, $JOB_OUTPUT_DIR/log, is uploaded to MinIO and served to
API clients, so only what these scripts choose to say may go into it.

Capturing the scripts' raw output instead would also capture config.py's import-time
log of app_config, which includes secret.yaml (MinIO keys, database password, hCaptcha
secret), and anything other modules log, such as the MinIO access key ID. So the file
handler goes on the script's own logger only, and records from other loggers never
reach it.
"""
import logging
import os
import sys

JOB_LOG_NAME = 'log'


def run_with_job_log(main, logger):
    """Run main(), also writing `logger`'s own records to $JOB_OUTPUT_DIR/log.

    An uncaught exception is logged here, with its traceback, and exits 1, so a crash
    leaves an explanation in the job log. SystemExit passes through unchanged.
    """
    output_dir = os.getenv('JOB_OUTPUT_DIR')
    if output_dir and os.path.isdir(output_dir):
        handler = logging.FileHandler(os.path.join(output_dir, JOB_LOG_NAME))
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)-8s %(message)s'))
        logger.addHandler(handler)
    try:
        main()
    except SystemExit:
        raise
    except Exception:
        logger.exception('Failed with an unexpected error')
        sys.exit(1)
