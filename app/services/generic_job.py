"""The generic job contract: tools that need no per-tool input or results code.

(Notification email still dispatches per tool in kubejob_service; mmli-backend#145.)

A tool opts in by declaring an `output_schema` in tools.yaml. Then:

  submission  the input, already validated against the tool's input schema, is written
              to MinIO as {job_id}/in/job.json alongside any uploaded files, and the
              image runs with no command, so its ENTRYPOINT reads job.json from
              $JOB_INPUT_DIR
  results     the container writes results.json to $JOB_OUTPUT_DIR; it is served as-is
              after validation against the output schema. Larger outputs are separate
              files, listed in results.json under `artifacts` and served by name from
              /v1/jobs/{id}/artifacts/{name}
  routes      /v1 only. The legacy create and upload routes take a client-chosen job_id,
              which would let anyone holding an id overwrite job.json, so both refuse
              these tools

moleculemaker/job-template implements the container side of this contract.
"""
import json
from typing import Any, List, Optional

from fastapi import HTTPException
from jsonschema import Draft202012Validator

from services import tool_registry
from services.minio_service import MinIOService

JOB_INPUT_FILENAME = 'job.json'
RESULTS_FILENAME = 'results.json'


class InvalidResults(Exception):
    """results.json exists but is not JSON or does not match the tool's output schema."""

    def __init__(self, errors: List[str]):
        super().__init__('; '.join(errors))
        self.errors = errors


def write_job_input(service: MinIOService, job_type: str, job_id: str, job_info: str) -> None:
    service.ensure_bucket_exists(job_type)
    if not service.upload_file(job_type, f'{job_id}/in/{JOB_INPUT_FILENAME}', job_info.encode('utf-8')):
        raise HTTPException(status_code=400, detail='Failed to upload file to MinIO')


def read_results(service: MinIOService, job_type: str, job_id: str) -> Optional[Any]:
    """Return the job's results.json, or None if the job wrote none.

    Raises InvalidResults when the file breaks the tool's published output schema, so a
    container that drifted from its schema fails loudly here rather than in a browser.
    """
    content = service.get_file(job_type, f'{job_id}/out/{RESULTS_FILENAME}')
    if content is None:
        return None

    try:
        results = json.loads(content)
    except ValueError as ex:
        raise InvalidResults([f'{RESULTS_FILENAME} is not valid JSON: {ex}'])

    schema = tool_registry.get_output_schema(job_type)
    errors = [
        f'{error.json_path}: {error.message}'
        for error in Draft202012Validator(schema).iter_errors(results)
    ]
    if errors:
        raise InvalidResults(errors)
    return results
