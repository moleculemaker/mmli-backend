"""Backend service for mutagenesis (primer design) jobs.

Two responsibilities, mirroring the CRISPR-COPIES / CLEAN service pattern:
  * build_mutagenesis_job_command - translate the frontend form JSON into the
    `python main.py ...` command run inside the job container (called from
    routers/job.py at job-creation time).
  * resultPostProcess - parse the job's three output CSVs into typed JSON for
    the frontend (called from routers/files.py when results are fetched).

Backend script (in the job image): main.py
  https://github.com/Zhao-Group/Primer_Design_and_Worklists (h_to_args_part2)
"""

import csv
import json
from io import StringIO
from os.path import basename
from typing import Any

from fastapi import HTTPException
from pydantic import ValidationError
from sqlmodel.ext.asyncio.session import AsyncSession

from config import get_logger
from models.mutagenesis_params import MutagenesisJobInfo
from services.minio_service import MinIOService

log = get_logger(__name__)

# Uploaded inputs are synced into ${JOB_INPUT_DIR}; results written to
# ${JOB_OUTPUT_DIR} are synced back to {job_id}/out/. We emit the literal env-var
# references and let the in-container shell expand them.
JOB_INPUT_DIR = "${JOB_INPUT_DIR}"
JOB_OUTPUT_DIR = "${JOB_OUTPUT_DIR}"

# Single combined results file we pass to main.py via -f, so resultPostProcess can
# read a predictable name regardless of main.py defaults.
RESULTS_FILE = "results.csv"


def _num(value: str) -> Any:
    """Coerce a numeric CSV value to int when whole, else float."""
    return float(value) if "." in value else int(value)


def _json_list(value: str) -> list:
    """Decode the JSON-encoded Alerts cell into a list; tolerant of bad/empty cells."""
    value = (value or "").strip()
    if not value:
        return []
    try:
        decoded = json.loads(value)
        return decoded if isinstance(decoded, list) else [decoded]
    except (ValueError, TypeError):
        return [value]


# Maps each results-CSV header to (jsonKey, value-converter). Mirrors the combined
# CSV emitted by build_results() in the method (Mutations..Well).
_RESULT_COLUMN_MAP = {
    "Mutations": ("mutations", str),
    "Assembly_Fragments": ("assemblyFragments", int),
    "Direction": ("direction", str),
    "Sequence": ("sequence", str),
    "Tm": ("tm", _num),
    "GC": ("gc", _num),
    "Length": ("length", int),
    "Alerts": ("alerts", _json_list),
    "Plate": ("plate", str),
    "Well": ("well", str),
}


class MutagenesisService:

    @staticmethod
    def build_mutagenesis_job_command(job_id, job_info) -> str:
        """Translate the mutagenesis form JSON into the job container command.

        Inputs (all uploaded files, materialized client-side; see spec Q1):
          orfFile       - ORF DNA sequence (.txt)                        (-orf)
          mutationList  - CSV with a `Mutations` header                  (-m)
          leftOverhang  - optional upstream flank (.txt)                 (-left)
          rightOverhang - optional downstream flank (.txt)               (-right)
          codonTableValue - NCBI codon translation table id (int)        (-c)
          tmMethod      - melting-temperature method                     (-tm)

        The codon-statistics file (-cod) is bundled in the job image, so we leave
        main.py's default (nuclear_codon_statistics.tsv) and omit the flag.
        """
        # Validate the (opaque) job_info against the typed contract → clean 422
        # rather than a broken `python main.py` invocation in the container.
        try:
            params = MutagenesisJobInfo.parse_obj(job_info)
        except ValidationError as e:
            raise HTTPException(status_code=422, detail=e.errors())

        args: list[tuple[str, Any]] = []
        args.append(("-orf", f"{JOB_INPUT_DIR}/{basename(params.orfFile.resolved_name())}"))
        args.append(("-m", f"{JOB_INPUT_DIR}/{basename(params.mutationList.resolved_name())}"))

        if params.leftOverhang and params.leftOverhang.resolved_name():
            args.append(("-left", f"{JOB_INPUT_DIR}/{basename(params.leftOverhang.resolved_name())}"))
        if params.rightOverhang and params.rightOverhang.resolved_name():
            args.append(("-right", f"{JOB_INPUT_DIR}/{basename(params.rightOverhang.resolved_name())}"))

        args.append(("-o", JOB_OUTPUT_DIR))
        args.append(("-f", RESULTS_FILE))
        args.append(("-c", str(params.codonTableValue)))
        args.append(("-tm", params.tmMethod))

        rendered = " ".join(f'{flag} "{value}"' for flag, value in args)
        inner = f"python main.py {rendered}"
        # `set -o pipefail` is REQUIRED: the command is piped to `tee`, so without it
        # the pipeline exit status is tee's (0) and a failed method reads as success
        # with no output (handoff Finding B). tee logs, list outputs on success, drop
        # an `error` sentinel + fail the container on any error.
        command = (
            "set -o pipefail; "
            f'((({inner} 2>&1 | tee {JOB_OUTPUT_DIR}/log) && ls -al {JOB_OUTPUT_DIR}/)'
            f' || (touch {JOB_OUTPUT_DIR}/error && false))'
        )
        return command

    @staticmethod
    def _parse_results_csv(bucket_name: str, job_id: str, file_name: str, service: MinIOService) -> list:
        result_path = f"{job_id}/out/{file_name}"
        csv_file = service.get_file(bucket_name, result_path)
        if csv_file is None:
            log.error(f"Mutagenesis result file not found: {result_path}")
            raise HTTPException(status_code=404, detail=f"404: Not Found - {result_path}")

        reader = csv.DictReader(StringIO(csv_file.decode("utf-8")))
        rows = []
        # Tolerant by design (handoff): skip unknown columns, null-out un-convertible
        # values, rather than raising and 500-ing the whole fetch.
        for i, row in enumerate(reader):
            converted = {"id": i + 1}
            for header, value in row.items():
                if header not in _RESULT_COLUMN_MAP:
                    continue
                json_key, convert = _RESULT_COLUMN_MAP[header]
                try:
                    converted[json_key] = convert(value)
                except (ValueError, TypeError):
                    log.warning(f"mutagenesis {job_id}: could not convert {header}={value!r}, nulling")
                    converted[json_key] = None
            rows.append(converted)
        return rows

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """Parse the single combined results CSV into the typed `primers` array."""
        return {
            "primers": MutagenesisService._parse_results_csv(
                bucket_name, job_id, RESULTS_FILE, service),
        }
