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
from io import StringIO
from os.path import basename
from typing import Any, Mapping

from fastapi import HTTPException
from sqlmodel.ext.asyncio.session import AsyncSession

from config import get_logger
from services.minio_service import MinIOService

log = get_logger(__name__)

# Uploaded inputs are synced into ${JOB_INPUT_DIR}; results written to
# ${JOB_OUTPUT_DIR} are synced back to {job_id}/out/. We emit the literal env-var
# references and let the in-container shell expand them.
JOB_INPUT_DIR = "${JOB_INPUT_DIR}"
JOB_OUTPUT_DIR = "${JOB_OUTPUT_DIR}"

# Canonical output filenames we pass to main.py via -f/-fwd/-rev, so
# resultPostProcess can read predictable names regardless of main.py defaults.
DESIGNED_PRIMERS_FILE = "designed_primers.csv"
FORWARD_PRIMERS_FILE = "forward_primers.csv"
REVERSE_PRIMERS_FILE = "reverse_primers.csv"


def _num(value: str) -> Any:
    """Coerce a numeric CSV value to int when whole, else float."""
    return float(value) if "." in value else int(value)


# Maps each primer-CSV header to (jsonKey, value-converter). The designed-primers
# file has the first five columns; forward/reverse files add "Well".
_PRIMER_COLUMN_MAP = {
    "Name": ("name", str),
    "Sequence": ("sequence", str),
    "Tm": ("tm", _num),
    "GC": ("gc", _num),
    "Length": ("length", int),
    "Well": ("well", str),
}


class MutagenesisService:

    @staticmethod
    def _get(form: Mapping[str, Any], path: str, default: Any = None) -> Any:
        cur: Any = form
        for key in path.split("."):
            if not isinstance(cur, Mapping) or key not in cur:
                return default
            cur = cur[key]
        return cur

    @staticmethod
    def _resolve_input_path(metadata: Any) -> str | None:
        """Resolve an uploaded file to its in-container ${JOB_INPUT_DIR} path."""
        if not isinstance(metadata, Mapping):
            return None
        name = metadata.get("filename") or metadata.get("name") or metadata.get("url")
        if isinstance(name, str) and name:
            return f"{JOB_INPUT_DIR}/{basename(name)}"
        return None

    @staticmethod
    def build_mutagenesis_job_command(job_id, job_info) -> str:
        """Translate the mutagenesis form JSON into the job container command.

        Inputs:
          mutationList  - uploaded CSV mapping mutations to well positions (-m)
          orfFile       - uploaded TXT of the ORF DNA sequence            (-orf)
          codonTableValue - NCBI codon translation table id (int)         (-c)
          tmMethod      - melting-temperature method                      (-tm)

        The codon-statistics file (-cod) is bundled in the job image, so we leave
        main.py's default (nuclear_codon_statistics.tsv) and omit the flag.
        """
        _get = MutagenesisService._get
        args: list[tuple[str, Any]] = []

        mutation_path = MutagenesisService._resolve_input_path(_get(job_info, "mutationList"))
        if not mutation_path:
            raise HTTPException(status_code=400, detail="mutationList is required")
        args.append(("-m", mutation_path))

        orf_path = MutagenesisService._resolve_input_path(_get(job_info, "orfFile"))
        if not orf_path:
            raise HTTPException(status_code=400, detail="orfFile is required")
        args.append(("-orf", orf_path))

        args.append(("-o", JOB_OUTPUT_DIR))
        args.append(("-f", DESIGNED_PRIMERS_FILE))
        args.append(("-fwd", FORWARD_PRIMERS_FILE))
        args.append(("-rev", REVERSE_PRIMERS_FILE))
        args.append(("-c", str(int(_get(job_info, "codonTableValue", 1)))))
        args.append(("-tm", str(_get(job_info, "tmMethod", "SantaLucia"))))

        rendered = " ".join(f'{flag} "{value}"' for flag, value in args)
        inner = f"python main.py {rendered}"
        # Mirror the CLEAN command wrapper: tee logs, list outputs on success,
        # and drop an `error` sentinel + fail the container on any error.
        command = (
            f'((({inner} 2>&1 | tee {JOB_OUTPUT_DIR}/log) && ls -al {JOB_OUTPUT_DIR}/)'
            f' || (touch {JOB_OUTPUT_DIR}/error && false))'
        )
        return command

    @staticmethod
    def _parse_primer_csv(bucket_name: str, job_id: str, file_name: str, service: MinIOService) -> list:
        result_path = f"{job_id}/out/{file_name}"
        csv_file = service.get_file(bucket_name, result_path)
        if csv_file is None:
            log.error(f"Mutagenesis result file not found: {result_path}")
            raise HTTPException(status_code=404, detail=f"404: Not Found - {result_path}")

        reader = csv.DictReader(StringIO(csv_file.decode("utf-8")))
        rows = []
        try:
            for i, row in enumerate(reader):
                converted = {}
                for header, value in row.items():
                    if header not in _PRIMER_COLUMN_MAP:
                        raise ValueError(f"Unknown CSV header: {header}")
                    json_key, convert = _PRIMER_COLUMN_MAP[header]
                    converted[json_key] = convert(value)
                converted["id"] = i + 1
                rows.append(converted)
        except ValueError as e:
            log.error(f"Failed to parse mutagenesis output {file_name} for job {job_id}: {e}")
            raise HTTPException(status_code=500, detail=f"Failed to parse mutagenesis output: {e}")
        return rows

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """Parse the three primer CSVs into typed JSON keyed by primer set."""
        return {
            "designedPrimers": MutagenesisService._parse_primer_csv(
                bucket_name, job_id, DESIGNED_PRIMERS_FILE, service),
            "forwardPrimers": MutagenesisService._parse_primer_csv(
                bucket_name, job_id, FORWARD_PRIMERS_FILE, service),
            "reversePrimers": MutagenesisService._parse_primer_csv(
                bucket_name, job_id, REVERSE_PRIMERS_FILE, service),
        }
