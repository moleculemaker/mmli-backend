"""Backend service for CRISPR-COPIES jobs.

Two responsibilities, mirroring the CLEAN service pattern:
  * build_crispr_copies_job_command - translate the frontend form JSON into the
    `python main.py ...` command run inside the job container (called from
    routers/job.py at job-creation time).
  * resultPostProcess - parse the job's output.csv into typed JSON for the
    frontend (called from routers/files.py when results are fetched).

Frontend form source: src/app/components/crispr-copies/crispr-copies.ts
Backend script (in the job image): src/scripts/main.py
"""

import csv
import json
from io import StringIO
from os.path import basename
from typing import Any, Mapping

from fastapi import HTTPException
from sqlmodel.ext.asyncio.session import AsyncSession

from config import get_logger
from services.minio_service import MinIOService

log = get_logger(__name__)

# Files the frontend uploads land in MinIO under {job_id}/in/<filename> and are
# synced by the container's prejob step into ${JOB_INPUT_DIR}. Results written to
# ${JOB_OUTPUT_DIR} are synced back to {job_id}/out/. We emit the literal env-var
# references and let the in-container shell expand them (same convention as the
# SimpleFold and CLEAN commands).
JOB_INPUT_DIR = "${JOB_INPUT_DIR}"
JOB_OUTPUT_DIR = "${JOB_OUTPUT_DIR}"

OUTPUT_FILE_NAME = "output.csv"


def _int_or_str(value: str) -> Any:
    """Chromosome may be numeric (1, 2, ...) or a name (e.g. 'chrX')."""
    return int(value) if value.isdigit() else value


# Maps each output.csv header to (jsonKey, value-converter). Mirrors
# convert_crispr_copies_output_csv_to_json.py.
_COLUMN_MAP = {
    "Guide Sequence": ("guideSequence", str),
    "PAM": ("pam", str),
    "Accession": ("accession", str),
    "Self-Complementarity": ("selfComplementarity", int),
    "Chromosome": ("chromosome", _int_or_str),
    "Strand": ("strand", str),
    "Location": ("location", int),
    "Chromosome Length": ("chromosomeLength", int),
    "Intergenic Size": ("intergenicSize", int),
    "Left Gene": ("leftGene", str),
    "Right Gene": ("rightGene", str),
    "Relative Orientation": ("relativeOrientation", str),
    "Gene Density": ("geneDensity", float),
    "Left HR": ("leftHR", str),
    "Right HR": ("rightHR", str),
    "On-target Score": ("onTargetScore", float),
    "Zone": ("zone", str),
}


class CRISPRCopiesService:

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
        """Resolve an uploaded file to its in-container path.

        The frontend sends a FileMetadata-shaped object per uploaded file. We
        only need the original filename - the file itself has already been
        uploaded to MinIO at {job_id}/in/<filename> and synced into
        ${JOB_INPUT_DIR}. Returns None when no file was provided.
        """
        if not isinstance(metadata, Mapping):
            return None
        name = metadata.get("filename") or metadata.get("name") or metadata.get("url")
        if isinstance(name, str) and name:
            return f"{JOB_INPUT_DIR}/{basename(name)}"
        return None

    @staticmethod
    def build_crispr_copies_job_command(job_id, job_info, num_threads: int = 1) -> str:
        """Translate the crispr-copies form JSON into the job container command.

        Ports build_main_args(): maps the nested form into main.py CLI flags.
        Every value is double-quoted - this is required for values containing
        spaces (e.g. --on_target "Doench et al. 2016") and still allows the
        shell to expand ${JOB_INPUT_DIR}/${JOB_OUTPUT_DIR} inside the quotes.
        """
        _get = CRISPRCopiesService._get
        args: list[tuple[str, Any]] = []

        genome_path = CRISPRCopiesService._resolve_input_path(_get(job_info, "organism.genome"))
        if not genome_path:
            raise HTTPException(status_code=400, detail="organism.genome is required")
        args.append(("--Genome", genome_path))

        annotations_path = CRISPRCopiesService._resolve_input_path(_get(job_info, "organism.annotations"))
        if not annotations_path:
            raise HTTPException(status_code=400, detail="organism.annotations is required")
        args.append(("--Gene_table", annotations_path))

        use_cas9_overrides = _get(job_info, "parameters.crisprSystem") == "CRISPR/Cas9"
        if use_cas9_overrides:
            pam = _get(job_info, "parameters.cas9Specific.pam3")
            guide_length = _get(job_info, "parameters.cas9Specific.sgRNALengthWithoutPAM")
            on_target = _get(job_info, "parameters.cas9Specific.efficiencyScore")
            # The cas9Specific.efficiencyScore default is "Doench et al. 2016 -
            # only for NGG PAM"; main.py's --on_target only accepts the bare
            # model name (everything before the " - " descriptor).
            if isinstance(on_target, str) and " - " in on_target:
                on_target = on_target.split(" - ", 1)[0]
        else:
            pam = _get(job_info, "parameters.guideRNA.pam")
            guide_length = _get(job_info, "parameters.guideRNA.guideLength")
            on_target = _get(job_info, "parameters.guideRNA.onTargetScore")

        args.append(("--PAM", str(pam)))
        args.append(("--Orientation", str(_get(job_info, "parameters.guideRNA.pamOrientation"))))
        args.append(("--Guide_Length", str(int(guide_length))))
        args.append(("--Seed_Length", str(int(_get(job_info, "parameters.guideRNA.seedLength")))))
        args.append(("--edit_dist", str(int(_get(job_info, "parameters.guideRNA.editDistance")))))

        gc_min = int(_get(job_info, "parameters.guideRNA.gcContentMin"))
        gc_max = int(_get(job_info, "parameters.guideRNA.gcContentMax"))
        args.append(("--GC_grna", f"{gc_min},{gc_max}"))

        args.append(("--dist_type", str(_get(job_info, "parameters.guideRNA.distanceType"))))
        args.append(("--polyG_grna", str(int(_get(job_info, "parameters.guideRNA.polyG", 0)))))
        args.append(("--polyT_grna", str(int(_get(job_info, "parameters.guideRNA.polyT", 0)))))

        re_grna = _get(job_info, "parameters.guideRNA.restrictionEnzyme")
        if re_grna:
            args.append(("--RE_grna", str(re_grna)))

        backbone = _get(job_info, "parameters.guideRNA.backboneSequence")
        if backbone:
            args.append(("--backbone_complementarity_check", str(backbone)))

        args.append(("--HR_Length", str(int(_get(job_info, "parameters.homologyArm.lengthBp")))))
        args.append(("--polyG_hr", str(int(_get(job_info, "parameters.homologyArm.polyG", 0)))))
        args.append(("--polyT_hr", str(int(_get(job_info, "parameters.homologyArm.polyT", 0)))))

        re_hr = _get(job_info, "parameters.homologyArm.restrictionEnzyme")
        if re_hr:
            args.append(("--RE_hr", str(re_hr)))

        args.append(("--intspace", str(int(_get(job_info, "parameters.intergenicRegion.distance")))))
        args.append((
            "--gene_density_len",
            str(int(_get(job_info, "parameters.intergenicRegion.distanceToMeasureGeneDenisty"))),
        ))
        args.append((
            "--distal_end_len",
            str(int(_get(job_info, "parameters.intergenicRegion.distalEndOfChromosomeToAvoid"))),
        ))

        blast_org = _get(job_info, "parameters.intergenicRegion.referenceOrganismIdentifier")
        if blast_org:
            args.append(("--blast_org", str(blast_org)))

        args.append(("--on_target", str(on_target)))

        protein_path = CRISPRCopiesService._resolve_input_path(_get(job_info, "organism.protein"))
        if protein_path:
            args.append(("--protein_file", protein_path))

        args.append(("--Output_file", f"{JOB_OUTPUT_DIR}/{OUTPUT_FILE_NAME}"))
        args.append(("--num_threads", str(int(num_threads))))

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
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """Parse output.csv into a typed JSON array.

        Ports convert_crispr_copies_output_csv_to_json.py: camelCase keys, typed
        values, and a 1-based `id` per row.
        """
        result_path = f"{job_id}/out/{OUTPUT_FILE_NAME}"
        csv_file = service.get_file(bucket_name, result_path)
        if csv_file is None:
            log.error(f"CRISPR-COPIES result file not found: {result_path}")
            raise HTTPException(status_code=404, detail=f"404: Not Found - {result_path}")

        reader = csv.DictReader(StringIO(csv_file.decode("utf-8")))
        results = []
        try:
            for i, row in enumerate(reader):
                converted = {}
                for header, value in row.items():
                    if header not in _COLUMN_MAP:
                        raise ValueError(f"Unknown CSV header: {header}")
                    json_key, convert = _COLUMN_MAP[header]
                    converted[json_key] = convert(value)
                converted["id"] = i + 1
                results.append(converted)
        except ValueError as e:
            log.error(f"Failed to parse CRISPR-COPIES output for job {job_id}: {e}")
            raise HTTPException(status_code=500, detail=f"Failed to parse CRISPR-COPIES output: {e}")

        return results
