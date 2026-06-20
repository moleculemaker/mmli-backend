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
from pydantic import ValidationError
from sqlmodel.ext.asyncio.session import AsyncSession

from config import get_logger
from models.crispr_copies_params import CrisprCopiesJobInfo
from services.minio_service import MinIOService
# Fixed filenames fetch_organism.py writes into $JOB_INPUT_DIR for the dropdown path.
# Imported (not redefined) so the builder and the fetcher can never drift apart.
from fetch_organism import (
    GENOME_NAME as FETCHED_GENOME_NAME,
    FEATURE_TABLE_NAME as FETCHED_FEATURE_TABLE_NAME,
    PROTEIN_NAME as FETCHED_PROTEIN_NAME,
)

log = get_logger(__name__)

# Files the frontend uploads land in MinIO under {job_id}/in/<filename> and are
# synced by the container's prejob step into ${JOB_INPUT_DIR}. Results written to
# ${JOB_OUTPUT_DIR} are synced back to {job_id}/out/. We emit the literal env-var
# references and let the in-container shell expand them (same convention as the
# SimpleFold and CLEAN commands).
JOB_INPUT_DIR = "${JOB_INPUT_DIR}"
JOB_OUTPUT_DIR = "${JOB_OUTPUT_DIR}"

OUTPUT_FILE_NAME = "output.csv"


# Maps each output.csv header to (jsonKey, value-converter). Mirrors
# convert_crispr_copies_output_csv_to_json.py.
_COLUMN_MAP = {
    "ID": ("id", int),                                  # method's 1-based row id (replaces the old auto-generated id)
    "Guide Sequence": ("guideSequence", str),
    "PAM": ("pam", str),
    "Accession": ("accession", str),
    "GC Content": ("gcContent", float),
    "Self-Complementarity": ("selfComplementarity", int),
    "Chromosome": ("chromosome", str),
    "Strand": ("strand", str),
    "Location": ("location", int),
    "Chromosome Length": ("chromosomeLength", int),
    "Distance, Closest off-target": ("closestOffTargetDistance", float),  # column NAME contains a comma
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
        # Validate the parameter contract up front so a malformed payload is rejected
        # with a clean 422 here, rather than crashing or silently misbehaving inside the
        # `python main.py` container (see models/crispr_copies_params.py).
        try:
            CrisprCopiesJobInfo.parse_obj(job_info)
        except ValidationError as e:
            raise HTTPException(
                status_code=422,
                detail={"message": "Invalid CRISPR-COPIES parameters", "errors": e.errors()},
            ) from e

        _get = CRISPRCopiesService._get
        args: list[tuple[str, Any]] = []

        # Organism inputs arrive by one of two mutually-exclusive paths (enforced by
        # CrisprCopiesJobInfo): a dropdown accession, which the init container's
        # fetch_organism.py resolves into $JOB_INPUT_DIR under the FETCHED_* names, or
        # uploaded files referenced by their MinIO-synced filenames.
        organism_id = _get(job_info, "organism.organismIdentifier.value")
        uploaded_genome = CRISPRCopiesService._resolve_input_path(_get(job_info, "organism.genome"))
        use_dropdown = bool(organism_id) and not uploaded_genome

        if use_dropdown:
            genome_path = f"{JOB_INPUT_DIR}/{FETCHED_GENOME_NAME}"
            annotations_path = f"{JOB_INPUT_DIR}/{FETCHED_FEATURE_TABLE_NAME}"
        else:
            genome_path = uploaded_genome
            if not genome_path:
                raise HTTPException(status_code=400, detail="organism.genome is required")
            annotations_path = CRISPRCopiesService._resolve_input_path(_get(job_info, "organism.annotations"))
            if not annotations_path:
                raise HTTPException(status_code=400, detail="organism.annotations is required")
        args.append(("--Genome", genome_path))
        args.append(("--Gene_table", annotations_path))

        # Single guide-parameter set for every CRISPR system. The PAM/length/score
        # appropriate to the selected system are set on guideRNA by the frontend's
        # crisprSystem coupling (the old parallel cas9Specific group has been folded away).
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
            str(int(_get(job_info, "parameters.intergenicRegion.distanceToMeasureGeneDensity"))),
        ))
        args.append((
            "--distal_end_len",
            str(int(_get(job_info, "parameters.intergenicRegion.distalEndOfChromosomeToAvoid"))),
        ))

        blast_org = _get(job_info, "parameters.intergenicRegion.referenceOrganismIdentifier")
        if blast_org:
            args.append(("--blast_org", str(blast_org)))

        args.append(("--on_target", str(on_target)))

        # Protein FASTA is optional. For uploads we know at build time whether one exists;
        # for the dropdown path NCBI may or may not provide one, so the flag is appended at
        # runtime only if fetch_organism.py actually produced a non-empty protein file.
        protein_runtime_suffix = ""
        if use_dropdown:
            fetched_protein = f"{JOB_INPUT_DIR}/{FETCHED_PROTEIN_NAME}"
            protein_runtime_suffix = (
                f' $( [ -s "{fetched_protein}" ] && printf -- \'--protein_file "%s"\' "{fetched_protein}" )'
            )
        else:
            protein_path = CRISPRCopiesService._resolve_input_path(_get(job_info, "organism.protein"))
            if protein_path:
                args.append(("--protein_file", protein_path))

        args.append(("--Output_file", f"{JOB_OUTPUT_DIR}/{OUTPUT_FILE_NAME}"))
        args.append(("--num_threads", str(int(num_threads))))

        rendered = " ".join(f'{flag} "{value}"' for flag, value in args)
        inner = f"python main.py {rendered}{protein_runtime_suffix}"
        # Mirror the CLEAN command wrapper: tee logs, list outputs on success,
        # and drop an `error` sentinel + fail the container on any error.
        # `set -o pipefail` is REQUIRED: without it the pipeline's exit status is
        # `tee`'s (always 0), so a non-zero `python main.py` (crash, argparse
        # rejection, or sys.exit(1)) would be masked, the `error` sentinel skipped,
        # and the job falsely reported Done with no output.csv. The job template runs
        # this under `bash -c` in an `&&` chain, so we join with `&&` (not `;`) to
        # keep that chain intact.
        command = (
            f'set -o pipefail && '
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
        for row in reader:
            converted = {}
            for header, value in row.items():
                mapping = _COLUMN_MAP.get(header)
                if mapping is None:
                    # Tolerate columns the method emits that we don't surface (forward-proof
                    # against new main.py output columns) rather than 500 the whole fetch.
                    log.warning(f"Unmapped CRISPR-COPIES output column '{header}' for job {job_id}; skipping it.")
                    continue
                json_key, convert = mapping
                try:
                    converted[json_key] = convert(value)
                except (ValueError, TypeError):
                    # e.g. 'NA' On-target Score when scoring is unavailable; null rather than crash.
                    converted[json_key] = None
            results.append(converted)

        return results
