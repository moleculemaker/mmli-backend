import csv
import io
import json
from typing import Any, Optional

from fastapi import HTTPException
from sqlmodel.ext.asyncio.session import AsyncSession

from config import get_logger
from services.minio_service import MinIOService

log = get_logger(__name__)


class EzSpecificityService:
    """Post-process ez-specificity job output into the shape the frontend expects.

    An ez-specificity *parent* job (bucket ``ez-specificity``) orchestrates two
    subjobs whose outputs are collected under ``{job_id}/out/``:

      * unidock   -> ``{job_id}/out/{unidock_id}/out/results.json`` (+ ``complexes/``)
      * inference -> ``{job_id}/out/{inference_id}/out/results.csv``

    The frontend result table consumes a flat JSON array, one row per docked
    (enzyme, substrate) pair, with ``enzyme``, ``substrate``, ``smiles``,
    ``complex`` (a directly-fetchable URL to the complex PDB) and ``ez_score``
    (the inference score). ``resultPostProcess`` joins the two subjob outputs
    into that shape; ``unidockResultPostProcess`` / ``inferenceResultPostProcess``
    expose the individual subjob outputs (bucket ``ezspec-unidock`` /
    ``ezspec-inference``) for debugging / reuse.
    """

    RESULTS_JSON_SUFFIX = "/results.json"
    RESULTS_CSV_SUFFIX = "/results.csv"

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession) -> list[dict[str, Any]]:
        object_names = EzSpecificityService._list_output_objects(bucket_name, job_id, service)

        unidock_results_path = EzSpecificityService._find_suffix(object_names, EzSpecificityService.RESULTS_JSON_SUFFIX)
        if unidock_results_path is None:
            raise HTTPException(status_code=404, detail=f"Docking results not found for job {job_id}")

        rows = EzSpecificityService._docking_rows_with_complexes(bucket_name, unidock_results_path, service)

        inference_results_path = EzSpecificityService._find_suffix(object_names, EzSpecificityService.RESULTS_CSV_SUFFIX)
        scores = EzSpecificityService._read_inference_scores(bucket_name, inference_results_path, service)

        for row in rows:
            key = EzSpecificityService._pair_key(row.get("enzyme_index"), row.get("substrate_index"))
            row["ez_score"] = scores.get(key) if key is not None else None
        return rows

    @staticmethod
    async def unidockResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession) -> list[dict[str, Any]]:
        """Docking-only output for a single ezspec-unidock subjob ({job_id}/out/results.json)."""
        object_names = EzSpecificityService._list_output_objects(bucket_name, job_id, service)
        results_path = EzSpecificityService._find_suffix(object_names, EzSpecificityService.RESULTS_JSON_SUFFIX)
        if results_path is None:
            raise HTTPException(status_code=404, detail=f"Docking results not found for job {job_id}")
        return EzSpecificityService._docking_rows_with_complexes(bucket_name, results_path, service)

    @staticmethod
    async def inferenceResultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession) -> list[dict[str, Any]]:
        """Score rows for a single ezspec-inference subjob ({job_id}/out/results.csv)."""
        object_names = EzSpecificityService._list_output_objects(bucket_name, job_id, service)
        results_path = EzSpecificityService._find_suffix(object_names, EzSpecificityService.RESULTS_CSV_SUFFIX)
        if results_path is None:
            raise HTTPException(status_code=404, detail=f"Inference results not found for job {job_id}")
        content = service.get_file(bucket_name, results_path)
        if content is None:
            raise HTTPException(status_code=404, detail=f"Unable to read inference results: {results_path}")
        return list(csv.DictReader(io.StringIO(content.decode("utf-8"))))

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _docking_rows_with_complexes(bucket_name: str, results_json_path: str, service: MinIOService) -> list[dict[str, Any]]:
        """Read a docking results.json and turn each ``complex`` filename into a URL.

        Complexes live alongside results.json: ``{.../out}/complexes/<complex>``.
        """
        docking_rows = EzSpecificityService._read_docking_results(bucket_name, results_json_path, service)
        complexes_prefix = results_json_path[: -len(EzSpecificityService.RESULTS_JSON_SUFFIX)]

        rows: list[dict[str, Any]] = []
        for entry in docking_rows:
            complex_name = entry.get("complex")
            complex_url = None
            if complex_name:
                complex_url = service.get_file_url(bucket_name, f"{complexes_prefix}/complexes/{complex_name}")
            rows.append({**entry, "complex": complex_url})
        return rows

    @staticmethod
    def _list_output_objects(bucket_name: str, job_id: str, service: MinIOService) -> list[str]:
        objects = service.list_files(bucket_name, f"{job_id}/out/", recursive=True)
        if objects is None:
            return []
        return [obj.object_name for obj in objects]

    @staticmethod
    def _find_suffix(object_names: list[str], suffix: str) -> Optional[str]:
        # Prefer the shallowest match so we don't pick up nested/duplicated files.
        matches = sorted((name for name in object_names if name.endswith(suffix)), key=lambda n: n.count("/"))
        return matches[0] if matches else None

    @staticmethod
    def _read_docking_results(bucket_name: str, object_name: str, service: MinIOService) -> list[dict[str, Any]]:
        content = service.get_file(bucket_name, object_name)
        if content is None:
            raise HTTPException(status_code=404, detail=f"Unable to read docking results: {object_name}")
        try:
            data = json.loads(content.decode("utf-8"))
        except (ValueError, UnicodeDecodeError) as err:
            log.error(f"Failed to parse docking results {object_name}: {err}")
            raise HTTPException(status_code=500, detail="Malformed docking results")
        if not isinstance(data, list):
            raise HTTPException(status_code=500, detail="Unexpected docking results format")
        return data

    @staticmethod
    def _read_inference_scores(bucket_name: str, object_name: Optional[str], service: MinIOService) -> dict[tuple[int, int], float]:
        """Map (enzyme_index, substrate_index) -> inference score.

        We deliberately join on the enzyme/substrate indices rather than the CSV's
        "Dock Index" column: the inference container overwrites that column with an
        internal value, so it no longer matches the docking output.
        """
        if object_name is None:
            return {}
        content = service.get_file(bucket_name, object_name)
        if content is None:
            return {}

        scores: dict[tuple[int, int], float] = {}
        reader = csv.DictReader(io.StringIO(content.decode("utf-8")))
        for row in reader:
            key = EzSpecificityService._pair_key(row.get("Enzyme Index"), row.get("Substrate Index"))
            score = EzSpecificityService._parse_float(row.get("score"))
            if key is not None and score is not None:
                scores[key] = score
        return scores

    @staticmethod
    def _pair_key(enzyme_index: Any, substrate_index: Any) -> Optional[tuple[int, int]]:
        enzyme = EzSpecificityService._parse_int(enzyme_index)
        substrate = EzSpecificityService._parse_int(substrate_index)
        if enzyme is None or substrate is None:
            return None
        return (enzyme, substrate)

    @staticmethod
    def _parse_int(value: Any) -> Optional[int]:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _parse_float(value: Any) -> Optional[float]:
        try:
            return float(value)
        except (TypeError, ValueError):
            return None
