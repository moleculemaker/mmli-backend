#!/bin/env python3
"""Init-container step: convert an uploaded GFF3 annotation file to feature_table.txt.

Runs after prejob.py (input sync) and fetch_organism.py. If the user uploaded a
`.gff`/`.gff3` annotation, the method (`main.py`) and the genome-viewer service can't
read it — both expect NCBI `feature_table.txt`. This converts the GFF3 in $JOB_INPUT_DIR
to a fixed-name feature table that `build_crispr_copies_job_command` points `--Gene_table`
at (and that the genome-viewer service's `*_feature_table.txt` glob also finds).

No-op when no GFF is present (dropdown/NCBI jobs, or feature_table.txt uploads).

UNVERIFIED in CI/sandbox — exercised on a real upload job on the cluster.
"""

import os
import sys
from glob import glob

from config import get_logger
from services.gff3_to_feature_table import convert_gff3_file

log = get_logger(__name__)

CONVERTED_NAME = "converted_feature_table.txt"


def main():
    input_dir = os.getenv("JOB_INPUT_DIR")
    if not input_dir or not os.path.isdir(input_dir):
        return  # nothing synced (e.g. no uploads); nothing to convert

    gffs = sorted(glob(os.path.join(input_dir, "*.gff3")) + glob(os.path.join(input_dir, "*.gff")))
    if not gffs:
        return  # dropdown path or a feature_table.txt upload — no conversion needed

    src = gffs[0]
    dest = os.path.join(input_dir, CONVERTED_NAME)
    try:
        n = convert_gff3_file(src, dest)
        log.info(f"Converted GFF3 {os.path.basename(src)} -> {CONVERTED_NAME} ({n} gene/CDS rows)")
    except Exception as ex:
        log.error(f"Failed to convert GFF3 {src}: {ex}")
        sys.exit(1)


if __name__ == "__main__":
    main()
