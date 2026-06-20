#!/bin/env python3
"""Init-container step: resolve a CRISPR-COPIES dropdown organism to input files.

Runs after prejob.py (which syncs uploaded inputs from MinIO). If the job used the
organism dropdown, the frontend sends an NCBI RefSeq assembly accession; the backend
passes it in as $ORGANISM_ACCESSION. This script lands three files in $JOB_INPUT_DIR
under fixed names that build_crispr_copies_job_command points main.py at:

    organism_genome.fna          (--Genome)         required
    organism_feature_table.txt   (--Gene_table)     required
    organism_protein.faa         (--protein_file)   optional (not all assemblies have one)

Fetch strategy: genome + protein via the NCBI `datasets` CLI; the RefSeq
`*_feature_table.txt` (the exact format main.py parses) is NOT a `datasets` include, so
it is pulled from the NCBI FTP path derived from the accession. Results are cached in
MinIO under `_organism_cache/<accession>/` so repeat jobs skip the network fetch.

If $ORGANISM_ACCESSION is unset/empty (the upload path), this is a no-op.
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from glob import glob

from config import get_logger
from services.minio_service import MinIOService

log = get_logger(__name__)

GENOME_NAME = "organism_genome.fna"
FEATURE_TABLE_NAME = "organism_feature_table.txt"
PROTEIN_NAME = "organism_protein.faa"
CACHE_PREFIX = "_organism_cache"


def _ftp_feature_table_url(full_accession: str) -> str:
    """Derive the RefSeq FTP feature-table URL from a full assembly dir name,
    e.g. 'GCF_000146045.2_R64' ->
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_feature_table.txt.gz
    """
    gcx = full_accession[:3]                          # GCF or GCA
    core = full_accession.split("_")[1].split(".")[0] # 9 digits, e.g. 000146045
    triplets = f"{core[0:3]}/{core[3:6]}/{core[6:9]}"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{gcx}/{triplets}/"
        f"{full_accession}/{full_accession}_feature_table.txt.gz"
    )


def _datasets_accession(full_accession: str) -> str:
    """The `datasets` CLI wants the bare assembly accession (GCF_000146045.2),
    not the full directory name that includes the assembly label (..._R64)."""
    parts = full_accession.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else full_accession


def _download_via_datasets(full_accession: str, workdir: str):
    """Returns (genome_path, protein_path_or_None) downloaded into workdir."""
    zip_path = os.path.join(workdir, "ncbi.zip")
    subprocess.run(
        ["datasets", "download", "genome", "accession", _datasets_accession(full_accession),
         "--include", "genome,protein", "--filename", zip_path],
        check=True,
    )
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(os.path.join(workdir, "unzipped"))
    data_root = os.path.join(workdir, "unzipped", "ncbi_dataset", "data")
    genome_matches = glob(os.path.join(data_root, "**", "*_genomic.fna"), recursive=True) \
        or glob(os.path.join(data_root, "**", "*.fna"), recursive=True)
    if not genome_matches:
        raise RuntimeError(f"datasets returned no genome .fna for {full_accession}")
    protein_matches = glob(os.path.join(data_root, "**", "*protein*.faa"), recursive=True) \
        or glob(os.path.join(data_root, "**", "*.faa"), recursive=True)
    return genome_matches[0], (protein_matches[0] if protein_matches else None)


def _download_feature_table(full_accession: str, dest_path: str):
    url = _ftp_feature_table_url(full_accession)
    log.info(f"Fetching feature table: {url}")
    with urllib.request.urlopen(url) as resp, open(dest_path + ".gz", "wb") as out:
        shutil.copyfileobj(resp, out)
    with gzip.open(dest_path + ".gz", "rb") as gz, open(dest_path, "wb") as out:
        shutil.copyfileobj(gz, out)
    os.remove(dest_path + ".gz")


def _try_cache_download(client, bucket: str, accession: str, input_dir: str) -> bool:
    """If genome+feature_table are cached for this accession, pull all cached files
    (incl. protein if present) into input_dir. Returns True on cache hit."""
    base = f"{CACHE_PREFIX}/{accession}"
    required = {GENOME_NAME, FEATURE_TABLE_NAME}
    try:
        for name in required:
            client.stat_object(bucket, f"{base}/{name}")
    except Exception:
        return False
    for name in (GENOME_NAME, FEATURE_TABLE_NAME, PROTEIN_NAME):
        try:
            client.fget_object(bucket, f"{base}/{name}", os.path.join(input_dir, name))
        except Exception:
            if name == PROTEIN_NAME:
                continue  # protein is optional
            raise
    log.info(f"Organism cache hit for {accession}")
    return True


def _cache_upload(client, bucket: str, accession: str, input_dir: str):
    base = f"{CACHE_PREFIX}/{accession}"
    for name in (GENOME_NAME, FEATURE_TABLE_NAME, PROTEIN_NAME):
        path = os.path.join(input_dir, name)
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            try:
                client.fput_object(bucket, f"{base}/{name}", path)
            except Exception as ex:
                log.warning(f"Failed to cache {name} for {accession}: {ex}")


def main():
    accession = (os.getenv("ORGANISM_ACCESSION") or "").strip()
    if not accession:
        return  # upload path; nothing to fetch

    input_dir = os.getenv("JOB_INPUT_DIR")
    bucket = os.getenv("JOB_TYPE")
    if not input_dir or not bucket:
        log.error("JOB_INPUT_DIR / JOB_TYPE not set; cannot fetch organism")
        sys.exit(1)
    os.makedirs(input_dir, exist_ok=True)

    svc = MinIOService()
    client = svc.client
    try:
        svc.ensure_bucket_exists(bucket)
    except Exception:
        pass

    if _try_cache_download(client, bucket, accession, input_dir):
        return

    log.info(f"Fetching organism {accession} from NCBI")
    with tempfile.TemporaryDirectory() as workdir:
        genome_src, protein_src = _download_via_datasets(accession, workdir)
        shutil.copyfile(genome_src, os.path.join(input_dir, GENOME_NAME))
        if protein_src:
            shutil.copyfile(protein_src, os.path.join(input_dir, PROTEIN_NAME))
        _download_feature_table(accession, os.path.join(input_dir, FEATURE_TABLE_NAME))

    _cache_upload(client, bucket, accession, input_dir)


if __name__ == "__main__":
    main()
