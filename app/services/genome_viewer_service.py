"""Data feeds for the frontend genome-viewer (gene track + protein/AA track).

Both read the job's inputs from the shared working volume the backend mounts
(`get_job_root_dir_from_id(job_type, job_id)/in`), where uploaded files are synced by
prejob.py and dropdown organisms are written by fetch_organism.py.

  gene_annotations()    -> {chromosome: [{name, start, end}]}   (B)
  protein_translation() -> AA string over a bp window, '-' for non-coding   (C)

Both target the NCBI RefSeq `*_feature_table.txt` format that the method consumes
(tab-delimited, columns incl. `# feature`, `chromosome`, `genomic_accession`, `start`,
`end`, `strand`, `symbol`/`name`/`locus_tag`, `product_accession`). NOTE: the frontend
currently accepts GFF/GFF3 uploads for annotations, which is NOT this format — uploaded-
GFF jobs will yield empty results here (and won't run in the method either); see the
annotations-format mismatch flagged in the spec. Dropdown/NCBI jobs use feature_table.txt.

UNVERIFIED in CI/sandbox — needs a real job on the cluster (shared volume + real files).
"""

import os
from glob import glob
from typing import Optional

import pandas as pd

from config import get_logger
from services.kubejob_service import get_job_root_dir_from_id

log = get_logger(__name__)

OUTPUT_FILE_NAME = "output.csv"
_GENE_WINDOW_PAD_BP = 50_000  # genes returned within this padding around the candidate-site span

# Standard genetic code (no biopython dependency in the backend).
_CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _job_input_dir(job_type: str, job_id: str) -> str:
    return os.path.join(get_job_root_dir_from_id(job_type, job_id), "in")


def _find_one(input_dir: str, preferred: str, *patterns: str) -> Optional[str]:
    """Return the preferred filename if present, else the first glob match."""
    p = os.path.join(input_dir, preferred)
    if os.path.isfile(p):
        return p
    for pat in patterns:
        hits = sorted(glob(os.path.join(input_dir, pat)))
        if hits:
            return hits[0]
    return None


def _read_feature_table(job_type: str, job_id: str) -> Optional[pd.DataFrame]:
    in_dir = _job_input_dir(job_type, job_id)
    path = _find_one(in_dir, "organism_feature_table.txt", "*_feature_table.txt", "*.txt")
    if not path:
        log.warning(f"No feature_table found for {job_type}/{job_id} in {in_dir}")
        return None
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as ex:
        log.warning(f"Could not parse feature table {path} as NCBI feature_table.txt: {ex}")
        return None
    if "# feature" not in df.columns:
        log.warning(f"{path} is not an NCBI feature_table.txt (no '# feature' column) — skipping")
        return None
    return df


def _site_windows(job_type: str, job_id: str) -> dict:
    """Per-chromosome (min,max) candidate-site bp from the job's output.csv, padded.
    Returns {} (meaning 'no filtering') if output is unavailable."""
    out = os.path.join(get_job_root_dir_from_id(job_type, job_id), "out", OUTPUT_FILE_NAME)
    if not os.path.isfile(out):
        return {}
    try:
        df = pd.read_csv(out)
    except Exception:
        return {}
    if "Chromosome" not in df.columns or "Location" not in df.columns or df.empty:
        return {}
    windows = {}
    loc = pd.to_numeric(df["Location"], errors="coerce")
    grouped = df.assign(_loc=loc).dropna(subset=["_loc"]).groupby("Chromosome")["_loc"]
    for chrom, series in grouped:
        windows[str(chrom)] = (max(0, int(series.min()) - _GENE_WINDOW_PAD_BP),
                               int(series.max()) + _GENE_WINDOW_PAD_BP)
    return windows


def gene_annotations(job_type: str, job_id: str) -> dict:
    """B: {chromosome: [{name, start, end}]}, keyed by the same `Chromosome` (#Name) the
    result rows use, filtered to the candidate-site window per chromosome."""
    df = _read_feature_table(job_type, job_id)
    if df is None:
        return {}
    genes = df[df["# feature"] == "gene"]
    windows = _site_windows(job_type, job_id)
    name_cols = [c for c in ("symbol", "name", "locus_tag") if c in genes.columns]
    chrom_col = "chromosome" if "chromosome" in genes.columns else "genomic_accession"

    out: dict = {}
    for _, g in genes.iterrows():
        try:
            start, end = int(g["start"]), int(g["end"])
        except (ValueError, KeyError, TypeError):
            continue
        chrom = str(g.get(chrom_col, "")).strip()
        if not chrom:
            continue
        w = windows.get(chrom)
        if w and (end < w[0] or start > w[1]):
            continue
        name = next((str(g[c]).strip() for c in name_cols if str(g.get(c, "")).strip()), "")
        out.setdefault(chrom, []).append({"name": name, "start": start, "end": end})
    for chrom in out:
        out[chrom].sort(key=lambda x: x["start"])
    return out


def _load_fasta_record(genome_path: str, accession: str) -> Optional[str]:
    """Return the sequence (uppercase) of the FASTA record whose header id matches
    `accession` (matched on the first whitespace-delimited token after '>')."""
    seq_parts: list[str] = []
    capturing = False
    with open(genome_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if capturing:
                    break
                header_id = line[1:].split()[0] if len(line) > 1 else ""
                capturing = header_id == accession
            elif capturing:
                seq_parts.append(line.strip())
    return "".join(seq_parts).upper() if seq_parts else None


def _translate(seq: str) -> str:
    return "".join(_CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def protein_translation(job_type: str, job_id: str, chromosome: str, start: int, end: int) -> str:
    """C: AA string across [start, end] (1-based bp) on `chromosome`, at codon (3 bp)
    resolution — coding codons show their amino acid, non-coding positions show '-'.
    The frontend samples this string across the plot width."""
    df = _read_feature_table(job_type, job_id)
    if df is None or start >= end:
        return ""
    chrom_col = "chromosome" if "chromosome" in df.columns else "genomic_accession"
    on_chrom = df[df[chrom_col].astype(str).str.strip() == str(chromosome)]
    if on_chrom.empty:
        return ""
    accession = str(on_chrom.iloc[0].get("genomic_accession", "")).strip()

    in_dir = _job_input_dir(job_type, job_id)
    genome_path = _find_one(in_dir, "organism_genome.fna", "*_genomic.fna", "*.fna", "*.fa", "*.fasta")
    if not genome_path or not accession:
        return ""
    chrom_seq = _load_fasta_record(genome_path, accession)
    if not chrom_seq:
        return ""

    # AA per codon position across the window; '-' outside CDS.
    n_codons = max(0, (end - start) // 3)
    aa = ["-"] * n_codons
    cds = on_chrom[on_chrom["# feature"] == "CDS"]
    for _, c in cds.iterrows():
        try:
            c_start, c_end = int(c["start"]), int(c["end"])
        except (ValueError, KeyError, TypeError):
            continue
        if c_end < start or c_start > end:
            continue
        sub = chrom_seq[c_start - 1:c_end]                       # 1-based inclusive
        if str(c.get("strand", "+")).strip() == "-":
            sub = sub.translate(_COMPLEMENT)[::-1]
        protein = _translate(sub)
        for k, residue in enumerate(protein):
            genomic_bp = (c_start + k * 3) if str(c.get("strand", "+")).strip() != "-" else (c_end - k * 3)
            idx = (genomic_bp - start) // 3
            if 0 <= idx < n_codons:
                aa[idx] = residue
    return "".join(aa)
