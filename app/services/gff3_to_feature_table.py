"""Convert an uploaded GFF3 annotation file to the NCBI RefSeq `feature_table.txt`
shape that CRISPR-COPIES (`main.py`) and the genome-viewer service consume.

The method reads the gene table as a tab-delimited table and selects these columns
(main.py:963-966), filtering rows to `assembly_unit == 'Primary Assembly'`:

    assembly_unit, # feature, class, chromosome, genomic_accession, start, end,
    strand, locus_tag, product_accession

GFF3 carries most of this directly; the friction points (documented in the spec):
  * `class == 'with_protein'`  — NOT a GFF field; DERIVED here (a gene is with_protein
     if it has a CDS descendant; CDS rows are with_protein by definition).
  * `chromosome` (friendly #Name) — GFF3 has no friendly name; we fall back to the seqid
     (accession), so result rows show the accession as the chromosome label.
  * `assembly_unit` — absent in GFF3; emitted as 'Primary Assembly' for every row so the
     method's filter keeps them.

Handles the common NCBI gene->mRNA->CDS and prokaryotic gene->CDS hierarchies, and
tolerates attribute-key dialects (Name/gene, locus_tag, protein_id, ID/Parent).
Pure functions — no I/O dependencies — so it is unit-testable in isolation.
"""

import csv
from typing import Optional

OUTPUT_COLUMNS = [
    "assembly_unit", "# feature", "class", "chromosome", "genomic_accession",
    "start", "end", "strand", "locus_tag", "product_accession", "symbol",
]


def _parse_attrs(field: str) -> dict:
    """Parse a GFF3 column-9 attribute string (key=value;key=value) into a dict.
    Tolerates URL-ish values and missing trailing semicolons."""
    attrs = {}
    for part in field.strip().split(";"):
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        attrs[k.strip()] = v.strip()
    return attrs


def _attr_any(attrs: dict, *keys: str) -> str:
    for k in keys:
        if attrs.get(k):
            return attrs[k]
    return ""


def convert_gff3_text(text: str) -> list[dict]:
    """Parse GFF3 text → list of feature_table-shaped row dicts (gene + CDS features)."""
    raw = []  # (type, seqid, start, end, strand, attrs)
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        seqid, _src, ftype, start, end, _score, strand, _phase, attr_field = cols[:9]
        raw.append((ftype, seqid, start, end, strand, _parse_attrs(attr_field)))

    # Index by feature ID and record parent links, to derive which genes have a CDS.
    id_to_type: dict[str, str] = {}
    id_to_parents: dict[str, list[str]] = {}
    for ftype, _seqid, _s, _e, _st, attrs in raw:
        fid = attrs.get("ID")
        if fid:
            id_to_type[fid] = ftype
            parent = attrs.get("Parent", "")
            if parent:
                id_to_parents[fid] = parent.split(",")

    def ancestor_gene(fid: str, _depth: int = 0) -> Optional[str]:
        """Walk Parent links up to the enclosing gene feature ID (depth-bounded)."""
        if _depth > 10:
            return None
        if id_to_type.get(fid) == "gene":
            return fid
        for p in id_to_parents.get(fid, []):
            g = ancestor_gene(p, _depth + 1)
            if g:
                return g
        return None

    genes_with_cds: set[str] = set()
    for ftype, _seqid, _s, _e, _st, attrs in raw:
        if ftype == "CDS":
            for p in attrs.get("Parent", "").split(",") if attrs.get("Parent") else []:
                g = ancestor_gene(p)
                if g:
                    genes_with_cds.add(g)

    rows: list[dict] = []
    for ftype, seqid, start, end, strand, attrs in raw:
        if ftype not in ("gene", "CDS"):
            continue
        if ftype == "gene":
            cls = "with_protein" if attrs.get("ID") in genes_with_cds else "without_protein"
        else:  # CDS
            cls = "with_protein"
        rows.append({
            "assembly_unit": "Primary Assembly",
            "# feature": ftype,
            "class": cls,
            "chromosome": seqid,            # no friendly name in GFF3 → use the accession
            "genomic_accession": seqid,
            "start": start,
            "end": end,
            "strand": strand,
            "locus_tag": _attr_any(attrs, "locus_tag"),
            "product_accession": _attr_any(attrs, "protein_id"),
            "symbol": _attr_any(attrs, "Name", "gene", "gene_name"),
        })
    return rows


def convert_gff3_file(src_path: str, dest_path: str) -> int:
    """Read a GFF3 file, write the feature_table TSV. Returns the row count."""
    with open(src_path) as fh:
        rows = convert_gff3_text(fh.read())
    with open(dest_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)
