"""FASTA measurement, kept free of HTTP and of rdkit so it can be imported anywhere."""


def _residue_count(lines: list) -> int:
    """Length of one record's body: whitespace removed, trailing stop codon dropped."""
    seq = ''.join(''.join(lines).split())
    return len(seq.rstrip('*'))


def count_fasta_residues(fasta: str) -> int:
    """Residue count of the longest record in a FASTA string.

    Longest record, not the total: the constraint being guarded is peak GPU memory
    during a fold, which is driven by the largest single chain, not by how many chains
    were submitted.

    Counts every non-whitespace body character, not just letters, on purpose. SimpleFold
    reads the file with Biopython's FASTA parser and hands `str(record.seq)` to the model
    unchanged, so gap characters, position digits and `;` lines after a header all reach
    the model as residue positions and all cost GPU memory. Counting only letters would
    under-count what the GPU sees. The one exception is a trailing '*': the frontend
    strips it when measuring, so it is dropped here for the same input to get the same
    answer on both sides.

    Tolerant of the shapes that actually arrive -- a bare sequence with no header,
    wrapped lines, blank lines -- because this runs before any tool-specific parsing.
    """
    longest = 0
    current = []
    for line in (fasta or '').splitlines():
        if line.startswith('>'):
            longest = max(longest, _residue_count(current))
            current = []
        else:
            current.append(line)
    return max(longest, _residue_count(current))
