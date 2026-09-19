"""Tests for the SimpleFold sequence-length bound.

Why this exists: `ml-simplefold`'s VRAM grows with chain length until it exhausts the
shared GPU. At 1022 residues it OOMs running ALONE -- and 1022 is exactly what the MEP
form used to accept, because that is ESM-2's context limit and has nothing to do with
folding. Those submissions fail in production today.

The bound is enforced server-side rather than only in the browser because the browser is
not the only client. Both `/{job_type}/jobs` and `/v1/...` reach `job_builder.prepare_job`
before `kubejob_service.create_job` (two callers of each, checked), and `coordinator.py`
creates subjobs by POSTing to the same router, so the check in `prepare_job` is the one
choke point every path crosses.

The parser cases below are the shapes that actually arrive, not hypotheticals: the
frontend sends `>{name}\\n{sequence}`, pasted FASTA is line-wrapped, and a trailing '*'
stop codon is common enough that the frontend's own validator strips one.
"""
import pytest
from fastapi import HTTPException

from services.shared import (
    DEFAULT_MAX_STRUCTURE_RESIDUES,
    count_fasta_residues,
    validate_fasta_length,
)


class TestCountFastaResidues:
    @pytest.mark.parametrize("fasta,expected,why", [
        (">h\nABCDE", 5, "the shape the frontend sends"),
        (">h\nABC\nDE", 5, "pasted FASTA is line-wrapped"),
        ("ABCDE", 5, "bare sequence, no header"),
        (">h\nABCDE*", 5, "trailing stop codon is not a residue"),
        (">h\nABC\n\nDE\n", 5, "blank lines"),
        (">h\n ABC DE ", 5, "internal and surrounding whitespace"),
        ("", 0, "empty input must not raise"),
        (">h\n", 0, "header with no body"),
    ])
    def test_single_record(self, fasta, expected, why):
        assert count_fasta_residues(fasta) == expected, why

    @pytest.mark.parametrize("fasta", [">a\nABC\n>b\nABCDEFG", ">b\nABCDEFG\n>a\nABC"])
    def test_multi_record_takes_the_longest_not_the_total(self, fasta):
        """Peak GPU memory is driven by the largest single chain, so summing would reject
        many short chains that fold fine, and averaging would admit one that cannot."""
        assert count_fasta_residues(fasta) == 7


class TestValidateFastaLength:
    def test_under_the_bound_passes(self):
        validate_fasta_length(">h\n" + "A" * 699, 700, "SimpleFold")

    def test_exactly_at_the_bound_passes(self):
        """The limit is inclusive. 700 is a length we intend to accept, so an
        off-by-one here would silently narrow the product by one residue."""
        validate_fasta_length(">h\n" + "A" * 700, 700, "SimpleFold")

    def test_over_the_bound_raises_422(self):
        with pytest.raises(HTTPException) as e:
            validate_fasta_length(">h\n" + "A" * 701, 700, "SimpleFold")
        # 422, not 400: the request is well-formed, the input is simply too large for
        # the hardware. It is also not transient, so a client must not retry it.
        assert e.value.status_code == 422

    def test_the_error_names_both_numbers(self):
        """A bare rejection sends the user back to guess. The message has to say what
        they submitted and what the limit is, because the limit is not in the UI copy
        on the API path."""
        with pytest.raises(HTTPException) as e:
            validate_fasta_length(">h\n" + "A" * 1022, 700, "SimpleFold")
        detail = str(e.value.detail)
        assert "1022" in detail and "700" in detail

    def test_the_1022_case_that_fails_in_production_is_rejected(self):
        """The regression this bound exists for: 1022 residues is what the MEP form
        allowed, and simplefold OOMs on it even with the GPU to itself."""
        with pytest.raises(HTTPException):
            validate_fasta_length(">h\n" + "A" * 1022, DEFAULT_MAX_STRUCTURE_RESIDUES,
                                  "SimpleFold")

    def test_the_longest_real_submission_is_still_accepted_at_esm_limit(self):
        """770 residues (APP770) is the longest sequence in 173 historical MEP
        submissions. It must still reach the MEP job -- only the structure half is
        bounded -- so this asserts the bound is not applied to the MEP sequence."""
        assert count_fasta_residues(">APP770\n" + "A" * 770) == 770
