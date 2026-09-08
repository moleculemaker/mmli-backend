"""Typed validation for the mutagenesis (SDM primer design) job parameters (Pydantic v1).

`POST /{job_type}/jobs` accepts `job_info` as an opaque JSON string, so this parameter
contract is NOT described by the OpenAPI spec. This module makes it explicit and validates
it server-side before a job is queued, returning a clean 422 instead of letting a malformed
payload reach (and crash, or silently misbehave in) the `python main.py` container.

Bounds/enums mirror the method's argparse + in-code checks in
`Primer_Design_and_Worklists/main.py` (the upstream source of truth). The frontend's
hand-written validators must match this module; neither is derived from the other.
See MUTAGENESIS_INTEGRATION_SPEC.md.
"""

from typing import Optional

from pydantic import BaseModel, Field, field_validator, validator

# Tm methods the method understands (`-tm`); SantaLucia (NN) is the default.
TM_METHODS = ("SantaLucia", "Wallace")

# NCBI genetic-code (codon table) ids supported by Biopython's
# CodonTable.unambiguous_dna_by_id. v1 defaults to 1 (Standard) and the FE does not
# expose a chooser, but we validate the value in case it is sent.
VALID_NCBI_CODON_TABLE_IDS = frozenset(
    {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26,
     27, 28, 29, 30, 31, 33}
)


class FileMetadataModel(BaseModel):
    """A frontend FileMetadata reference. Only the name is needed downstream; the file
    itself is already in MinIO. Extra fields are ignored (Pydantic v1 default)."""
    filename: Optional[str] = None
    name: Optional[str] = None
    url: Optional[str] = None

    def resolved_name(self) -> Optional[str]:
        return self.filename or self.name or self.url


class MutagenesisJobInfo(BaseModel):
    """Top-level mutagenesis `job_info` payload. Unknown keys (e.g. `email`) are ignored
    by Pydantic v1's default behavior.

    Inputs are uploaded files (materialized client-side from the form's text/upload
    fields, see Q1 in the spec):
      orfFile        - required; ORF DNA sequence (.txt)            -> -orf
      mutationList   - required; CSV with a `Mutations` header      -> -m
      leftOverhang   - optional; upstream flank sequence (.txt)     -> -left
      rightOverhang  - optional; downstream flank sequence (.txt)   -> -right
    """
    orfFile: FileMetadataModel
    mutationList: FileMetadataModel
    leftOverhang: Optional[FileMetadataModel] = None
    rightOverhang: Optional[FileMetadataModel] = None
    codonTableValue: int = Field(1)
    tmMethod: str = Field("SantaLucia")

    # field_validator rather than the deprecated validator shim: this is the one
    # validator here that needs the field's name, and v2 supplies that through a
    # ValidationInfo argument which the shim does not accept. The message is unchanged.
    @field_validator("orfFile", "mutationList")
    @classmethod
    def _required_file_has_name(cls, v, info):
        if not v or not v.resolved_name():
            raise ValueError(f"{info.field_name} must reference an uploaded file (filename required)")
        return v

    @validator("codonTableValue")
    def _valid_codon_table(cls, v):
        if v not in VALID_NCBI_CODON_TABLE_IDS:
            raise ValueError(
                f"codonTableValue must be a valid NCBI genetic-code id "
                f"(one of {sorted(VALID_NCBI_CODON_TABLE_IDS)})"
            )
        return v

    @validator("tmMethod")
    def _valid_tm_method(cls, v):
        if v not in TM_METHODS:
            raise ValueError(f"tmMethod must be one of {TM_METHODS}")
        return v
