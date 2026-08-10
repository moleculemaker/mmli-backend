"""Typed validation for the CRISPR-COPIES job parameters (Pydantic v2).

`POST /{job_type}/jobs` accepts `job_info` as an opaque JSON string, so the CRISPR-COPIES
parameter contract is NOT described by the OpenAPI spec. This module makes that contract
explicit and validates it server-side before a job is queued, returning a clean 422 instead
of letting a malformed payload reach (and crash, or silently misbehave in) the
`python main.py` container.

Bounds/enums mirror the method's argparse + in-code checks in COPIES `code/main.py`
(the upstream source of truth). The frontend's hand-written validators must match this
module; neither is derived from the other (see CLAUDE.md / the integration spec).
"""

import re
from typing import Optional

# Using Pydantic V2 field_validator, model_validator, and ValidationInfo
from pydantic import BaseModel, Field, field_validator, model_validator, ValidationInfo
from models.deg_organisms import DEG_ORGANISMS

# --- shared vocab (mirrors main.py) ----------------------------------------- #

# IUPAC nucleotide alphabet (incl. ambiguity codes) for PAM / restriction enzymes.
_IUPAC_RE = re.compile(r"[ACGTRYSWKMBDHVN]+", re.IGNORECASE)

# Unambiguous bases only, for the backbone self-complementarity sequence.
_ACGT_RE = re.compile(r"[ACGT]+", re.IGNORECASE)

ON_TARGET_MODELS = (
    "Doench et al. 2016",
    "CROPSR",
    "DeepGuide (Cas9)",
    "DeepGuide (Cas12a)",
    "sgRNA_ecoli (Cas9)",
    "sgRNA_ecoli (eSpCas9)",
)

# Models that only make sense for a given enzyme (used for system-coupling checks).
_CAS9_ONLY_MODELS = {"DeepGuide (Cas9)", "sgRNA_ecoli (Cas9)", "sgRNA_ecoli (eSpCas9)"}
_CAS12A_ONLY_MODELS = {"DeepGuide (Cas12a)"}

PAM_ORIENTATIONS = ("3prime", "5prime")
DISTANCE_TYPES = ("hamming", "levenshtein")
CRISPR_SYSTEMS = ("CRISPR/Cas9", "CRISPR/Cas12a")  # Cas13d intentionally unsupported


# --- leaf models -------------------------------------------------------------

class FileMetadataModel(BaseModel):
    """A frontend FileMetadata reference. Only the name is needed downstream; the file
    itself is already in MinIO. Extra fields are ignored (Pydantic v2 default)."""
    filename: Optional[str] = None
    name: Optional[str] = None
    url: Optional[str] = None

    def resolved_name(self) -> Optional[str]:
        return self.filename or self.name or self.url


class OrganismIdentifierModel(BaseModel):
    label: Optional[str] = None
    value: Optional[str] = None


class OrganismModel(BaseModel):
    organismIdentifier: Optional[OrganismIdentifierModel] = None
    genome: Optional[FileMetadataModel] = None
    annotations: Optional[FileMetadataModel] = None
    protein: Optional[FileMetadataModel] = None

    # In V2, root_validator(skip_on_failure=True) becomes model_validator(mode="after")
    # By default, "after" validators are only executed if individual fields pass validation.
    @model_validator(mode="after")
    def exactly_one_source(self) -> "OrganismModel":
        ident = self.organismIdentifier
        has_dropdown = bool(ident and ident.value)
        has_upload = bool(self.genome) and bool(self.annotations)

        if has_dropdown and has_upload:
            raise ValueError(
                "Provide an organism EITHER from the dropdown OR as uploaded "
                "genome+annotations files, not both."
            )
        if not has_dropdown and not has_upload:
            raise ValueError(
                "An organism is required: select one from the dropdown, or upload "
                "both a genome and an annotations file."
            )
        return self


class GuideRNAModel(BaseModel):
    pam: str = Field(..., min_length=1)
    pamOrientation: str
    onTargetScore: str
    guideLength: int = Field(..., ge=10, le=40)
    seedLength: int = Field(..., ge=0, le=27)
    editDistance: int = Field(..., ge=0, le=20)
    gcContentMin: int = Field(..., ge=0, le=100)
    gcContentMax: int = Field(..., ge=0, le=100)
    distanceType: str
    backboneSequence: Optional[str] = None
    restrictionEnzyme: Optional[str] = None
    polyG: int = Field(0, ge=0, le=10)
    polyT: int = Field(0, ge=0, le=10)

    @field_validator("pamOrientation")
    @classmethod
    def _orientation(cls, v: str) -> str:
        if v not in PAM_ORIENTATIONS:
            raise ValueError(f"pamOrientation must be one of {PAM_ORIENTATIONS}")
        return v

    @field_validator("distanceType")
    @classmethod
    def _distance_type(cls, v: str) -> str:
        if v not in DISTANCE_TYPES:
            raise ValueError(f"distanceType must be one of {DISTANCE_TYPES}")
        return v

    @field_validator("onTargetScore")
    @classmethod
    def _on_target(cls, v: str) -> str:
        if v not in ON_TARGET_MODELS:
            raise ValueError(f"onTargetScore must be one of {ON_TARGET_MODELS}")
        return v

    @field_validator("pam")
    @classmethod
    def _pam_iupac(cls, v: str) -> str:
        if not _IUPAC_RE.fullmatch(v.strip()):
            raise ValueError("PAM may only contain IUPAC nucleotide codes (ACGTRYSWKMBDHVN)")
        return v

    @field_validator("restrictionEnzyme")
    @classmethod
    def _re_iupac(cls, v: Optional[str]) -> Optional[str]:
        if v and not _IUPAC_RE.fullmatch(v.strip()):
            raise ValueError("restrictionEnzyme may only contain IUPAC nucleotide codes (ACGTRYSWKMBDHVN)")
        return v

    @field_validator("backboneSequence")
    @classmethod
    def _backbone_acgt(cls, v: Optional[str]) -> Optional[str]:
        if v and not _ACGT_RE.fullmatch(v.strip()):
            raise ValueError("backboneSequence may only contain unambiguous bases (ACGT)")
        return v

    @model_validator(mode="after")
    def _cross_field(self) -> "GuideRNAModel":
        if self.seedLength is not None and self.guideLength is not None and self.seedLength >= self.guideLength:
            raise ValueError("seedLength must be less than guideLength")
        if self.gcContentMin is not None and self.gcContentMax is not None and self.gcContentMin >= self.gcContentMax:
            raise ValueError("gcContentMin must be less than gcContentMax")
        return self


class HomologyArmModel(BaseModel):
    lengthBp: int = Field(..., ge=5, le=1000)
    restrictionEnzyme: Optional[str] = None
    polyG: int = Field(0, ge=0, le=10)
    polyT: int = Field(0, ge=0, le=10)

    @field_validator("restrictionEnzyme")
    @classmethod
    def _re_iupac(cls, v: Optional[str]) -> Optional[str]:
        if v and not _IUPAC_RE.fullmatch(v.strip()):
            raise ValueError("restrictionEnzyme may only contain IUPAC nucleotide codes (ACGTRYSWKMBDHVN)")
        return v


class IntergenicRegionModel(BaseModel):
    distance: int = Field(..., ge=1)  # --intspace
    distanceToMeasureGeneDensity: int = Field(..., ge=1)  # --gene_density_len
    distalEndOfChromosomeToAvoid: int = Field(..., ge=0)  # --distal_end_len
    referenceOrganismIdentifier: Optional[str] = None  # --blast_org

    @field_validator("referenceOrganismIdentifier")
    @classmethod
    def _deg_member(cls, v: Optional[str]) -> Optional[str]:
        if v in (None, ""):
            return v
        if v not in DEG_ORGANISMS:
            raise ValueError(
                "referenceOrganismIdentifier must be one of the curated DEG organisms "
                "(see deg.csv); essential-gene annotation cannot run for other organisms."
                )
        return v


class ParametersModel(BaseModel):
    crisprSystem: str
    guideRNA: GuideRNAModel
    homologyArm: HomologyArmModel
    intergenicRegion: IntergenicRegionModel

    @field_validator("crisprSystem")
    @classmethod
    def _system(cls, v: str) -> str:
        if v not in CRISPR_SYSTEMS:
            raise ValueError(
                f"crisprSystem must be one of {CRISPR_SYSTEMS} "
                "(CRISPR/Cas13d is not supported by the method)"
            )
        return v

    @model_validator(mode="after")
    def _system_coupling(self) -> "ParametersModel":
        system = self.crisprSystem
        guide = self.guideRNA
        arm = self.homologyArm
        intergenic = self.intergenicRegion

        if not (system and guide and arm and intergenic):
            return self

        use_cas9 = system == "CRISPR/Cas9"

        # Scoring model must match the enzyme
        if use_cas9 and guide.onTargetScore in _CAS12A_ONLY_MODELS:
            raise ValueError(f"on-target model '{guide.onTargetScore}' is Cas12a-specific; not valid for a Cas9 job")
        if (not use_cas9) and guide.onTargetScore in _CAS9_ONLY_MODELS:
            raise ValueError(f"on-target model '{guide.onTargetScore}' is Cas9-specific; not valid for a {system} job")

        # Mirror main.py:985-988 — for NGG/3prime, the distal-end cutoff must exceed the homology-arm length
        if guide.pam == "NGG" and guide.pamOrientation == "3prime":
            if intergenic.distalEndOfChromosomeToAvoid < arm.lengthBp:
                raise ValueError(
                    "distalEndOfChromosomeToAvoid must be >= homology-arm length when "
                    "PAM is NGG and orientation is 3prime"
                )
        return self


class CrisprCopiesJobInfo(BaseModel):
    """Top-level CRISPR-COPIES `job_info` payload. Unknown keys (e.g. `email`) are ignored
    by Pydantic v2's default behavior."""
    organism: OrganismModel
    parameters: ParametersModel
