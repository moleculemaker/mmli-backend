from enum import Enum
from typing import Optional

from sqlmodel import SQLModel, Field, Relationship
from pydantic import BaseModel

from ..enums import JobType, JobStatus

# CREATE TABLE IF NOT EXISTS `job`(
#     `id` int NOT NULL AUTO_INCREMENT,
#     `job_id` varchar(32) NOT NULL,
#     `run_id` varchar(64) NOT NULL,
#     `user_id` varchar(50) NOT NULL,
#     `command` text NOT NULL DEFAULT '',
#     `type` varchar(50) NOT NULL,
#     `phase` varchar(50) NOT NULL,
#     `time_created` datetime NOT NULL DEFAULT 0,
#     `time_start` datetime NOT NULL DEFAULT 0,
#     `time_end` datetime NOT NULL DEFAULT 0,
#     `user_agent` varchar(256) NOT NULL DEFAULT '',
#     `email` varchar(128) NOT NULL DEFAULT '',
#     `job_info` text NOT NULL DEFAULT '',
#     `deleted` boolean NOT NULL DEFAULT 0,
#     `queue_position` int,
#     PRIMARY KEY (`id`), UNIQUE KEY `id` (`id`), UNIQUE KEY `job_id` (`job_id`)
# )


# What is required for any Job (base)
class JobBase(SQLModel):
    # Job configuration
    job_info: Optional[str] = None
    email: Optional[str] = None

    job_id: Optional[str] = Field(default=None, nullable=False, primary_key=True)
    run_id: Optional[str] = None


# What is stored in the database for each Job
class Job(JobBase, table=True):
    # Job metadata
    #queue_position: int = Field(default=0, nullable=False)
    phase: JobStatus = Field(default=JobStatus.QUEUED, nullable=False)
    type: JobType = Field(default=None, nullable=False)
    image: str = Field(default=None, nullable=True)
    command: Optional[str] = Field(default=None, nullable=True)

    # The image reference as configured, e.g. "moleculemaker/novostoic", is frequently
    # untagged or points at a mutable tag, so it does not identify what actually ran.
    # This is the immutable digest reported by the pod once it starts, which is what a
    # result has to cite to be reproducible. Null for jobs that predate this column and
    # for any job whose pod was reaped before the watcher could read it.
    image_digest: Optional[str] = Field(default=None, nullable=True)

    # Set on subjobs created by a parent job's coordinator. Until now this relationship
    # existed only as keys inside the job_info JSON blob, so it could not be queried and
    # was invisible to anything reading the schema.
    parent_job_id: Optional[str] = Field(default=None, nullable=True, index=True)

    # Job timestamps
    time_created: int = Field(default=None, nullable=False)
    time_start: Optional[int] = Field(default=0, nullable=False)
    time_end: Optional[int] = Field(default=0, nullable=False)
    deleted: int = Field(default=0, nullable=False)

    # User metadata
    # owner = relationship("User", back_populates="jobs")
    # owner: int = Relationship(link_model="User", back_populates="jobs")
    user_agent: str = Field(default=None, nullable=False)

    # Which API the submission arrived through: legacy, v1 or mcp. Recorded so that
    # adoption of the versioned API can be reported without inferring it from user
    # agent strings, which are neither stable nor trustworthy.
    client_surface: Optional[str] = Field(default=None, nullable=True, index=True)

    # Origin header, when the caller sent one. Distinguishes our own frontends from a
    # script or a notebook; scripts send no Origin at all, which is itself the signal.
    client_origin: Optional[str] = Field(default=None, nullable=True)

    # Pseudonymous, salted derivation of the client address, used only to approximate a
    # distinct-user count for submissions that carry no email. Never the address itself.
    # Null when no salt is configured, which is the default. See services/analytics.py.
    client_fingerprint: Optional[str] = Field(default=None, nullable=True, index=True)


# Anything additional that is passed to the API to create a new Job
class JobCreate(JobBase):
    pass


# Anything additional that is passed to the API to create a new Job
class JobUpdate(SQLModel):
    # Immutable metadata
    job_id: str 
    run_id: int

    # Updatable properties
    time_start: Optional[int] = None
    time_end: Optional[int] = None
    job_info: Optional[str] = None
    email: Optional[str] = None
    image: Optional[str] = None
    command: Optional[str] = None
    phase: Optional[JobStatus] = None

class IdempotencyKey(SQLModel, table=True):
    """Maps a client-supplied Idempotency-Key to the job it created.

    Submission is not naturally idempotent: every POST starts a container. A script
    whose request times out mid-flight currently has no safe way to retry - it either
    abandons a job that may be running or starts a second one. Replaying the same key
    returns the original job instead.

    Scoped by job type so two tools cannot collide on the same key.
    """
    __tablename__ = 'idempotency_key'

    key: str = Field(primary_key=True)
    job_type: str = Field(primary_key=True)
    job_id: str = Field(index=True, nullable=False)
    time_created: int = Field(nullable=False)


class FlaggedMolecule(SQLModel, table=True):
    smile: str = Field(default=None, primary_key=True)
    job_id: str = Field(default=None, primary_key=True, foreign_key="job.job_id")
    # For future: which doc contains the molecule within the job when multiple docs are allowed
    doc_id: str = Field(default=None)
    time_created: Optional[int] = Field(default=None, nullable=False)

class FlaggedMoleculeDelete(BaseModel):
    smile: str 
    job_id: str

class ChemicalIdentifier(SQLModel, table=True):
    # change table name to chemical_identifier
    __tablename__ = 'chemical_identifier'

    id: int = Field(default=None, primary_key=True)
    metanetx_id: Optional[str] = Field(default=None, nullable=True, index=True)
    inchi: Optional[str] = Field(default=None, nullable=True)
    inchi_key: Optional[str] = Field(default=None, nullable=True)
    name: Optional[str] = Field(default=None, nullable=True, index=True)
    smiles: Optional[str] = Field(default=None, nullable=True, index=True)
    reference: Optional[str] = Field(default=None, nullable=True)
    kegg_id: Optional[str] = Field(default=None, nullable=True, index=True)
    formula: Optional[str] = Field(default=None, nullable=True)
    is_cofactor: Optional[bool] = Field(default=False, nullable=True)

class SavedMolecule(SQLModel, table=True):
    __tablename__ = 'saved_molecule'
    
    id: Optional[int] = Field(default=None, primary_key=True)
    email: Optional[str] = Field(default=None, index=True, nullable=True)
    job_id: str = Field(default=None, index=True)
    molecule_id: str = Field(default=None, index=True)
    time_created: Optional[int] = Field(default=None)

class SavedMoleculeDelete(BaseModel):
    id: int
