from enum import Enum
from typing import Optional

from sqlalchemy.types import TypeDecorator
from sqlmodel import SQLModel, Field, Relationship
from sqlmodel.sql.sqltypes import AutoString
from pydantic import BaseModel

from ..enums import JobType, JobStatus


class EnumValueString(TypeDecorator):
    """Store an Enum by VALUE in a VARCHAR column, and read it back as the Enum.

    Needed because the two ways this schema gets built disagree about enum columns, and
    only one of them is what production runs.

    SQLModel 0.0.8 mapped a str-Enum field to AutoString, so every migration in this repo
    created `job.phase` and `job.type` as VARCHAR holding enum VALUES ('ml-simplefold').
    SQLModel 0.0.9 changed that mapping to a native sa.Enum, so the upgrade to 0.0.22
    silently redeclared both columns as Postgres types 'jobstatus' and 'jobtype' -- which
    no migration creates and which exist in no deployed database. Every query filtering on
    them rendered '$1::jobtype' and failed with UndefinedObjectError, which is every job
    submission, via the duplicate check in routers/job.py.

    Switching to a native enum for real is not an option: sa.Enum persists NAMES, so its
    labels are 'ML_SIMPLEFOLD' while every existing row -- and every frontend,
    coordinator.py, and the published /v1 input schema -- uses 'ml-simplefold'. Migrating
    would have to rewrite every row and every client.

    Plain AutoString fixes the outage but reads back a bare str, which leaves the model
    holding a str where it declares an enum. Pydantic 2 warns on every serialization, and
    the test suite -- which builds its schema from the metadata, where these were real
    enums -- has always seen enum members. Converting on the way out keeps the column as
    the VARCHAR the database has, and the attribute as the enum the code expects.
    """
    impl = AutoString
    cache_ok = True

    def __init__(self, enum_class, *args, **kwargs):
        self.enum_class = enum_class
        super().__init__(*args, **kwargs)

    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        return value.value if isinstance(value, Enum) else str(value)

    def process_result_value(self, value, dialect):
        if value is None:
            return None
        try:
            return self.enum_class(value)
        except ValueError:
            # A stored value that is no longer a member -- a tool renamed or retired since
            # the row was written. Those rows have to stay readable, so hand back the raw
            # string rather than failing every query that happens to touch one.
            return value

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
    # VARCHAR-backed, not native enums. See EnumValueString above: declaring these as bare
    # JobStatus/JobType makes SQLModel 0.0.22 ask for Postgres types that no migration
    # creates, which took down every job submission. These columns carry no validation of
    # their own -- SQLModel skips Pydantic validation on table models -- so the input guard
    # is the `job_type in JobTypes` check in routers/job.py, which is deliberate and
    # commented there.
    phase: JobStatus = Field(default=JobStatus.QUEUED, nullable=False,
                             sa_type=EnumValueString(JobStatus))
    type: JobType = Field(default=None, nullable=False,
                          sa_type=EnumValueString(JobType))
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
