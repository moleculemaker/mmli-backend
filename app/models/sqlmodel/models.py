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

    # Job timestamps
    time_created: int = Field(default=None, nullable=False)
    time_start: Optional[int] = Field(default=0, nullable=False)
    time_end: Optional[int] = Field(default=0, nullable=False)
    deleted: int = Field(default=0, nullable=False)

    # User metadata
    # owner = relationship("User", back_populates="jobs")
    # owner: int = Relationship(link_model="User", back_populates="jobs")
    user_agent: str = Field(default=None, nullable=False)


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
