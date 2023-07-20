from enum import Enum
from typing import Optional

from sqlmodel import SQLModel, Field, Relationship

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

    job_id: Optional[str]
    run_id: Optional[str]


# What is stored in the database for each Job
class Job(JobBase, table=True):
    id: int = Field(default=None, nullable=False, primary_key=True)

    # Job metadata
    #queue_position: int = Field(default=0, nullable=False)
    phase: JobStatus = Field(default=None, nullable=False)
    type: JobType = Field(default=None, nullable=False)
    image: str = Field(default=None, nullable=False)
    command: Optional[str] = Field(default=None, nullable=True)

    # Job timestamps
    time_created: int = Field(default=0, nullable=False)
    time_start: int = Field(default=0, nullable=False)
    time_end: int = Field(default=0, nullable=False)
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
    id: int
    job_id: str
    run_id: str

    # Updatable properties
    time_start: Optional[int] = None
    time_end: Optional[int] = None
    job_info: Optional[str] = None
    email: Optional[str] = None
    image: Optional[str] = None
    command: Optional[str] = None
    phase: Optional[JobStatus] = None
