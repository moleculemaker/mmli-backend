import os
import time

from dotenv import load_dotenv
from sqlmodel import SQLModel
from sqlmodel.ext.asyncio.session import AsyncSession

from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine

from config import get_logger, SQLALCHEMY_DATABASE_URL, SERVER_LOGLEVEL
from models.enums import JobStatus
from models.sqlmodel.models import Job

# By default, use SQLite (for local development)
load_dotenv()


# app_secrets['db']['url'] = "postgresql://user:password@postgresserver/db"
def create_db_engine():
    # SQLAlchemy 2.0 builds async engines through create_async_engine. The old
    # AsyncEngine(create_engine(...)) wrapper is gone, and sqlmodel no longer re-exports
    # AsyncEngine at all. `future=True` is dropped because 2.0 behavior is the only
    # behavior now and passing it is an error.
    return create_async_engine(
        SQLALCHEMY_DATABASE_URL,
        echo=(True if SERVER_LOGLEVEL == 'DEBUG' else False),
    )


engine = create_db_engine()


log = get_logger(__name__)


async def init_db():
    async with engine.begin() as conn:
        # await conn.run_sync(SQLModel.metadata.drop_all)
        await conn.run_sync(SQLModel.metadata.create_all)


async def get_session() -> AsyncSession:
    async_session = async_sessionmaker(
        engine, class_=AsyncSession, expire_on_commit=False
    )
    async with async_session() as session:
        yield session


async def update_job_phase(session, job_id, phase: JobStatus):
        db_job: Job = await session.get(Job, job_id)
        db_job.phase = phase
        if phase == JobStatus.PROCESSING:
            db_job.time_start = int(time.time())
        else:
            db_job.time_end = int(time.time())
        session.add(db_job)
        await session.commit()
        log.info('Updated job phase: %s -> %s' % (job_id, phase))

