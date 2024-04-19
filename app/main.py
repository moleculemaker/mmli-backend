import logging
from contextlib import asynccontextmanager
from alembic.config import Config
from alembic import command


import uvicorn
from fastapi import FastAPI

from config import app_config, get_logger
from routers import chemscraper, job, files, novostoic, somn
from fastapi.middleware.cors import CORSMiddleware

from services.kubejob_service import KubeEventWatcher

from models.sqlmodel.models import Job
from models.sqlmodel.db import init_db

log = get_logger(__name__)


def run_migrations():
    alembic_cfg = Config("../alembic.ini")
    command.upgrade(alembic_cfg, "head")


watcher = KubeEventWatcher()


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger = get_logger('main:LifecycleEvents')
    logger.info("Starting up...")
    logger.info("Running alembic upgrade head...")
    run_migrations()
    logger.info("Starting KubeWatcher...")
    watcher.run()
    yield
    logger.info("Shutting down...")
    watcher.close()


app = FastAPI(lifespan=lifespan)

app.include_router(files.router)
app.include_router(job.router)
app.include_router(chemscraper.router)
app.include_router(novostoic.router)
app.include_router(somn.router)

origins = [
    "http://test.mydomain.com",
    "http://localhost:4200",
    "https://chemscraper.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://chemscraper.frontend.mmli1.ncsa.illinois.edu"
    "https://somn.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://somn.frontend.mmli1.ncsa.illinois.edu",
    "https://novostoic.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://novostoic.frontend.mmli1.ncsa.illinois.edu",
    # "http://another.allowed-origin.com", # Add more origins if needed
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
