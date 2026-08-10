#!/bin/env python3
import asyncio
import logging
from contextlib import asynccontextmanager


from fastapi import FastAPI

from config import app_config, get_logger
from routers import chemscraper, job, files, somn, novostoic, molli, shared, reactionminer
from routers.v1 import create_v1_app
from fastapi.middleware.cors import CORSMiddleware

from services.kubejob_service import KubeEventWatcher

from models.sqlmodel.models import Job
from models.sqlmodel.db import init_db

log = get_logger(__name__)


watcher = KubeEventWatcher()
#asyncio.run(init_db())


@asynccontextmanager
async def lifespan(app_: FastAPI):
    global watcher
    global log
    log.info("Starting up...")
    # KubeEventWatcher starts its own daemon thread in __init__, so there is nothing to start
    # here. Do NOT call watcher.run() - it is the thread body (an endless watch loop), and
    # calling it inline would block the event loop and hang startup.
    log.info(f"KubeWatcher running: {watcher.is_alive()}")
    yield
    log.info("Shutting down...")
    watcher.close()


app = FastAPI(lifespan=lifespan)

# Mounted BEFORE the routers below, and that order is load-bearing. Starlette matches
# routes in registration order, and the legacy router declares
# /{job_type}/jobs/{job_id}, which happily matches /v1/jobs/<id> with job_type="v1".
# Registering the legacy routes first therefore shadows most of the versioned API and
# answers it with "Invalid job type: v1".
#
# Mounted rather than included as a router because /v1 needs its own CORS policy (any
# origin, no credentials) and its own error format (RFC 9457 problem+json); a
# sub-application is what gives it an independent middleware stack and exception
# handlers without disturbing the legacy surface. Its OpenAPI document is at
# /v1/openapi.json and its docs at /v1/docs.
#
# Kept as a module-level name, not an inline expression: a mounted sub-application has
# its own dependency-override registry, so tests need a handle on it to stub MinIO.
v1_app = create_v1_app()
app.mount("/v1", v1_app)

app.include_router(files.router)
app.include_router(job.router)
app.include_router(chemscraper.router)
app.include_router(reactionminer.router)
app.include_router(novostoic.router)
app.include_router(somn.router)
app.include_router(molli.router)
app.include_router(shared.router)

origins = [
    "http://test.mydomain.com",
    "http://localhost:4200",
    "http://127.0.0.1:4200",
    "https://chemscraper.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://chemscraper.frontend.mmli2.ncsa.illinois.edu",
    "https://novostoic.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://novostoic.frontend.mmli2.ncsa.illinois.edu",
    "https://somn.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://somn.frontend.mmli2.ncsa.illinois.edu",
    "https://frontend.staging.openenzymedb.mmli2.ncsa.illinois.edu",
    "https://openenzymedb.platform.moleculemaker.org",
    "https://reactionminer.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://reactionminer.frontend.mmli2.ncsa.illinois.edu",
    "https://ezspecificity.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://ezspecificity.frontend.mmli2.ncsa.illinois.edu",
    "https://crispr-copies.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://crispr-copies.frontend.mmli2.ncsa.illinois.edu",
    "https://mutagenesis.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://mutagenesis.frontend.mmli2.ncsa.illinois.edu",
    # "http://another.allowed-origin.com", # Add more origins if needed
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)



