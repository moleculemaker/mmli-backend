#!/bin/env python3
import asyncio
import logging
from contextlib import asynccontextmanager


from fastapi import FastAPI

from config import app_config, get_logger
from routers import chemscraper, job, files, somn, novostoic, molli, shared, reactionminer
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
    # KubeEventWatcher starts its own daemon thread in __init__, so there is nothing to
    # start here. Do NOT call watcher.run(): it is the thread body, an endless watch
    # loop, and calling it inline blocks the event loop so startup never completes.
    #
    # This mattered only in theory until now. FastAPI 0.89 discarded the lifespan
    # argument entirely, so this handler never ran; the upgrade in this commit makes it
    # live, which turns a dormant hazard into a hang on boot. The open PR that fixes the
    # stranded-job bug removes the same call, so expect a small conflict if both land.
    log.info(f"KubeWatcher running: {watcher.is_alive()}")
    yield
    log.info("Shutting down...")
    watcher.close()


app = FastAPI(lifespan=lifespan)

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
    "https://chemscraper.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://chemscraper.frontend.mmli1.ncsa.illinois.edu",
    "https://chemscraper.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://chemscraper.frontend.mmli2.ncsa.illinois.edu",
    "https://novostoic.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://novostoic.frontend.mmli1.ncsa.illinois.edu",
    "https://novostoic.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://novostoic.frontend.mmli2.ncsa.illinois.edu",
    "https://somn.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://somn.frontend.mmli1.ncsa.illinois.edu",
    "https://somn.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://somn.frontend.mmli2.ncsa.illinois.edu",
    "https://frontend.staging.openenzymedb.mmli1.ncsa.illinois.edu",
    "https://frontend.staging.openenzymedb.mmli2.ncsa.illinois.edu",
    "https://openenzymedb.platform.moleculemaker.org",
    "https://reactionminer.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://reactionminer.frontend.mmli1.ncsa.illinois.edu",
    "https://reactionminer.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://reactionminer.frontend.mmli2.ncsa.illinois.edu",
    "https://ezspecificity.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://ezspecificity.frontend.mmli1.ncsa.illinois.edu",
    "https://ezspecificity.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://ezspecificity.frontend.mmli2.ncsa.illinois.edu",
    "https://crispr-copies.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://crispr-copies.frontend.mmli2.ncsa.illinois.edu",
    "https://mutagenesis.frontend.staging.mmli2.ncsa.illinois.edu",
    "https://mutagenesis.frontend.mmli2.ncsa.illinois.edu",
    # The hostnames the charts/config actually serve CRISPR-Copies and Mutagenesis at.
    # These tools use the `frontend.staging.<tool>.<cluster>` shape on staging (see
    # chart/values.staging.yaml, chart/values.mmli2.staging.yaml, app/cfg/config.yaml)
    # and ibiofoundry.illinois.edu on prod (chart/values.prod.yaml,
    # chart/values.mmli2.prod.yaml) - NOT the `<tool>.frontend.<cluster>` shape above.
    # Without these, both frontends are CORS-blocked in every environment.
    "https://frontend.staging.crispr-copies.mmli1.ncsa.illinois.edu",
    "https://frontend.staging.crispr-copies.mmli2.ncsa.illinois.edu",
    "https://crispr-copies.platform.ibiofoundry.illinois.edu",
    "https://frontend.staging.mutagenesis.mmli1.ncsa.illinois.edu",
    "https://frontend.staging.mutagenesis.mmli2.ncsa.illinois.edu",
    "https://mutagenesis.platform.ibiofoundry.illinois.edu",
    # "http://another.allowed-origin.com", # Add more origins if needed
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


