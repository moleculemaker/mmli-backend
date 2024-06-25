#!/bin/env python3
import asyncio
import logging
from contextlib import asynccontextmanager


from fastapi import FastAPI

from config import app_config, get_logger
from routers import chemscraper, job, files, somn, novostoic, molli
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
    log.info("Starting KubeWatcher...")
    watcher.run()
    yield
    log.info("Shutting down...")
    watcher.close()


app = FastAPI(lifespan=lifespan)

app.include_router(files.router)
app.include_router(job.router)
app.include_router(chemscraper.router)
app.include_router(novostoic.router)
app.include_router(somn.router)
app.include_router(molli.router)

origins = [
    "http://test.mydomain.com",
    "http://localhost:4200",
    "http://127.0.0.1:4200",
    "https://chemscraper.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://chemscraper.frontend.mmli1.ncsa.illinois.edu",
    "https://novostoic.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://novostoic.frontend.mmli1.ncsa.illinois.edu",
    "https://somn.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://somn.frontend.mmli1.ncsa.illinois.edu",
    # "http://another.allowed-origin.com", # Add more origins if needed
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


