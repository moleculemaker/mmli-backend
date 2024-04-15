import logging
from contextlib import asynccontextmanager
from alembic.config import Config
from alembic import command


import uvicorn
from fastapi import FastAPI

from config import app_config
from routers import chemscraper, job, files
from fastapi.middleware.cors import CORSMiddleware

from models.sqlmodel.models import Job
from models.sqlmodel.db import init_db

log = logging.getLogger(__name__)

app = FastAPI()

app.include_router(files.router)
app.include_router(job.router)
app.include_router(chemscraper.router)

origins = [
    "http://test.mydomain.com",
    "http://localhost:4200",
    "https://chemscraper.frontend.staging.mmli1.ncsa.illinois.edu",
    "https://chemscraper.frontend.mmli1.ncsa.illinois.edu"
    # "http://another.allowed-origin.com", # Add more origins if needed
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


def run_migrations():
    alembic_cfg = Config("alembic.ini")
    command.upgrade(alembic_cfg, "head")


@asynccontextmanager
async def lifespan(app_: FastAPI):
    log.info("Starting up...")
    log.info("Running alembic upgrade head...")
    run_migrations()
    yield
    log.info("Shutting down...")


if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=app_config['server']['port'], reload=True)
