import os

import uvicorn
from fastapi import FastAPI

from routers import chemscraper, job, files
from fastapi.middleware.cors import CORSMiddleware

from models.sqlmodel.models import Job
from models.sqlmodel.db import init_db

app = FastAPI()

DEBUG = 't' in str.lower(os.getenv("DEBUG", "true"))


@app.on_event("startup")
async def on_startup():
    await init_db()


app.include_router(files.router)
app.include_router(job.router)
app.include_router(chemscraper.router)

origins = [
    "http://test.mydomain.com",
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


if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8080, reload=True)
