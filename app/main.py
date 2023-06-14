from fastapi import FastAPI
from routers import chemscraper

app = FastAPI()

app.include_router(chemscraper.router)
