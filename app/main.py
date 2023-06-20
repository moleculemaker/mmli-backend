from fastapi import FastAPI
from routers import chemscraper
from fastapi.middleware.cors import CORSMiddleware


app = FastAPI()

app.include_router(chemscraper.router)

origins = [
    "http://test.mydomain.com",
    # "http://another.allowed-origin.com", # Add more origins if needed
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
