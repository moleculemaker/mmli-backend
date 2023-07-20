import os

from dotenv import load_dotenv
from sqlmodel import create_engine, SQLModel, Session


# SQLALCHEMY_DATABASE_URL = "postgresql://user:password@postgresserver/db"
# By default, use SQLite (for local development)
load_dotenv()
DATABASE_URL = os.getenv("SQLALCHEMY_DATABASE_URL", "sqlite:///./sql_app.db")

engine = create_engine(DATABASE_URL, echo=True)


def init_db():
    SQLModel.metadata.create_all(engine)


def get_session():
    with Session(engine) as session:
        yield session
