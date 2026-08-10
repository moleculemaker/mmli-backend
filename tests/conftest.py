"""Test harness for the mmli-backend FastAPI app.

Importing the application has three side effects that make it untestable as-is, so all
three are neutralized here BEFORE any app module is imported:

1. `app/config.py` reads `cfg/config.yaml` relative to the working directory at import
   time, and derives the database URL from a secrets file that does not exist in CI.
2. `models/sqlmodel/db.py` builds a global engine from `SQLALCHEMY_DATABASE_URL` at
   import time, so the URL has to be in the environment before the import happens.
3. `main.py` constructs a `KubeEventWatcher()` at module scope, whose `__init__` opens a
   real database connection and starts a thread running an endless Kubernetes watch loop.

Order matters in this file. Environment setup and the watcher patch must both happen
above the application imports, which is why those imports are not at the top.
"""
import os
import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
APP_DIR = REPO_ROOT / "app"

# (1) + (2): configure the app before anything imports `config`.
#
# Defaults to a throwaway SQLite file so the suite needs no services. Set
# TEST_DATABASE_URL to an async Postgres URL to run the same tests against the driver
# production actually uses:
#
#   TEST_DATABASE_URL=postgresql+asyncpg://postgres:pw@host:5432/mmli pytest
#
# Worth doing when changing anything in the ORM layer: SQLite exercises neither asyncpg
# nor Postgres type handling, so it cannot catch a whole class of driver-level problem.
_db_file = REPO_ROOT / ".pytest-characterization.db"
_async_url = os.getenv("TEST_DATABASE_URL", f"sqlite+aiosqlite:///{_db_file}")

os.environ["CONFIG_FILEPATH"] = str(APP_DIR / "cfg" / "config.yaml")
os.environ["SECRET_FILEPATH"] = str(APP_DIR / "cfg" / "does-not-exist.yaml")
os.environ["SQLALCHEMY_DATABASE_URL"] = _async_url

# The schema reset below runs through a synchronous driver (see fresh_database), so the
# async URL needs its driver swapped for the blocking equivalent.
_sync_url = _async_url.replace("+aiosqlite", "").replace("+asyncpg", "+psycopg2")
os.environ.setdefault("MINIO_SERVER", "localhost:9000")
os.environ.setdefault("MINIO_ACCESS_KEY", "test")
os.environ.setdefault("MINIO_SECRET_KEY", "test")

if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

# (3): stub the watcher before `main` instantiates one.
from services import kubejob_service  # noqa: E402


def _noop(self, *args, **kwargs):
    return None


kubejob_service.KubeEventWatcher.__init__ = _noop
kubejob_service.KubeEventWatcher.run = _noop
kubejob_service.KubeEventWatcher.close = _noop
kubejob_service.KubeEventWatcher.is_alive = lambda self: False

import logging  # noqa: E402

import pytest  # noqa: E402
from fastapi.testclient import TestClient  # noqa: E402
from sqlalchemy.ext.asyncio import create_async_engine  # noqa: E402
from sqlalchemy.pool import NullPool  # noqa: E402
from sqlmodel import SQLModel, create_engine  # noqa: E402

import main  # noqa: E402
from models.sqlmodel import db as db_module  # noqa: E402
from services.minio_service import MinIOService  # noqa: E402

# Rebuild the app's engine without connection pooling, for tests only.
#
# TestClient runs the ASGI app on a fresh event loop per instantiation, and this suite
# builds one per test. asyncpg connections are bound to the loop that opened them, so a
# pooled connection handed to a later test raises "attached to a different loop". The
# production process has a single long-lived loop and is unaffected, so this is a
# property of the harness rather than of the application. NullPool opens and closes a
# connection per checkout, which removes the cross-loop reuse entirely.
db_module.engine = create_async_engine(_async_url, poolclass=NullPool)

# db.py builds its engine with echo=True. SQLAlchemy implements echo by setting the
# logger's level when the engine is constructed, so this has to run after that import
# rather than before it, or the engine simply raises the level back to INFO and every
# test emits a full SQL log.
logging.getLogger("sqlalchemy.engine").setLevel(logging.WARNING)


class FakeMinIO:
    """In-memory stand-in for MinIOService.

    Only the methods the routers actually call are implemented; anything else should
    fail loudly rather than silently succeed, so unimplemented calls raise AttributeError
    in the normal way.
    """

    def __init__(self):
        # {bucket: {object_name: bytes}}
        self.objects = {}

    # -- helpers for tests -------------------------------------------------
    def put(self, bucket, object_name, data: bytes):
        self.objects.setdefault(bucket, {})[object_name.lstrip("/")] = data

    # -- MinIOService surface ----------------------------------------------
    def ensure_bucket_exists(self, bucket_name):
        self.objects.setdefault(bucket_name, {})
        return True

    def upload_file(self, bucket_name, object_name, data):
        self.put(bucket_name, object_name, data)
        return True

    def get_file(self, bucket_name, object_name):
        return self.objects.get(bucket_name, {}).get(object_name.lstrip("/"))

    def get_file_urls(self, bucket_name, prefix):
        names = [n for n in self.objects.get(bucket_name, {}) if n.startswith(prefix.lstrip("/"))]
        return [f"http://minio.test/{bucket_name}/{n}" for n in names] or None

    def list_files(self, bucket_name, prefix="", recursive=True):
        return [n for n in self.objects.get(bucket_name, {}) if n.startswith(prefix.lstrip("/"))]


@pytest.fixture(autouse=True)
def fresh_database():
    """Recreate every table before each test so ordering cannot leak state.

    Deliberately uses a synchronous engine over the same file rather than the app's
    async one. TestClient runs the ASGI app on its own event loop in a worker thread,
    and aiosqlite connections are bound to the loop that opened them, so driving the
    async engine from the test's loop risks cross-loop errors that surface as flakes.
    DDL through a sync connection sidesteps that entirely.
    """
    sync_engine = create_engine(_sync_url)
    SQLModel.metadata.drop_all(sync_engine)
    SQLModel.metadata.create_all(sync_engine)
    sync_engine.dispose()
    yield


@pytest.fixture
def sync_db_url():
    """Blocking URL for the test database.

    KubeEventWatcher uses a synchronous engine, so tests that drive it directly need
    this rather than the async URL the app is configured with.
    """
    return _sync_url


@pytest.fixture
def fake_minio():
    return FakeMinIO()


@pytest.fixture
def created_k8s_jobs(monkeypatch):
    """Record calls to kubejob_service.create_job instead of talking to Kubernetes.

    Returns the list of recorded kwargs so a test can assert on what would have been
    submitted to the cluster.
    """
    calls = []

    def _record(**kwargs):
        calls.append(kwargs)
        return {"metadata": {"name": f"mmli-job-{kwargs.get('job_type')}-{kwargs.get('job_id')}"}}

    monkeypatch.setattr(kubejob_service, "create_job", _record)
    return calls


@pytest.fixture
def deleted_k8s_jobs(monkeypatch):
    """Record calls to kubejob_service.delete_job instead of talking to Kubernetes."""
    calls = []

    def _record(job_type=None, job_id=None):
        calls.append({"job_type": job_type, "job_id": job_id})

    monkeypatch.setattr(kubejob_service, "delete_job", _record)
    return calls


@pytest.fixture
def client(fake_minio, created_k8s_jobs, deleted_k8s_jobs):
    main.app.dependency_overrides[MinIOService] = lambda: fake_minio
    with TestClient(main.app) as test_client:
        test_client.minio = fake_minio
        test_client.k8s_jobs = created_k8s_jobs
        yield test_client
    main.app.dependency_overrides.clear()
