"""Characterization tests for application startup.

These record a latent hazard so that it fails loudly rather than in production.
"""
import inspect
import pathlib

from fastapi import FastAPI

import main


def test_fastapi_silently_ignores_the_lifespan_argument():
    """main.py passes `lifespan=` to FastAPI, but this version does not support it.

    FastAPI gained lifespan support in 0.93.0. On the pinned 0.89.0 the argument is
    swallowed by **extra and discarded, so `main.lifespan` never executes -- neither its
    startup half nor its `watcher.close()` shutdown half.
    """
    assert "lifespan" not in inspect.signature(FastAPI.__init__).parameters


def test_the_custom_lifespan_is_not_wired_into_the_app():
    """Corollary of the above, asserted against the app object rather than the library."""
    lifespan_context = main.app.router.lifespan_context

    # Starlette's default, not main.lifespan
    assert getattr(lifespan_context, "__name__", None) != "lifespan"


def test_the_inert_lifespan_hides_a_call_that_would_hang_startup():
    """UPGRADE HAZARD, recorded here so that whoever upgrades FastAPI meets it.

    `main.lifespan` calls `watcher.run()`, which is the KubeEventWatcher thread body: an
    unbroken `while True` with no break or return. The watcher is already running -- its
    constructor starts the thread -- so if this lifespan ever became live it would run a
    second copy of that loop synchronously on the event loop and startup would never
    complete.

    It is harmless today only because FastAPI 0.89 discards the lifespan entirely. The
    moment FastAPI is upgraded past 0.93 this becomes a hang on boot, so the call has to
    be removed as part of that upgrade rather than after it.
    """
    assert "watcher.run()" in inspect.getsource(main.lifespan)


def test_the_watcher_thread_is_started_by_its_constructor():
    """Why the inert lifespan has gone unnoticed: startup does not depend on it.

    Asserted against the module's source text because conftest replaces the real
    __init__ with a stub, so inspecting the live class would only describe the stub.
    """
    from services import kubejob_service

    source = pathlib.Path(kubejob_service.__file__).read_text()

    assert "self.thread.start()" in source
