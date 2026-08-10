"""Characterization tests for application startup.

These record an unusual situation: the application's lifespan handler is written, is
correct, and has never executed.
"""
import inspect
import pathlib

from fastapi import FastAPI

import main


def test_fastapi_silently_ignores_the_lifespan_argument():
    """main.py passes `lifespan=` to FastAPI, but this version does not support it.

    FastAPI gained lifespan support in 0.93.0. On the pinned 0.89.0 the argument is
    swallowed by **extra and discarded, so neither the startup half of `main.lifespan`
    nor its `watcher.close()` shutdown half has ever run.

    The KubeEventWatcher works regardless, because its constructor starts its own
    daemon thread. Nothing depends on the lifespan today.
    """
    assert "lifespan" not in inspect.signature(FastAPI.__init__).parameters


def test_the_custom_lifespan_is_not_wired_into_the_app():
    """Corollary of the above, asserted against the app object rather than the library."""
    lifespan_context = main.app.router.lifespan_context

    # Starlette's default, not main.lifespan
    assert getattr(lifespan_context, "__name__", None) != "lifespan"


def test_lifespan_does_not_call_the_watcher_thread_body():
    """Regression guard on the fix that makes the lifespan safe to activate.

    `watcher.run()` is the KubeEventWatcher thread body: an unbroken `while True` with
    no break or return. Calling it from the lifespan would run a second copy of the
    loop synchronously on the event loop, and startup would never complete.

    That call was removed, but because the lifespan is inert (see above) the fix is
    currently dormant -- no test or deployment exercises it. It becomes load-bearing
    the moment FastAPI is upgraded past 0.93 and the lifespan goes live. This test
    exists so the call cannot be reintroduced in the meantime, when doing so would
    appear harmless.

    Asserted against source text because the lifespan is never executed, so there is
    no runtime behavior to observe.
    """
    lifespan_source = inspect.getsource(main.lifespan)

    assert "watcher.run()" not in lifespan_source.replace("Do NOT call watcher.run()", "")


def test_the_watcher_thread_is_started_by_its_constructor():
    """Why the inert lifespan has gone unnoticed: startup does not depend on it.

    Asserted against the module's source text because conftest replaces the real
    __init__ with a stub, so inspecting the live class would only describe the stub.
    """
    from services import kubejob_service

    source = pathlib.Path(kubejob_service.__file__).read_text()

    assert "self.thread.start()" in source
