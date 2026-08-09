"""Characterization tests for application startup.

Before the framework upgrade, FastAPI 0.89 did not accept a `lifespan` argument -- it
gained that in 0.93.0 -- so the one main.py passes was swallowed by **extra and
silently discarded. The application's lifespan handler had never run, in either its
startup or its shutdown half.

It runs now. These tests record that transition and guard the hazard it activates.
"""
import inspect
import pathlib

from fastapi import FastAPI
from fastapi.testclient import TestClient

import main


def test_fastapi_supports_lifespan():
    """Was: asserted the opposite, because 0.89 silently discarded the argument."""
    assert "lifespan" in inspect.signature(FastAPI.__init__).parameters


def test_the_custom_lifespan_actually_executes(monkeypatch):
    """Was: asserted the opposite, since 0.89 discarded the handler entirely.

    Checked by observing a side effect rather than by inspecting the callable's name:
    FastAPI wraps user lifespans (the attribute is `merged_lifespan`), so the name says
    nothing about whether ours is among them. main.lifespan calls watcher.close() on
    the way out, so spying on that proves both halves ran.
    """
    shutdown_ran = []
    monkeypatch.setattr(main.watcher, "close", lambda: shutdown_ran.append(True))

    with TestClient(main.app) as client:
        assert client.get("/openapi.json").status_code == 200
        assert shutdown_ran == []  # startup done, shutdown not yet

    assert shutdown_ran == [True]


def test_lifespan_does_not_call_the_watcher_thread_body():
    """The hazard this upgrade activates.

    `watcher.run()` is the KubeEventWatcher thread body: an unbroken `while True` with
    no break or return. The watcher is already running, started by its own constructor.
    Calling run() from the lifespan runs a second copy of that loop synchronously on
    the event loop, and startup never completes.

    Under the previous FastAPI this was inert and unreachable. It is now live, so this
    assertion is the thing standing between the service and a permanent hang on boot.

    Asserted against source text on purpose: conftest replaces the watcher with a stub,
    so test_the_custom_lifespan_actually_executes above would pass even if the call
    were reintroduced. A runtime check here would hang the suite rather than fail it.
    """
    lifespan_source = inspect.getsource(main.lifespan)

    assert "watcher.run()" not in lifespan_source.replace("Do NOT call watcher.run()", "")


def test_the_watcher_thread_is_started_by_its_constructor():
    """Why the lifespan being inert went unnoticed for so long: startup never needed it.

    Asserted against the module's source text because conftest replaces the real
    __init__ with a stub, so inspecting the live class would only describe the stub.
    """
    from services import kubejob_service

    source = pathlib.Path(kubejob_service.__file__).read_text()

    assert "self.thread.start()" in source
