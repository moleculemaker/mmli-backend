"""RFC 9457 problem details for the versioned API.

The legacy API returns FastAPI's `{"detail": ...}`, which tells a client that something
went wrong but not what kind of thing, and never says which field was at fault. That is
workable for a browser app whose developer can read the source; it is not workable for a
script or an agent, which has to decide what to do next from the response alone.

Every /v1 error is `application/problem+json` with a stable `type` URI, so a client can
branch on the failure kind without string-matching prose. Validation failures also carry
JSON Pointers into the submitted document, so "what exactly is wrong" is machine-readable.

These handlers are registered on the /v1 sub-application only. The legacy surface keeps
its existing shape.
"""
from typing import Any, Dict, List, Optional

from fastapi import Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from starlette.exceptions import HTTPException as StarletteHTTPException

PROBLEM_CONTENT_TYPE = 'application/problem+json'

# Problem type URIs are identifiers, not necessarily fetchable documents. Keeping them
# under a path this API owns means they can be made resolvable later without changing
# what clients already match on.
TYPE_BASE = 'https://mmli.fastapi.mmli2.ncsa.illinois.edu/v1/problems'

INVALID_INPUT = f'{TYPE_BASE}/invalid-input'
UNKNOWN_TOOL = f'{TYPE_BASE}/unknown-tool'
UNKNOWN_JOB = f'{TYPE_BASE}/unknown-job'
JOB_NOT_FINISHED = f'{TYPE_BASE}/job-not-finished'
JOB_FAILED = f'{TYPE_BASE}/job-failed'
NOT_CANCELABLE = f'{TYPE_BASE}/not-cancelable'
IDEMPOTENCY_CONFLICT = f'{TYPE_BASE}/idempotency-conflict'
UPSTREAM_FAILURE = f'{TYPE_BASE}/upstream-failure'
ABOUT_BLANK = 'about:blank'


def problem(
    status_code: int,
    title: str,
    detail: Optional[str] = None,
    type_uri: str = ABOUT_BLANK,
    errors: Optional[List[Dict[str, Any]]] = None,
    headers: Optional[Dict[str, str]] = None,
    **extra: Any,
) -> JSONResponse:
    """Build an RFC 9457 response."""
    body: Dict[str, Any] = {'type': type_uri, 'title': title, 'status': status_code}
    if detail:
        body['detail'] = detail
    if errors:
        body['errors'] = errors
    body.update(extra)
    return JSONResponse(
        status_code=status_code, content=body, media_type=PROBLEM_CONTENT_TYPE,
        headers=headers or {},
    )


class ProblemException(Exception):
    """Raised by /v1 handlers to produce a problem+json response.

    Deliberately not an HTTPException subclass: FastAPI's own handler would catch that
    and render it in the legacy shape.
    """

    def __init__(self, status_code: int, title: str, detail: Optional[str] = None,
                 type_uri: str = ABOUT_BLANK, errors: Optional[List[Dict[str, Any]]] = None,
                 headers: Optional[Dict[str, str]] = None, **extra: Any):
        super().__init__(title)
        self.status_code = status_code
        self.title = title
        self.detail = detail
        self.type_uri = type_uri
        self.errors = errors
        self.headers = headers
        self.extra = extra


async def problem_exception_handler(request: Request, exc: ProblemException) -> JSONResponse:
    return problem(exc.status_code, exc.title, exc.detail, exc.type_uri,
                   exc.errors, exc.headers, **exc.extra)


async def http_exception_handler(request: Request, exc: StarletteHTTPException) -> JSONResponse:
    """Render anything raised as a plain HTTPException inside /v1 as a problem.

    Mostly this catches 404s for unrouted paths and errors raised by shared code that
    predates this API, so that /v1 never emits two different error shapes.
    """
    return problem(exc.status_code, str(exc.detail), type_uri=ABOUT_BLANK)


async def validation_exception_handler(request: Request, exc: RequestValidationError) -> JSONResponse:
    """Convert FastAPI's request-validation errors into pointer-carrying problems."""
    errors = []
    for error in exc.errors():
        # loc is ("body", "field", 0, ...) -- drop the source segment and render the
        # rest as a JSON Pointer into the submitted document.
        location = [str(part) for part in error.get('loc', [])[1:]]
        errors.append({
            'pointer': '/' + '/'.join(location) if location else '',
            'detail': error.get('msg', 'invalid value'),
        })
    return problem(
        422, 'Request failed validation', type_uri=INVALID_INPUT, errors=errors,
    )
