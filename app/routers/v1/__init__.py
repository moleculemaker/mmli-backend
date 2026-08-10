"""The versioned, self-describing API.

Mounted as a sub-application rather than added as a router, because /v1 needs its own
middleware and error handling and the legacy surface must keep its current behavior
exactly. A sub-application gets an independent middleware stack, its own exception
handlers, and its own OpenAPI document at /v1/openapi.json.

The two differ deliberately in two ways:

  CORS      the legacy app allows a fixed list of frontend origins with credentials.
            /v1 allows any origin without credentials, which is the combination the
            CORS spec permits and which lets a notebook or another institution's page
            call it. Since /v1 has no cookie authentication, there is no credential to
            expose by doing so.
  errors    the legacy app returns {"detail": ...}; /v1 returns RFC 9457 problem+json.
"""
from typing import Any, Dict

from fastapi import FastAPI
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.openapi.utils import get_openapi
from starlette.exceptions import HTTPException as StarletteHTTPException

from routers.v1 import jobs, service_info, tools
from routers.v1.problems import (
    ProblemException, http_exception_handler, problem_exception_handler,
    validation_exception_handler,
)

DESCRIPTION = """
A self-describing API for running the AlphaSynthesis computational tools.

Start at `/v1/tools` to see what is available, `/v1/tools/{tool}/input-schema` for what
a tool accepts, then `POST /v1/tools/{tool}/jobs` to run it. Job responses carry a
`links` object, so clients need never construct a URL.

No authentication. A job is reachable by anyone holding its id, which is an unguessable
UUIDv4; see `/v1/service-info` for the access model.
"""


def create_v1_app() -> FastAPI:
    app = FastAPI(
        title='AlphaSynthesis Tools API',
        version='1.0.0',
        description=DESCRIPTION,
        docs_url='/docs',
        openapi_url='/openapi.json',
    )

    # '*' with credentials is forbidden by the CORS spec, and rightly: any origin being
    # able to send a user's cookies is exactly the attack the rule prevents. Splitting
    # this from the legacy app's credentialed allowlist is what makes '*' safe here.
    app.add_middleware(
        CORSMiddleware,
        allow_origins=['*'],
        allow_credentials=False,
        allow_methods=['GET', 'POST', 'OPTIONS'],
        allow_headers=['*'],
        expose_headers=['Retry-After'],
    )

    app.add_exception_handler(ProblemException, problem_exception_handler)
    app.add_exception_handler(StarletteHTTPException, http_exception_handler)
    app.add_exception_handler(RequestValidationError, validation_exception_handler)

    app.include_router(service_info.router)
    app.include_router(tools.router)
    app.include_router(jobs.router)
    app.openapi = lambda: _openapi(app)
    return app


# Description shown for the JSON body on submission. The real schema is per tool and
# lives at /v1/tools/{tool}/input-schema, so nothing more specific can be stated here.
JSON_BODY_SCHEMA: Dict[str, Any] = {
    'title': 'Tool input document',
    'description': (
        "The tool's input document, exactly as published at "
        "/v1/tools/{tool}/input-schema. It is validated against that schema before the "
        "job starts, and a rejection names the offending field. No wrapper object: the "
        "body is the document itself.\n\n"
        "Most tools take an object; a few, such as somn, take an array. The example "
        "below is for novostoic-optstoic -- fetch the tool's own schema for anything else."
    ),
    # oneOf rather than an unconstrained schema: without a type, Swagger UI renders the
    # example as the bare string "string", which suggests the body is a JSON string
    # rather than a document. This cannot be narrowed further because the real schema is
    # per tool.
    'oneOf': [
        {'type': 'object', 'title': 'Object input'},
        {'type': 'array', 'title': 'Array input', 'items': {'type': 'object'}},
    ],
    'example': {'primary_precursor': 'MNXM1137670', 'target_molecule': 'MNXM26'},
}


def _openapi(app: FastAPI) -> Dict[str, Any]:
    """Generate the OpenAPI document, then fix how submission bodies are presented.

    Submission accepts either a JSON document or multipart with files, but only the
    multipart form is declared through function parameters, so that is all FastAPI
    infers. Left alone, the docs offer multipart as the only -- and therefore the
    default -- option, which makes the simple case look like the complicated one: a
    reader sees an `inputs` string field and file uploads, with no indication that
    posting the input document as JSON is both supported and easier.

    So the JSON alternative is added here, and listed first. Ordering is the whole point:
    Swagger UI selects whichever content type comes first, and dict insertion order
    survives into the served JSON. `openapi_extra` on the route cannot do this -- FastAPI
    deep-merges it into the generated operation, so a key added that way is appended
    after the multipart entry rather than placed before it.
    """
    if app.openapi_schema:
        return app.openapi_schema

    schema = get_openapi(
        title=app.title, version=app.version, description=app.description,
        routes=app.routes,
        # Required, and easy to lose by overriding openapi(): this app is mounted at
        # /v1, and FastAPI's own /openapi.json route records that by inserting the mount
        # path into app.servers just before calling this. Dropping it makes Swagger UI
        # issue every "Try it out" against the unmounted path, so the docs answer 404
        # for every request a reader makes.
        servers=app.servers,
    )

    for operations in schema.get('paths', {}).values():
        for operation in operations.values():
            content = (operation.get('requestBody') or {}).get('content')
            if not content or 'multipart/form-data' not in content:
                continue
            content.setdefault('application/json', {'schema': JSON_BODY_SCHEMA})
            operation['requestBody']['content'] = {
                key: content[key]
                for key in sorted(content, key=lambda k: k != 'application/json')
            }

    app.openapi_schema = schema
    return schema
