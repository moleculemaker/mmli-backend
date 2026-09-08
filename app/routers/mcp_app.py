"""Model Context Protocol server.

Exposes the same tools to AI agents that /v1 exposes to scripts. An agent connected
here sees one typed tool per scientific tool, with that tool's real JSON Schema as its
inputSchema, so it can construct a valid call from the schema rather than from prose.

Implemented as an adapter over /v1 rather than as a second implementation. Every handler
below issues an in-process request against the versioned API through an ASGI transport,
so validation, error shapes, idempotency and result semantics are by construction the
same ones an HTTP client gets. The alternative -- calling the services directly -- would
have produced a second code path that drifts from the first the moment either changes.

Tool surface:

  submit_<tool>       one per tool, carrying that tool's published input schema
  get_job_status      poll a job
  get_job_results     fetch results once finished
  list_job_artifacts  raw output files
  cancel_job          stop a running job
  describe_tool       full metadata for one tool, including citation and license

Rate limiting is enforced here rather than at the ingress. MCP multiplexes every call
over a single HTTP request stream, so a path-and-method rule at the edge cannot see
individual tool invocations; from outside, one submission and a hundred look identical.
"""
import json
import time
from contextlib import asynccontextmanager
from collections import defaultdict, deque
from contextvars import ContextVar
from typing import Any, Dict, List, Optional, Tuple

import httpx
import mcp.types as types
from starlette.responses import Response
from mcp.server.lowlevel import Server
from mcp.server.streamable_http_manager import StreamableHTTPSessionManager

from config import get_logger
from services import tool_registry
from services.analytics import SURFACE_MCP, SURFACE_HEADER, SURFACE_TOKEN, SURFACE_TOKEN_HEADER

log = get_logger(__name__)

SUBMIT_PREFIX = 'submit_'

# Mirrors the ingress limit applied to /v1 submissions. In-process and therefore
# per-replica: with N replicas the effective ceiling is N times this. Shared state would
# be needed to do better, which this service does not have. Documented rather than
# quietly approximated.
SUBMIT_LIMIT = 10
SUBMIT_WINDOW_SECONDS = 60

_submit_history: Dict[str, deque] = defaultdict(deque)

# Set per request from the ASGI scope, read by the rate limiter. A context variable
# rather than an argument because MCP tool handlers are invoked by the SDK and never see
# the underlying request.
_current_client: ContextVar[str] = ContextVar('mcp_client_key', default='unknown')


def client_key_from_scope(scope) -> str:
    """Identify the calling client well enough to give it its own budget.

    Honors the proxy header the ingress sets, falling back to the peer address. Both are
    caller-controlled to some degree, so this bounds accidental starvation between
    clients rather than resisting a determined evader -- the same standing this has at
    the ingress, where a forged X-Forwarded-For evades the limit too.
    """
    headers = {k.lower(): v for k, v in (scope.get('headers') or [])}
    forwarded = headers.get(b'x-forwarded-for')
    if forwarded:
        return forwarded.decode('latin-1').split(',')[0].strip()
    client = scope.get('client')
    return client[0] if client else 'unknown'


def _rate_limited(client_key: str) -> bool:
    """Sliding-window check for submissions from one client."""
    now = time.monotonic()
    history = _submit_history[client_key]
    while history and now - history[0] > SUBMIT_WINDOW_SECONDS:
        history.popleft()
    if len(history) >= SUBMIT_LIMIT:
        return True
    history.append(now)
    return False


def _submit_tool_name(identifier: str) -> str:
    """MCP tool names allow a restricted character set; job types use hyphens."""
    return SUBMIT_PREFIX + identifier.replace('-', '_')


def _identifier_from_submit_tool(name: str) -> Optional[str]:
    if not name.startswith(SUBMIT_PREFIX):
        return None
    stem = name[len(SUBMIT_PREFIX):]
    for identifier in tool_registry.list_tools():
        if identifier.replace('-', '_') == stem:
            return identifier
    return None


def _submit_description(identifier: str, entry: Dict[str, Any]) -> str:
    """What an agent reads when deciding whether this is the tool it wants."""
    parts = [entry.get('description') or entry.get('summary') or identifier]
    if entry.get('input_files'):
        parts.append(
            f"Requires files to be uploaded before submission ({entry['input_files']}), "
            f"which this interface cannot do; use the HTTP API for this tool."
        )
    parts.append(
        'Returns a job_id immediately; the job runs asynchronously. Poll '
        'get_job_status until it reports "completed", then call get_job_results.'
    )
    return ' '.join(parts)


def _job_id_schema(description: str) -> Dict[str, Any]:
    return {
        '$schema': 'https://json-schema.org/draft/2020-12/schema',
        'type': 'object',
        'required': ['job_id'],
        'additionalProperties': False,
        'properties': {'job_id': {'type': 'string', 'description': description}},
    }


def build_tool_list() -> List[types.Tool]:
    """One MCP tool per submittable scientific tool, plus the job-lifecycle tools."""
    tools: List[types.Tool] = []

    for identifier, entry in sorted(tool_registry.list_tools().items()):
        schema = entry.get('input_schema')
        if schema is None:
            # Nothing to describe, so an agent could not construct a valid call.
            continue

        # MCP requires an object at the top level of inputSchema. Not every published
        # schema is one: somn's is an array, and cleandb-mepesm's is an anyOf that
        # accepts either an object or a bare string, so it has no top-level type at all.
        # Those are nested under a named property rather than reshaped, so the published
        # schema stays the authority on what the tool accepts.
        if schema.get('type') == 'object':
            input_schema = dict(schema)
            input_schema.pop('$schema', None)
        else:
            input_schema = {
                'type': 'object',
                'required': ['inputs'],
                'properties': {
                    'inputs': {**{k: v for k, v in schema.items() if k != '$schema'},
                               'description': 'The tool input document.'},
                },
            }

        tools.append(types.Tool(
            name=_submit_tool_name(identifier),
            title=entry.get('name') or identifier,
            description=_submit_description(identifier, entry),
            inputSchema=input_schema,
        ))

    tools.extend([
        types.Tool(
            name='get_job_status',
            title='Get job status',
            description=(
                'Current status of a job: queued, processing, completed, error or '
                'canceled. Also returns provenance, including the image digest that '
                'actually ran once the cluster reports it.'
            ),
            inputSchema=_job_id_schema('Job id returned by a submit_* tool.'),
        ),
        types.Tool(
            name='get_job_results',
            title='Get job results',
            description=(
                'Results of a finished job. Reports an error while the job is still '
                'running, and a distinct one if the job failed or produced no output, '
                'so "not ready" is never confused with "nothing to show".'
            ),
            inputSchema=_job_id_schema('Job id of a completed job.'),
        ),
        types.Tool(
            name='list_job_artifacts',
            title='List job output files',
            description=(
                'Raw files a job wrote. Logs and status markers are excluded unless '
                'include_logs is true.'
            ),
            inputSchema={
                '$schema': 'https://json-schema.org/draft/2020-12/schema',
                'type': 'object',
                'required': ['job_id'],
                'additionalProperties': False,
                'properties': {
                    'job_id': {'type': 'string', 'description': 'Job id to inspect.'},
                    'include_logs': {
                        'type': 'boolean', 'default': False,
                        'description': 'Include stdout/stderr logs and status markers.',
                    },
                },
            },
        ),
        types.Tool(
            name='cancel_job',
            title='Cancel a job',
            description='Stop a queued or running job. Jobs that already finished cannot be canceled.',
            inputSchema=_job_id_schema('Job id to cancel.'),
        ),
        types.Tool(
            name='describe_tool',
            title='Describe a tool',
            description=(
                'Full metadata for one tool: what it does, its source repository, the '
                'image it runs, and its license and citation where those are known.'
            ),
            inputSchema={
                '$schema': 'https://json-schema.org/draft/2020-12/schema',
                'type': 'object',
                'required': ['tool'],
                'additionalProperties': False,
                'properties': {
                    'tool': {
                        'type': 'string',
                        'enum': sorted(tool_registry.list_tools()),
                        'description': 'Tool identifier.',
                    },
                },
            },
        ),
    ])
    return tools


def _text(payload: Any) -> List[types.ContentBlock]:
    return [types.TextContent(type='text', text=json.dumps(payload, indent=2))]


def create_mcp_app(v1_app: Any, client_key_provider=None):
    """Build the MCP ASGI application.

    v1_app is the versioned FastAPI sub-application. Requests are issued against it
    in-process, so this never leaves the pod and needs no knowledge of the external URL.
    """
    server = Server('alphasynthesis-tools')

    # Default to the per-request client key. Previously this fell back to the constant
    # 'mcp', which gave every caller in the world a single shared budget: one busy agent
    # could exhaust the allowance for all of them.
    if client_key_provider is None:
        client_key_provider = _current_client.get

    def _client() -> httpx.AsyncClient:
        return httpx.AsyncClient(
            transport=httpx.ASGITransport(app=v1_app),
            base_url='http://mcp.internal',
            timeout=60.0,
            # Declare the surface so submissions are recorded as agent traffic rather
            # than as whatever user agent this client happens to send -- which, without
            # this, is httpx's default and indistinguishable from an external script.
            # The token proves the claim came from in-process code.
            headers={SURFACE_HEADER: SURFACE_MCP, SURFACE_TOKEN_HEADER: SURFACE_TOKEN},
        )

    @server.list_tools()
    async def list_tools() -> List[types.Tool]:
        return build_tool_list()

    @server.call_tool()
    async def call_tool(name: str, arguments: Dict[str, Any]) -> List[types.ContentBlock]:
        arguments = arguments or {}

        identifier = _identifier_from_submit_tool(name)
        if identifier is not None:
            return await _submit(identifier, arguments)

        if name == 'describe_tool':
            return await _get(f"/tools/{arguments['tool']}")
        if name == 'get_job_status':
            return await _get(f"/jobs/{arguments['job_id']}")
        if name == 'get_job_results':
            return await _get(f"/jobs/{arguments['job_id']}/results")
        if name == 'list_job_artifacts':
            suffix = '?include=logs' if arguments.get('include_logs') else ''
            return await _get(f"/jobs/{arguments['job_id']}/artifacts{suffix}")
        if name == 'cancel_job':
            return await _post(f"/jobs/{arguments['job_id']}/cancel", None)

        raise ValueError(f'Unknown tool: {name}')

    async def _submit(identifier: str, arguments: Dict[str, Any]) -> List[types.ContentBlock]:
        key = client_key_provider() if client_key_provider else 'mcp'
        if _rate_limited(key):
            raise ValueError(
                f'Rate limit exceeded: at most {SUBMIT_LIMIT} submissions per '
                f'{SUBMIT_WINDOW_SECONDS} seconds. Poll get_job_status on the jobs you '
                f'have already started rather than resubmitting.'
            )

        entry = tool_registry.get_tool(identifier) or {}
        schema = entry.get('input_schema') or {}
        # Undo the nesting applied above to schemas that are not objects at the root.
        if schema.get('type') == 'object':
            payload = arguments
        else:
            payload = arguments.get('inputs')

        return await _post(f'/tools/{identifier}/jobs', payload)

    async def _get(path: str) -> List[types.ContentBlock]:
        async with _client() as client:
            response = await client.get(path)
        return _render(response)

    async def _post(path: str, payload: Any) -> List[types.ContentBlock]:
        async with _client() as client:
            response = await client.post(path, json=payload)
        return _render(response)

    def _render(response: httpx.Response) -> List[types.ContentBlock]:
        try:
            body = response.json()
        except ValueError:
            body = {'raw': response.text}

        if response.is_success:
            return _text(body)

        # Surface the problem document rather than a bare status code: it names the
        # failure kind and, for validation errors, points at the offending field, which
        # is what an agent needs in order to correct its own call.
        raise ValueError(json.dumps(body))

    return MCPMount(server)


class MCPMount:
    """Owns the session manager's lifecycle and the ASGI entry point.

    A StreamableHTTPSessionManager can only be entered once, so one is created per
    lifespan rather than once per process. Without that, anything which starts the
    application a second time in the same interpreter -- a test suite, most obviously --
    fails on the second start with an error that says nothing about the real cause.
    """

    def __init__(self, server: Server):
        self._server = server
        self._manager: Optional[StreamableHTTPSessionManager] = None

    @asynccontextmanager
    async def run(self):
        self._manager = StreamableHTTPSessionManager(
            app=self._server,
            # Stateless: no session affinity, so this survives running behind a load
            # balancer with several replicas. json_response avoids requiring SSE.
            stateless=True,
            json_response=True,
        )
        try:
            async with self._manager.run():
                yield
        finally:
            self._manager = None

    async def __call__(self, scope, receive, send):
        """ASGI entry point.

        Implemented on the instance so the same object can be registered both as an
        exact route and as a mount; see main.py for why both are needed. Also the only
        place with access to the request, so it is where the caller's identity is
        recorded for the rate limiter.
        """
        token = _current_client.set(client_key_from_scope(scope))
        try:
            await self.asgi_app(scope, receive, send)
        finally:
            _current_client.reset(token)

    async def asgi_app(self, scope, receive, send):
        if self._manager is None:
            # Reached only if something routes here outside the application lifespan.
            response = Response('MCP server is not running', status_code=503)
            await response(scope, receive, send)
            return
        await self._manager.handle_request(scope, receive, send)
