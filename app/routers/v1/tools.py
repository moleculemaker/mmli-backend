"""Tool discovery: what exists, what it does, and what it accepts.

These are the endpoints that make the API self-describing. A client that knows only the
base URL can list the tools, read what each one is for, fetch a schema describing its
inputs, and validate a payload before spending a submission against the rate limit.

Descriptors are schema.org JSON-LD so that a registry crawler can ingest them without
bespoke parsing. Fields nobody has supplied are emitted as null rather than omitted:
an absent key is ambiguous between "not applicable" and "nobody has told us", and for
license and citation the second is the truth today.
"""
from typing import Any, Dict

from fastapi import APIRouter, Request

from config import app_config, get_logger
from routers.v1.problems import UNKNOWN_TOOL, ProblemException
from services import tool_registry

router = APIRouter()
log = get_logger(__name__)

SCHEMA_ORG = 'https://schema.org'


def _base_url(request: Request) -> str:
    """External base URL for /v1, used to build absolute ids and links."""
    return str(request.base_url).rstrip('/') + '/v1'


def _descriptor(identifier: str, entry: Dict[str, Any], base: str) -> Dict[str, Any]:
    """Render one registry entry as a schema.org SoftwareApplication."""
    runtime = app_config['kubernetes_jobs'].get(identifier, {})
    has_schema = entry.get('input_schema') is not None

    return {
        '@context': SCHEMA_ORG,
        '@type': 'SoftwareApplication',
        '@id': f'{base}/tools/{identifier}',
        'identifier': identifier,
        'name': entry.get('name') or identifier,
        'description': entry.get('description'),
        'abstract': entry.get('summary'),
        'applicationCategory': 'Computational Science',
        'codeRepository': entry.get('repository'),
        # The image reference as configured. What a given run actually used is on that
        # job's provenance block, because this can be an untagged or mutable reference.
        'softwareVersion': runtime.get('image'),
        'license': entry.get('license'),
        'citation': entry.get('citation'),
        # EDAM terms are the interoperability vocabulary. Null until an SME supplies them.
        'edam': entry.get('edam', {'operation': None, 'topic': None}),
        'x-execution': entry.get('execution'),
        'x-input-files': entry.get('input_files'),
        # How much to trust the input schema, since they were derived several ways.
        'x-schema-status': entry.get('schema_status'),
        'inputSchema': f'{base}/tools/{identifier}/input-schema' if has_schema else None,
        'links': {
            'self': f'{base}/tools/{identifier}',
            'input_schema': f'{base}/tools/{identifier}/input-schema' if has_schema else None,
            'submit': f'{base}/tools/{identifier}/jobs',
        },
    }


@router.get('/tools', tags=['Discovery'], summary='List every available tool')
async def list_tools(request: Request) -> Dict[str, Any]:
    """Return the catalog.

    Subjob types created by a parent job's coordinator are omitted: they cannot be
    submitted directly, so listing them as runnable would be misleading.
    """
    base = _base_url(request)
    registry = tool_registry.list_tools()
    return {
        '@context': SCHEMA_ORG,
        '@type': 'ItemList',
        'numberOfItems': len(registry),
        'itemListElement': [
            _descriptor(identifier, entry, base)
            for identifier, entry in sorted(registry.items())
        ],
    }


@router.get('/tools/{tool}', tags=['Discovery'], summary='Describe one tool')
async def get_tool(tool: str, request: Request) -> Dict[str, Any]:
    entry = tool_registry.get_tool(tool)
    if entry is None:
        raise ProblemException(
            404, 'Unknown tool', f'No tool is registered with identifier "{tool}".',
            type_uri=UNKNOWN_TOOL,
            known_tools=sorted(tool_registry.list_tools()),
        )
    return _descriptor(tool, entry, _base_url(request))


@router.get('/tools/{tool}/input-schema', tags=['Discovery'],
            summary='JSON Schema describing this tool\'s inputs')
async def get_input_schema(tool: str, request: Request) -> Dict[str, Any]:
    """Return the tool's JSON Schema 2020-12 document.

    Self-contained by construction: every $ref resolves within the document, so it can
    be handed straight to a validator or embedded in a tool definition.
    """
    entry = tool_registry.get_tool(tool)
    if entry is None:
        raise ProblemException(
            404, 'Unknown tool', f'No tool is registered with identifier "{tool}".',
            type_uri=UNKNOWN_TOOL,
        )

    schema = entry.get('input_schema')
    if schema is None:
        raise ProblemException(
            404, 'Tool has no input schema',
            f'"{tool}" takes no job_info. See the tool descriptor for what it does need.',
            type_uri=UNKNOWN_TOOL,
        )

    # $id is applied at serve time rather than stored, because the correct value depends
    # on the host this API is reached at.
    return {'$id': f'{_base_url(request)}/tools/{tool}/input-schema', **schema}
