"""Service description.

`/v1/service-info` follows the GA4GH service-info shape, which is a small, widely
implemented spec for "what is this service and who runs it". Using it means a client
that already speaks to other life-science services can identify this one without
bespoke code.

The access model is stated explicitly. FAIR does not require that data be open; it
requires that the conditions for getting at it be stated. Saying "no authentication,
jobs are reachable by whoever holds the id" is a real answer, and a more useful one
than silence.
"""
from typing import Any, Dict

from fastapi import APIRouter, Request

from config import app_config, get_logger
from services import tool_registry

router = APIRouter()
log = get_logger(__name__)


@router.get('/service-info', tags=['Discovery'], summary='Describe this service')
async def service_info(request: Request) -> Dict[str, Any]:
    base = str(request.base_url).rstrip('/') + '/v1'
    return {
        'id': 'org.moleculemaker.alphasynthesis',
        'name': 'AlphaSynthesis Tools API',
        'type': {
            'group': 'org.moleculemaker',
            'artifact': 'tools',
            'version': '1.0.0',
        },
        'description': (
            'Runs computational chemistry and enzymology tools as asynchronous jobs. '
            'Tools, their input schemas and their results are discoverable from this API.'
        ),
        'organization': {
            'name': 'Molecule Maker Lab Institute',
            'url': 'https://moleculemaker.org',
        },
        'contactUrl': 'https://github.com/moleculemaker/mmli-backend/issues',
        'documentationUrl': f'{base}/docs',
        'environment': app_config.get('server', {}).get('hostName', 'unknown'),
        'version': '1.0.0',

        # Not part of the GA4GH spec; stated because a client needs to know it and
        # because leaving it unsaid would be the least honest option.
        'auth': {
            'type': 'none',
            'model': 'capability-url',
            'description': (
                'No authentication. A job is readable and cancelable by anyone who has '
                'its job_id, which is a server-assigned UUIDv4 and is never listed. '
                'Treat a job_id as a secret. There is no endpoint that enumerates jobs.'
            ),
        },
        'rateLimit': {
            'description': (
                'Job submission is rate limited per client IP at the ingress. Discovery '
                'and job polling are not limited.'
            ),
        },
        'toolCount': len(tool_registry.list_tools()),
        'links': {
            'self': f'{base}/service-info',
            'tools': f'{base}/tools',
            'openapi': f'{base}/openapi.json',
            'docs': f'{base}/docs',
        },
    }
