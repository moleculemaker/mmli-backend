"""Caller attribution recorded at submission time.

Deliberately small. Nothing here logs requests, and nothing runs on a read: the two
questions this exists to answer -- how much each tool is used, and whether the versioned
API is being adopted outside our own frontends -- are both about submissions, and every
submission is already a durable row. Recording a little more on that row is enough, and
avoids a request log that a single polling client could dominate.

What is stored, and what is not:

  stored      which API the call arrived through, the Origin header if present, the
              user agent, and a salted derivation of the client address
  not stored  the client address itself, and nothing at all about reads
"""
import hashlib
import hmac
import secrets
import time
from typing import Optional

from config import ANALYTICS_SALT, get_logger

log = get_logger(__name__)

SURFACE_LEGACY = 'legacy'
SURFACE_V1 = 'v1'
SURFACE_MCP = 'mcp'

# Header the MCP adapter uses to declare which surface a submission came from. It is
# accompanied by a token minted at startup, so a caller who guesses the header name
# cannot mislabel their traffic as agent traffic.
SURFACE_HEADER = 'X-MMLI-Surface'
SURFACE_TOKEN_HEADER = 'X-MMLI-Surface-Token'

# Minted per process and never leaves it. The MCP adapter reaches /v1 over an in-process
# ASGI transport and needs some way to say "this submission came from an agent"; without
# a token, any external caller could send the same header and have their traffic counted
# as agent traffic, quietly corrupting the adoption numbers this exists to produce.
SURFACE_TOKEN = secrets.token_urlsafe(32)


def _salt_period(when: Optional[float] = None) -> str:
    """Rotation period for the fingerprint salt: the calendar year.

    Rotating at all limits how long any one pseudonym remains linkable. Rotating
    annually rather than monthly is a deliberate trade: it keeps a unique-researcher
    count computable across a full reporting year, at the cost of a pseudonym that
    persists for that year.
    """
    return time.strftime('%Y', time.gmtime(when if when is not None else time.time()))


def client_fingerprint(ip: Optional[str], when: Optional[float] = None) -> Optional[str]:
    """Derive a pseudonymous, period-scoped identifier for a client address.

    Returns None when no salt is configured, which is the default. That is deliberate:
    a fingerprint derived from an empty or guessable salt is a reversible encoding of
    the address rather than a pseudonym, so it is better to record nothing.
    """
    if not ANALYTICS_SALT or not ip:
        return None

    message = f'{_salt_period(when)}:{ip}'.encode('utf-8')
    digest = hmac.new(ANALYTICS_SALT.encode('utf-8'), message, hashlib.sha256).hexdigest()
    # Truncated: enough to distinguish clients within a period, short enough that the
    # stored value is not a full-strength handle on anything.
    return digest[:16]


def client_ip(request) -> Optional[str]:
    """Best-effort client address, honoring the proxy header the ingress sets.

    Only ever used to derive a fingerprint; the address is not stored. X-Forwarded-For
    is caller-controlled and therefore spoofable, which is acceptable here -- a client
    that forges it makes its own traffic harder to count, and nothing security-relevant
    depends on the value.
    """
    forwarded = request.headers.get('x-forwarded-for')
    if forwarded:
        return forwarded.split(',')[0].strip()
    return request.client.host if request.client else None


def attribution(request, surface: str, surface_token: Optional[str] = None) -> dict:
    """Collect the attribution fields for a submission.

    `surface` is what the handler believes it is; a request may override it to 'mcp' by
    presenting the process-local token, which only the in-process adapter holds.
    """
    declared = request.headers.get(SURFACE_HEADER)
    presented = request.headers.get(SURFACE_TOKEN_HEADER)
    if declared and surface_token and hmac.compare_digest(presented or '', surface_token):
        surface = declared

    return {
        'client_surface': surface,
        'client_origin': request.headers.get('origin'),
        'user_agent': request.headers.get('user-agent', ''),
        'client_fingerprint': client_fingerprint(client_ip(request)),
    }
