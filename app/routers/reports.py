"""Usage reporting.

Aggregate-only, by design. These endpoints answer "how much is each tool used" and "is
the versioned API being adopted", and nothing else: no endpoint here returns a job row,
an email address, or a fingerprint. That is a deliberate constraint rather than an
oversight -- an endpoint that returns rows becomes, sooner or later, the way somebody
exports the user table.

Access is restricted to members of a configured OIDC group. The reports describe who is
using the service in aggregate, which is not something to publish alongside an
unauthenticated API.
"""
import time
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from sqlalchemy import func, select
from sqlmodel.ext.asyncio.session import AsyncSession

from config import REPORTING_GROUP, get_logger
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job
from services.userinfo_service import validate_auth_cookie

router = APIRouter()
log = get_logger(__name__)


async def require_reporting_group(request: Request) -> Dict[str, Any]:
    """Authorize by OIDC group membership.

    `validate_auth_cookie` accepts a `required_scopes` argument and has always ignored
    it, so authentication here has to be paired with an explicit membership check;
    resolving a user proves only that somebody is logged in.
    """
    if not REPORTING_GROUP:
        # Refuse rather than fall open. An unset group would otherwise mean "any
        # authenticated user", which is not what an operator who left it unset intended.
        raise HTTPException(
            status_code=503,
            detail='Reporting is not configured: no reporting group has been set.',
        )

    user = validate_auth_cookie(request)
    if not user:
        raise HTTPException(status_code=401, detail='Authentication required.')

    if REPORTING_GROUP not in (user.get('groups') or []):
        raise HTTPException(status_code=403, detail='Not a member of the reporting group.')

    return user


def _bounds(from_: Optional[str], to: Optional[str]) -> tuple:
    """Parse ISO dates into an epoch range, defaulting to the last 365 days."""
    def _parse(value: str) -> int:
        try:
            return int(time.mktime(time.strptime(value, '%Y-%m-%d')))
        except ValueError:
            raise HTTPException(
                status_code=400,
                detail=f'Invalid date "{value}". Expected YYYY-MM-DD.',
            )

    now = int(time.time())
    start = _parse(from_) if from_ else now - 365 * 86400
    end = _parse(to) if to else now
    if start > end:
        raise HTTPException(status_code=400, detail='"from" is after "to".')
    return start, end


@router.get('/internal/reports/usage', tags=['Reports'],
            summary='Aggregate usage counts by tool and by API surface')
async def usage_report(
    request: Request,
    from_: Optional[str] = Query(default=None, alias='from',
                                 description='Inclusive start date, YYYY-MM-DD.'),
    to: Optional[str] = Query(default=None, description='Inclusive end date, YYYY-MM-DD.'),
    _user: Dict[str, Any] = Depends(require_reporting_group),
    db: AsyncSession = Depends(get_session),
) -> Dict[str, Any]:
    """Return usage aggregates for a period.

    Counts reported per tool:

      jobs               submissions in the period
      identified_users   distinct email addresses, where a submitter supplied one
      distinct_clients   distinct pseudonymous fingerprints, an approximation of how
                         many separate anonymous callers there were. Null-safe: absent
                         when no salt is configured, which is the default.

    identified_users and distinct_clients measure different populations and must not be
    added together. A submitter who supplies an email is counted in the first; one who
    does not may be counted in the second, if fingerprinting is enabled.
    """
    start, end = _bounds(from_, to)
    period = (Job.time_created >= start, Job.time_created <= end)

    by_tool = (await db.execute(
        select(
            Job.type,
            func.count().label('jobs'),
            func.count(func.distinct(Job.email)).label('identified_users'),
            func.count(func.distinct(Job.client_fingerprint)).label('distinct_clients'),
        ).where(*period).group_by(Job.type)
    )).all()

    by_surface = (await db.execute(
        select(Job.client_surface, func.count())
        .where(*period).group_by(Job.client_surface)
    )).all()

    by_phase = (await db.execute(
        select(Job.phase, func.count()).where(*period).group_by(Job.phase)
    )).all()

    total = sum(row.jobs for row in by_tool)

    return {
        'period': {
            'from': time.strftime('%Y-%m-%d', time.gmtime(start)),
            'to': time.strftime('%Y-%m-%d', time.gmtime(end)),
        },
        'total_jobs': total,
        'by_tool': [
            {
                'tool': row.type,
                'jobs': row.jobs,
                'identified_users': row.identified_users,
                'distinct_clients': row.distinct_clients,
            }
            for row in sorted(by_tool, key=lambda r: r.jobs, reverse=True)
        ],
        # 'legacy' covers everything submitted before the versioned API existed, which
        # is why historical periods report entirely as legacy.
        'by_surface': {(surface or 'unknown'): count for surface, count in by_surface},
        'by_status': {str(phase): count for phase, count in by_phase},
    }
