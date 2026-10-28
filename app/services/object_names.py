"""Checks for caller-chosen parts of a MinIO object name, `{job_id}/in/{filename}`.

The routes that check whether a job_id is taken match it exactly, but MinIO sees the
whole key: a job_id of `victim/in/x` writes under another job's in/ with no row found,
and `..` segments can be normalized away on the way (urllib3 does this client-side).
So neither part may carry a path separator or be a dot segment.
"""
import re

# Every id in use fits: server-generated uuid4 hex (32), crispr-copies' crypto.randomUUID()
# hex, and its fallback `${Date.now()}${Math.random()}`, which has a '.'.
_JOB_ID = re.compile(r'[A-Za-z0-9._-]{1,64}')


def is_valid_job_id(job_id: str) -> bool:
    return bool(_JOB_ID.fullmatch(job_id)) and job_id not in ('.', '..')


def is_plain_file_name(name: str) -> bool:
    return name not in ('', '.', '..') and '/' not in name and '\\' not in name
