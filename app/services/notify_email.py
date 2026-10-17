"""Validation for the address a job's completion/failure notification is sent to.

The address is stored as submitted and later used as both the SMTP envelope recipient and
the `To:` header, by every job surface (legacy, /v1 and the MCP server that goes through
it, and the deprecated chemscraper analyze route). All of them call this before the row is
written, so nothing unvalidated reaches the column.

email_validator rejects what matters here: display names, comma lists, quoted local parts,
control characters (including NUL and DEL), and anything over 254 characters.
"""
from typing import Optional

from email_validator import EmailNotValidError, validate_email


class InvalidNotifyEmail(ValueError):
    """The submitted address cannot be used. The message says why, for the client."""


def normalize_notify_email(value: Optional[str]) -> Optional[str]:
    """Return the address to store, or None when the submitter gave none.

    Blank counts as none: frontends send an empty string when the field is left empty,
    and every reader of Job.email already treats a falsy value as "do not notify".
    """
    if value is None or not value.strip():
        return None
    try:
        validated = validate_email(
            value.strip(),
            # A DNS lookup on every submission, for a check a typo'd domain would fail
            # anyway at send time.
            check_deliverability=False,
            # The SMTP relay is not known to support SMTPUTF8, so a non-ASCII local part
            # could never be delivered. Internationalized domains are still accepted and
            # stored in their ASCII (punycode) form.
            allow_smtputf8=False,
        )
    except EmailNotValidError as ex:
        raise InvalidNotifyEmail(str(ex)) from ex
    # Domain lowercased; the local part keeps its case, which is the mailbox's to interpret.
    return validated.ascii_email
