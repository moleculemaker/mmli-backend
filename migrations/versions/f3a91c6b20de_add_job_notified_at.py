"""add job.notified_at, the durable notification marker

Notification idempotency was a `{job_id}/email-sent` object in MinIO, read through
MinIOService.check_file_exists. That gate fails OPEN: it returns False for every S3Error,
and a stat_object against a missing bucket is a bodyless 404 that minio 7.1.17
synthesizes as NoSuchKey - the same code as a genuinely absent marker - while
mark_email_as_sent cannot lay one down either, because put_object does not create
buckets. A reconcile pass that revisits every terminal Job every 600s for 12h turns that
into ~72 copies of the same mail.

job.notified_at replaces it. The guarantee is at-least-once rather than exactly-once -
the send precedes the record of it, in its own transaction - but the duplicate window
shrinks from "as long as MinIO is confused" to "a database failure in the moment after a
successful send".

BACKFILL. Null means "not yet notified", so leaving existing rows null would mail every
user whose job finished inside the current ttlSecondsAfterFinished window the moment this
deploys. Rows already in a terminal phase are therefore stamped, using the job's own end
time (or its creation time where time_end was never written, which is every job on the
Kubernetes path). Those values are NOT observed notification times and should not be read
as such - they are a suppression backfill, chosen over null because the alternative is a
mail storm and over a sentinel because the column is typed as a timestamp. Rows that are
not terminal keep null and notify normally when they finish.

The MinIO markers are left in place, unread, because removing the objects is not this
migration's business. They are NOT a general rollback safety net, though: they cover only
rows notified before this deploy. A rollback to the immediately preceding revision is
safe regardless, because that code gates on the phase transition and never reads them. A
rollback further back, to the marker-gated code, re-notifies every job this code emailed
that is still inside the 12h TTL, because those jobs were never given a marker.

Revision ID: f3a91c6b20de
Revises: c2f81a05d7e6
"""
import sqlalchemy as sa
from alembic import op


revision = 'f3a91c6b20de'
down_revision = 'c2f81a05d7e6'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column('job', sa.Column('notified_at', sa.Integer(), nullable=True))
    op.execute(
        """
        UPDATE job
           SET notified_at = COALESCE(NULLIF(time_end, 0), time_created)
         WHERE phase IN ('completed', 'error')
           AND notified_at IS NULL
        """
    )


def downgrade() -> None:
    op.drop_column('job', 'notified_at')
