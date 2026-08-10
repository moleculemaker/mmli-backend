"""add job provenance columns and the idempotency key table

Adds:
  job.image_digest   the immutable image digest a job actually ran, captured from the
                     pod. The configured image reference is often untagged or points at
                     a mutable tag, so it does not identify what produced a result.
  job.parent_job_id  the parent of a coordinator-created subjob. This relationship
                     previously existed only as keys inside the job_info JSON blob, so
                     it could not be queried.
  idempotency_key    maps a client-supplied Idempotency-Key to the job it created, so a
                     retried submission returns the original job rather than starting a
                     second container.

Both columns are nullable with no backfill. Existing rows genuinely do not have this
information: their pods are long gone, and the parent link for historical subjobs would
have to be reconstructed by parsing job_info, which would produce guesses rather than
facts. Null is the honest value.

Revision ID: a1c4e7f92b03
Revises: bdfd696c3c79
"""
import sqlalchemy as sa
import sqlmodel
from alembic import op


revision = 'a1c4e7f92b03'
down_revision = 'bdfd696c3c79'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column('job', sa.Column('image_digest', sqlmodel.sql.sqltypes.AutoString(), nullable=True))
    op.add_column('job', sa.Column('parent_job_id', sqlmodel.sql.sqltypes.AutoString(), nullable=True))
    op.create_index(op.f('ix_job_parent_job_id'), 'job', ['parent_job_id'], unique=False)

    op.create_table(
        'idempotency_key',
        sa.Column('key', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.Column('job_type', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.Column('job_id', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.Column('time_created', sa.Integer(), nullable=False),
        # Composite key: two tools may legitimately be handed the same key by different
        # clients, and one must not shadow the other.
        sa.PrimaryKeyConstraint('key', 'job_type'),
    )
    op.create_index(op.f('ix_idempotency_key_job_id'), 'idempotency_key', ['job_id'], unique=False)


def downgrade() -> None:
    op.drop_index(op.f('ix_idempotency_key_job_id'), table_name='idempotency_key')
    op.drop_table('idempotency_key')
    op.drop_index(op.f('ix_job_parent_job_id'), table_name='job')
    op.drop_column('job', 'parent_job_id')
    op.drop_column('job', 'image_digest')
