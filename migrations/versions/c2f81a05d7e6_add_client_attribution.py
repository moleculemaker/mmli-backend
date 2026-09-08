"""add caller attribution columns to job

Records which API a submission arrived through, the Origin header if one was sent, and
a pseudonymous derivation of the client address. Together these answer how much each
tool is used and whether the versioned API is being adopted outside our own frontends.

The backfill is a statement of fact rather than an inference: /v1 and the MCP server did
not exist before this migration, so every row that predates it necessarily arrived
through the legacy API. Populating it means the first usage report can cover prior years
and show the legacy-to-/v1 shift from a real baseline, instead of starting at zero on
deploy day with a large unattributed bucket to explain.

client_fingerprint stays null for historical rows. Unlike the surface, it cannot be
reconstructed -- the addresses were never recorded, which is the correct outcome.

Revision ID: c2f81a05d7e6
Revises: a1c4e7f92b03
"""
import sqlalchemy as sa
import sqlmodel
from alembic import op


revision = 'c2f81a05d7e6'
down_revision = 'a1c4e7f92b03'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column('job', sa.Column('client_surface', sqlmodel.sql.sqltypes.AutoString(), nullable=True))
    op.add_column('job', sa.Column('client_origin', sqlmodel.sql.sqltypes.AutoString(), nullable=True))
    op.add_column('job', sa.Column('client_fingerprint', sqlmodel.sql.sqltypes.AutoString(), nullable=True))

    op.create_index(op.f('ix_job_client_surface'), 'job', ['client_surface'], unique=False)
    op.create_index(op.f('ix_job_client_fingerprint'), 'job', ['client_fingerprint'], unique=False)

    # Every existing row came through the legacy API, because nothing else existed.
    op.execute("UPDATE job SET client_surface = 'legacy' WHERE client_surface IS NULL")


def downgrade() -> None:
    op.drop_index(op.f('ix_job_client_fingerprint'), table_name='job')
    op.drop_index(op.f('ix_job_client_surface'), table_name='job')
    op.drop_column('job', 'client_fingerprint')
    op.drop_column('job', 'client_origin')
    op.drop_column('job', 'client_surface')
