"""add saved_molecule table

Revision ID: 58fe2962368d
Revises: 30b240622d34
Create Date: 2024-06-21 11:22:55.825619

"""
from alembic import op
import sqlalchemy as sa
import sqlmodel


# revision identifiers, used by Alembic.
revision = '58fe2962368d'
down_revision = '30b240622d34'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table('saved_molecule',
    sa.Column('id', sa.Integer, primary_key=True, autoincrement=True),
    sa.Column('email', sqlmodel.sql.sqltypes.AutoString(), nullable=True),
    sa.Column('job_id', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('molecule_id', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('time_created', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.ForeignKeyConstraint(['job_id'], ['job.job_id']),
    sa.UniqueConstraint('email', 'job_id', 'molecule_id', name='unique_email_job_molecule')
    )
    op.create_index('idx_email_job', 'saved_molecule', ['email', 'job_id'])
    op.create_index('idx_email_job_molecule', 'saved_molecule', ['email', 'job_id', 'molecule_id'], unique=True)
    pass


def downgrade() -> None:
    op.drop_index('idx_email_job', table_name='saved_molecule')
    op.drop_index('idx_email_job_molecule', table_name='saved_molecule')
    op.drop_table('saved_molecule')
    pass
