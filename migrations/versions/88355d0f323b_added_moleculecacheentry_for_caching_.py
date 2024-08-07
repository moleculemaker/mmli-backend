"""added moleculecacheentry for caching molecules, modified job schema, added flaggedmolecule for saving flagged molecules

Revision ID: 88355d0f323b
Revises: d775ee615d7b
Create Date: 2023-10-25 19:19:18.903467

"""
from alembic import op
import sqlalchemy as sa
import sqlmodel


# revision identifiers, used by Alembic.
revision = '88355d0f323b'
down_revision = 'd775ee615d7b'
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.alter_column('job', 'job_id',
               existing_type=sa.VARCHAR(),
               nullable=False)
    op.alter_column('job', 'image',
               existing_type=sa.VARCHAR(),
               nullable=True)
    op.drop_column('job', 'id')
    op.create_primary_key('job_id_pk', 'job', ['job_id'])
    op.create_table('moleculecacheentry',
    sa.Column('smile', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('pub_chem_id', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('name', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('molecular_formula', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('molecular_weight', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('chemical_safety', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('description', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.PrimaryKeyConstraint('smile')
    )
    op.create_table('flaggedmolecule',
    sa.Column('smile', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('job_id', sqlmodel.sql.sqltypes.AutoString(), nullable=False),
    sa.Column('doc_id', sqlmodel.sql.sqltypes.AutoString(), nullable=True),
    sa.Column('time_created', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['job_id'], ['job.job_id'], ),
    sa.ForeignKeyConstraint(['smile'], ['moleculecacheentry.smile'], ),
    sa.PrimaryKeyConstraint('smile', 'job_id')
    )
    
    # ### end Alembic commands ###


def downgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.add_column('job', sa.Column('id', sa.INTEGER(), autoincrement=True, nullable=False))
    op.alter_column('job', 'image',
               existing_type=sa.VARCHAR(),
               nullable=False)
    op.alter_column('job', 'job_id',
               existing_type=sa.VARCHAR(),
               nullable=True)
    op.drop_table('flaggedmolecule')
    op.drop_table('moleculecacheentry')
    # ### end Alembic commands ###
