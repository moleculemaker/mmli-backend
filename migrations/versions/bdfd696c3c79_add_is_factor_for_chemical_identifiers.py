"""add is_factor for chemical_identifiers

Revision ID: bdfd696c3c79
Revises: 58fe2962368d
Create Date: 2024-07-28 15:16:59.193857

"""
from alembic import op
import sqlalchemy as sa
import sqlmodel


# revision identifiers, used by Alembic.
revision = 'bdfd696c3c79'
down_revision = '58fe2962368d'
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.add_column('chemical_identifier', sa.Column('is_cofactor', sa.Boolean(), nullable=True))
    # ### end Alembic commands ###


def downgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_column('chemical_identifier', 'is_cofactor')
    # ### end Alembic commands ###
