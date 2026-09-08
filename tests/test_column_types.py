"""Pins the SQL types of the enum-valued columns on `job`.

These assertions exist because the rest of the suite structurally cannot make them. Every
other test builds its schema from SQLModel.metadata.create_all(), so it tests against
whatever the models currently declare. A deployed database is built by Alembic instead,
and the two are not the same artifact -- so a change in how SQLModel maps a Python type is
invisible to a suite that only ever asks the models what they think.

That is not hypothetical. Upgrading SQLModel 0.0.8 -> 0.0.22 changed str-Enum fields from
AutoString to a native sa.Enum. Every migration in this repo created `job.phase` and
`job.type` as VARCHAR, and no migration creates a Postgres enum type, so the upgrade left
the models describing columns that exist in no database. Every job submission 500'd with
`type "jobtype" does not exist` -- while the full suite passed.

A type assertion is the cheap half of the guard: it needs no database, so it runs
everywhere, and it fails the moment a dependency bump moves the mapping again.
"""
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

from models.enums import JobStatus, JobType
from models.sqlmodel.models import Job

PG = postgresql.dialect()


class TestEnumColumnsAreStoredAsStrings:
    def test_type_is_not_a_native_enum(self):
        assert not isinstance(Job.__table__.c.type.type, sa.Enum), (
            "job.type is declared as a native SQL enum, but every migration creates it as "
            "VARCHAR and no migration creates the type. Queries filtering on it fail with "
            "UndefinedObjectError. Keep it on EnumValueString."
        )

    def test_phase_is_not_a_native_enum(self):
        assert not isinstance(Job.__table__.c.phase.type, sa.Enum), (
            "job.phase is declared as a native SQL enum; see test_type_is_not_a_native_enum."
        )

    def test_type_compiles_to_varchar(self):
        assert "VARCHAR" in str(Job.__table__.c.type.type.compile(dialect=PG))

    def test_phase_compiles_to_varchar(self):
        assert "VARCHAR" in str(Job.__table__.c.phase.type.compile(dialect=PG))


class TestTheStoredFormIsTheValueNotTheName:
    """sa.Enum persists 'ML_SIMPLEFOLD'; everything outside Python uses 'ml-simplefold' --
    the frontends, coordinator.py, and the published /v1 input schema. Storing names would
    not merely be a different encoding: it would orphan every row already written and every
    client already deployed."""

    def test_a_job_type_is_written_as_its_value(self):
        write = Job.__table__.c.type.type.bind_processor(PG)
        assert write(JobType.ML_SIMPLEFOLD) == "ml-simplefold"

    def test_a_phase_is_written_as_its_value(self):
        write = Job.__table__.c.phase.type.bind_processor(PG)
        assert write(JobStatus.QUEUED) == "queued"

    def test_a_plain_string_is_accepted_as_written(self):
        """Routers compare and assign raw path strings, so a str must bind unchanged."""
        write = Job.__table__.c.type.type.bind_processor(PG)
        assert write("cleandb-mepesm") == "cleandb-mepesm"

    def test_null_stays_null(self):
        assert Job.__table__.c.type.type.bind_processor(PG)(None) is None


class TestReadsComeBackAsTheEnum:
    """The attribute is declared JobType, so it should hold one. Plain AutoString hands
    back a bare str, which makes Pydantic 2 warn on every serialization and diverges from
    what the metadata-built schema the suite uses has always produced."""

    def test_a_stored_value_becomes_its_enum_member(self):
        read = Job.__table__.c.type.type.result_processor(PG, None)
        assert read("ml-simplefold") is JobType.ML_SIMPLEFOLD

    def test_a_stored_phase_becomes_its_enum_member(self):
        read = Job.__table__.c.phase.type.result_processor(PG, None)
        assert read("queued") is JobStatus.QUEUED

    def test_an_unrecognized_stored_value_is_returned_rather_than_raising(self):
        """A tool renamed or retired since the row was written. Those rows have to stay
        readable -- failing here would break every query that touched one, including the
        usage report, which aggregates over all of history."""
        read = Job.__table__.c.type.type.result_processor(PG, None)
        assert read("a-tool-we-since-removed") == "a-tool-we-since-removed"

    def test_null_stays_null(self):
        assert Job.__table__.c.type.type.result_processor(PG, None)(None) is None


class TestRoundTrip:
    def test_every_job_type_survives_a_write_and_a_read(self):
        write = Job.__table__.c.type.type.bind_processor(PG)
        read = Job.__table__.c.type.type.result_processor(PG, None)
        for member in JobType:
            assert read(write(member)) is member

    def test_every_phase_survives_a_write_and_a_read(self):
        write = Job.__table__.c.phase.type.bind_processor(PG)
        read = Job.__table__.c.phase.type.result_processor(PG, None)
        for member in JobStatus:
            assert read(write(member)) is member
