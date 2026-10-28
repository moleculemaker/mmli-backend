"""Tests for how /v1 validates against JSON Schema (#147, #149).

JSON Schema's `pattern` is an ECMA-262 regular expression, where `$` matches only at the
very end of the string. jsonschema evaluates it with Python's `re.search`, where `$` also
matches just before a final newline -- so the API accepted "a\\n" for `^[a-z]+$` while the
browser's RegExp rejected it.
"""
import pytest

from services import generic_job
from services.json_schema import SchemaValidator, ecma_pattern


def _valid(pattern, value):
    return SchemaValidator({"type": "string", "pattern": pattern}).is_valid(value)


class TestDollarMatchesOnlyAtTheEnd:
    def test_a_trailing_newline_fails_an_anchored_pattern(self):
        assert _valid("^[a-z]+$", "a") is True
        assert _valid("^[a-z]+$", "a\n") is False

    def test_dollar_inside_a_character_class_is_a_literal(self):
        assert _valid("^[a$]+$", "a$a") is True
        assert _valid("^[a$]+$", "a$a\n") is False

    def test_an_escaped_dollar_is_a_literal(self):
        assert _valid(r"^cost \$5$", "cost $5") is True
        assert _valid(r"^cost \$5$", "cost $5\n") is False

    def test_dollar_inside_a_group(self):
        assert _valid("^(a|b$)", "b") is True
        assert _valid("^(a|b$)", "b\n") is False
        assert _valid("^(a|b$)", "a\n") is True  # the `a` branch is unanchored

    def test_a_pattern_without_dollar_is_unchanged(self):
        assert ecma_pattern("^[a-z]+") == "^[a-z]+"
        assert _valid("^[a-z]+", "a\n") is True

    def test_an_escaped_backslash_before_dollar_leaves_the_dollar_an_anchor(self):
        """`\\\\$` is a literal backslash, then the end anchor."""
        assert _valid(r"^a\\$", "a\\") is True
        assert _valid(r"^a\\$", "a\\\n") is False

    def test_an_escaped_bracket_does_not_open_a_class(self):
        assert _valid(r"^\[$", "[") is True
        assert _valid(r"^\[$", "[\n") is False


class TestSubmissionRejectsATrailingNewline:
    """The reproduction from #147, through the API, against a schema that ships today."""

    def test_clean_rejects_a_sequence_ending_in_a_newline(self, client):
        payload = {"input_fasta": [{"header": "seq1", "sequence": "MEDIPDTSRPPLKYVK\n", "DNA_sequence": ""}]}

        resp = client.post("/v1/tools/clean/jobs", json=payload)

        assert resp.status_code == 422
        assert any(e["pointer"] == "/input_fasta/0/sequence" for e in resp.json()["errors"])

    def test_clean_still_accepts_the_sequence_without_it(self, client):
        payload = {"input_fasta": [{"header": "seq1", "sequence": "MEDIPDTSRPPLKYVK", "DNA_sequence": ""}]}

        assert client.post("/v1/tools/clean/jobs", json=payload).status_code == 201


class TestRootPointer:
    """#149: RFC 6901 names the whole document "", not "/" (the empty-string key)."""

    def test_a_missing_top_level_property_points_at_the_root(self, client):
        resp = client.post("/v1/tools/novostoic-optstoic/jobs", json={"primary_precursor": "MNXM1"})

        assert resp.status_code == 422
        assert [e["pointer"] for e in resp.json()["errors"]] == [""]

    def test_nested_pointers_are_unchanged(self, client):
        payload = {"input_fasta": [{"header": "seq1", "sequence": "X1", "DNA_sequence": ""}]}

        resp = client.post("/v1/tools/clean/jobs", json=payload)

        assert [e["pointer"] for e in resp.json()["errors"]] == ["/input_fasta/0/sequence"]


class TestResultsUseTheSameValidator:
    def test_results_validation_applies_ecma_dollar(self, monkeypatch):
        schema = {"type": "object", "properties": {"id": {"type": "string", "pattern": "^[a-z]+$"}}}
        monkeypatch.setattr(generic_job.tool_registry, "get_output_schema", lambda job_type: schema)

        class _Store:
            def get_file(self, bucket, name):
                return b'{"id": "a\\n"}'

        with pytest.raises(generic_job.InvalidResults) as raised:
            generic_job.read_results(_Store(), "t", "j1")
        assert raised.value.errors[0]["pointer"] == "/id"
