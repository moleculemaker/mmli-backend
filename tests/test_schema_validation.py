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


class TestTheRuleHoldsUnderARef:
    """jsonschema re-picks the validator class for every subschema (descend -> evolve ->
    validator_for), and a subschema carrying `$schema` -- the root, reached through
    `$ref: "#"` -- picked the stock Draft202012Validator, dropping the fix."""

    SCHEMA = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"id": {"type": "string", "pattern": "^[a-z]+$"}, "child": {"$ref": "#"}},
    }

    def test_a_recursive_ref_still_rejects_a_trailing_newline(self):
        assert SchemaValidator(self.SCHEMA).is_valid({"child": {"id": "a"}}) is True
        assert SchemaValidator(self.SCHEMA).is_valid({"child": {"id": "a\n"}}) is False

    def test_the_stock_validator_is_not_changed(self):
        """The MCP server validates with jsonschema in the same process."""
        from jsonschema import Draft202012Validator, validators
        assert validators.validator_for(self.SCHEMA) is Draft202012Validator
        assert Draft202012Validator(self.SCHEMA).is_valid({"id": "a\n"}) is True

    def test_the_evolve_copy_matches_the_installed_jsonschema(self):
        """services/json_schema.py's _evolve copies jsonschema 4.23.0's evolve. Upgrading
        jsonschema means comparing it with the new evolve, then updating this version."""
        from importlib.metadata import version
        assert version("jsonschema") == "4.23.0"


class TestTheRewriteFollowsPythonsParse:
    def test_a_bracket_first_in_a_negated_class_is_a_literal(self):
        """Python reads `[^]$]` as a class of `]` and `$`. Treating the `]` as closing
        it put `\\Z` inside the class, which `re` rejects (bad escape)."""
        assert ecma_pattern("^[^]$]+$") == r"^[^]$]+\Z"
        assert _valid("^[^]$]+$", "ab") is True
        assert _valid("^[^]$]+$", "a]") is False

    def test_a_bracket_first_in_a_class_is_a_literal(self):
        assert ecma_pattern("^[]$]+$") == r"^[]$]+\Z"
        assert _valid("^[]$]+$", "]$") is True
        assert _valid("^[]$]+$", "]$\n") is False

    def test_flags_that_change_neither_dollar_nor_classes_are_still_rewritten(self):
        assert ecma_pattern("(?i)^a$") == r"(?i)^a\Z"

    @pytest.mark.parametrize("pattern", ["(?m:^a$)", "(?m)^a$", "(?im)^a$", "(?x) a # [comment\n $"])
    def test_a_pattern_with_multiline_or_verbose_flags_is_left_alone(self, pattern):
        """`(?m)` makes `$` a line anchor and `(?x)` makes `[` inside a comment a
        literal; the scanner would read both wrong. ECMA has no `(?m)`/`(?x)` syntax,
        so no frontend-valid schema uses them; shipped schemas are checked in
        test_tool_registry."""
        assert ecma_pattern(pattern) == pattern


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
