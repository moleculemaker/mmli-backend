"""Tests for the tool registry.

Most of these guard the registry against drifting away from the code: a new job type
with no metadata, a schema file that stops parsing, a `$ref` that reaches outside its
own file. Cheap to run, and they fail at the moment the mistake is made rather than
when a client hits a missing entry.
"""
import json

import pytest
from jsonschema import Draft202012Validator

from models.enums import JobTypes
from services import tool_registry


REQUIRED_FIELDS = ["name", "summary", "description", "execution", "schema_status"]
SME_FIELDS = ["license", "citation"]


class TestRegistryCoverage:
    def test_every_job_type_is_registered(self):
        """A new job type without metadata would be invisible to discovery."""
        missing = tool_registry.unregistered_job_types()

        assert missing == set(), f"job types missing from tools.yaml: {sorted(missing)}"

    def test_registry_describes_only_real_job_types(self):
        for identifier in tool_registry.get_registry():
            assert identifier in JobTypes

    def test_defaults_is_not_registered(self):
        """`defaults` computes pi in Perl. It is a harness fixture, not a tool."""
        assert tool_registry.get_tool("defaults") is None

    def test_internal_subjobs_are_hidden_from_listings(self):
        listed = tool_registry.list_tools()

        assert "ezspec-unidock" not in listed
        assert "ezspec-inference" not in listed
        assert "ez-specificity" in listed

    def test_internal_subjobs_are_still_retrievable_by_id(self):
        assert tool_registry.get_tool("ezspec-unidock") is not None
        assert tool_registry.list_tools(include_internal=True).get("ezspec-unidock")


class TestRegistryContent:
    @pytest.mark.parametrize("identifier", sorted(tool_registry.get_registry()))
    def test_required_fields_are_present_and_non_empty(self, identifier):
        tool = tool_registry.get_tool(identifier)

        for field in REQUIRED_FIELDS:
            assert tool.get(field), f"{identifier}: {field} is missing or empty"

    @pytest.mark.parametrize("identifier", sorted(tool_registry.get_registry()))
    def test_execution_is_a_known_value(self, identifier):
        assert tool_registry.get_tool(identifier)["execution"] in {
            "kubernetes-job", "parent-job", "subjob",
        }

    @pytest.mark.parametrize("identifier", sorted(tool_registry.get_registry()))
    def test_schema_status_is_a_known_value(self, identifier):
        assert tool_registry.get_tool(identifier)["schema_status"] in {
            "verified", "from-container", "unverified",
        }

    @pytest.mark.parametrize("identifier", sorted(tool_registry.get_registry()))
    def test_sme_fields_are_present_as_keys_even_when_unknown(self, identifier):
        """Published as null rather than omitted.

        An absent key is ambiguous -- it could mean "no license" or "nobody has said".
        An explicit null says the second thing, which is the truth today.
        """
        tool = tool_registry.get_tool(identifier)

        for field in SME_FIELDS:
            assert field in tool, f"{identifier}: {field} key is missing entirely"
        assert "edam" in tool
        assert set(tool["edam"]) == {"operation", "topic"}


class TestSchemas:
    @pytest.mark.parametrize("identifier", sorted(tool_registry.get_registry()))
    def test_schema_is_valid_json_schema_2020_12(self, identifier):
        schema = tool_registry.get_input_schema(identifier)
        if schema is None:
            return

        assert schema.get("$schema") == "https://json-schema.org/draft/2020-12/schema"
        Draft202012Validator.check_schema(schema)

    @pytest.mark.parametrize("identifier", sorted(tool_registry.get_registry()))
    def test_schema_is_self_contained(self, identifier):
        """Every $ref must resolve inside its own document.

        MCP tool definitions carry a single inline inputSchema, and registry crawlers
        fetch one file. A cross-file $ref would force a bundling step on every consumer.
        """
        schema = tool_registry.get_input_schema(identifier)
        if schema is None:
            return

        def refs(node):
            if isinstance(node, dict):
                for key, value in node.items():
                    if key == "$ref":
                        yield value
                    else:
                        yield from refs(value)
            elif isinstance(node, list):
                for item in node:
                    yield from refs(item)

        for ref in refs(schema):
            assert ref.startswith("#/"), f"{identifier}: non-local $ref {ref}"

    @pytest.mark.parametrize("identifier", sorted(tool_registry.get_registry()))
    def test_schema_documents_its_properties(self, identifier):
        """An agent picks fields by reading descriptions, so undocumented ones are dead weight."""
        schema = tool_registry.get_input_schema(identifier)
        if schema is None:
            return

        def check(node, path="#"):
            if not isinstance(node, dict):
                return
            for name, subschema in (node.get("properties") or {}).items():
                assert isinstance(subschema, dict)
                assert subschema.get("description") or subschema.get("$ref"), \
                    f"{identifier}: {path}/{name} has no description"
                check(subschema, f"{path}/{name}")
            for key in ("items", "additionalProperties"):
                if isinstance(node.get(key), dict):
                    check(node[key], f"{path}/{key}")
            for key in ("anyOf", "oneOf", "allOf"):
                for i, sub in enumerate(node.get(key) or []):
                    check(sub, f"{path}/{key}/{i}")
            for name, sub in (node.get("$defs") or {}).items():
                check(sub, f"#/$defs/{name}")

        check(schema)


class TestSchemasAcceptKnownGoodInput:
    """Sanity-check the schemas against payloads taken from this repo or the containers.

    A schema that rejects a real request is worse than no schema, so the examples here
    are copied from container fixtures and from docstrings in job.py rather than
    invented.
    """

    CASES = {
        # from the docstring in routers/job.py
        "aceretro": {"smiles": "O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO"},
        # from novostoic-container in/input.optstoic.json
        "novostoic-optstoic": {"primary_precursor": "MNXM1137670", "target_molecule": "MNXM26"},
        # from novostoic-container in/input.enzrank.json
        "novostoic-enzrank": {"enzyme_sequence": "A0A4P8WFA8:MTKRVLVTGG", "primary_precursor": "C00149"},
        # from novostoic-container in/input.dgpredictor.json
        "novostoic-dgpredictor": {
            "ph": 7, "ionic_strength": 0.3,
            "reactions": [{
                "type": "keggId",
                "reaction_keggid": "C01745 + C00004 <=> N00001 + C00003 + C00001",
                "molecule_number": "N00001",
                "molecule_inchi_or_smiles": "InChI=1S/C14H12O",
            }],
        },
        # from cheminfo-service input/job.json
        "oed-cheminfo": {
            "query_smiles": "COC1=C(C=CC(=C1)CCN)O",
            "smiles_list": ["COc1ccc2c(c1OC)C(O)O[C@@H]2"],
            "algorithms": ["tanimoto", "fragment", "mcs"],
            "config": {"tanimoto": {"fptype": "rdkit"}},
        },
        # from the docstring in routers/job.py; the three kinetic tools share a schema
        "oed-dlkcat": {"input_pairs": [
            {"name": "example", "sequence": "MEDIPDTSRPPLKYVK", "type": "FASTA",
             "smiles": "OC1=CC=C(C[C@@H](C(O)=O)N)C=C1"},
        ]},
        "oed-unikp": {"input_pairs": [
            {"name": "example", "sequence": "MEDIPDTSRPPLKYVK", "type": "FASTA",
             "smiles": "OC1=CC=C(C[C@@H](C(O)=O)N)C=C1"},
        ]},
        "oed-catpred": {"input_pairs": [
            {"name": "example", "sequence": "MEDIPDTSRPPLKYVK", "type": "FASTA",
             "smiles": "OC1=CC=C(C[C@@H](C(O)=O)N)C=C1"},
        ]},
        "cleandb-mepesm": {"sequence": "MEDIPDTSRPPLKYVK"},
        "chemscraper": {"input_file": "paper.pdf"},
        "molli": {"CORES_FILE_NAME": "cores.cdxml", "SUBS_FILE_NAME": "subs.cdxml"},
        "ml-simplefold": {"fasta": ">seq\nMEDIPDTSRPPLKYVK"},
        "clean": {"input_fasta": [
            {"header": "seq1", "sequence": "MEDIPDTSRPPLKYVK", "DNA_sequence": ""},
        ]},
        "somn": [{
            "reactant_pair_name": "pair1",
            "el": "CC(C)c1ccc(Br)cc1", "el_name": "el1", "el_input_type": "smi", "el_idx": "0",
            "nuc": "NCc1ccccc1", "nuc_name": "nuc1", "nuc_input_type": "smi", "nuc_idx": "0",
        }],
        "ez-specificity": {
            "enzymes": [{"filename": "enzyme.pdb"}],
            "substrates": ["CCO"],
        },
        "novostoic-pathways": {
            "substrate": {"amount": 1, "molecule": "MNXM732866"},
            "product": {"amount": 1, "molecule": "MNXM5188"},
            "reactants": [{"amount": 1, "molecule": "MNXM10"}],
            "products": [{"amount": 1, "molecule": "MNXM8"}],
            "max_steps": 2,
            "iterations": 10,
        },
        "reactionminer": {},
    }

    @pytest.mark.parametrize("identifier,payload", sorted(CASES.items()))
    def test_known_good_payload_validates(self, identifier, payload):
        schema = tool_registry.get_input_schema(identifier)

        Draft202012Validator(schema).validate(payload)

    def test_every_tool_with_a_schema_has_a_case(self):
        """Otherwise a schema could be wrong in a way nothing notices."""
        with_schema = {
            k for k, v in tool_registry.list_tools().items() if v.get("input_schema") is not None
        }

        assert with_schema - set(self.CASES) == set()


class TestSchemasRejectBadInput:
    @pytest.mark.parametrize("identifier,payload", [
        ("aceretro", {}),                                    # missing smiles
        ("aceretro", {"smiles": ""}),                        # empty smiles
        ("chemscraper", {"input_file": "x.pdf", "extra": 1}),  # unknown property
        ("somn", []),                                        # empty array
        ("oed-cheminfo", {"query_smiles": "C", "smiles_list": ["C"],
                          "algorithms": ["nope"]}),          # unknown algorithm
        ("ez-specificity", {"enzymes": [], "substrates": ["CCO"]}),   # too few enzymes
        ("ez-specificity", {"enzymes": [{"filename": "a.pdb"}],
                            "substrates": ["C"] * 11}),      # too many substrates
    ])
    def test_invalid_payload_is_rejected(self, identifier, payload):
        schema = tool_registry.get_input_schema(identifier)

        with pytest.raises(Exception):
            Draft202012Validator(schema).validate(payload)
