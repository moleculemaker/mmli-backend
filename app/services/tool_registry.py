"""Loads the tool registry: scientific metadata plus per-tool input schemas.

The registry is data, not code. `app/cfg/tools.yaml` describes each job type for a
human or an agent choosing between them, and `app/schemas/*.input.json` describes what
each one accepts. Neither knows anything about Kubernetes; runtime facts stay in
config.yaml.

Schemas are deliberately self-contained -- every `$ref` resolves inside its own file --
so a schema can be handed to a client, an MCP tool definition, or a registry crawler
without a bundling step.

`$id` is not stored in the files, because the correct value depends on the hostname the
API is served from. Callers that need one should set it when serving.
"""
import json
import os
import pathlib
from typing import Any, Dict, Optional

import yaml

from config import get_logger
from models.enums import JobTypes

log = get_logger(__name__)

_APP_DIR = pathlib.Path(__file__).resolve().parent.parent
TOOLS_FILEPATH = pathlib.Path(os.getenv('TOOLS_FILEPATH', _APP_DIR / 'cfg' / 'tools.yaml'))
SCHEMAS_DIR = pathlib.Path(os.getenv('SCHEMAS_DIR', _APP_DIR / 'schemas'))


class ToolRegistryError(RuntimeError):
    """Raised when the registry on disk is inconsistent with the application."""


def _load_yaml() -> Dict[str, Any]:
    with open(TOOLS_FILEPATH, 'r') as handle:
        document = yaml.safe_load(handle) or {}
    tools = document.get('tools')
    if not isinstance(tools, dict):
        raise ToolRegistryError(f'{TOOLS_FILEPATH} must contain a "tools" mapping')
    return tools


def _load_schema(filename: str) -> Dict[str, Any]:
    path = SCHEMAS_DIR / filename
    if not path.is_file():
        raise ToolRegistryError(f'input_schema file not found: {path}')
    with open(path, 'r') as handle:
        return json.load(handle)


def _build() -> Dict[str, Dict[str, Any]]:
    registry: Dict[str, Dict[str, Any]] = {}

    for identifier, entry in _load_yaml().items():
        entry = dict(entry or {})

        # A registry entry for a job type the application does not know about would be
        # undiscoverable and un-runnable, so treat it as a mistake rather than ignore it.
        if identifier not in JobTypes:
            raise ToolRegistryError(
                f'tools.yaml describes "{identifier}", which is not a JobType. '
                f'Add it to models/enums.py or remove the entry.'
            )

        schema_filename = entry.get('input_schema')
        entry['input_schema'] = _load_schema(schema_filename) if schema_filename else None
        entry['identifier'] = identifier
        entry.setdefault('internal', False)
        registry[identifier] = entry

    return registry


_REGISTRY: Optional[Dict[str, Dict[str, Any]]] = None


def get_registry() -> Dict[str, Dict[str, Any]]:
    """Return every registered tool, keyed by job type. Loaded once, then cached."""
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = _build()
        log.info(f'Loaded tool registry: {len(_REGISTRY)} tools from {TOOLS_FILEPATH}')
    return _REGISTRY


def get_tool(identifier: str) -> Optional[Dict[str, Any]]:
    """Return one tool's registry entry, or None if it is not registered."""
    return get_registry().get(identifier)


def list_tools(include_internal: bool = False) -> Dict[str, Dict[str, Any]]:
    """Return registered tools, omitting internal subjob types by default.

    Subjob types exist only so a parent job's coordinator can create them. Listing them
    as things to run would be misleading.
    """
    registry = get_registry()
    if include_internal:
        return dict(registry)
    return {k: v for k, v in registry.items() if not v.get('internal')}


def get_input_schema(identifier: str) -> Optional[Dict[str, Any]]:
    """Return a tool's JSON Schema, or None when it accepts no job_info."""
    tool = get_tool(identifier)
    return tool.get('input_schema') if tool else None


def unregistered_job_types() -> set:
    """Job types with no registry entry.

    Used by the test suite to make sure a newly added job type does not silently miss
    its metadata. `defaults` is excluded deliberately: it runs a Perl one-liner that
    computes pi and exists as a harness fixture, not as a scientific tool.
    """
    return {t for t in JobTypes if t != 'defaults'} - set(get_registry())
