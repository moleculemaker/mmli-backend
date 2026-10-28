"""The JSON Schema validator /v1 uses for tool inputs and results.

JSON Schema's `pattern` is an ECMA-262 regular expression, but jsonschema evaluates it
with Python's `re.search`. The two disagree on `$`: in ECMA-262 (without the m flag) it
matches only at the very end of the string, while in Python it also matches just before
a final newline. So `^[a-z]+$` accepted "a\\n" here while the browser's RegExp, and so
every frontend, rejected it (#147).

Rather than fix each schema, the validator rewrites `$` to Python's `\\Z` (end of string
only) when it evaluates a pattern. Not handled, and kept out of shipped schemas by
tests/test_tool_registry.py where that matters:
  * patterns with an inline `m` or `x` flag are left as written (ECMA has no such syntax);
  * `patternProperties`, which jsonschema also matches in additionalProperties and
    unevaluatedProperties, so overriding the one keyword would make them disagree;
  * other ECMA/Python differences: `.` matching `\\r`, `\\d` matching non-ASCII digits.
"""
import functools
import re

from attrs import fields
from jsonschema import Draft202012Validator, ValidationError, validators

# An inline flag group that sets m (`$` becomes a line anchor) or x (whitespace and
# `#` comments, inside which `[` is a literal) -- either breaks the scan below.
INLINE_FLAGS = re.compile(r'\(\?[a-zA-Z-]+[:)]')
_SCAN_BREAKING_FLAGS = re.compile(r'\(\?[a-zA-Z]*[mx][a-zA-Z-]*[:)]')


@functools.lru_cache(maxsize=None)
def ecma_pattern(pattern: str) -> str:
    """Rewrite every unescaped `$` outside a character class to `\\Z`.

    Classes are read the way Python's `re` reads them, since `re` evaluates the result:
    a `]` first in a class (`[]a]`, `[^]a]`) is a literal, not the end.
    """
    if _SCAN_BREAKING_FLAGS.search(pattern):
        return pattern
    out = []
    in_class = False
    i = 0
    while i < len(pattern):
        char = pattern[i]
        if char == '\\':
            out.append(pattern[i:i + 2])
            i += 2
            continue
        if in_class:
            if char == ']':
                in_class = False
        elif char == '[':
            in_class = True
            out.append(char)
            i += 1
            if pattern.startswith('^', i):
                out.append('^')
                i += 1
            if pattern.startswith(']', i):
                out.append(']')
                i += 1
            continue
        elif char == '$':
            char = r'\Z'
        out.append(char)
        i += 1
    return ''.join(out)


def _pattern(validator, pattern, instance, schema):
    # Same as jsonschema._keywords.pattern (4.23.0), with the pattern rewritten
    if validator.is_type(instance, 'string') and not re.search(ecma_pattern(pattern), instance):
        yield ValidationError(f'{instance!r} does not match {pattern!r}')


SchemaValidator = validators.extend(Draft202012Validator, {'pattern': _pattern})


def _evolve(self, **changes):
    """The stock evolve (jsonschema 4.23.0), keeping this class for 2020-12 subschemas.

    jsonschema re-picks the class for every subschema it descends into, by the
    subschema's `$schema`. Every root schema here declares 2020-12, so a `$ref` to a
    root (`"#"`) got the stock Draft202012Validator and lost the `$` rule. Replacing
    evolve on this class only, rather than registering the class for the 2020-12 URI,
    leaves validator_for -- and the MCP server, which validates with jsonschema in the
    same process -- unchanged.
    """
    cls = self.__class__
    schema = changes.setdefault('schema', self.schema)
    new_validator = validators.validator_for(schema, default=cls)
    if new_validator is Draft202012Validator:
        new_validator = cls
    for field in fields(cls):
        if field.init and field.alias not in changes:
            changes[field.alias] = getattr(self, field.name)
    return new_validator(**changes)


SchemaValidator.evolve = _evolve
