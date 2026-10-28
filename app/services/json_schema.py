"""The JSON Schema validator /v1 uses for tool inputs and results.

JSON Schema's `pattern` is an ECMA-262 regular expression, but jsonschema evaluates it
with Python's `re.search`. The two disagree on `$`: in ECMA-262 (without the m flag) it
matches only at the very end of the string, while in Python it also matches just before
a final newline. So `^[a-z]+$` accepted "a\\n" here while the browser's RegExp, and so
every frontend, rejected it (#147).

Rather than fix each schema, the validator rewrites `$` to Python's `\\Z` (end of string
only) when it evaluates a pattern. Other ECMA/Python differences (`\\d` matching
non-ASCII digits in Python, for one) are not handled.
"""
import functools
import re

from jsonschema import Draft202012Validator, ValidationError, validators


@functools.lru_cache(maxsize=None)
def ecma_pattern(pattern: str) -> str:
    """Rewrite every unescaped `$` outside a character class to `\\Z`."""
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
