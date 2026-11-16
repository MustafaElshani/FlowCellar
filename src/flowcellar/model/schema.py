"""The bundle schema that ships inside the package, and validation against it.

Validation needs the ``jsonschema`` package (extra ``validate``); everything else here does not.
"""
from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources

from flowcellar import SCHEMA_VERSION
from flowcellar.util.errors import MissingExtraError

SCHEMA_VERSIONS = ("0.1.1",)          # every published schema version ships inside the package


def schema_text(version: str = SCHEMA_VERSION) -> str:
    """The JSON Schema as text: what ``flowcellar schema`` prints."""
    if version not in SCHEMA_VERSIONS:
        raise KeyError(f"no bundle schema {version} in this package (have: {', '.join(SCHEMA_VERSIONS)})")
    return resources.files("flowcellar.model").joinpath("schemas", f"bundle-{version}.schema.json").read_text("utf-8")


@lru_cache(maxsize=4)
def load_schema(version: str = SCHEMA_VERSION) -> dict:
    return json.loads(schema_text(version))


@lru_cache(maxsize=4)
def _validator(version: str):
    try:
        from jsonschema import Draft202012Validator
    except ImportError:
        raise MissingExtraError("jsonschema", "validate", "Schema validation") from None
    return Draft202012Validator(load_schema(version))


def schema_errors(bundle: dict, *, limit: int = 50) -> list[str]:
    """Human-readable schema errors, each with the JSON pointer of the place; empty when the bundle is valid."""
    out = []
    for e in sorted(_validator(SCHEMA_VERSION).iter_errors(bundle), key=lambda e: list(e.absolute_path)):
        where = "/" + "/".join(str(p) for p in e.absolute_path)
        out.append(f"{where}: {e.message[:300]}")
        if len(out) >= limit:
            break
    return out
