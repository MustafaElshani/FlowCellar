"""flowcellar schema: print the bundle schema."""
from __future__ import annotations

import click

from flowcellar.util.errors import UsageError


@click.command("schema", short_help="Print the JSON Schema of a bundle.")
@click.option("--schema-version", metavar="X.Y.Z", help="A published schema version. [default: the current one]")
@click.option("--list", "list_", is_flag=True, help="List the schema versions in this package.")
@click.pass_context
def cmd(ctx, schema_version, list_):
    """Prints the schema this FlowCellar writes (flowcellar.model.schema.schema_text)."""
    from flowcellar import SCHEMA_VERSION
    from flowcellar.model import schema
    if list_:
        for v in schema.SCHEMA_VERSIONS:
            click.echo(v + ("  (current)" if v == SCHEMA_VERSION else ""))
        return
    try:
        click.echo(schema.schema_text(schema_version or SCHEMA_VERSION), nl=False)
    except KeyError:
        raise UsageError(f"no bundle schema {schema_version} in this package",
                         hint="flowcellar schema --list shows the versions") from None
