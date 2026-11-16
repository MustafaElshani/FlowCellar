"""flowcellar serve: hand the hub out over HTTP on this computer."""
from __future__ import annotations

import os

import click


@click.command("serve", short_help="Show a library in the browser, served from this computer.")
@click.argument("library", type=click.Path(file_okay=False))
@click.option("--host", default=None, help="Address to listen on. [default: 127.0.0.1, or [serve] host]")
@click.option("--port", type=click.IntRange(1, 65535), default=None, help="[default: 8737, or [serve] port]")
@click.option("--allow-lan", is_flag=True, help="Accept requests from other computers (needed with a --host other than 127.0.0.1).")
@click.option("--open", "open_", is_flag=True, help="Open the hub in the default browser.")
@click.pass_context
def cmd(ctx, library, host, port, allow_lan, open_):
    """Serve the hub of LIBRARY: the same files 'flowcellar hub' writes, read-only. Stops with Ctrl-C.

    Only this computer can reach it unless --allow-lan is given. There is no editing and no password.
    """
    from flowcellar.library.store import LibraryStore
    from flowcellar.report import api as report_api
    from flowcellar.serve import HubServer
    from flowcellar.util.errors import InputError
    store = LibraryStore.open(library, config_file=(ctx.obj or {}).get("config"))
    host = host or store.config["serve"]["host"]
    port = port or int(store.config["serve"]["port"])
    if report_api.serve_file(store, "/index.html") is None:
        raise InputError(f"{store.path}: the hub has not been written yet", hint=f"flowcellar hub {store.path}")
    server = HubServer(store, host, port, allow_lan=allow_lan, quiet=(ctx.obj or {}).get("verbosity", 0) < 1)
    click.echo(f"serving {store.name} at {server.url}   (read-only; Ctrl-C stops it)", err=True)
    if allow_lan:
        click.echo("  --allow-lan: other computers on the network can read this library; there is no password", err=True)
    elif os.environ.get("SSH_CONNECTION") or os.environ.get("SSH_TTY"):
        # a session on a server: the address above is the server's own. Say how to reach it from the person's computer.
        click.echo(f"  to look at it from your own computer, forward the port:  ssh -N -L {server.server_address[1]}:127.0.0.1:"
                   f"{server.server_address[1]} USER@THIS_SERVER   then open http://127.0.0.1:{server.server_address[1]}/", err=True)
    if open_:
        import webbrowser
        webbrowser.open(server.url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        click.echo("stopped", err=True)
    finally:
        server.server_close()
