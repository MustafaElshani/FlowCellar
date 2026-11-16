"""flowcellar serve: a read-only server for the hub on this computer."""
from __future__ import annotations

import http.client
import os
import socket

import pytest

from flowcellar.library.store import LibraryStore
from flowcellar.serve import HubServer
from flowcellar.serve.server import SECURITY_HEADERS, clean_path, host_allowed
from flowcellar.util.errors import InputError, UsageError


@pytest.fixture
def server(scanned, lib_dir):
    srv = HubServer(LibraryStore.open(lib_dir), "127.0.0.1", 0, quiet=True)      # port 0: any free port
    thread = srv.start_in_thread()
    yield srv
    srv.shutdown()
    srv.server_close()
    thread.join(5)


def ask(server, method: str, path: str, host: str | None = "default", body: bytes | None = None):
    port = server.server_address[1]
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    headers = {} if host is None else {"Host": f"127.0.0.1:{port}" if host == "default" else host}
    conn.putrequest(method, path, skip_host=True, skip_accept_encoding=True)
    for k, v in headers.items():
        conn.putheader(k, v)
    if body is not None:
        conn.putheader("Content-Length", str(len(body)))
    conn.endheaders(body)
    resp = conn.getresponse()
    data = resp.read()
    conn.close()
    return resp.status, dict(resp.getheaders()), data


def test_the_four_addresses_and_their_headers(server, tree, lib_dir, parts):
    key = tree.runs["a1"].run_key
    hub = os.path.join(lib_dir, "hub")
    for path, ctype, want in (("/", "text/html; charset=utf-8", open(os.path.join(hub, "index.html"), "rb").read()),
                              ("/index.html", "text/html; charset=utf-8", open(os.path.join(hub, "index.html"), "rb").read()),
                              ("/index.js", "text/javascript; charset=utf-8", open(os.path.join(hub, "index.js"), "rb").read()),
                              (f"/data/{key}.js", "text/javascript; charset=utf-8", None)):
        status, headers, body = ask(server, "GET", path)
        assert status == 200 and headers["Content-Type"] == ctype, path
        if want is not None:
            assert body == want
        else:
            assert body.startswith(f'__FCB_SIDECAR__("{key}","'.encode())
        for k, v in SECURITY_HEADERS.items():
            assert headers[k] == v, (path, k)
        assert headers["Content-Security-Policy"] == parts.report.CSP
        assert headers["Content-Length"] == str(len(body))
        h_status, h_headers, h_body = ask(server, "HEAD", path)
        assert h_status == 200 and h_body == b"" and h_headers["Content-Length"] == str(len(body))
    # the viewer adds a counter to fetch a fresh index while a scan runs: the query string is ignored
    assert ask(server, "GET", "/index.js?t=17")[0] == 200


def test_everything_else_is_not_found_and_nothing_can_be_written(server, lib_dir):
    for path in ("/flowcellar.toml", "/LIBRARY", "/bundles/", "/data/", "/data/../flowcellar.toml", "/data/%2e%2e/flowcellar.toml",
                 "/..%2f..%2fetc/passwd", "/index.html/", "/data//index.html", "/data/TST00099_ffffffff.js", "/hub/index.html",
                 "/index.js\\..\\x", "/" + "a" * 400):
        status, headers, body = ask(server, "GET", path)
        assert status == 404, path
        assert b"flowcellar.toml" not in body and b"[[roots]]" not in body
        assert headers["X-Content-Type-Options"] == "nosniff" and headers["Cache-Control"] == "no-store"
    before = sorted(os.listdir(lib_dir))
    for method in ("POST", "PUT", "DELETE", "PATCH", "OPTIONS"):
        status, headers, body = ask(server, method, "/index.js", body=b'{"x": 1}' if method in ("POST", "PUT", "PATCH") else None)
        assert status == 405 and headers["Allow"] == "GET, HEAD", method
    assert sorted(os.listdir(lib_dir)) == before


def test_only_this_computer_is_answered(server):
    port = server.server_address[1]
    # a loopback name with any port: a port forward may use another local port (ssh -L 9000:127.0.0.1:8737)
    for host in (f"localhost:{port}", f"127.0.0.1:{port}", f"[::1]:{port}", f"LOCALHOST:{port}", "127.0.0.1", f"127.0.0.1:{port + 1}",
                 "localhost:9000"):
        assert ask(server, "GET", "/index.js", host=host)[0] == 200, host
    for host in ("evil.example", f"evil.example:{port}", f"127.0.0.1.evil.example:{port}", "localhost.evil.example", "127.0.0.1:x",
                 "127.0.0.1@evil.example", None):
        status, headers, body = ask(server, "GET", "/index.js", host=host)
        assert status == 403 and b"__FC_INDEX__" not in body, host
        # an answer, with the address to use: never a closed connection
        assert f"http://127.0.0.1:{port}/".encode() in body and b"ssh -N -L" in body, body
    assert server.server_address[0] == "127.0.0.1"
    assert host_allowed("192.168.1.5:8737", 8737, allow_lan=True) and not host_allowed("192.168.1.5:8737", 8737, allow_lan=False)


def test_listening_for_other_computers_must_be_asked_for(scanned, lib_dir):
    store = LibraryStore.open(lib_dir)
    with pytest.raises(UsageError) as e:
        HubServer(store, "0.0.0.0", 0)
    assert "--allow-lan" in e.value.hint
    srv = HubServer(store, "0.0.0.0", 0, allow_lan=True, quiet=True)
    port = srv.server_address[1]
    srv.server_close()
    busy = socket.socket()
    busy.bind(("127.0.0.1", 0))
    busy.listen(1)
    try:
        with pytest.raises(InputError) as e:
            HubServer(store, "127.0.0.1", busy.getsockname()[1])
        assert e.value.exit_code == 4 and "--port" in e.value.hint
    finally:
        busy.close()
    assert port > 0


def test_path_cleaning():
    assert clean_path("/") == "/" and clean_path("/index.js?t=3") == "/index.js" and clean_path("/data/TST00001_0a1b2c3d.js") == "/data/TST00001_0a1b2c3d.js"
    for bad in ("", "index.js", "/a b", "/%2e%2e/x", "/../x", "/a/../b", "//x", "/a\\b", "/\x00", "http://x/index.js", "/é"):
        assert clean_path(bad) is None, bad


def test_serve_command_refuses_a_library_without_a_hub(parts, tree, lib_dir, capsys):
    from library_helpers import run_scan

    from flowcellar.cli import main
    run_scan(lib_dir, tree.roots["arch"], no_hub=True)
    with pytest.raises(SystemExit) as e:
        main(["serve", lib_dir, "--port", "1"])
    assert e.value.code in (1, 4) and "hub" in capsys.readouterr().err
