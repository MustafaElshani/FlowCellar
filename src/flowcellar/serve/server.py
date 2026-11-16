"""``flowcellar serve``: hand the hub of a library out over HTTP on this computer.

A read-only static server from the standard library. It knows four addresses and asks the viewer's Python
(``report.api.serve_file``) for each; nothing a request says is ever joined to a path on disk here.

* binds 127.0.0.1 unless told otherwise;
* answers only requests whose ``Host`` is this computer by a loopback name (127.0.0.1, localhost, [::1]), with
  any port: the name is what stops another web page from reaching it through a made-up host name, and a port
  forward may well use another local port. Any other ``Host`` gets a 403 that names the address to use, unless
  ``allow_lan``;
* ``GET`` and ``HEAD`` only: anything else is 405; no editing, no API, no directory listing;
* every answer carries ``nosniff``, ``no-referrer``, ``same-origin`` opener policy, ``no-store`` and a ban on
  being shown inside another page.
"""
from __future__ import annotations

import re
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from flowcellar import __version__
from flowcellar.report import api as report_api
from flowcellar.util.errors import FlowCellarError, InputError, UsageError

LOOPBACK_HOSTS = ("127.0.0.1", "localhost", "[::1]")
LOOPBACK_BIND = ("127.0.0.1", "localhost", "::1")
SECURITY_HEADERS = {
    "X-Content-Type-Options": "nosniff",
    "Referrer-Policy": "no-referrer",
    "Cross-Origin-Opener-Policy": "same-origin",
    "Cache-Control": "no-store",
    "X-Frame-Options": "DENY",                 # stands in for frame-ancestors 'none', which a page policy cannot carry
}
_PATH_OK = re.compile(r"^/[A-Za-z0-9._/-]{0,160}$")


def clean_path(raw: str) -> str | None:
    """The path of a request line when it is one this server may be asked for; None otherwise. The query string
    is dropped (the viewer adds one to fetch a fresh index). Nothing encoded, no ``..``, no backslash."""
    path = raw.split("?", 1)[0].split("#", 1)[0]
    if not _PATH_OK.match(path) or ".." in path or "//" in path:
        return None
    return path


def host_allowed(host_header: str | None, port: int, allow_lan: bool) -> bool:
    """True for a loopback name with any port (``127.0.0.1:9000`` behind ``ssh -L 9000:127.0.0.1:8737``). ``port``
    is the port listened on; it is not compared, the name is the protection."""
    if allow_lan:
        return True
    if not host_header:
        return False
    h = host_header.strip().lower()
    m = re.fullmatch(r"(\[[0-9a-f:]+\]|[a-z0-9.-]+)(?::(\d{1,5}))?", h)
    return bool(m) and m.group(1) in LOOPBACK_HOSTS


class _Handler(BaseHTTPRequestHandler):
    server_version = "FlowCellar/" + __version__
    sys_version = ""
    protocol_version = "HTTP/1.1"

    def _send(self, status: int, body: bytes, content_type: str, extra: dict | None = None, head: bool = False) -> None:
        self.send_response(status)
        headers = {**SECURITY_HEADERS, **(extra or {})}
        headers["Content-Type"] = content_type
        headers["Content-Length"] = str(len(body))
        for k, v in headers.items():
            self.send_header(k, v)
        self.end_headers()
        if not head:
            self.wfile.write(body)

    def _plain(self, status: int, text: str, head: bool = False, extra: dict | None = None) -> None:
        self._send(status, (text + "\n").encode("utf-8"), "text/plain; charset=utf-8", extra, head)

    def _answer(self, head: bool) -> None:
        srv = self.server
        if not host_allowed(self.headers.get("Host"), srv.server_address[1], srv.allow_lan):
            port = srv.server_address[1]
            self._plain(403, "This FlowCellar hub answers only requests addressed to this computer by the name 127.0.0.1 or "
                             f"localhost. Open http://127.0.0.1:{port}/ on the computer it runs on; from another computer forward "
                             f"the port first (ssh -N -L {port}:127.0.0.1:{port} USER@SERVER) and open the same address there.", head)
            return
        path = clean_path(self.path)
        if path is None:
            self._plain(404, "Not found.", head)
            return
        try:
            served = report_api.serve_file(srv.library, path)
        except FlowCellarError as e:
            self._plain(500, f"FlowCellar could not answer: {e.message}", head)
            return
        except Exception:  # noqa: BLE001  never show a traceback to a browser
            self._plain(500, "FlowCellar could not answer.", head)
            return
        if served is None:
            self._plain(404, "Not found.", head)
            return
        self._send(served.status, served.body, served.content_type, {**served.headers, **SECURITY_HEADERS}, head)

    def do_GET(self) -> None:  # noqa: N802
        self._answer(head=False)

    def do_HEAD(self) -> None:  # noqa: N802
        self._answer(head=True)

    def _refuse(self) -> None:
        length = self.headers.get("Content-Length")
        if length and length.isdigit():
            self.rfile.read(min(int(length), 1 << 16))       # drain a small body so the connection stays usable
        self._plain(405, "This hub is read-only: GET and HEAD only.", extra={"Allow": "GET, HEAD"})
        self.close_connection = True

    do_POST = do_PUT = do_DELETE = do_PATCH = do_OPTIONS = do_TRACE = do_CONNECT = _refuse

    def log_message(self, fmt: str, *args) -> None:
        if self.server.quiet:
            return
        from flowcellar.util import log
        log.event("serve.request", "%s %s" % (self.address_string(), fmt % args), level=10)


class HubServer(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True

    def __init__(self, library, host: str = "127.0.0.1", port: int = 8737, *, allow_lan: bool = False, quiet: bool = False):
        if host not in LOOPBACK_BIND and not allow_lan:
            raise UsageError(f"--host {host} makes the hub reachable from other computers",
                             hint="add --allow-lan if that is what you want; there is no password")
        self.library = library
        self.allow_lan = allow_lan
        self.quiet = quiet
        if ":" in host:
            import socket
            self.address_family = socket.AF_INET6
        try:
            super().__init__((host, port), _Handler)
        except OSError as e:
            raise InputError(f"cannot listen on {host}:{port} ({e.strerror or e})",
                             hint="another program uses the port; try --port with another number") from None

    @property
    def url(self) -> str:
        host, port = self.server_address[0], self.server_address[1]
        shown = "127.0.0.1" if host in ("0.0.0.0", "127.0.0.1") else f"[{host}]" if ":" in host else host
        return f"http://{shown}:{port}/"

    def start_in_thread(self) -> threading.Thread:
        t = threading.Thread(target=self.serve_forever, name="flowcellar-serve", daemon=True)
        t.start()
        return t
