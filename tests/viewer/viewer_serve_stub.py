"""A minimal read-only HTTP server around ``flowcellar.report.api.serve_file`` for the viewer's browser checks.

    python tests/viewer/viewer_serve_stub.py LIBRARY_DIR [PORT]

It answers GET and HEAD with exactly what ``serve_file`` returns (body, type, headers) and 404 otherwise. It is a
test double for ``flowcellar serve``, which calls the same function; it binds
127.0.0.1 only and prints the address it listens on.
"""
from __future__ import annotations

import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from viewer_fake_library import FakeLibrary  # noqa: E402

from flowcellar.report import api  # noqa: E402


def handler(library):
    class Handler(BaseHTTPRequestHandler):
        def _answer(self, body: bool):
            served = api.serve_file(library, self.path)
            if served is None:
                self.send_error(404)
                return
            self.send_response(served.status)
            self.send_header("Content-Type", served.content_type)
            self.send_header("Content-Length", str(len(served.body)))
            for k, v in served.headers.items():
                self.send_header(k, v)
            self.end_headers()
            if body:
                self.wfile.write(served.body)

        def do_GET(self):          # noqa: N802
            self._answer(True)

        def do_HEAD(self):         # noqa: N802
            self._answer(False)

        def log_message(self, *args):
            pass

    return Handler


def serve(library_dir: str, port: int = 0) -> ThreadingHTTPServer:
    return ThreadingHTTPServer(("127.0.0.1", port), handler(FakeLibrary(library_dir)))


if __name__ == "__main__":
    srv = serve(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    print(f"http://127.0.0.1:{srv.server_address[1]}/", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
