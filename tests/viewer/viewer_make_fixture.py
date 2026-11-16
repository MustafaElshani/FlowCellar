"""Everything the browser checks open, written through the VIEWER's Python from the synthetic library.

    python tests/viewer/viewer_make_fixture.py OUT_DIR

    OUT_DIR/library/            the synthetic library (viewer_synth) with hub/ and reports/
    OUT_DIR/library/hub/index.next.js   a newer index, as a running scan would write it (for the Reload check)
    OUT_DIR/empty/              a library without runs, with its hub
    OUT_DIR/noindex/index.html  a hub page with no index beside it
    OUT_DIR/report/<run_key>.flowcellar.html   reports written from a bundle alone (flowcellar run -o)
"""
from __future__ import annotations

import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import viewer_synth  # noqa: E402
from viewer_fake_library import FakeLibrary  # noqa: E402

from flowcellar.report import api  # noqa: E402
from flowcellar.report import html as page  # noqa: E402
from flowcellar.util import atomic  # noqa: E402


def make(out: str) -> dict:
    out = os.path.abspath(out)
    shutil.rmtree(out, ignore_errors=True)
    lib_dir = os.path.join(out, "library")
    made = viewer_synth.make_library(lib_dir)
    lib = FakeLibrary(lib_dir)
    api.write_hub(lib, reports=True)
    scanning = viewer_synth.make_library(os.path.join(out, "_scanning"), scanning=True)
    scanning["index"]["written"] = "2026-11-13T12:30:00Z"
    with open(os.path.join(lib_dir, "hub", "index.next.js"), "wb") as fh:
        fh.write(page.index_sidecar(atomic.gzip_bytes(atomic.dumps(scanning["index"]).encode("utf-8"))))
    shutil.rmtree(os.path.join(out, "_scanning"))
    empty_dir = os.path.join(out, "empty")
    viewer_synth.make_library(empty_dir, empty=True, name="Empty demo library")
    api.write_hub(FakeLibrary(empty_dir))
    os.makedirs(os.path.join(out, "noindex"))
    atomic.atomic_write_text(os.path.join(out, "noindex", "index.html"), page.hub_page())
    reports = {}
    os.makedirs(os.path.join(out, "report"))
    for spec, bundle in zip(made["specs"], made["bundles"]):
        if spec.key in ("TST00001_0a1b2c3d", "TST00009_0e1f2a3b", "TST00006_7b8c9d0e", "TST00004_3d4e5f6a"):
            reports[spec.key] = api.write_report(bundle, os.path.join(out, "report", spec.key + ".flowcellar.html"))
    return {"out": out, "library": lib_dir, "hub": os.path.join(lib_dir, "hub"), "empty": os.path.join(empty_dir, "hub"),
            "noindex": os.path.join(out, "noindex"), "reports": reports, "runs": [s.key for s in made["specs"]]}


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    print(json.dumps(make(sys.argv[1])))
