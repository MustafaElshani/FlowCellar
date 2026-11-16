"""The viewer build: ``node viewer/build.mjs`` turns viewer/ and the package catalogue into viewer_dist/.

These tests need Node and the packages of viewer/node_modules; they are skipped where those are missing (the
package itself never needs Node: it ships the built files).
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
VIEWER = REPO / "viewer"
pytestmark = pytest.mark.viewer


def _node():
    node = shutil.which("node")
    if not node or not (VIEWER / "node_modules" / "esbuild").is_dir():
        pytest.skip("Node or viewer/node_modules is missing (run 'npm ci' in viewer/)")
    return node


def _build(out: Path, *flags: str) -> dict:
    r = subprocess.run([_node(), str(VIEWER / "build.mjs"), "--out", str(out), *flags], cwd=VIEWER, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    return {p.name: p.read_bytes() for p in sorted(out.iterdir())}


def test_two_builds_give_identical_files(tmp_path):
    one, two = _build(tmp_path / "a"), _build(tmp_path / "b")
    assert sorted(one) == ["favicon.svg", "manifest.json", "template.html", "viewer.core.js", "viewer.css"]
    assert one == two
    manifest = json.loads(one["manifest.json"])
    assert "built" not in manifest                                  # nothing in the output depends on the time
    for name in ("viewer.core.js", "viewer.css", "favicon.svg", "template.html"):
        assert manifest["files"][name]["sha256"] == hashlib.sha256(one[name]).hexdigest()
        assert manifest["files"][name]["bytes"] == len(one[name])


def test_the_built_viewer_in_the_package_is_up_to_date(tmp_path):
    dist = REPO / "src" / "flowcellar" / "viewer_dist"
    if not (dist / "manifest.json").exists():
        pytest.skip("the viewer is not built: run 'node viewer/build.mjs'")
    have = json.loads((dist / "manifest.json").read_text())
    fresh = json.loads(_build(tmp_path / "fresh", *([] if have.get("minified", True) else ["--dev"]))["manifest.json"])
    assert {k: v["sha256"] for k, v in have["files"].items()} == {k: v["sha256"] for k, v in fresh["files"].items()}, \
        "src/flowcellar/viewer_dist is older than viewer/: run 'node viewer/build.mjs'"


def test_labels_come_from_the_package_catalogue(tmp_path):
    _build(tmp_path / "x")
    gen = (VIEWER / "src" / "gen" / "catalogue.js").read_text("utf-8")
    cat = json.loads((REPO / "src" / "flowcellar" / "catalogue" / "metrics_nanopore.json").read_text("utf-8"))
    data = json.loads(gen[gen.index("export default ") + len("export default "): gen.rindex(";")])
    assert data["version"] == cat["catalogue_version"]
    labelled = {m["id"]: m["display"]["label"] for m in cat["metrics"] + cat.get("derived", []) if (m.get("display") or {}).get("label")}
    assert labelled and {k: v["label"] for k, v in data["metrics"].items() if "label" in v} == labelled
    assert [p["id"] for p in data["panels"]] == [p["id"] for p in cat["panels"]]


def test_the_build_reads_only_product_files():
    text = (VIEWER / "build.mjs").read_text("utf-8")
    code = "\n".join(line.split("//", 1)[0] for line in text.splitlines())
    assert not re.search(r"""['"`/](design|prototype|local)['"`/]""", code)
    reads = re.findall(r"path\.join\((here|repo)((?:, '[^']+')+)\)", code)
    tops = {(base, args.split("'")[1]) for base, args in reads}
    assert tops <= {("here", "src"), ("here", "styles"), ("here", "brand"), ("repo", "src")}, tops


def test_the_viewer_sources_build_no_markup_from_text():
    pat = re.compile(r"innerHTML|outerHTML|insertAdjacentHTML|dangerouslySetInnerHTML|document\.write|createContextualFragment|DOMParser|srcdoc")
    hits = [str(p.relative_to(REPO)) for p in (VIEWER / "src").rglob("*.js*") if pat.search(p.read_text("utf-8"))]
    assert hits == []


def test_the_viewer_sources_do_not_call_the_product_a_prototype():
    hits = []
    for p in list((VIEWER / "src").rglob("*.js*")) + list((VIEWER / "styles").glob("*.css")) + [VIEWER / "build.mjs", VIEWER / "package.json"]:
        if "gen" in p.parts:
            continue
        for i, line in enumerate(p.read_text("utf-8").splitlines(), 1):
            if re.search(r"\bprototype\b", line, re.I) and ".prototype" not in line:
                hits.append(f"{p.relative_to(REPO)}:{i}")
    assert hits == []


def test_storage_keys_carry_the_library_id():
    state = (VIEWER / "src" / "core" / "state.js").read_text("utf-8")
    assert "'fc:' + (libraryId || '_') + ':'" in state
    for p in (VIEWER / "src").rglob("*.js*"):
        assert "fc:prototype" not in p.read_text("utf-8"), p
        if p.name != "state.js":
            assert "localStorage" not in p.read_text("utf-8") and "sessionStorage" not in p.read_text("utf-8"), p
