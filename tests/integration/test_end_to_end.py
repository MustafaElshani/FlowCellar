"""FlowCellar end to end, through the command line, on synthetic data only.

The journey a lab makes: ``scan`` three folders into a new library (run health and read data), open the hub,
``scan`` again, ``serve`` it, and ``run`` on one folder. Every step starts the real console program in its own
process, so nothing here can pass because of a test double.
"""
from __future__ import annotations

import base64
import gzip
import hashlib
import html
import http.client
import json
import os
import re
import socket
import subprocess
import sys
import time

import pytest

from flowcellar.model import schema as schema_model
from flowcellar.testing import readdata as rd
from flowcellar.testing import trees

pytestmark = [pytest.mark.slow]

CONFIG = """\
[library]
name = "demo library"

[[inputs]]
role = "dorado.bam"
glob = "**/*_basecalled/bam_pass/{sample}/*.bam"
[[inputs]]
role = "align.genome_bam"
glob = "**/analysis/alignments/{sample}_genome_sorted.bam"
[[inputs]]
role = "pychopper.stats"
glob = "**/analysis/pychopper_reports/{sample}_pychopper_stats.tsv"
[[inputs]]
role = "quant.salmon"
glob = "**/analysis/counts/salmon/{sample}_salmon/aux_info/meta_info.json"
[[inputs]]
role = "basecaller.sample_sheet"
sniff = "csv_columns:barcode,alias,flow_cell_id|position_id"

[[link.rules]]
name = "experiment-folder"
roles = ["pychopper.stats", "quant.salmon"]
via = "sample"
within = "experiment"
"""

# Runs the console program with every file it opens written to a list (an audit hook sees each ``open`` of this
# process and of the children it forks). Directory listings and ``stat`` calls are not opens.
AUDITED = """\
import os, sys
fd = os.open(sys.argv[1], os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
def hook(event, args):
    if event == "open" and isinstance(args[0], (str, bytes)):
        p = os.fsdecode(args[0])
        os.write(fd, (os.path.abspath(p) + "\\n").encode())
sys.addaudithook(hook)
from flowcellar.cli import main
main(sys.argv[2:])
"""


def flowcellar(*args: str, timeout: float = 600) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    return subprocess.run([sys.executable, "-m", "flowcellar", *args], capture_output=True, text=True, timeout=timeout, env=env)


def audited(log: str, *args: str, timeout: float = 600) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    return subprocess.run([sys.executable, "-c", AUDITED, log, *args], capture_output=True, text=True, timeout=timeout, env=env)


def opened_under(log: str, folders: list[str]) -> list[str]:
    if not os.path.exists(log):
        return []
    with open(log) as fh:
        paths = [line.rstrip("\n") for line in fh]
    return sorted({p for p in paths if any(p.startswith(f.rstrip("/") + "/") for f in folders)})


def payload(text: str, call: str) -> list[str]:
    """The string arguments of the one call a sidecar file holds: ``__FC_INDEX__("...");``."""
    m = re.fullmatch(re.escape(call) + r"\((.*)\);\n", text, re.S)
    assert m, f"not exactly one {call}(...) call"
    return json.loads("[" + m.group(1) + "]")


def other_origins(page: str) -> list[str]:
    """Every absolute address in a page, apart from the XML namespace names the script and the logo carry."""
    return [u for u in re.findall(r"https?://[^\s\"'<>)]+", page) if not u.startswith("http://www.w3.org/")]


def unpack(b64: str) -> dict:
    return json.loads(gzip.decompress(base64.b64decode(b64)))


@pytest.fixture(scope="module")
def journey(tmp_path_factory):
    try:
        import pysam  # noqa: F401
    except ImportError:
        pytest.skip("needs pysam (the 'bam' extra) to write the fake read data")
    from flowcellar.report import api as report_api
    from flowcellar.util.errors import FlowCellarError
    try:
        report_api.viewer_manifest()
    except FlowCellarError as e:
        pytest.skip(f"the viewer is not built: {e.message}")
    base = tmp_path_factory.mktemp("journey")
    tree = trees.make_three_roots(base / "data", reads=1500, hours=8.0, aged=False)
    a3 = tree.runs["a3"]
    project = rd.make_project(os.path.join(tree.roots["store"], "group_one", "alpha_project"), rd.run_ref(a3),
                              reads_per_sample=300, unclassified=60)
    trees.age_tree(tree.base)                       # nothing is "still being copied"
    lib = str(base / "library")
    os.makedirs(lib)
    with open(os.path.join(lib, "flowcellar.toml"), "w") as fh:
        fh.write(CONFIG)
    t0 = time.monotonic()
    first = flowcellar("scan", *tree.root_list(), "-L", lib)
    return {"base": str(base), "tree": tree, "project": project, "lib": lib, "first": first, "seconds": time.monotonic() - t0}


def index_of(lib: str) -> dict:
    with gzip.open(os.path.join(lib, "index", "index.json.gz")) as fh:
        return json.load(fh)


def column(idx: dict, name: str) -> list:
    col = idx["runs"]["columns"][name]
    return [col["dict"][c] for c in col["codes"]] if isinstance(col, dict) else col


def test_scan_of_three_roots_makes_a_library(journey):
    first, lib, tree = journey["first"], journey["lib"], journey["tree"]
    assert first.returncode == 0, first.stderr
    err = first.stderr
    assert "created library" in err and "(3 roots" in err
    assert "found 6 runs (5 with MinKNOW reports, 1 with raw data only)" in err
    assert "hub ready:" in err and err.index("hub ready:") < err.index("done:")
    assert re.search(r"done: 6 runs · .* · read data: \d+ of \d+ files read", err), err
    for line in err.splitlines():
        assert "Traceback" not in line
    assert first.stdout == "", "events go to stderr"
    with open(os.path.join(lib, "flowcellar.toml")) as fh:
        toml = fh.read()
    for root in tree.root_list():
        assert f'path = "{root}"' in toml, "the roots are remembered"
    idx = index_of(lib)
    assert idx["format"] == "flowcellar.index" and idx["runs"]["n"] == 6
    ids = column(idx, "id")
    by_id = {k: i for i, k in enumerate(ids)}
    for run in tree.runs.values():
        i = by_id[run.run_key]
        assert column(idx, "reads")[i] == run.expected["reads"], "the read count comes from the full pass"
        assert column(idx, "detail_t1")[i] == "full"
    assert sorted(column(idx, "detail_t1"), key=str) == ["full"] * 5 + ["none"]
    projects = {p["id"]: p for p in idx["projects"]}
    assert projects["alpha_project"]["runs"] == 2 and projects["alpha_project"]["roots"] == ["root_seq"]
    assert projects["group_one"]["runs"] == 1, "by default a project is the first folder under its root"
    assert projects["_none"]["runs"] == 1, "a run directly under a root has no project"
    assert [r["id"] for r in idx["roots"]] == ["root_seq", "root_store", "root_arch"]
    assert idx["progress"]["state"] == "idle" and idx["library"]["counts"]["runs"] == 6
    assert idx["library"]["path"] == lib


def test_read_data_is_exact_and_linked_to_its_run(journey):
    lib, tree, project = journey["lib"], journey["tree"], journey["project"]
    idx = index_of(lib)
    i = column(idx, "id").index(tree.runs["a3"].run_key)
    assert column(idx, "read_data")[i] in ("complete", "partial")
    assert column(idx, "t2_reads")[i] == project.expected["basecalled"], "every read of the Dorado files, none sampled"
    assert column(idx, "t2_samples")[i] == len(project.aliases)
    others = [v for j, v in enumerate(column(idx, "read_data")) if j != i]
    assert set(others) == {"absent"}
    names = {s["name"] for s in idx["samples"] if s["run"] == tree.runs["a3"].run_key}
    assert set(project.aliases) <= names


def test_the_hub_and_every_sidecar_are_written(journey):
    lib = journey["lib"]
    hub = os.path.join(lib, "hub")
    with open(os.path.join(hub, "index.html"), encoding="utf-8") as fh:
        page = fh.read()
    assert page.lower().startswith("<!doctype html") and '"mode": "hub"' in page.replace('"mode":"hub"', '"mode": "hub"')
    assert other_origins(page) == [], "the page names no other origin"
    with open(os.path.join(hub, "index.js"), encoding="utf-8") as fh:
        (b64,) = payload(fh.read(), "__FC_INDEX__")
    served = unpack(b64)
    assert served["runs"]["n"] == 6 and served == index_of(lib)
    ids = column(served, "id")
    assert sorted(os.listdir(os.path.join(hub, "data"))) == sorted(k + ".js" for k in ids), "one sidecar per run, no other file"
    for key in ids:
        with open(os.path.join(hub, "data", key + ".js"), encoding="utf-8") as fh:
            got_key, b64, sha = payload(fh.read(), "__FCB_SIDECAR__")
        assert got_key == key and hashlib.sha256(base64.b64decode(b64)).hexdigest() == sha
        bundle = unpack(b64)
        assert bundle["bundle_id"] == key and "status" in bundle, "a sidecar holds the composed bundle"
        assert schema_model.schema_errors(bundle) == [], key


def test_every_bundle_of_the_library_validates(journey):
    out = flowcellar("validate", journey["lib"])
    assert out.returncode == 0, out.stdout + out.stderr


def test_the_second_scan_reads_no_input(journey):
    lib, tree = journey["lib"], journey["tree"]
    keep = ("bundles", "status", "cards", "partials")
    def stamps():
        return {os.path.join(d, n): os.stat(os.path.join(d, n)).st_mtime_ns
                for sub in keep for d, _, names in os.walk(os.path.join(lib, sub)) for n in names}
    before = stamps()
    log = os.path.join(journey["base"], "opens_second_scan.txt")
    t0 = time.monotonic()
    second = audited(log, "scan", "-L", lib)
    seconds = time.monotonic() - t0
    assert second.returncode == 0, second.stderr
    assert "nothing new" in second.stderr
    assert opened_under(log, tree.root_list()) == [], "a scan with nothing new opens no file under any root"
    assert opened_under(log, [lib]), "the audit does see the files a scan opens (its own library)"
    assert stamps() == before, "no bundle, status file, card or partial is rewritten"
    logs = sorted(os.listdir(os.path.join(lib, "logs")))
    with open(os.path.join(lib, "logs", logs[-1])) as fh:
        done = [json.loads(x) for x in fh if '"scan.done"' in x][-1]
    for what in ("identify", "quick", "full_pass", "read_pass", "link", "assemble", "judge", "bundles_written"):
        assert done["counts"][what] == 0, what
    assert seconds < max(10.0, journey["seconds"]), "and it is quick"


def test_the_audit_sees_a_read_of_a_run_folder(journey, tmp_path):
    """The control of the test above: when a run folder IS read, the same audit reports it."""
    run = journey["tree"].runs["b1"]
    log = str(tmp_path / "opens_run.txt")
    out = audited(log, "run", run.run_dir, "-o", str(tmp_path / "qc"), "--quick")
    assert out.returncode == 0, out.stderr
    assert opened_under(log, [run.run_dir]) != []


def test_run_writes_a_bundle_and_one_offline_report(journey, tmp_path):
    run = journey["tree"].runs["a1"]
    out_dir = str(tmp_path / "qc")
    out = flowcellar("run", run.run_dir, "-o", out_dir)
    assert out.returncode == 0, out.stderr
    lines = out.stderr.splitlines()
    assert lines[0].startswith(run.run_key), "the first line names the run"
    assert any("quick report" in x for x in lines) and any("full report" in x for x in lines)
    assert any(x.strip().startswith("Run health") for x in lines) and any(x.strip().startswith("Read data") for x in lines)
    report = os.path.join(out_dir, run.run_key + ".flowcellar.html")
    bundle_file = os.path.join(out_dir, run.run_key + ".fcb.json.gz")
    assert sorted(n for n in os.listdir(out_dir) if not n.startswith(".")) == sorted(os.path.basename(p) for p in (report, bundle_file))
    with gzip.open(bundle_file) as fh:
        bundle = json.load(fh)
    assert bundle["bundle_id"] == run.run_key and bundle["status"]["verdict"] in ("pass", "warn", "fail", "incomplete")
    assert bundle["modules"]["yield"]["summary"]["yield_reads"] == run.expected["reads"]
    assert flowcellar("validate", bundle_file).returncode == 0
    with open(report, encoding="utf-8") as fh:
        page = fh.read()
    assert page.lower().startswith("<!doctype html") and other_origins(page) == []
    boot = json.loads(re.search(r'<script type="application/json" id="fc-boot">(.*?)</script>', page, re.S).group(1))
    assert boot["mode"] == "report" and boot["run"] == run.run_key
    m = re.search(r'<script type="application/vnd\.flowcellar\.bundle\+gzip" id="fcb:([^"]+)" data-sha256="([0-9a-f]{64})">([A-Za-z0-9+/=]+)</script>', page)
    assert m and m.group(1) == run.run_key
    assert hashlib.sha256(base64.b64decode(m.group(3))).hexdigest() == m.group(2)
    embedded = unpack(m.group(3))
    assert embedded["bundle_id"] == run.run_key and embedded["status"]["verdict"] == bundle["status"]["verdict"]
    policy = html.unescape(re.search(r'<meta http-equiv="Content-Security-Policy" content="([^"]*)"', page).group(1))
    assert "connect-src 'none'" in policy and "default-src 'none'" in policy and "'self'" not in policy


def test_run_with_fail_on_trips_the_gate(journey, tmp_path):
    """A fake run is far below the fixed output levels, so its verdict is at least a warning: exit code 3."""
    run = journey["tree"].runs["a1"]
    out = flowcellar("run", run.run_dir, "-o", str(tmp_path / "gate"), "--quick", "--no-report", "--fail-on", "incomplete")
    with gzip.open(os.path.join(str(tmp_path / "gate"), run.run_key + ".fcb.json.gz")) as fh:
        verdict = json.load(fh)["status"]["verdict"]
    assert out.returncode == (3 if verdict != "pass" else 0), out.stderr


def test_run_on_something_else_is_an_input_error(journey, tmp_path):
    out = flowcellar("run", str(tmp_path), "-o", str(tmp_path / "x"))
    assert out.returncode == 4 and "Traceback" not in out.stderr


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _get(port: int, path: str, host: str | None = None) -> tuple[int, dict, bytes]:
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    try:
        conn.request("GET", path, headers={"Host": host} if host else {})
        resp = conn.getresponse()
        return resp.status, dict(resp.getheaders()), resp.read()
    finally:
        conn.close()


def test_serve_answers_on_this_computer_only(journey):
    lib, tree = journey["lib"], journey["tree"]
    port = _free_port()
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    proc = subprocess.Popen([sys.executable, "-m", "flowcellar", "serve", lib, "--port", str(port)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    try:
        deadline = time.monotonic() + 30
        while True:
            try:
                status, headers, body = _get(port, "/")
                break
            except OSError:
                assert proc.poll() is None, proc.stderr.read()
                assert time.monotonic() < deadline, "the server did not start"
                time.sleep(0.1)
        assert status == 200 and body[:15].lower() == b"<!doctype html>" and headers["Content-Type"].startswith("text/html")
        with open(os.path.join(lib, "hub", "index.html"), "rb") as fh:
            assert body == fh.read(), "the same page the static hub holds"
        for name in ("Content-Security-Policy", "X-Content-Type-Options", "Referrer-Policy", "Cache-Control"):
            assert name in headers, name
        assert "connect-src 'none'" in headers["Content-Security-Policy"]
        status, headers, body = _get(port, "/index.js")
        assert status == 200 and headers["Content-Type"].startswith("text/javascript")
        assert unpack(payload(body.decode(), "__FC_INDEX__")[0])["runs"]["n"] == 6
        key = tree.runs["a3"].run_key
        status, _, body = _get(port, f"/data/{key}.js")
        assert status == 200 and unpack(payload(body.decode(), "__FCB_SIDECAR__")[1])["bundle_id"] == key
        for path in ("/flowcellar.toml", "/data/../flowcellar.toml", "/bundles/" + key + ".fcb.json.gz", "/data/", "/hub/index.html"):
            assert _get(port, path)[0] == 404, path
        assert _get(port, "/", host="example.org")[0] in (400, 403), "a request for another host name is refused"
        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
        conn.request("POST", "/", body=b"x")
        assert conn.getresponse().status == 405
        conn.close()
    finally:
        proc.terminate()
        try:
            proc.wait(10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(10)
    assert proc.returncode is not None, "the server is stopped"


def test_a_new_root_is_added_later_and_ls_and_show_read_the_library(journey, tmp_path):
    lib = journey["lib"]
    from flowcellar.testing.fixtures import RunSpec, make_run
    extra_root = tmp_path / "root_new"
    new = make_run(extra_root / "gamma_project", RunSpec(reads=300, hours=2.0, flow_cell_id="TST00021", experiment="gamma_project",
                                                         sample="pool_g", start="2026-05-02T09:00:00Z", seed=31))
    trees.age_tree(str(extra_root))
    out = flowcellar("scan", str(extra_root), "-L", lib)
    assert out.returncode == 0, out.stderr
    assert "added root" in out.stderr and "found 1 new run" in out.stderr
    idx = index_of(lib)
    assert idx["runs"]["n"] == 7 and len(idx["roots"]) == 4
    assert any(p["id"] == "gamma_project" for p in idx["projects"])
    assert os.path.isfile(os.path.join(lib, "hub", "data", new.run_key + ".js"))
    ls = flowcellar("ls", "-L", lib, "--format", "tsv")
    assert ls.returncode == 0 and new.run_key in ls.stdout and len(ls.stdout.strip().splitlines()) == 8
    show = flowcellar("show", new.run_key, "-L", lib, "--status")
    assert show.returncode == 0 and "run_stop_cause" in show.stdout
    assert flowcellar("doctor", "-L", lib).returncode == 0


def test_a_project_under_two_roots_is_one_project(journey):
    """The second root keeps its projects one level down. Saying so in flowcellar.toml regroups the runs without
    reading or computing anything again."""
    lib = journey["lib"]
    path = os.path.join(lib, "flowcellar.toml")
    with open(path) as fh:
        text = fh.read()
    head, sep, tail = text.partition('id = "root_store"')
    assert sep and "project_depth = 1" in tail
    with open(path, "w") as fh:
        fh.write(head + sep + tail.replace("project_depth = 1", "project_depth = 2", 1))
    bundles = os.path.join(lib, "bundles")
    before = {n: os.stat(os.path.join(bundles, n)).st_mtime_ns for n in os.listdir(bundles)}
    out = flowcellar("scan", "-L", lib)
    assert out.returncode == 0, out.stderr
    assert {n: os.stat(os.path.join(bundles, n)).st_mtime_ns for n in os.listdir(bundles)} == before
    projects = {p["id"]: p for p in index_of(lib)["projects"]}
    assert "group_one" not in projects
    assert projects["alpha_project"]["runs"] == 3 and projects["alpha_project"]["roots"] == ["root_seq", "root_store"]
    assert projects["alpha_project"]["samples"] >= 3, "the samples of its read data are counted with the project"
