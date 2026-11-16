#!/usr/bin/env python3
"""Try an installed FlowCellar from the outside, on invented data: the check a wheel must pass before a release.

    /path/to/environment/bin/python tools/smoke_test.py [--work DIR] [--installed] [--allow-network]

Run it with the Python of the environment under test. It writes fake run folders (and, when pysam is there,
fake read data) with the generators the package ships, then starts the real console program for every step:

    flowcellar --version, doctor, schema
    flowcellar run <one run folder> -o qc/           a bundle and one offline report
    flowcellar scan <three folders> -L library       a new library with its hub
    flowcellar scan -L library                       a second scan
    flowcellar validate <bundle>, validate library, ls, show, discover, hub --reports, serve

Every child process runs with a guard that refuses and records any attempt to reach another computer (a
connection or a name lookup that is not this computer's own loopback address), so the test also shows that
nothing needs the network at run time. ``--installed`` also demands that the package is imported from the
environment and not from a source tree. Exit code 0 when every step holds.

Standard library only; nothing of this file is part of the package.
"""
from __future__ import annotations

import argparse
import http.client
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

GUARD = '''\
import os, sys
_LOG = os.environ.get("FC_SMOKE_NETLOG")
_LOCAL = ("127.0.0.1", "::1", "localhost", "", "0.0.0.0", None)
def _note(kind, target, refused):
    if _LOG:
        fd = os.open(_LOG, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        os.write(fd, ("%s\\t%s\\t%s\\n" % ("refused" if refused else "local", kind, target)).encode())
        os.close(fd)
def _hook(event, args):
    if event in ("socket.connect", "socket.sendto", "socket.sendmsg"):
        address = args[1] if len(args) > 1 else None
        if isinstance(address, tuple) and address:
            local = address[0] in _LOCAL
            _note(event, address[0], not local)
            if not local:
                raise OSError("smoke test: FlowCellar tried to reach %r" % (address[0],))
    elif event in ("socket.getaddrinfo", "socket.gethostbyname", "socket.gethostbyaddr"):
        host = args[0]
        host = host.decode() if isinstance(host, bytes) else host
        local = host in _LOCAL
        _note(event, host, not local)
        if not local:
            raise OSError("smoke test: FlowCellar tried to look up %r" % (host,))
sys.addaudithook(_hook)
'''

CONFIG = """\
[library]
name = "smoke test library"

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


class Smoke:
    def __init__(self, work: Path, guard: bool):
        self.work = work
        self.failed: list[str] = []
        self.steps = 0
        self.env = {k: v for k, v in os.environ.items() if not k.startswith(("FLOWCELLAR_", "PYTHON"))}
        self.netlog = work / "network.log"
        if guard:
            guard_dir = work / "guard"
            guard_dir.mkdir(parents=True)
            (guard_dir / "sitecustomize.py").write_text(GUARD)
            self.env.update(PYTHONPATH=str(guard_dir), FC_SMOKE_NETLOG=str(self.netlog))
        self.exe = [sys.executable, "-m", "flowcellar"]
        script = Path(sys.executable).parent / "flowcellar"
        if script.exists():
            self.exe = [str(script)]

    def run(self, *args: str, timeout: float = 900) -> subprocess.CompletedProcess:
        return subprocess.run([*self.exe, *args], cwd=self.work, env=self.env, capture_output=True, text=True, timeout=timeout)

    def python(self, code: str, *args: str) -> subprocess.CompletedProcess:
        return subprocess.run([sys.executable, "-c", code, *args], cwd=self.work, env=self.env, capture_output=True, text=True, timeout=900)

    def check(self, what: str, ok: bool, detail: str = "") -> bool:
        self.steps += 1
        print(f"  {'ok  ' if ok else 'FAIL'}  {what}" + (f"\n        {detail.strip()[-1500:]}" if detail and not ok else ""), flush=True)
        if not ok:
            self.failed.append(what)
        return ok


def other_origins(page: str) -> list[str]:
    return [u for u in re.findall(r"https?://[^\s\"'<>)]+", page) if not u.startswith("http://www.w3.org/")]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work", type=Path, help="folder to work in (default: a temporary one, removed afterwards)")
    ap.add_argument("--installed", action="store_true", help="the package must be imported from this environment, not from a source tree")
    ap.add_argument("--allow-network", action="store_true", help="run without the guard that refuses the network")
    ap.add_argument("--expect-version", help="the version the program must report")
    args = ap.parse_args(argv)

    tmp = None
    if args.work:
        work = args.work.resolve()
        if work.exists() and any(work.iterdir()):
            raise SystemExit(f"{work} is not empty")
        work.mkdir(parents=True, exist_ok=True)
    else:
        tmp = tempfile.mkdtemp(prefix="flowcellar-smoke-")
        work = Path(tmp)
    s = Smoke(work, guard=not args.allow_network)
    print(f"smoke test with {sys.executable} in {work}", flush=True)
    try:
        # ---- the installation ------------------------------------------------------------------------
        out = s.python("import flowcellar, sys; print(flowcellar.__file__); print(flowcellar.__version__); print(sys.prefix)")
        where, version, prefix = (out.stdout.splitlines() + ["", "", ""])[:3]
        s.check("the package can be imported", out.returncode == 0, out.stderr)
        if args.installed:
            s.check(f"it is imported from this environment ({os.path.relpath(where, prefix) if where.startswith(prefix) else where})",
                    where.startswith(prefix + os.sep))
        out = s.run("--version")
        s.check(f"flowcellar --version  ->  {out.stdout.strip()}", out.returncode == 0 and out.stdout.strip() == f"flowcellar {version}", out.stderr)
        if args.expect_version:
            s.check(f"the version is {args.expect_version}", version == args.expect_version, version)
        out = s.run("doctor")
        viewer_line = next((ln.strip() for ln in out.stdout.splitlines() if "viewer" in ln), "")
        s.check(f"flowcellar doctor  ->  {viewer_line}", out.returncode == 0 and viewer_line.startswith("ok"), out.stdout + out.stderr)
        out = s.run("doctor", "--format", "json")
        try:
            packages = json.loads(out.stdout)["packages"]
        except (ValueError, KeyError):
            packages = {}
        s.check("doctor --format json names the packages: " + ", ".join(f"{k} {v}" for k, v in packages.items() if v),
                all(packages.get(k) for k in ("numpy", "pyarrow", "click")), out.stdout + out.stderr)
        have_pysam = bool(packages.get("pysam"))
        have_jsonschema = bool(packages.get("jsonschema"))

        def validate(what: str, target: str) -> None:
            out = s.run("validate", target)
            if have_jsonschema:
                s.check(f"flowcellar validate {what}", out.returncode == 0, out.stdout + out.stderr)
            else:
                s.check(f"flowcellar validate {what}: without jsonschema it stops and names the extra to install",
                        out.returncode == 4 and "flowcellar[validate]" in out.stdout + out.stderr, out.stdout + out.stderr)
        out = s.run("schema")
        try:
            schema_id = json.loads(out.stdout).get("$id")
        except ValueError:
            schema_id = None
        s.check(f"flowcellar schema  ->  {schema_id}", out.returncode == 0 and bool(schema_id), out.stderr)
        out = s.python("from flowcellar.report import api; from flowcellar import catalogue; "
                       "print(sorted(api.viewer_manifest()['files'])); print(len(catalogue.metrics()))")
        s.check("the built viewer and the catalogue are inside the package: " + " · ".join(out.stdout.splitlines()),
                out.returncode == 0 and "viewer.core.js" in out.stdout, out.stderr)

        # ---- invented data -----------------------------------------------------------------------------
        make = (
            "import json, os, sys\n"
            "from flowcellar.testing import trees\n"
            "tree = trees.make_three_roots('data', reads=1200, hours=8.0, aged=False)\n"
            "out = {'roots': tree.root_list(), 'runs': {k: [r.run_key, r.run_dir, r.expected['reads']] for k, r in tree.runs.items()}}\n"
            "if sys.argv[1] == 'reads':\n"
            "    from flowcellar.testing import readdata as rd\n"
            "    a3 = tree.runs['a3']\n"
            "    p = rd.make_project(os.path.join(tree.roots['store'], 'group_one', 'alpha_project'), rd.run_ref(a3), reads_per_sample=250, unclassified=40)\n"
            "    out['read_data'] = {'run': a3.run_key, 'basecalled': p.expected['basecalled'], 'samples': len(p.aliases)}\n"
            "trees.age_tree(tree.base)\n"
            "print(json.dumps(out))\n"
        )
        out = s.python(make, "reads" if have_pysam else "none")
        if not s.check("invented run folders are written" + (" with read data" if have_pysam else " (no pysam: without read data)"),
                       out.returncode == 0, out.stderr):
            return 1
        data = json.loads(out.stdout.splitlines()[-1])
        runs = data["runs"]

        # ---- one run -------------------------------------------------------------------------------------
        key, run_dir, reads = runs["a1"] if "a1" in runs else next(iter(runs.values()))
        out = s.run("run", run_dir, "-o", "qc")
        made = sorted(n for n in os.listdir(work / "qc") if not n.startswith(".")) if (work / "qc").is_dir() else []
        s.check(f"flowcellar run  ->  {', '.join(made)}", out.returncode == 0 and made == [f"{key}.fcb.json.gz", f"{key}.flowcellar.html"],
                out.stdout + out.stderr)
        report = work / "qc" / f"{key}.flowcellar.html"
        if report.exists():
            page = report.read_text("utf-8")
            s.check(f"the report is one self-contained page ({len(page) / 1024:,.0f} kB) that names no other origin",
                    page.lower().startswith("<!doctype html") and 'id="fc-boot"' in page and "Content-Security-Policy" in page
                    and other_origins(page) == [] and key in page, str(other_origins(page)[:3]))
        validate("<bundle>", str(work / "qc" / f"{key}.fcb.json.gz"))

        # ---- a library ---------------------------------------------------------------------------------
        lib = work / "library"
        lib.mkdir()
        (lib / "flowcellar.toml").write_text(CONFIG)
        t0 = time.monotonic()
        out = s.run("scan", *data["roots"], "-L", str(lib))
        first_seconds = time.monotonic() - t0
        done = next((ln.strip() for ln in out.stderr.splitlines() if ln.strip().startswith("done:")), "")
        s.check(f"flowcellar scan of {len(data['roots'])} folders ({first_seconds:.0f} s)  ->  {done}",
                out.returncode == 0 and f"done: {len(runs) + 1} runs" in done and "hub ready:" in out.stderr and "Traceback" not in out.stderr,
                out.stderr)
        hub = lib / "hub"
        sidecars = sorted(os.listdir(hub / "data")) if (hub / "data").is_dir() else []
        s.check(f"the hub is written: index.html, index.js and {len(sidecars)} run files",
                (hub / "index.html").is_file() and (hub / "index.js").is_file() and len(sidecars) == len(runs) + 1)
        if (hub / "index.html").is_file():
            s.check("the hub page names no other origin", other_origins((hub / "index.html").read_text("utf-8")) == [])
        validate("<library>", str(lib))
        out = s.run("ls", "-L", str(lib), "--format", "json", "--columns", "id,reads,detail_t1,read_data,t2_reads,t2_samples,verdict")
        try:
            rows = {r["id"]: r for r in json.loads(out.stdout)}
        except (ValueError, KeyError, TypeError):
            rows = {}
        s.check(f"flowcellar ls  ->  {len(rows)} runs", out.returncode == 0 and len(rows) == len(runs) + 1, out.stdout[:500] + out.stderr)
        exact = {k: rows.get(k, {}).get("reads") for k, _dir, _n in runs.values()}
        s.check(f"every run was counted in full and the read counts are exact ({', '.join(str(v) for v in exact.values())})",
                all(exact[k] == n and rows[k]["detail_t1"] == "full" for k, _dir, n in runs.values()), json.dumps(exact))
        if "read_data" in data:
            row = rows.get(data["read_data"]["run"], {})
            s.check(f"read data is linked to its run and every read is counted ({row.get('t2_reads')} reads, {row.get('t2_samples')} samples, "
                    f"state {row.get('read_data')})",
                    row.get("t2_reads") == data["read_data"]["basecalled"] and row.get("t2_samples") == data["read_data"]["samples"]
                    and row.get("read_data") in ("complete", "partial"), json.dumps(row))
        out = s.run("show", key, "-L", str(lib), "--status")
        s.check("flowcellar show --status lists the checks", out.returncode == 0 and "run_health" in out.stdout and "fc_occupancy_first_hour" in out.stdout,
                out.stdout + out.stderr)
        t0 = time.monotonic()
        out = s.run("scan", "-L", str(lib))
        s.check(f"a second scan with nothing new ({time.monotonic() - t0:.1f} s)", out.returncode == 0 and "done:" in out.stderr, out.stderr)
        out = s.run("discover", "-L", str(lib), "--format", "json")
        s.check("flowcellar discover", out.returncode == 0 and out.stdout.lstrip().startswith(("{", "[")), out.stderr)
        out = s.run("hub", str(lib), "--reports")
        reports = sorted(os.listdir(lib / "reports")) if (lib / "reports").is_dir() else []
        s.check(f"flowcellar hub --reports  ->  {len(reports)} reports", out.returncode == 0 and len(reports) >= len(runs), out.stderr)

        # ---- serve: this computer only -------------------------------------------------------------------
        with socket.socket() as probe:
            probe.bind(("127.0.0.1", 0))
            port = probe.getsockname()[1]
        server = subprocess.Popen([*s.exe, "serve", str(lib), "--port", str(port)], cwd=work, env=s.env,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            answers = {}
            deadline = time.monotonic() + 30
            while time.monotonic() < deadline and not answers:
                try:
                    for path in ("/", "/index.js", f"/data/{key}.js", "/flowcellar.toml"):
                        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
                        conn.request("GET", path)
                        r = conn.getresponse()
                        answers[path] = (r.status, r.getheader("Content-Security-Policy") or "", len(r.read()))
                        conn.close()
                except OSError:
                    answers = {}
                    time.sleep(0.3)
            ok = bool(answers) and [answers[p][0] for p in ("/", "/index.js", f"/data/{key}.js", "/flowcellar.toml")] == [200, 200, 200, 404] \
                and "default-src 'none'" in answers["/"][1]
            s.check("flowcellar serve answers /, /index.js and a run file on 127.0.0.1, and nothing else", ok, str(answers))
        finally:
            server.terminate()
            try:
                server.wait(timeout=10)
            except subprocess.TimeoutExpired:
                server.kill()

        # ---- the network ---------------------------------------------------------------------------------
        if not args.allow_network:
            lines = s.netlog.read_text().splitlines() if s.netlog.exists() else []
            refused = [ln for ln in lines if ln.startswith("refused")]
            s.check(f"no step tried to reach another computer ({len(lines) - len(refused)} loopback events of serve, {len(refused)} refused)",
                    not refused, "\n".join(refused[:10]))
    finally:
        if tmp and not s.failed:
            shutil.rmtree(tmp, ignore_errors=True)
    print(f"{s.steps - len(s.failed)} of {s.steps} steps passed" + (f"; FAILED: {'; '.join(s.failed)}" if s.failed else ""), flush=True)
    return 1 if s.failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
