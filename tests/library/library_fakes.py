"""Stand-ins for the other three parts, so that the library is tested on its own.

The library calls ENGINE, READ DATA and VIEWER only through their ``api`` modules. These doubles implement
those contracts in the simplest way that still lets a test see what the library asked for and when: every call
is counted, every pass over a large input is recorded in ``reads``, and the history each run was judged against
is kept. They read the fake run folders of ``flowcellar.testing.fixtures`` and invented read-data files
("fake BAMs": text files with one read per line).
"""
from __future__ import annotations

import base64
import fnmatch
import os
import time
from collections import Counter

from flowcellar.contracts import Candidate, InputRule, LinkedInput, PartialRef, ReadData, RunIdentity, Served
from flowcellar.model import bundle as bundle_model
from flowcellar.model import partial as partial_model
from flowcellar.model import status as status_model
from flowcellar.model.runfolder import TIER1_ROLES
from flowcellar.model.runkey import make_run_key, run_key_from_folder_name
from flowcellar.util import atomic
from flowcellar.util.errors import InputError
from flowcellar.util.fingerprint import fingerprint, hash_obj
from flowcellar.util.timeutil import iso_z, now_z

SUMMARY = "minknow.sequencing_summary"
PORE_SCAN = "minknow.pore_scan"


class FakeEngine:
    def __init__(self):
        self.calls: Counter = Counter()
        self.reads: list[str] = []                 # large inputs a pass really read, in order
        self.judged: list[tuple[str, list[str]]] = []   # (run key, the history keys it was given)
        self.fail_identify: set[str] = set()       # run folders (basename) that raise
        self.fail_pass: set[str] = set()           # run keys whose pass raises
        self.warn_flow_cells: set[str] = set()     # flow cells whose runs get a Warning
        self.on_pass = None                        # called before a pass reads (tests move the clock here)
        self.stale_after_s = 6 * 3600

    # ---------------------------------------------------------------------------------------------- identify
    def identify(self, run_dir: str) -> RunIdentity:
        self.calls["identify"] += 1
        run_dir = os.path.abspath(run_dir)
        if os.path.basename(run_dir) in self.fail_identify:
            raise InputError(f"{run_dir}: cannot be read (test)")
        try:
            names = sorted(os.listdir(run_dir))
        except OSError as e:
            raise InputError(f"{run_dir}: {e.strerror}") from None
        files, sizes, newest = {}, {}, 0.0
        for role, pat in TIER1_ROLES.items():
            sub, _, pat_name = pat.rpartition("/")
            folder = os.path.join(run_dir, sub) if sub else run_dir
            pool = sorted(os.listdir(folder)) if os.path.isdir(folder) else []
            hit = [n for n in pool if fnmatch.fnmatchcase(n, pat_name)]
            if hit:
                p = os.path.join(folder, hit[0])
                files[role], sizes[role] = p, os.path.getsize(p)
        for n in names:
            p = os.path.join(run_dir, n)
            if os.path.isfile(p):
                newest = max(newest, os.path.getmtime(p))
        final = files.get("minknow.final_summary")
        sheet = files.get("minknow.sample_sheet")
        if final:
            kv = dict(line.strip().split("=", 1) for line in open(final) if "=" in line)
            fc, prot, started = kv["flow_cell_id"], kv["protocol_run_id"], iso_z(kv.get("started"), 0)
            facts = {"experiment_id": kv.get("protocol_group_id"), "sample_id": kv.get("sample_id")}
            state = "finished"
        elif sheet:
            head, row = [line.strip().split(",") for line in open(sheet)][:2]
            d = dict(zip(head, row))
            fc, prot, started = d["flow_cell_id"], d["protocol_run_id"], None
            facts = {"experiment_id": d.get("experiment_id"), "sample_id": d.get("sample_id"), "kit": d.get("kit")}
            state = "running" if time.time() - newest < self.stale_after_s else "interrupted"
        else:
            key = run_key_from_folder_name(run_dir)
            if not key or not any(n.startswith(("pod5", "fast5")) or n == "other_reports" for n in names):
                raise InputError(f"{run_dir}: not a run folder")
            fc, short = key.split("_")
            return RunIdentity(run_key=key, flow_cell_id=fc, run_short=short, run_dir=run_dir, kind="run_raw_only",
                               run_state="raw_only", newest_mtime=newest)
        if sheet and "kit" not in facts:
            head, row = [line.strip().split(",") for line in open(sheet)][:2]
            facts["kit"] = dict(zip(head, row)).get("kit")
        return RunIdentity(run_key=make_run_key(fc, prot), flow_cell_id=fc, run_short=prot[:8], run_dir=run_dir, kind="run",
                           run_state=state, protocol_run_id=prot, started=started, files=files, sizes=sizes,
                           newest_mtime=newest, facts=facts)

    # ---------------------------------------------------------------------------------------------- bundles
    def _envelope(self, run: RunIdentity, ctx) -> dict:
        ident = {"run_key": run.run_key, "flow_cell": {"id": run.flow_cell_id},
                 "protocol_run": {"id": run.protocol_run_id, "short_id": run.run_short, "started": run.started,
                                  "experiment_id": run.facts.get("experiment_id"), "sample_id": run.facts.get("sample_id")}}
        b = bundle_model.new_bundle(ident)
        b["provenance"]["run_dir"] = {"path": run.run_dir, "root": ctx.root_id, "relpath": ctx.relpath}
        b["provenance"]["library"] = {"id": ctx.library_id, "name": ctx.library_name} if ctx.library_id else None
        b["x_fake"] = {"kind": run.kind, "run_state": run.run_state, "kit": run.facts.get("kit"), "reads": None}
        return b

    def quick_bundle(self, run: RunIdentity, ctx=None) -> dict:
        from flowcellar.contracts import RunContext
        self.calls["quick_bundle"] += 1
        b = self._envelope(run, ctx or RunContext())
        if run.kind == "run":
            b["tiers"]["run_health"] = {"state": "live" if run.run_state == "running" else "complete", "detail": "quick"}
        return b

    def full_pass(self, run: RunIdentity, partials_dir: str, opts=None) -> list[PartialRef]:
        self.calls["full_pass"] += 1
        if run.run_key in self.fail_pass:
            raise RuntimeError("the worker died (test)")
        refs = []
        for role in (SUMMARY, PORE_SCAN):
            path = run.files.get(role)
            if not path:
                continue
            fp = fingerprint(path)
            pid = partial_model.partial_id(fp, role, {"fake": "1"})
            hit = partial_model.find_partial(partials_dir, pid)
            if hit:
                refs.append(PartialRef(pid, hit, role, path, "ok", True, True))
                continue
            if self.on_pass:
                self.on_pass(run, role)
            self.reads.append(path)
            with open(path) as fh:
                n = sum(1 for _ in fh) - 1
            p = partial_model.new_partial(role=role, input_path=path, fp=fp, size=os.path.getsize(path), mtime=None,
                                          modules={"fake": "1"})
            p["coverage"].update(records=n, records_total=n)
            p["routing"]["runs"] = {run.run_key: n}
            p["fragments"] = [partial_model.fragment(run.run_key, "yield", facts={"records": n, "role": role})]
            out = partial_model.write_partial(partials_dir, p)
            refs.append(PartialRef(pid, out, role, path, "ok", True, False, 0.01))
        return refs

    def assemble(self, run: RunIdentity, partials=(), read_data=None, ctx=None) -> dict:
        from flowcellar.contracts import RunContext
        self.calls["assemble"] += 1
        b = self.quick_bundle(run, ctx or RunContext())
        self.calls["quick_bundle"] -= 1
        ids = []
        for ref in partials:
            p = partial_model.read_partial(ref.path)
            ids.append(ref.partial_id)
            if ref.role == SUMMARY and ref.complete:
                b["x_fake"]["reads"] = p["fragments"][0]["facts"]["records"]
                b["tiers"]["run_health"] = {"state": "complete", "detail": "full"}
        if read_data is not None:
            b["modules"].update(read_data.modules)
            b["samples"] = list(read_data.samples)
            b["tiers"]["read_level"] = dict(read_data.read_level)
            b["provenance"]["inputs"] = list(read_data.inputs)
            ids.extend(read_data.partial_ids)
            b["x_fake"]["t2_reads"] = sum((m.get("x_reads") or 0) for m in read_data.modules.values())
            b["x_fake"]["waiting"] = len(read_data.waiting)
        b["provenance"]["partials"] = sorted(ids)
        return b

    # ---------------------------------------------------------------------------------------------- judging
    @staticmethod
    def _before(history, started):
        return sorted(c["columns"]["id"] for c in history
                      if c["columns"].get("started") and started and c["columns"]["started"] < started)

    def judge(self, bundle: dict, history, *, declared=None, profile=None, baseline_groups=(), flow_cell=None,
              in_library: bool = True) -> dict:
        self.calls["judge"] += 1
        key = bundle["bundle_id"]
        x = bundle["x_fake"]
        started = bundle["identity"]["protocol_run"].get("started")
        self.judged.append((key, sorted(c["columns"]["id"] for c in history)))
        earlier = self._before(history, started)
        fc = bundle["identity"]["flow_cell"]["id"]
        same_fc = [k for k in earlier if k.startswith(fc + "_")]
        detail = bundle["tiers"]["run_health"].get("detail")
        rules, cause, provisional, pending = [], None, False, []
        if x["kind"] != "run":
            health, cause, headline = "incomplete", "raw_only", "Only raw data was found for this run."
        elif fc in self.warn_flow_cells:
            health, headline = "warn", "Occupancy in the first hour: 61 %. The guide is 70 % or more."
            rules.append({"id": "fc_occupancy_first_hour", "role": "supporting", "tier": "run_health", "severity": "warn",
                          "outcome": "triggered", "message": headline})
        elif detail == "quick":
            health, headline, provisional, pending = "pass", "The flow cell ran well so far.", True, ["dyn_speed_median"]
        else:
            health, headline = "pass", f"The flow cell ran well ({x['reads']} reads)."
        data = "pass" if bundle["tiers"]["read_level"]["state"] != "absent" else "not_available"
        status = {"verdict": status_model.combine_verdicts(health, data), "headline": headline,
                  "by_tier": {"run_health": health, "read_level": data}, "rules": rules, "cause": cause,
                  "modifiers": [{"code": "reused_flow_cell", "message": "The flow cell was used before."}] if same_fc else [],
                  "provisional": provisional, "pending": pending, "evaluated": now_z(),
                  "judged": {"at": now_z(), "as_of": started, "history_hash": hash_obj(earlier, 16), "library": in_library}}
        rel = bundle_model.default_relations()
        rel["flow_cell_use"].update(index=len(same_fc) + 1, source="library" if in_library else "unknown", prior_runs=len(same_fc))
        return status_model.new_status_file(key, bundle.get("content_sha256"), status, rel)

    def card(self, bundle: dict, status_file: dict | None = None) -> dict:
        self.calls["card"] += 1
        x = bundle["x_fake"]
        pr = bundle["identity"]["protocol_run"]
        st = (status_file or {}).get("status") or {}
        why = status_model.reasons(st) if st else []
        detail = bundle["tiers"]["run_health"].get("detail")
        reads = x.get("reads")
        cols = {
            "id": bundle["bundle_id"], "fc": bundle["identity"]["flow_cell"]["id"],
            "use": ((status_file or {}).get("relations") or {}).get("flow_cell_use", {}).get("index", 1),
            "started": pr.get("started"), "hours": 2.0 if x["kind"] == "run" else None,
            "run_state": x["run_state"], "kit": x.get("kit"), "assay_group": "cdna" if (x.get("kit") or "").startswith("SQK-PC") else "gdna",
            "product": "FLO-MIN114", "device_type": "MINION_MK1D", "experiment": pr.get("experiment_id"),
            "run_name": pr.get("sample_id"),
            "verdict": st.get("verdict"), "health": (st.get("by_tier") or {}).get("run_health"),
            "data": (st.get("by_tier") or {}).get("read_level"), "provisional": bool(st.get("provisional")),
            "pending": len(st.get("pending") or []), "checks_done": 3, "cause": st.get("cause"),
            "modifiers": [m["code"] for m in st.get("modifiers") or []], "headline": st.get("headline"), "detail": None,
            "reasons": [r["message"] for r in why], "reason_levels": [r["severity"] for r in why],
            "reason_checks": [r["check"] for r in why], "reason_tiers": [r["tier"] for r in why],
            "reads": reads, "gb": round(reads * 900 / 1e9, 6) if reads else (0.0001 if x["kind"] == "run" else None),
            "pores": 1200 if x["kind"] == "run" else None, "occ_first": 0.9 if x["kind"] == "run" else None,
            "detail_t1": detail if x["kind"] == "run" else "none",
            "history_hash": (st.get("judged") or {}).get("history_hash"), "limits": "fake-1",
        }
        return {"columns": cols, "facts": {"started": pr.get("started"), "reads": reads}}

    def baselines_today(self, cards) -> dict:
        self.calls["baselines_today"] += 1
        return {"groups": []}

    def history_hash(self, card, history) -> str | None:
        return hash_obj(self._before(history, card["columns"].get("started")), 16)

    def install(self, monkeypatch) -> FakeEngine:
        from flowcellar.engine import api
        for name in ("identify", "quick_bundle", "full_pass", "assemble", "judge", "card", "baselines_today", "history_hash"):
            monkeypatch.setattr(api, name, getattr(self, name))
        return self


class FakeReadData:
    """Read data told by file name: ``<run_key>__<sample>.fakebam`` holds one read per line (its length)."""

    def __init__(self):
        self.calls: Counter = Counter()
        self.reads: list[str] = []
        self.fail: set[str] = set()                # basenames whose pass raises
        self.on_pass = None

    def default_input_rules(self) -> list[InputRule]:
        return []

    def sniff(self, path: str) -> str | None:
        self.calls["sniff"] += 1
        with open(path, errors="replace") as fh:
            head = fh.readline().strip().split(",")
        return "basecaller.sample_sheet" if {"barcode", "alias", "flow_cell_id"} <= set(head) else None

    def candidates_for_run(self, run, paths=(), explicit=None, sample_sheet=None) -> list[Candidate]:
        return []

    def link(self, candidates, runs, *, config=None, probes_dir=None) -> list[LinkedInput]:
        self.calls["link"] += 1
        known = {r.run_key for r in runs}
        out = []
        for c in candidates:
            name = os.path.basename(c.path)
            key = name.split("__", 1)[0]
            fp = fingerprint(c.path)
            big = name.endswith(".fakebam")
            li = LinkedInput(candidate=c, fingerprint=fp, run_keys=[key] if key in known else [],
                             link={"method": "bam_header_rg", "confidence": "certain", "evidence": "named after its run (test)"},
                             sample=c.sample, variant=c.variant, needs=["reads"] if big else ["collected"],
                             counted_for={key: ["reads"] if big else ["collected"]} if key in known else {},
                             partial_id=partial_model.partial_id(fp, c.role, {"fakerd": "1"}) if big else None)
            out.append(li)
        return out

    def read_pass(self, inp: LinkedInput, partials_dir: str, known_runs, opts=None) -> PartialRef | None:
        self.calls["read_pass"] += 1
        c = inp.candidate
        if inp.partial_id is None:
            return None
        if os.path.basename(c.path) in self.fail:
            raise RuntimeError("truncated file (test)")
        hit = partial_model.find_partial(partials_dir, inp.partial_id)
        if hit:
            return PartialRef(inp.partial_id, hit, c.role, c.path, "ok", True, True)
        if self.on_pass:
            self.on_pass(inp)
        self.reads.append(c.path)
        with open(c.path) as fh:
            lengths = [int(x) for x in fh.read().split()]
        p = partial_model.new_partial(role=c.role, input_path=c.path, fp=inp.fingerprint, size=c.size, mtime=None,
                                      modules={"fakerd": "1"}, sample=inp.sample)
        p["coverage"].update(records=len(lengths), records_total=len(lengths))
        p["fragments"] = [partial_model.fragment(inp.run_keys[0], "len", sample=inp.sample,
                                                 facts={"reads": len(lengths), "bases": sum(lengths)})]
        return PartialRef(inp.partial_id, partial_model.write_partial(partials_dir, p), c.role, c.path, "ok", True, False, 0.01)

    def read_data(self, run, linked, partials_dir, *, tier1=None, config=None) -> ReadData:
        self.calls["read_data"] += 1
        rd = ReadData()
        reads = bases = 0
        for i, li in enumerate(linked, 1):
            hit = partial_model.find_partial(partials_dir, li.partial_id) if li.partial_id else None
            state = "ok"
            if li.partial_id and hit is None:
                rd.waiting.append(li.candidate.path)
                state = "skipped"
            elif hit:
                f = partial_model.read_partial(hit)["fragments"][0]["facts"]
                reads, bases = reads + f["reads"], bases + f["bases"]
                rd.partial_ids.append(li.partial_id)
            rd.inputs.append({"id": f"rd{i}", "role": li.candidate.role, "path": li.candidate.path, "state": state})
            if li.sample:
                rd.samples.append({"key": f"s{len(rd.samples) + 1}", "name": li.sample, "alias": li.sample,
                                   "sample_uid": li.sample, "barcode": None, "project": None})
        if rd.partial_ids:
            rd.modules["len"] = {"x_reads": reads, "x_bases": bases}
        rd.read_level = {"state": "absent" if not linked else "partial" if rd.waiting or not rd.partial_ids else "complete"}
        return rd

    def read_checks(self, bundle) -> list[dict]:
        return []

    def card_columns(self, bundle) -> dict:
        x = bundle.get("x_fake") or {}
        return {"read_data": bundle["tiers"]["read_level"]["state"], "t2_reads": x.get("t2_reads"),
                "t2_samples": len(bundle.get("samples") or []) or None, "t2_files_waiting": x.get("waiting") or 0}

    def sample_rows(self, bundle) -> list[dict]:
        return [{"sample_uid": s["sample_uid"], "name": s["name"], "alias": s["alias"], "barcode": s["barcode"],
                 "project": s["project"], "run": bundle["bundle_id"], "key": s["key"], "reads": None, "bases": None,
                 "n50": None, "q_med": None} for s in bundle.get("samples") or []]

    def install(self, monkeypatch) -> FakeReadData:
        from flowcellar.readdata import api
        for name in ("default_input_rules", "sniff", "candidates_for_run", "link", "read_pass", "read_data", "read_checks",
                     "card_columns", "sample_rows"):
            monkeypatch.setattr(api, name, getattr(self, name))
        return self


class FakeReport:
    """The hub writer and the file server of the viewer, as the real ones behave."""
    CSP = "default-src 'none'; script-src 'sha256-fake' 'self'"

    def __init__(self):
        self.calls: Counter = Counter()
        self.hub_writes: list[list[str] | None] = []

    def write_hub(self, library, *, sidecars: bool = True, reports: bool = False, only=None) -> dict:
        self.calls["write_hub"] += 1
        self.hub_writes.append(None if only is None else list(only))
        hub = library.hub_dir()
        os.makedirs(os.path.join(hub, "data"), exist_ok=True)
        page = os.path.join(hub, "index.html")
        if not os.path.isfile(page):
            atomic.atomic_write_text(page, "<!doctype html><title>fake hub</title>")
        atomic.atomic_write_text(os.path.join(hub, "index.js"),
                                 '__FC_INDEX__("' + base64.b64encode(library.index_gzip()).decode() + '");\n')
        keys = library.run_keys()
        written, removed = [], []
        for k in (keys if only is None else [k for k in only if k in keys]):
            path = os.path.join(hub, "data", k + ".js")
            if only is None and os.path.isfile(path) and os.path.getmtime(path) >= os.path.getmtime(
                    os.path.join(library.path, "bundles", k + ".fcb.json.gz")):
                continue
            if sidecars:
                atomic.atomic_write_text(path, self._sidecar(library, k))
                written.append(k)
        for n in os.listdir(os.path.join(hub, "data")):
            if n.endswith(".js") and n[:-3] not in keys:
                os.unlink(os.path.join(hub, "data", n))
                removed.append(n[:-3])
        return {"index_html": page, "written": written, "removed": removed}

    @staticmethod
    def _sidecar(library, k: str) -> str:
        gz = library.composed_gzip(k)
        return f'__FCB_SIDECAR__("{k}","{base64.b64encode(gz).decode()}","{hash_obj(k)}");\n'

    def serve_file(self, library, url_path: str) -> Served | None:
        self.calls["serve_file"] += 1
        hub = library.hub_dir()
        headers = {"Content-Security-Policy": self.CSP}
        if url_path in ("/", "/index.html"):
            return Served(atomic.read_bytes(os.path.join(hub, "index.html")), "text/html; charset=utf-8", headers)
        if url_path == "/index.js":
            return Served(atomic.read_bytes(os.path.join(hub, "index.js")), "text/javascript; charset=utf-8", headers)
        if url_path.startswith("/data/") and url_path.endswith(".js"):
            k = url_path[len("/data/"):-3]
            if k in library.run_keys():
                return Served(self._sidecar(library, k).encode(), "text/javascript; charset=utf-8", headers)
        return None

    def install(self, monkeypatch) -> FakeReport:
        from flowcellar.report import api
        monkeypatch.setattr(api, "write_hub", self.write_hub)
        monkeypatch.setattr(api, "serve_file", self.serve_file)
        return self
