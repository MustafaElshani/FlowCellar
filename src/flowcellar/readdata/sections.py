"""The read-level sections of one run, from its partials and its small collected files.

Nothing large is read here: a BAM contributes only through the partial its pass wrote. The work is selection and
bookkeeping: which input is the source of record of which need, which sample a fragment belongs to, which reads a
section describes (its *population*), and what is not known.
"""
from __future__ import annotations

import copy
import os

from flowcellar.aggregates import build, derive, merge
from flowcellar.contracts import LinkedInput, ReadData, RunIdentity
from flowcellar.model import axes
from flowcellar.model import partial as partial_model
from flowcellar.readdata import fate, link, record, sections_align, tables
from flowcellar.readdata.bam import model_tier
from flowcellar.readdata.sample_sheet import read_sample_sheet
from flowcellar.util.errors import FlowCellarError
from flowcellar.util.text import clean_text
from flowcellar.util.timeutil import from_epoch

VERSION = "0.1.2"                 # 0.1.1: reads of a file named otherwise go to the sample of the sheet; a source version is named. 0.1.2: samples nothing names any more are dropped
STANDARD_MAX_SAMPLES = 12
MAX_SAMPLES = 384
UNCLASSIFIED = "unclassified"
SHORT_READ_EDGE = float(axes.EDGES["read_length_log40"][52])        # 199.5 b: with whole lengths, "under 200 bases"
COARSE = {"read_length_log40": "read_length_log20", "qscore_0p25": "qscore_1", "identity_q_0p1": "identity_q_0p5"}
_ROLE_ORDER = ("basecaller.sample_sheet", "dorado.bam", "dorado.summary", "align.genome_bam", "align.transcriptome_bam",
               "pychopper.stats", "quant.salmon", "quant.oarfish", "quant.kallisto")
_FAIL_PARTS = ("bam_fail", "fastq_fail", "fail")


# ------------------------------------------------------------------------------------------ samples
class Samples:
    """The sample table of a bundle while it is built."""

    def __init__(self, tier1: dict | None, config: dict | None):
        self.config = config or {}
        self.rows: list[dict] = []
        self.unclassified: dict | None = None
        for s in (tier1 or {}).get("samples") or []:
            if s.get("key") == UNCLASSIFIED:
                self.unclassified = dict(s)
            else:
                self.rows.append(dict(s))
        self.warnings: list[dict] = []
        #: file sample name -> key of the sheet's sample its reads belong to on this run (``learn_file``)
        self.file_alias: dict[str, str] = {}
        #: keys this assembly gave to a file, a fragment or a row of a sample sheet. The rows a bundle already holds
        #: are the starting point (so that a sample keeps its key); a row that nothing names any more, and that no
        #: sample sheet lists, is left out by ``finish`` instead of staying in the bundle for ever.
        self.touched: set[str] = set()

    def in_sheet(self, label: str | None) -> dict | None:
        """The sample of this run's sample sheet that a name or a read label means: its alias (after the library's
        alias map) or its barcode."""
        name = link.canonical(label, self.config)
        if not name:
            return None
        for s in self.rows:
            if s.get("expected") and (name in (s.get("name"), s.get("alias"), link.canonical(s.get("alias"), self.config)) or name == s.get("barcode")):
                return s
        return None

    def learn_file(self, file_sample: str | None, labels) -> str | None:
        """A file named after a sample the sheet of this run does not list, whose reads carry exactly one label the
        sheet does list: the file is that sample on this run. This is a file that holds one sample from two flow
        cells, each under the alias of its own sample sheet (the identifiers inside a file beat its name).

        One file name per sample of the sheet: a second file name whose reads carry the same label keeps its own
        name, so that two candidate files are never added up under one sample."""
        if not file_sample or file_sample == record.ALL or link.is_unclassified(file_sample, self.config):
            return None
        if file_sample in self.file_alias:
            return self.file_alias[file_sample]
        if self.in_sheet(file_sample) is not None:
            return None
        hits = {s["key"] for s in (self.in_sheet(x) for x in labels or []) if s is not None}
        if len(hits) != 1:
            return None
        key = next(iter(hits))
        if key in self.file_alias.values():
            return None
        self.file_alias[file_sample] = key
        self.touched.add(key)
        s = next(x for x in self.rows if x["key"] == key)
        self.warnings.append({"level": "info", "code": "file_sample_is_sheet_alias",
                              "message": f"The files named '{file_sample}' hold the reads of this flow cell under the name '{s.get('alias') or s.get('name')}' "
                                         "of its sample sheet; they are shown as that sample."})
        return key

    def _find(self, *, name=None, barcode=None) -> dict | None:
        for s in self.rows:
            if name and name in (s.get("name"), s.get("alias"), link.canonical(s.get("alias"), self.config)):
                return s
        if barcode:
            for s in self.rows:
                if s.get("barcode") == barcode:
                    return s
        return None

    def _new(self, **kw) -> dict:
        used = {s["key"] for s in self.rows}
        n = len(self.rows) + 1
        while f"s{n}" in used:
            n += 1
        s = {"key": f"s{n}", "barcode": None, "alias": None, "name": None, "sample_uid": None, "type": None, "expected": False,
             "source": "none", "project": None}
        s.update(kw)
        self.rows.append(s)
        return s

    def from_sheet(self, row: dict, project: str | None) -> None:
        alias = row.get("alias") or row.get("sample_id")
        name = link.canonical(alias, self.config) or row.get("barcode")
        s = self._find(name=name, barcode=row.get("barcode"))
        if s is None:
            s = self._new()
        elif s.get("alias") and alias and link.canonical(s["alias"], self.config) != name and s.get("source") == "basecaller_sample_sheet":
            self.warnings.append({"level": "warning", "code": "sample_sheets_disagree",
                                  "message": f"Two sample sheets name barcode {row.get('barcode')} differently ({s['alias']} and {alias}); the first one is used."})
            return
        s.update(barcode=row.get("barcode") or s.get("barcode"), alias=alias, name=name, expected=True, source="basecaller_sample_sheet",
                 type=row.get("type") or s.get("type") or "test_sample")
        self.touched.add(s["key"])
        if project and not s.get("project"):
            s["project"] = project

    def key_for(self, label: str | None, file_sample: str | None, *, barcoded: bool, source: str, project: str | None) -> str:
        """The sample key of a fragment: the sample of its file when the file is one sample's, else the label its
        records carry, else "no barcode found" on a barcoded run and the run's only sample otherwise."""
        name = None
        if file_sample and file_sample != record.ALL:
            if link.is_unclassified(file_sample, self.config):
                return self._unclassified()
            # the file's name is not in the sample sheet but the label on its reads is: the sheet's sample
            key = self.learn_file(file_sample, [label] if label else [])
            if key is not None:
                s = next(x for x in self.rows if x["key"] == key)
                if project and not s.get("project"):
                    s["project"] = project
                return key
            name = file_sample
        elif label:
            name = link.canonical(label, self.config)
            if link.is_unclassified(name, self.config):
                return self._unclassified()
        if name is None:
            if barcoded or len(self.rows) > 1:
                return self._unclassified()
            if not self.rows:
                self._new(name="sample", source="none")
            self.touched.add(self.rows[0]["key"])
            return self.rows[0]["key"]
        s = self._find(name=name, barcode=name if name.startswith("barcode") else None)
        if s is None:
            s = self._new(name=name, alias=label or name, source=source, barcode=name if name.startswith("barcode") else None)
            if any(x.get("expected") for x in self.rows):
                self.warnings.append({"level": "warning", "code": "sample_not_in_sheet",
                                      "message": f"The sample name '{name}' is in the data but not in the sample sheet; it is shown as its own sample."})
        if project and not s.get("project"):
            s["project"] = project
        self.touched.add(s["key"])
        return s["key"]

    def _unclassified(self) -> str:
        if self.unclassified is None:
            self.unclassified = {"key": UNCLASSIFIED, "barcode": None, "alias": UNCLASSIFIED, "name": "No barcode found", "sample_uid": None,
                                 "type": None, "expected": False, "source": "none", "project": None}
        return UNCLASSIFIED

    def finish(self) -> list[dict]:
        # a row carried over from the bundle as it was, that no sample sheet lists and no file or read names today
        self.rows = [s for s in self.rows if s["key"] in self.touched or s.get("source") == "basecaller_sample_sheet"]
        for s in self.rows:
            s.setdefault("name", s.get("alias"))
            if not s.get("sample_uid") and s.get("name"):
                s["sample_uid"] = f"{s['project']}/{s['name']}" if s.get("project") else s["name"]
            for k in ("alias", "name", "sample_uid", "barcode", "project", "type"):
                if s.get(k) is not None:
                    s[k] = clean_text(s[k], 256)
        return self.rows + ([self.unclassified] if self.unclassified else [])

    def order(self) -> list[str]:
        return [s["key"] for s in self.rows] + [UNCLASSIFIED]

    def name(self, key: str) -> str:
        if key == UNCLASSIFIED:
            return "no barcode found"
        return next((s.get("name") or key for s in self.rows if s["key"] == key), key)

    def expected_keys(self) -> list[str]:
        exp = [s["key"] for s in self.rows if s.get("expected")]
        return exp or [s["key"] for s in self.rows]


def project_of(li: LinkedInput, config: dict | None) -> str | None:
    """The project a file lies in: the path component named by its root's ``project_depth``."""
    c = li.candidate
    if not c.root_id or not c.relpath:
        return None
    root = next((r for r in (config or {}).get("roots") or [] if r.get("id") == c.root_id), None)
    if not root:
        return None
    depth = int(root.get("project_depth", 1) or 0)
    parts = [p for p in c.relpath.split("/") if p]
    name = parts[depth - 1] if 0 < depth < len(parts) else None
    # a folder name that an [[inputs]] pattern spells out (a pipeline's "analysis") is not a project
    from flowcellar.config.loader import not_project_names
    return None if name is not None and name.lower() in not_project_names(config) else name


# -------------------------------------------------------------------------------------------- cells
def _coarsen_hist(h: dict, axis: str) -> dict:
    if h["axis"] == axis:
        return h
    f = axes.coarsen_factor(h["axis"], axis)
    n = axes.coarsen(build.hist_full(h, "n"), f)
    w = axes.coarsen(build.hist_full(h, "w"), f) if "w" in h else None
    return build.hist(axis, n, w, w_name=h.get("w_name"), missing=h.get("missing"), vmin=h.get("min"), vmax=h.get("max"))


def _filter_of(li: LinkedInput) -> str:
    """pass or fail, from where the basecaller put the file (``bam_fail/``) or how it named it (``_fail_``)."""
    parts = li.candidate.path.split(os.sep)
    if any(p in _FAIL_PARTS for p in parts[:-1]):
        return "fail"
    from flowcellar.model.runkey import parse_minknow_name
    p = parse_minknow_name(parts[-1])
    return "fail" if p and p.status == "fail" else "pass"


def _mtime(ns: int | None) -> str | None:
    return from_epoch(ns / 1e9) if ns else None


def _weak(li: LinkedInput) -> bool:
    return li.link.get("confidence") in ("medium", "low") and not li.link.get("confirmed")


class _Cube:
    """Cells of one section keyed by (sample key, filter), merged as fragments arrive."""

    def __init__(self):
        self.cells: dict[tuple, dict] = {}

    def add(self, key: tuple, cell: dict) -> None:
        cur = self.cells.get(key)
        self.cells[key] = copy.deepcopy(cell) if cur is None else merge.merge_cell(cur, cell)

    def listed(self, order: list[str]) -> list[dict]:
        pos = {k: i for i, k in enumerate(order)}
        out = []
        for (s, f) in sorted(self.cells, key=lambda k: (pos.get(k[0], len(pos)), k[1] != "pass", k)):
            out.append({"at": [s, f], **self.cells[(s, f)]})
        return out


def _total(cells: list[dict], *, pass_only: bool) -> dict | None:
    where = {"filter": ["pass"]} if pass_only and any(c["at"][1] == "pass" for c in cells) else None
    rolled = merge.rollup(cells, ["sample", "filter"], [], where=where)
    return rolled[0] if rolled else None


def _num(x, digits: int = 0):
    if x is None:
        return None
    return int(round(x)) if digits == 0 else round(float(x), digits)


def len_summary(cell: dict | None) -> dict:
    h = ((cell or {}).get("hists") or {}).get("len_hist")
    if not h or not derive.total(h):
        return {"len_n50": None, "len_median": None, "len_mean": None, "len_short_read_frac": None}
    short = derive.frac_above(h, SHORT_READ_EDGE)
    return {"len_n50": _num(derive.nx(h, 0.5)), "len_median": _num(derive.median(h)), "len_mean": _num(derive.mean(h), 1),
            "len_short_read_frac": None if short is None else round(1.0 - short, 6)}


def qual_median(cell: dict | None) -> float | None:
    h = ((cell or {}).get("hists") or {}).get("qual_hist")
    return _num(derive.median(h), 2) if h and derive.total(h) else None


def aln_summary(cell: dict | None) -> dict:
    sc = (cell or {}).get("scalars") or {}
    h = ((cell or {}).get("hists") or {}).get("aln_identity_hist")
    out = {"aln_mapped_frac": None, "aln_identity_median": None, "aln_supplementary_frac": None}
    if sc.get("reads"):
        out["aln_mapped_frac"] = round(sc.get("mapped_reads", 0) / sc["reads"], 6)
    if sc.get("mapped_reads"):
        out["aln_supplementary_frac"] = round(sc.get("reads_with_supplementary", 0) / sc["mapped_reads"], 6)
    if h and derive.total(h):
        q = derive.median(h)
        out["aln_identity_median"] = round(1.0 - 10.0 ** (-q / 10.0), 6)
    return out


# --------------------------------------------------------------------------------------- populations
def population(role: str, facts_list: list[dict], *, barcoded: bool, has_unclassified: bool, has_fail: bool, labelled: bool,
               unmapped_reads: int) -> dict:
    """Which reads a source holds, from its role, its @PG lines and its read names."""
    ev: list[str] = []
    dor: dict = {}
    pych_names = names = 0
    pych_prog = full_length = False
    for f in facts_list:
        d = f.get("dorado") or {}
        for k, v in d.items():
            dor.setdefault(k, v)
        pych_names += f.get("pychopper_names", 0)
        names += f.get("names_checked", 0)
        pych_prog = pych_prog or bool(f.get("pychopper_program"))
        full_length = full_length or bool((f.get("aligner") or {}).get("input_full_length"))
    min_q = dor.get("min_qscore")
    if min_q is not None:
        ev.append(f"@PG dorado command line holds --min-qscore {min_q:g}")
    both = True if dor.get("barcode_both_ends") else None
    if both:
        ev.append("@PG dorado command line holds --barcode-both-ends")
    trimmed = None
    if "no_trim" in dor:
        trimmed = not dor["no_trim"]
        ev.append("@PG dorado command line holds --no-trim" if dor["no_trim"] else "@PG dorado basecaller ran without --no-trim")
    pych = pych_prog or full_length or (names > 0 and pych_names * 10 >= names * 9)
    if role in ("dorado.bam", "dorado.summary"):
        if barcoded and labelled and not has_unclassified:
            stage = "classified"
            ev.append("the files hold the reads given to samples; no file of reads without a barcode is linked")
        elif has_fail or min_q is None:
            stage = "basecalled_all"
            ev.append("passed and failed reads are both present" if has_fail else "no quality filter is named in the header")
        else:
            stage = "basecalled_pass"
    elif pych:
        stage, trimmed = "full_length_trimmed", True
        if names and pych_names:
            ev.append(f"{pych_names:,} of the first {names:,} read names carry the segment coordinates a cDNA pre-processor writes")
        if full_length:
            ev.append("the aligner's input is a full-length FASTQ (@PG)")
        if pych_prog:
            ev.append("a cDNA pre-processor is named in the header (@PG)")
    elif unmapped_reads == 0:
        stage = "aligned_primary"
        ev.append("the file holds no unmapped read")
    elif barcoded and labelled:
        stage = "classified"
        ev.append("the file also holds unmapped reads, and its reads carry a sample name")
    else:
        stage = "basecalled_pass" if min_q is not None else "basecalled_all"
        ev.append("the file also holds unmapped reads: it is taken to hold every read given to the aligner")
    ev.append(f"input role {role}")
    return {"stage": stage, "min_qscore": min_q, "trimmed": trimmed, "both_end_barcodes": both, "source_role": role, "evidence": ev}


# ------------------------------------------------------------------------------------------ building
class _Input:
    """One linked input of this run while the bundle is built."""

    def __init__(self, n: int, li: LinkedInput, run_key: str):
        self.id = f"rd{n}"
        self.li = li
        self.role = li.candidate.role
        self.counted = list(li.counted_for.get(run_key, []))
        self.partial: dict | None = None
        self.waiting = False
        self.records = None
        self.samples: set[str] = set()
        self.used_by: set[str] = set()
        self.notes: list[dict] = []
        self.state = "ok"
        self.coverage: dict | None = None
        self.read_error: str | None = None


def _fragments(p: dict, run: RunIdentity, sole_run: bool) -> list[dict]:
    """The fragments of a partial that belong to this run: routed to its key, routed to its run id when the file
    gives no flow cell, and the unrouted ones when the file is linked to this run only."""
    short = "?_" + run.run_short[:8]
    return [f for f in p.get("fragments", []) if f.get("run_key") == run.run_key or f.get("run_key") == short
            or (f.get("run_key") is None and sole_run)]


def build_read_data(run: RunIdentity, linked, partials_dir: str | None, *, tier1: dict | None = None, config: dict | None = None) -> ReadData:
    config = config or {}
    K = run.run_key
    mine = [li for li in linked if K in li.run_keys]
    if not mine:
        return ReadData()
    mine.sort(key=lambda li: (_ROLE_ORDER.index(li.candidate.role) if li.candidate.role in _ROLE_ORDER else 99, li.candidate.path))
    ins = [_Input(i + 1, li, K) for i, li in enumerate(mine)]
    rd = ReadData()
    warnings: list[dict] = []
    samples = Samples(tier1, config)
    assay = ((tier1 or {}).get("identity") or {}).get("assay") or {}

    # ---- sample sheets first: they say which samples to expect
    for x in ins:
        if x.role != "basecaller.sample_sheet":
            continue
        if _weak(x.li):
            x.notes.append({"level": "info", "code": "weak_link", "message": "Linked by a rule, not by what the file contains: its sample names are shown but samples are not expected from it."})
            continue
        try:
            sheet = read_sample_sheet(x.li.candidate.path)
        except (OSError, FlowCellarError) as e:
            x.state, x.read_error = "error", str(e)
            continue
        rows = [r for r in sheet["rows"] if not r.get("flow_cell_id") or r["flow_cell_id"] == run.flow_cell_id]
        proj = project_of(x.li, config)
        for r in rows:
            samples.from_sheet(r, proj)
        x.records = len(rows)
        x.used_by.update(["funnel.samples"])
    barcoded = bool(assay.get("barcoded")) or any(s.get("barcode") for s in samples.rows) or len(samples.rows) > 1
    # a file named after a sample this sheet does not list, whose reads carry one label it does list (the labels the
    # link probe saw): the files of that name are the sheet's sample here, before and after they are read
    for x in ins:
        if x.role in link.BAM_ROLES and not _weak(x.li):
            samples.learn_file(x.li.sample, x.li.link.get("labels") or [])
    # every file that is one sample's names that sample, counted or not: a sample in conflict is listed, without data
    for x in ins:
        if x.li.sample and x.li.sample != record.ALL and x.role != "basecaller.sample_sheet" and not _weak(x.li):
            x.samples.add(samples.key_for(None, x.li.sample, barcoded=barcoded, source="file_name", project=project_of(x.li, config)))
    barcoded = barcoded or len(samples.rows) > 1

    # ---- partials of the inputs that are the source of record of something
    for x in ins:
        if x.role not in link.PASS_ROLES or not x.counted:
            continue
        path = partial_model.find_partial(partials_dir, x.li.partial_id) if (partials_dir and x.li.partial_id) else None
        if not path:
            x.waiting = True
            x.state = "skipped"
            x.notes.append({"level": "info", "code": "waiting", "message": "Not read yet: the pass over this file has not run."})
            rd.waiting.append(x.li.candidate.path)
            continue
        try:
            x.partial = partial_model.read_partial(path)
        except (OSError, ValueError, FlowCellarError) as e:
            x.waiting = True
            x.state = "error"
            x.read_error = f"its cached result cannot be read ({e})"
            rd.waiting.append(x.li.candidate.path)
            continue
        rd.partial_ids.append(x.partial["partial_id"])
        st = x.partial.get("state")
        x.state = {"ok": "ok", "truncated": "truncated", "partially_read": "partially_read"}.get(st, "error")
        x.coverage = x.partial.get("coverage")
        for e in x.partial.get("errors") or []:
            x.notes.append({"level": "warning", "code": clean_text(e.get("code") or "note", 60), "message": clean_text(e.get("message") or "", 600)})

    # ---- reads and genome: cells per sample and filter
    cubes = {"len": _Cube(), "qual": _Cube(), "aln": _Cube()}
    col = sections_align.Collector(cubes, _Cube)       # over time, poly(A), depth, the alignment to transcripts
    src = {"reads": [], "genome": [], "transcriptome": []}     # the inputs behind each need, for populations and notes
    unrouted_notes: list[str] = []
    labelled = {"reads": False, "genome": False}
    for x in ins:
        if x.partial is None:
            continue
        p = x.partial
        sole = len(x.li.run_keys) == 1
        frags = _fragments(p, run, sole)
        proj = project_of(x.li, config)
        frag_filter = _filter_of(x.li)
        reads_here = 0
        for f in frags:
            mod = f.get("module")
            need = sections_align.need_of(mod, x.counted)
            if need not in x.counted:
                continue
            if f.get("sample"):
                labelled[need] = True
            skey = samples.key_for(f.get("sample"), x.li.sample, barcoded=barcoded, project=proj,
                                   source="file_name" if x.li.sample and x.li.sample != record.ALL else "bam_read_group")
            x.samples.add(skey)
            for cell in f.get("cells") or []:
                flt = cell["at"][0] if f.get("dims") == ["filter"] and cell.get("at") else frag_filter
                body = {k: v for k, v in cell.items() if k != "at"}
                cubes[col.route(mod, need)].add((skey, flt), col.take(f, need, skey, flt, body))
                if mod == "len":
                    reads_here += (body.get("scalars") or {}).get("reads", 0)
                elif mod == "aln" and "reads" not in x.counted:
                    reads_here += (body.get("scalars") or {}).get("reads", 0)
        x.records = reads_here
        for need in ("reads", "genome", "transcriptome"):
            if need in x.counted:
                src[need].append(x)
        routing = p.get("routing") or {}
        if not sole and routing.get("unrouted"):
            unrouted_notes.append(f"{os.path.basename(x.li.candidate.path)}: {routing['unrouted']:,} reads name no run and were given to none")
        others = {k: v for k, v in (routing.get("runs") or {}).items() if k not in (K, "?_" + run.run_short[:8])}
        if others:
            x.notes.append({"level": "info", "code": "shared_file_routed", "value": int(reads_here),
                            "message": f"This file also holds reads of other runs ({', '.join(f'{k} {v:,}' for k, v in sorted(others.items())[:4])}); "
                                       f"{reads_here:,} reads were routed to this run by the identifiers on each read."})

    order = samples.order()
    _thin_stray_samples(cubes, samples, warnings)
    n_named = len({k[0] for k in cubes["len"].cells if k[0] != UNCLASSIFIED} | {k[0] for k in cubes["aln"].cells if k[0] != UNCLASSIFIED})
    detail = "standard" if n_named <= STANDARD_MAX_SAMPLES else "coarse"
    conflicts = _conflicts(ins, K)

    def facts(xs):
        return [x.partial.get("facts") or {} for x in xs if x.partial]

    def state_of(xs, need: str, cells: list[dict]) -> tuple[str, str | None, dict]:
        """(state, partial_reason, coverage) of a section from the inputs behind it."""
        waiting = [x for x in ins if x.waiting and need in x.counted]
        part = [x for x in xs if x.state in ("partially_read", "truncated", "error")]
        cov = {"kind": "full", "records": int(sum((c.get("scalars") or {}).get("reads", 0) for c in cells)), "records_total": None, "fraction": 1.0}
        if waiting and not part:
            # files of this section have not been read: the share is the share of their bytes that was, never "full"
            every = [x for x in ins if need in x.counted and x.role in link.PASS_ROLES]
            total = sum(int(x.li.candidate.size or 0) for x in every)
            read = sum(int(x.li.candidate.size or 0) for x in every if not x.waiting)
            cov.update(kind="byte_ranges", fraction=round(read / total, 6) if total else None)
        if part:
            fr = [x.coverage.get("fraction") for x in part if x.coverage and x.coverage.get("fraction") is not None]
            cov.update(kind="truncated" if any(x.state == "truncated" for x in part) else "byte_ranges", fraction=min(fr) if fr else None)
            return "partial", "partly_read", cov
        if waiting:
            return "partial", "waiting", cov
        have = {c["at"][0] for c in cells}
        if conflicts.get(need) or any(k not in have for k in samples.expected_keys()):
            return "partial", "samples_missing", cov
        return "ok", None, cov

    def section(mod: str, need: str, cells: list[dict], pop: dict | None, summary: dict, params: dict, notes: list[dict]) -> dict:
        xs = src[need]
        state, reason, cov = state_of(xs, need, cells)
        if not cells:
            state, reason = ("partial", reason or "waiting") if any(x.waiting and need in x.counted for x in ins) else ("skipped", None)
        # a source that is not the primary version of its file (a second basecalling of the same reads) is named
        primary_variant = (config.get("variants") or {}).get("primary")
        other = sorted({x.li.variant for x in xs if x.li.variant and x.li.variant != primary_variant})
        if other:
            notes = notes + [{"level": "info", "code": "source_variant", "value": len(other),
                              "message": clean_text(f"Read from the files of version '{', '.join(other)}'"
                                                    + (f", not from the primary version '{primary_variant}'" if primary_variant else "") + ".", 300)}]
        sec = {"version": VERSION, "tier": "read_level", "state": state, "inputs": [x.id for x in xs] + [x.id for x in ins if x.waiting and need in x.counted and x not in xs],
               "dims": ["sample", "filter"], "cells": cells, "summary": summary, "params": params, "notes": notes, "detail": detail,
               "population": pop, "coverage": cov}
        if state == "partial":
            sec["partial_reason"] = reason
        for x in xs:
            x.used_by.add(mod)
        return sec

    def conflict_notes(need: str) -> list[dict]:
        return [{"level": "warning", "code": "two_candidate_sources", "message": clean_text("Not assessed: " + m, 1200)} for m in conflicts.get(need, [])]

    def weak_notes(xs) -> list[dict]:
        weak = [x for x in xs if _weak(x.li)]
        if not weak:
            return []
        return [{"level": "warning", "code": "weak_link",
                 "message": f"{len(weak)} of the {len(xs)} files behind this section are linked by a rule, not by what they contain "
                            f"({os.path.basename(weak[0].li.candidate.path)}{', ...' if len(weak) > 1 else ''}). The numbers are shown but no check is assessed on them."}]

    # ---- yield.basecalled, len, qual
    reads_pop = None
    if src["reads"] or conflicts.get("reads") or any(x.waiting and "reads" in x.counted for x in ins):
        len_cells_all = cubes["len"].listed(order)
        roles = sorted({x.role for x in src["reads"]}, key=lambda r: record.READS_PREFERENCE.index(r))
        role = max(roles, key=lambda r: sum(1 for x in src["reads"] if x.role == r)) if roles else None
        tot_all = _total(len_cells_all, pass_only=False)
        has_fail = any(c["at"][1] == "fail" for c in len_cells_all)
        has_uncl = any(c["at"][0] == UNCLASSIFIED and (c.get("scalars") or {}).get("reads") for c in len_cells_all)
        unmapped = sum((c.get("scalars") or {}).get("unmapped_reads", 0) for c in cubes["aln"].listed(order)) if role and role.startswith("align.") else 0
        if role:
            reads_pop = population(role, facts(src["reads"]), barcoded=barcoded, has_unclassified=has_uncl, has_fail=has_fail,
                                   labelled=labelled["reads"] or any(x.li.sample and x.li.sample != record.ALL for x in src["reads"]),
                                   unmapped_reads=unmapped)
        notes = conflict_notes("reads") + weak_notes(src["reads"])
        if len(roles) > 1:
            notes.append({"level": "warning", "code": "mixed_sources",
                          "message": "The samples of this run were read from different kinds of file (" + ", ".join(roles) + "), so their numbers describe different sets of reads; no check is assessed on them together."})
        for t in unrouted_notes:
            notes.append({"level": "info", "code": "unrouted_reads", "message": t})
        params_common = {"summary_filter": "pass" if has_fail else "all", "quantile_rule": "read from the stored histogram by rank inside a bin",
                         "mixed_sources": len(roles) > 1}
        # len
        len_cells, qual_cells, yb_cells = [], [], []
        qual_all = {tuple(c["at"]): c for c in cubes["qual"].listed(order)}
        grid_total = None
        for c in len_cells_all:
            sc = c.get("scalars") or {}
            h = (c.get("hists") or {}).get("len_hist")
            lc = {"at": c["at"], "scalars": {"reads": sc.get("reads", 0), "bases": sc.get("bases", 0)}}
            if h:
                lc["hists"] = {"len_hist": h if detail == "standard" else _coarsen_hist(h, COARSE[h["axis"]])}
            len_cells.append(lc)
            yb_cells.append({"at": c["at"], "scalars": {"reads": sc.get("reads", 0), "bases": sc.get("bases", 0),
                                                        "split_children": sc.get("split_children", 0), "split_parents": sc.get("split_parents", 0)}})
            q = qual_all.get(tuple(c["at"]))
            if q:
                qh = (q.get("hists") or {}).get("qual_hist")
                qc = {"at": c["at"], "scalars": {"reads": sc.get("reads", 0), "quality_from_base_qualities": (q.get("scalars") or {}).get("quality_from_base_qualities", 0)}}
                if qh:
                    qc["hists"] = {"qual_hist": qh if detail == "standard" else _coarsen_hist(qh, COARSE[qh["axis"]])}
                g = (q.get("hists2d") or {}).get("qual_vs_length")
                if g:
                    if detail == "standard":
                        qc["hists2d"] = {"qual_vs_length": g}
                    grid_total = g if grid_total is None else merge.merge_hist2d(grid_total, g)
                qual_cells.append(qc)
        tot = _total(len_cells_all, pass_only=True)
        ls = len_summary(tot)
        sec_len = section("len", "reads", len_cells, reads_pop, ls, dict(params_common, short_read_below_bases=200), list(notes))
        rd.modules["len"] = sec_len
        qtot = _total(list(qual_all.values()), pass_only=True)
        sec_q = section("qual", "reads", qual_cells, reads_pop, {"qual_median": qual_median(qtot)},
                        dict(params_common, quality="mean quality of a read as the basecaller reports it (qs tag); from the base qualities when a read has no qs tag"),
                        list(notes))
        if detail == "coarse" and grid_total is not None:
            sec_q["run_level"] = {"hists2d": {"qual_vs_length": grid_total}}
        rd.modules["qual"] = sec_q
        if reads_pop and reads_pop["stage"] in ("basecalled_all", "basecalled_pass", "classified"):
            ysum = {"yield_basecalled_reads": (tot_all.get("scalars") or {}).get("reads") if tot_all else None,
                    "yield_basecalled_bases": (tot_all.get("scalars") or {}).get("bases") if tot_all else None,
                    "yield_pass_frac": None}
            if reads_pop["stage"] == "basecalled_all" and has_fail and tot_all and tot:
                ysum["yield_pass_frac"] = round((tot.get("scalars") or {}).get("reads", 0) / max((tot_all.get("scalars") or {}).get("reads", 0), 1), 6)
            rd.modules["yield.basecalled"] = section("yield.basecalled", "reads", yb_cells, reads_pop, ysum, dict(params_common), list(notes))

    # ---- aln
    aln_pop = None
    if src["genome"] or conflicts.get("genome") or any(x.waiting and "genome" in x.counted for x in ins):
        raw = cubes["aln"].listed(order)
        aln_cells = []
        for c in raw:
            hs = c.get("hists") or {}
            cell = {"at": c["at"], "scalars": dict(c.get("scalars") or {})}
            if hs:
                cell["hists"] = {name: (h if detail == "standard" or h["axis"] not in COARSE else _coarsen_hist(h, COARSE[h["axis"]]))
                                 for name, h in hs.items()}
            aln_cells.append(cell)
        tot = _total(raw, pass_only=True)
        if src["genome"]:
            unmapped = (tot.get("scalars") or {}).get("unmapped_reads", 0) if tot else 0
            aln_pop = population("align.genome_bam", facts(src["genome"]), barcoded=barcoded, has_unclassified=False, has_fail=False,
                                 labelled=labelled["genome"] or any(x.li.sample and x.li.sample != record.ALL for x in src["genome"]),
                                 unmapped_reads=unmapped)
        notes = conflict_notes("genome") + weak_notes(src["genome"])
        rd.modules["aln"] = section("aln", "genome", aln_cells, aln_pop, aln_summary(tot),
                                    {"identity_rule": "1 - NM / (M + I + D) from the NM tag and the CIGAR, stored as -10 log10(1 - identity); identity 1 is the overflow cell",
                                     "reads_counted": "primary records; secondary and supplementary records are counted, not read as reads"}, notes)

    # ---- what the same pass gives beside these: strand, CIGAR, grids, depth, over time, poly(A), transcripts
    sections_align.add_sections(col, rd, section=section, samples=samples, order=order, detail=detail, tier1=tier1, src=src, ins=ins,
                                conflicts=conflicts, reads_pop=reads_pop, aln_pop=aln_pop, barcoded=barcoded, labelled=labelled)

    # ---- cdna and the collected numbers of read fate
    collected = _collect(ins, samples, config, barcoded)
    if collected["cdna_cells"] or conflicts.get("pychopper.stats"):
        cd = collected["cdna_cells"]
        pos = {k: i for i, k in enumerate(order)}
        cells = [{"at": [k], **cd[k]} for k in sorted(cd, key=lambda k: pos.get(k, 1 << 30))]
        num = sum(c["scalars"]["primers_found"] + c["scalars"]["rescued_reads"] for c in cells)
        den = sum(c["scalars"]["pass_reads"] for c in cells)
        missing = [k for k in samples.expected_keys() if k not in cd]
        sec = {"version": VERSION, "tier": "assay", "state": "ok" if cells and not missing and not conflicts.get("pychopper.stats") else "partial" if cells else "skipped",
               "inputs": [x.id for x in ins if x.role == "pychopper.stats" and "pychopper.stats" in x.counted], "dims": ["sample"], "cells": cells,
               "summary": {"cdna_full_length_frac": round(num / den, 6) if den else None, "cdna_samples": len(cells),
                           "cdna_rescued_frac": round(sum(c["scalars"]["rescued_reads"] for c in cells) / den, 6) if den else None},
               "params": {"full_length_rule": "(Primers_found + rescued reads) / PassReads",
                          "written_rule": "reads written = Primers_found + Rescue - LenFail; a fused read that was cut apart counts once per piece"},
               "notes": [{"level": "info", "code": "source_is_report", "message": "Collected from each sample's pychopper statistics file; exact."}]
               + collected["cdna_notes"] + [{"level": "warning", "code": "two_candidate_sources", "message": clean_text("Not assessed: " + m, 1200)} for m in conflicts.get("pychopper.stats", [])]}
        if sec["state"] == "partial":
            sec["partial_reason"] = "samples_missing"
        rd.modules["cdna"] = sec

    # ---- library preparation tables: every category of the pychopper tables, and the annotated read tables
    from flowcellar.readdata import libprep
    libprep.add_sections(rd, ins, samples, config=config, barcoded=barcoded, conflicts=conflicts)

    # ---- read fate
    weak_ids = {x.id for x in ins if _weak(x.li)}

    def per_sample(need: str) -> dict:
        """Which samples have their files of this need read whole, and which are in a two-candidate conflict."""
        read, unread, conflict = set(), set(), set()
        for x in ins:
            if need in ((x.li.link.get("conflicts") or {}).get(K) or {}):
                conflict.update(x.samples or {None})
            if need in x.counted and x.role in link.PASS_ROLES:
                (read if (not x.waiting and x.state == "ok") else unread).update(x.samples)
        return {"read": read - unread, "conflict": conflict}
    funnel = fate.build(run=run, tier1=tier1, samples=samples, len_cells=cubes["len"].listed(order) if not any(x.id in weak_ids for x in src["reads"]) else [],
                        aln_cells=cubes["aln"].listed(order) if not any(x.id in weak_ids for x in src["genome"]) else [],
                        reads_pop=reads_pop, aln_pop=aln_pop, collected=collected, assay_group=assay.get("group"), barcoded=barcoded,
                        reads_inputs=[x.id for x in src["reads"]], genome_inputs=[x.id for x in src["genome"]],
                        reads_complete=not any((x.waiting or x.state != "ok") and "reads" in x.counted for x in ins) and not conflicts.get("reads"),
                        genome_complete=not any((x.waiting or x.state != "ok") and "genome" in x.counted for x in ins) and not conflicts.get("genome"),
                        conflicts=conflicts, reads_roles=sorted({x.role for x in ins if "reads" in x.counted}),
                        reads_waiting=any(x.waiting and "reads" in x.counted for x in ins),
                        genome_waiting=any(x.waiting and "genome" in x.counted for x in ins),
                        sample_state={"reads": per_sample("reads"), "genome": per_sample("genome")})
    rd.modules.update(funnel["sections"])
    warnings.extend(funnel["warnings"])
    for x in ins:
        if x.id in funnel["used"]:
            x.used_by.update(funnel["used"][x.id])

    # ---- the tier, the inputs, basecalling
    rd.samples = samples.finish()
    warnings.extend(samples.warnings)
    for msgs in conflicts.values():
        for m in msgs:
            warnings.append({"level": "warning", "code": "two_candidate_sources", "message": clean_text("Not assessed: " + m, 1200)})
    for x in ins:
        if x.read_error:
            warnings.append({"level": "warning", "code": "input_not_read", "input": x.id, "message": clean_text(f"{os.path.basename(x.li.candidate.path)}: {x.read_error}", 600)})
    seq = fate.sequenced_reads(tier1)
    reads_cells = (rd.modules.get("len") or {}).get("cells") or []
    covered = sum((c.get("scalars") or {}).get("reads", 0) for c in reads_cells)
    with_data = {c["at"][0] for sec in ("len", "aln", "cdna", "funnel.samples") for c in (rd.modules.get(sec) or {}).get("cells", [])} - {UNCLASSIFIED}
    expected = samples.expected_keys()
    any_section = any((rd.modules.get(m) or {}).get("cells") for m in ("len", "aln", "cdna", "funnel.samples"))
    complete = bool(any_section) and not rd.waiting and all(k in with_data for k in expected) and not any(conflicts.values()) \
        and all((rd.modules.get(m) or {}).get("state", "ok") == "ok" for m in ("len", "qual", "aln", "cdna"))
    live = any(x.li.link.get("method") == "run_dir" and x.role in ("dorado.bam", "dorado.summary") for x in src["reads"])
    missing_roles = []
    if not src["reads"] and not any("reads" in x.counted for x in ins):
        missing_roles.append("dorado.bam")
    if not src["genome"] and not any("genome" in x.counted for x in ins):
        missing_roles.append("align.genome_bam")
    pop = reads_pop or aln_pop
    cov_parts = [(rd.modules.get(m) or {}).get("coverage") for m in ("len", "aln") if rd.modules.get(m)]
    rd.read_level = {"state": "complete" if complete else "partial" if (any_section or rd.waiting or rd.modules) else "absent",
                     "origin": "minknow_live_basecalling" if live else "offline_basecalling",
                     "samples_expected": len(expected), "samples_with_data": len(with_data),
                     "reads_expected": seq, "reads_covered": int(covered) if reads_cells else None,
                     "missing_roles": missing_roles, "population": pop,
                     "coverage": {"kind": "full" if all((c or {}).get("kind") == "full" for c in cov_parts) else next(c["kind"] for c in cov_parts if c.get("kind") != "full"),
                                  "records": int(covered) if reads_cells else None, "records_total": seq,
                                  "fraction": (round(min(covered / seq, 1.0), 6) if (seq and reads_cells) else None)}}
    rd.basecalling = _basecalling(ins, src["reads"])
    rd.inputs = [_input_entry(x, K) for x in ins]
    rd.warnings = warnings
    rd.partial_ids = sorted(set(rd.partial_ids))
    rd.waiting = sorted(set(rd.waiting))
    return rd


def _thin_stray_samples(cubes: dict[str, _Cube], samples: Samples, warnings: list[dict]) -> None:
    """A sample that no sample sheet lists and that holds under 0.1 % of the reads keeps its counts and loses its
    histograms; so does every sample beyond the 384 a bundle stores in full. Nothing is dropped from a total."""
    expected = {s["key"] for s in samples.rows if s.get("expected")}
    reads: dict[str, int] = {}
    for mod in ("len", "aln"):
        for (skey, _), cell in cubes[mod].cells.items():
            if mod == "len" or skey not in reads:
                reads[skey] = reads.get(skey, 0) + (cell.get("scalars") or {}).get("reads", 0)
    total = sum(reads.values())
    named = sorted((k for k in reads if k != UNCLASSIFIED), key=lambda k: -reads[k])
    thin = {k for k in named[MAX_SAMPLES:]}
    if expected and total:
        thin |= {k for k in named if k not in expected and reads[k] < 0.001 * total}
    if not thin:
        return
    for cube in cubes.values():
        for (skey, flt), cell in list(cube.cells.items()):
            if skey in thin:
                cube.cells[(skey, flt)] = {"scalars": cell.get("scalars") or {}}
    warnings.append({"level": "info", "code": "samples_counts_only", "value": len(thin),
                     "message": f"{len(thin)} sample name{'s' if len(thin) != 1 else ''} found in the data but not in the sample sheet hold"
                                f"{'s' if len(thin) == 1 else ''} under 0.1 % of the reads (or lie beyond the {MAX_SAMPLES} samples a run stores in full): "
                                "their reads and bases are counted, their distributions are not stored."})


def _conflicts(ins: list[_Input], run_key: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for x in ins:
        for need, msg in ((x.li.link.get("conflicts") or {}).get(run_key) or {}).items():
            if msg not in out.setdefault(need, []):
                out[need].append(msg)
    return out


def _collect(ins: list[_Input], samples: Samples, config: dict, barcoded: bool) -> dict:
    """Read the small result files that are counted for this run: pychopper cells and quantification counts."""
    cdna: dict[str, dict] = {}
    notes: list[dict] = []
    quant: dict[str, dict[str, dict]] = {}          # sample key -> role -> {counted, processed, tool, input, scope}
    pych: dict[str, dict] = {}
    for x in ins:
        if x.role not in ("pychopper.stats",) + tables.QUANT_PREFERENCE or x.role not in x.counted:
            continue
        li = x.li
        try:
            if x.role == "pychopper.stats":
                d = tables.read_pychopper(li.candidate.path)
            else:
                d = tables.read_quant(x.role, li.candidate.path)
        except (OSError, FlowCellarError) as e:
            x.state, x.read_error = "error", str(e)
            continue
        skey = samples.key_for(None, li.sample, barcoded=barcoded, source="file_name", project=project_of(li, config))
        x.samples.add(skey)
        several = len(li.run_keys) > 1
        if several:
            x.notes.append({"level": "info", "code": "covers_several_runs", "value": len(li.run_keys),
                            "message": f"This file has no per-read identifiers and covers the sample on {len(li.run_keys)} flow cells "
                                       f"({', '.join(sorted(li.run_keys))}); it is attached unsplit and must be counted once."})
        if x.role == "pychopper.stats":
            sc = d["scalars"]
            cell = {"scalars": dict(sc),
                    "counts": {"strand": build.counts(list(d["strand"]), list(d["strand"].values()), role="categories")}}
            if d["primer_configs"]:
                t = build.counts([k for k, _ in d["primer_configs"]], [v for _, v in d["primer_configs"]], role="categories")
                t["other"] = int(d["primer_configs_other"])
                cell["counts"]["cdna_primer_config_table"] = t
            cdna[skey] = cell
            pych[skey] = {"scalars": sc, "input": x.id, "several": several, "weak": _weak(li), "consistent": d["consistent"]}
            x.records = sc["input"]
            x.used_by.update(["cdna"])
            if not d["consistent"]:
                notes.append({"level": "warning", "code": "pychopper_layout", "input": x.id,
                              "message": f"{samples.name(skey)}: the classes of the pychopper table do not add up to PassReads; its numbers are shown as written."})
            if several:
                notes.append({"level": "info", "code": "covers_several_runs", "input": x.id,
                              "message": f"{samples.name(skey)}: one pychopper analysis for {len(li.run_keys)} flow cells; their reads were merged before pychopper ran."})
        else:
            quant.setdefault(skey, {})[x.role] = {**d, "input": x.id, "several": several, "weak": _weak(li)}
            x.records = d.get("processed")
    return {"cdna_cells": cdna, "cdna_notes": notes, "quant": quant, "pychopper": pych}


def _basecalling(ins: list[_Input], reads_src: list[_Input]) -> list[dict]:
    out: dict[tuple, dict] = {}
    for x in ins:
        if x.partial is None or x.role not in link.BAM_ROLES:
            continue
        f = x.partial.get("facts") or {}
        d = f.get("dorado") or {}
        models = f.get("models") or []
        if not d and not models:
            continue
        for model in (models or [None]):
            key = (model, x.li.variant)
            e = out.get(key)
            if e is None:
                e = out[key] = {"basecaller": "dorado", "version": d.get("version"), "model": model, "modbase_models": [],
                                "live": x.li.link.get("method") == "run_dir", "min_qscore": d.get("min_qscore"), "inputs": [],
                                "model_tier": model_tier(model), "trim": (not d["no_trim"]) if "no_trim" in d else None,
                                "barcode_both_ends": True if d.get("barcode_both_ends") else None,
                                "poly_a_estimation": d.get("estimate_poly_a"), "variant": x.li.variant, "primary": False}
            e["inputs"].append(x.id)
            if x in reads_src:
                e["primary"] = True
            for k in ("version", "min_qscore"):                  # an aligned file knows the model only; the basecalled file the rest
                if e.get(k) is None and d.get(k) is not None:
                    e[k] = d[k]
            if e["trim"] is None and "no_trim" in d:
                e["trim"] = not d["no_trim"]
            if d.get("barcode_both_ends"):
                e["barcode_both_ends"] = True
            if e["poly_a_estimation"] is None and d.get("estimate_poly_a") is not None:
                e["poly_a_estimation"] = d["estimate_poly_a"]
    return sorted(out.values(), key=lambda e: (not e["primary"], str(e["model"]), str(e["variant"])))


def _input_entry(x: _Input, run_key: str) -> dict:
    li, c = x.li, x.li.candidate
    lk = li.link
    not_counted = None
    conflicts = (lk.get("conflicts") or {}).get(run_key) or {}
    if not x.counted:
        not_counted = next(iter(conflicts.values()), None) or li.not_counted_reason or lk.get("not_counted") or "another input is the source of record"
    elif conflicts:
        not_counted = next(iter(conflicts.values()))
    state = x.state if x.counted else "skipped"
    e = {"id": x.id, "role": c.role, "path": clean_text(c.path), "root": c.root_id, "relpath": c.relpath, "size": c.size,
         "mtime": _mtime(c.mtime_ns),
         "fingerprint": li.fingerprint if li.fingerprint.get("method") in ("sha256", "quick-v1") else None,
         "state": state, "records": x.records,
         "scope": li.scope, "link": {"method": lk.get("method"), "confidence": lk.get("confidence"), "evidence": clean_text(lk.get("evidence") or "", 1000),
                                     "rule": lk.get("rule"), "confirmed": bool(lk.get("confirmed"))},
         "samples": sorted(x.samples), "used_by": sorted(x.used_by), "notes": x.notes, "alt_paths": [clean_text(p) for p in li.alt_paths],
         "needs": sorted({record.bundle_need(n) for n in li.needs}), "counted_for": sorted({record.bundle_need(n) for n in x.counted}),
         "variant": li.variant, "not_counted_reason": clean_text(not_counted, 1500) if not_counted else None}
    others = [k for k in li.run_keys if k != run_key]
    if others:
        e["shared_with"] = others
    if x.coverage:
        e["coverage"] = {k: x.coverage.get(k) for k in ("kind", "records", "records_total", "fraction")}
    return e
