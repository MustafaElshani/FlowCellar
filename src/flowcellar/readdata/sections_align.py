"""The alignment, over-time, poly(A) and transcriptome sections of a run, from the fragments ``alnpass`` wrote.

``sections.build_read_data`` hands every fragment to a ``Collector`` while it fills its cubes, then calls
``add_sections``. Nothing large is read here; the work is adding up clock tables and per-transcript counts of the
files of one sample, turning clock bins into bins counted from the run's time origin, and writing the sections:

==================  =====================  =========================================================================
section             need                   holds
==================  =====================  =========================================================================
``aln`` (extended)  ``genome``             strand, CIGAR totals, mismatches, aligned length, soft-clipped share, three
                                           grids, reads and aligned bases per contig, aligned length and identity over
                                           time
``aln.windows``     ``genome``             aligned bases per window along the genome (depth = bases / window length)
``dyn.basecalled``  ``reads``              output per ten minutes, length, quality and measured speed over time
``rna.polya``       ``reads``              poly(A) tail length from the ``pt`` tag
``rna.tx``          ``transcriptome``      alignment to the transcripts: status, strand, identity, reads per
                                           transcript, the saturation curve
``rna.coverage``    ``transcriptome``      coverage along the transcript body per transcript length class
==================  =====================  =========================================================================

**Time.** The pass bins start times on the clock (it does not know when a run started). Here bin 0 becomes the
clock step that holds the run's time origin, so a bin of the stored series begins ``params.clock.offset_*`` seconds
(zero or less, more than minus one step) from a whole multiple of the step on the run-time axis. Nothing is
interpolated and no read changes bin. A read stamped before bin 0 is counted in bin 0.

**Per sample and added up.** Up to twelve samples (``detail: standard``) every cell carries its own grids, depth
windows and over-time curves. The quartiles per hour are worked out here from the full tables of each cell, because
quartiles cannot be added; the tables themselves are stored once, for the passed reads of all samples, under
``run_level``, so the viewer can draw "all samples added up" from exact counts.
"""
from __future__ import annotations

import datetime as dt

import numpy as np

from flowcellar.aggregates import build, merge
from flowcellar.model import axes
from flowcellar.readdata import saturation
from flowcellar.readdata.alnacc import GRID_ID, GRID_LEN, GRID_Q, SPEED_AXIS, TimeTable, dense_counts
from flowcellar.readdata.reference import OTHER_CONTIGS, TX_CLASSES

VERSION = "0.1.1"                 # 0.1.1: whole-number read qualities give the value itself as the quartile of an hour
EXTRA_CUBES = ("time", "polya", "tx", "aln@tx")
UNCLASSIFIED = "unclassified"
# name in the facts of a fragment -> (value axis, curve id, unit of the curve, decimals)
READ_TABLES = {"length": (GRID_LEN, "len_over_time", "bases", 0), "quality": (GRID_Q, "qual_over_time", "phred", 2),
               "speed": (SPEED_AXIS, "dyn_speed_over_time", "bases/s", 1)}
ALN_TABLES = {"aligned_length": (GRID_LEN, "aln_aligned_length_over_time", "bases", 0),
              "identity": (GRID_ID, "aln_identity_over_time", "phred", 2)}
COARSE_GRID = {"aln_identity_vs_aligned": (GRID_LEN, "identity_q_1"), "aln_q_calibration": (GRID_Q, "identity_q_1")}
QUARTILES = (0.25, 0.5, 0.75)
ALN_TX_SCALARS = ("reads", "mapped_reads", "unmapped_reads", "secondary_records", "supplementary_records", "reads_with_supplementary",
                  "forward_reads", "reverse_reads", "alignment_columns", "edit_distance", "perfect_reads")


def need_of(mod: str, counted) -> str:
    """The need a fragment serves: an ``aln`` fragment is the genome alignment, or the alignment to the transcripts
    when its file is the source of record of that; ``tx`` is always the transcripts; everything else the reads."""
    if mod == "aln":
        return "transcriptome" if "transcriptome" in counted and "genome" not in counted else "genome"
    if mod == "tx":
        return "transcriptome"
    return "reads"


def _epoch(iso: str | None) -> float | None:
    if not iso:
        return None
    try:
        text = iso.replace("Z", "+00:00")
        if "." in text:                                      # MinKNOW writes nanoseconds; Python reads microseconds
            head, tail = text.split(".", 1)
            frac = tail[:tail.index("+")] if "+" in tail else tail
            text = head + "." + frac[:6] + "+00:00"
        return dt.datetime.fromisoformat(text).timestamp()
    except (ValueError, TypeError):
        return None


def _iso(epoch: float) -> str:
    return dt.datetime.fromtimestamp(epoch, dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def table_quantiles(a: np.ndarray, axis: str, qs=QUARTILES, whole: bool = False) -> np.ndarray:
    """Quantiles of every row of a table of counts over the cells of ``axis`` (under, bins..., over), by rank inside
    a bin, exactly as ``aggregates.derive.quantile`` reads one histogram. NaN for an empty row. Shape (rows, len(qs)).

    ``whole``: every value of a bin is the bin's lower edge (whole-number quality tags on a 1 Q axis): the quantile is
    that value, and nothing is read from a rank inside the bin."""
    e = axes.EDGES[axis]
    log = axes.AXES[axis]["type"] == "log10"
    c = a.astype(np.float64)
    cum = np.cumsum(c, axis=1)
    tot = cum[:, -1]
    out = np.full((len(a), len(qs)), np.nan)
    rows = np.flatnonzero(tot > 0)
    if not len(rows):
        return out
    for k, q in enumerate(qs):
        target = q * tot[rows]
        b = (cum[rows] < target[:, None]).sum(axis=1)                      # first cell whose cumulative count reaches the rank
        b = np.minimum(b, c.shape[1] - 1)
        before = np.where(b > 0, cum[rows, np.maximum(b - 1, 0)], 0.0)
        inside = c[rows, b]
        frac = np.where(inside > 0, (target - before) / np.maximum(inside, 1.0), 0.0)
        lo = e[np.clip(b - 1, 0, len(e) - 1)]
        hi = e[np.clip(b, 0, len(e) - 1)]
        val = lo + 0.0 * frac if whole else lo * (hi / lo) ** frac if log else lo + (hi - lo) * frac
        val = np.where(b <= 0, e[0], np.where(b >= len(e), e[-1], val))
        out[rows, k] = val
    return out


class _Tables:
    """Clock tables of one cell, added up over the files of its sample."""

    def __init__(self, spec: dict):
        self.t = {name: TimeTable(axes.n_cells(spec[name][0])) for name in spec}

    def add(self, facts: dict | None) -> None:
        for name, table in self.t.items():
            d = (facts or {}).get(name)
            if d:
                table.merge(TimeTable.from_json(table.width, d))


class Collector:
    """What the fragments of a run carry beside their cells, gathered while ``build_read_data`` fills its cubes."""

    def __init__(self, cubes: dict, cube_cls):
        for name in EXTRA_CUBES:
            cubes[name] = cube_cls()
        self.cubes = cubes
        self.reference: dict[str, dict] = {}                 # need -> the reference of the first fragment seen
        self.other_reference: dict[str, set] = {"genome": set(), "transcriptome": set()}
        self.read_tables: dict[tuple, _Tables] = {}          # (sample key, filter) -> length, quality, speed by hour
        self.output: dict[tuple, TimeTable] = {}             # (sample key, filter) -> reads and bases per ten minutes
        self.aln_tables: dict[tuple, _Tables] = {}
        self.tx_counts: dict[tuple, dict] = {}               # (sample key, filter) -> transcripts, genes (arrays or None)
        self.tx_untold: set[tuple] = set()                   # cells with a file whose per-transcript counts were not kept
        # the large members of the genome cells are added up here, in arrays, and never enter the cube
        self.grids: dict[tuple, dict[str, dict]] = {}        # (sample key, filter) -> name -> {x, y, m, outside, missing}
        self.windows: dict[tuple, np.ndarray] = {}           # (sample key, filter) -> aligned bases per window
        self.contigs: dict[tuple, dict] = {}                 # (sample key, filter) -> the contig table

    def _heavy(self, cell: tuple, body: dict, keep: bool) -> dict:
        """Take the grids, the depth windows and the contig table out of an ``aln`` cell (and add them up when
        ``keep``); what is left goes into the cube."""
        if keep:
            for name, h in (body.get("hists2d") or {}).items():
                g = self.grids.setdefault(cell, {}).get(name)
                if g is None:
                    g = self.grids[cell][name] = {"x": h["x"], "y": h["y"], "m": np.zeros((axes.n_bins(h["x"]), axes.n_bins(h["y"])), np.int64),
                                                  "outside": 0, "missing": 0}
                np.add.at(g["m"], (np.asarray(h.get("i", []), np.int64), np.asarray(h.get("j", []), np.int64)), np.asarray(h["n"], np.int64))
                g["outside"] += int(h.get("outside", 0))
                g["missing"] += int(h.get("missing", 0))
            s = (body.get("series") or {}).get("aln_depth_windows")
            if s:
                v = np.asarray(s["values"]["aligned_bases"], np.int64)
                cur = self.windows.get(cell)
                if cur is None or len(cur) != len(v):
                    self.windows[cell] = v.copy()
                else:
                    cur += v
            t = (body.get("counts") or {}).get("aln_contig_table")
            if t:
                self.contigs[cell] = t if cell not in self.contigs else merge.merge_counts(self.contigs[cell], t)
                self.contigs[cell]["w_name"] = "aligned_bases"       # the merge of two tables does not carry the name over
        return {k: v for k, v in body.items() if k not in ("hists2d", "series", "counts")}

    @staticmethod
    def route(mod: str, need: str) -> str:
        """The cube a fragment's cell goes to."""
        return "aln@tx" if mod == "aln" and need == "transcriptome" else mod

    def take(self, f: dict, need: str, skey: str, flt: str, body: dict) -> dict:
        """Gather the facts of one fragment and return the cell that goes into its cube."""
        mod = f.get("module")
        facts = f.get("facts") or {}
        cell = (skey, flt)
        if mod == "aln":
            ref = facts.get("reference")
            same = True
            if ref and ref.get("kind"):
                first = self.reference.setdefault(need, ref)
                if first.get("id") != ref.get("id"):         # another reference: its windows and contigs cannot be added
                    self.other_reference[need].add(skey)
                    same = False
            if need == "genome" and facts.get("time"):
                self.aln_tables.setdefault(cell, _Tables(ALN_TABLES)).add(facts["time"])
            if need == "genome" and same:
                body = self._heavy(cell, body, True)
                if (body.get("hists") or {}).get("aln_identity_hist") is None:      # no mapped read in this file
                    pass
            else:
                body = self._heavy(cell, body, False)
        elif mod == "time":
            t = facts.get("time") or {}
            self.read_tables.setdefault(cell, _Tables(READ_TABLES)).add(t)
            if t.get("output"):
                self.output.setdefault(cell, TimeTable(2)).merge(TimeTable.from_json(2, t["output"]))
        elif mod == "tx":
            got = self.tx_counts.setdefault(cell, {"transcripts": None, "genes": None, "ref": (facts.get("reference") or {}).get("id")})
            if not facts.get("transcripts") or (facts.get("reference") or {}).get("id") != got["ref"]:
                self.tx_untold.add(cell)
            else:
                for kind in ("transcripts", "genes"):
                    if facts.get(kind):
                        a = dense_counts(facts[kind])
                        got[kind] = a if got[kind] is None else got[kind] + a
        return body


# ------------------------------------------------------------------------------------------- the clock
class Clock:
    """Clock bins -> bins counted from the step of the clock that holds the run's time origin."""

    def __init__(self, tier1: dict | None, tables: list[TimeTable], outputs: list[TimeTable]):
        origin = _epoch((((tier1 or {}).get("identity") or {}).get("protocol_run") or {}).get("time_origin"))
        self.known = origin is not None
        if origin is None:                                   # no run folder behind the reads: count from the first read
            firsts = [t.trimmed()[0] for t in outputs]
            firsts = [f * 600.0 for f in firsts if f is not None] + [f * 3600.0 for f in (t.trimmed()[0] for t in tables) if f is not None]
            origin = min(firsts) if firsts else 0.0
        self.origin = origin
        self.bin0_600 = int(origin // 600)
        self.bin0_3600 = int(origin // 3600)

    def params(self) -> dict:
        return {"time_origin": _iso(self.origin) if self.known else None,
                "origin_known": self.known,
                "bin0_600s": _iso(self.bin0_600 * 600.0), "offset_600s_s": round(self.bin0_600 * 600.0 - self.origin, 3),
                "bin0_3600s": _iso(self.bin0_3600 * 3600.0), "offset_3600s_s": round(self.bin0_3600 * 3600.0 - self.origin, 3),
                "rule": "start times are binned on the UTC clock; bin 0 is the step that holds the run's time origin, so on "
                        "the run-time axis bin k begins offset + k x step seconds after the origin; a read stamped before "
                        "bin 0 is counted in bin 0"}

    def rows(self, table: TimeTable, step: int) -> np.ndarray:
        """The rows of a table from bin 0 on (earlier rows are added to row 0). Shape (bins, width); (0, width) when empty."""
        lo, a = table.trimmed()
        if lo is None:
            return np.zeros((0, table.width), np.int64)
        first = lo - (self.bin0_600 if step == 600 else self.bin0_3600)
        if first >= 0:
            out = np.zeros((first + len(a), table.width), np.int64)
            out[first:] = a
            return out
        early = min(-first, len(a))
        out = np.zeros((max(len(a) + first, 1), table.width), np.int64)
        out[0] += a[:early].sum(axis=0)
        rest = a[early:]
        out[:len(rest)] += rest
        return out


def _curve(rows: np.ndarray, axis: str, unit: str, digits: int, whole: bool = False) -> dict | None:
    """Reads and the quartiles of one quantity per clock hour, for the hours that hold a read."""
    if not len(rows):
        return None
    n = rows.sum(axis=1)
    keep = np.flatnonzero(n > 0)
    if not len(keep):
        return None
    q = table_quantiles(rows[keep], axis, whole=whole)

    def col(k):
        return [int(round(v)) if digits == 0 else round(float(v), digits) for v in q[:, k]]
    return {"x": [int(v) for v in keep], "y": {"reads": [int(v) for v in n[keep]], "q25": col(0), "median": col(1), "q75": col(2)},
            "x_unit": "bin of clock_3600s", "y_unit": unit, "x_scale": "linear"}


def whole_number_qualities(qual_cells: list[dict]) -> bool:
    """True when every read quality of these cells is a whole number (files whose ``qs`` tag is an integer): on the
    0.25 Q axis of the quality histogram only the bins that begin at a whole number hold a read."""
    seen = 0
    for c in qual_cells:
        h = (c.get("hists") or {}).get("qual_hist")
        if not h or h.get("axis") != "qscore_0p25":
            if h and build.hist_full(h, "n").sum() > 0:
                return False                                  # stored on another axis: nothing can be said
            continue
        n = np.asarray(build.hist_full(h, "n"))[1:-1]         # without the cells under and over the axis
        if n.reshape(-1, 4)[:, 1:].sum() > 0:
            return False
        seen += int(n.sum())
    return seen >= 100


def _grid(rows: np.ndarray, axis: str) -> dict:
    """A table by clock hour as a stored 2-D histogram on clock_3600s x axis."""
    inner = rows[:, 1:-1]
    return build.hist2d("clock_3600s", axis, inner, outside=int(rows.sum() - inner.sum()))


def _sum_rows(parts: list[np.ndarray], width: int) -> np.ndarray:
    n = max((len(p) for p in parts), default=0)
    out = np.zeros((n, width), np.int64)
    for p in parts:
        out[:len(p)] += p
    return out


def _coarsen_grid(h: dict, x: str, y: str) -> dict:
    """A stored 2-D histogram re-binned exactly onto coarser axes of the same families."""
    if (h["x"], h["y"]) == (x, y):
        return h
    fx, fy = axes.coarsen_factor(h["x"], x), axes.coarsen_factor(h["y"], y)
    m = np.zeros((axes.n_bins(x), axes.n_bins(y)), np.int64)
    np.add.at(m, (np.asarray(h.get("i", []), np.int64) // fx, np.asarray(h.get("j", []), np.int64) // fy), np.asarray(h["n"], np.int64))
    return build.hist2d(x, y, m, outside=int(h.get("outside", 0)), missing=h.get("missing"))


def _rate(num, den, digits: int = 6):
    return None if not den else round(num / den, digits)


# ---------------------------------------------------------------------------------------- the sections
def add_sections(col: Collector, rd, *, section, samples, order: list[str], detail: str, tier1: dict | None, src: dict,
                 ins: list, conflicts: dict, reads_pop: dict | None, aln_pop: dict | None, barcoded: bool, labelled: dict) -> None:
    """Extend ``rd.modules["aln"]`` and add the sections of the table above. ``section`` is the builder of
    ``build_read_data`` (it fills state, coverage, inputs and the partial reason from the inputs behind a need)."""
    standard = detail == "standard"
    clock = Clock(tier1, [t for g in list(col.read_tables.values()) + list(col.aln_tables.values()) for t in g.t.values()],
                  list(col.output.values()))
    pos = {k: i for i, k in enumerate(order)}

    def by_order(keys):
        return sorted(keys, key=lambda k: (pos.get(k[0], len(pos)), k[1] != "pass", k))

    def waiting(need: str) -> bool:
        return any(x.waiting and need in x.counted for x in ins)

    def conflict_notes(need: str) -> list[dict]:
        return [{"level": "warning", "code": "two_candidate_sources", "message": "Not assessed: " + m[:1150]} for m in conflicts.get(need, [])]

    def over_time(tables: dict, spec: dict, cells: list[dict], whole: tuple = ()) -> dict:
        """Curves into the cells (standard detail); returns the run-level grids of the passed and failed reads.
        ``whole``: the tables whose values are whole numbers (see ``table_quantiles``)."""
        by_filter: dict[str, dict[str, list]] = {}
        have = {tuple(c["at"]): c for c in cells}
        for key in by_order(tables):
            cell = have.get(key)
            if cell is None:
                continue
            for name, (axis, cid, unit, digits) in spec.items():
                rows = clock.rows(tables[key].t[name], 3600)
                if not len(rows):
                    continue
                by_filter.setdefault(key[1], {}).setdefault(name, []).append(rows)
                if standard:
                    cv = _curve(rows, axis, unit, digits, whole=name in whole)
                    if cv:
                        cell.setdefault("curves", {})[cid] = cv
        grids = {}
        for flt, per in sorted(by_filter.items()):
            for name, parts in per.items():
                axis, cid = spec[name][0], spec[name][1]
                grids[f"{cid}_{flt}"] = _grid(_sum_rows(parts, axes.n_cells(axis)), axis)
        return grids

    # ---- aln: the genome alignment
    aln = rd.modules.get("aln")
    ref = col.reference.get("genome")
    if aln is not None and aln.get("cells"):
        pooled: dict[str, dict] = {}
        for cell in aln["cells"]:
            key = tuple(cell["at"])
            if not cell.get("hists"):                        # a stray sample keeps its counts only
                continue
            for name, g in (col.grids.get(key) or {}).items():
                if standard:
                    h = build.hist2d(g["x"], g["y"], g["m"], outside=g["outside"], missing=g["missing"])
                    cell.setdefault("hists2d", {})[name] = _coarsen_grid(h, *COARSE_GRID[name]) if name in COARSE_GRID else h
                elif key[1] == "pass":
                    p = pooled.setdefault(name, {"x": g["x"], "y": g["y"], "m": np.zeros_like(g["m"]), "outside": 0, "missing": 0})
                    p["m"] += g["m"]
                    p["outside"] += g["outside"]
                    p["missing"] += g["missing"]
            if key in col.contigs:
                cell.setdefault("counts", {})["aln_contig_table"] = col.contigs[key]
        run_level = dict(aln.get("run_level") or {})
        for name, g in pooled.items():
            h = build.hist2d(g["x"], g["y"], g["m"], outside=g["outside"], missing=g["missing"])
            run_level.setdefault("hists2d", {})[name] = _coarsen_grid(h, *COARSE_GRID[name]) if name in COARSE_GRID else h
        if run_level:
            aln["run_level"] = run_level
        light = [{"at": c["at"], "scalars": c.get("scalars") or {}, "counts": c.get("counts") or {}} for c in aln["cells"]]
        tot = merge.rollup(light, ["sample", "filter"], [], where={"filter": ["pass"]} if any(c["at"][1] == "pass" for c in light) else None)
        sc = (tot[0].get("scalars") if tot else None) or {}
        mapped = sc.get("mapped_reads") or 0
        contig = ((tot[0].get("counts") or {}).get("aln_contig_table") if tot else None) or {}
        mito_name = (ref or {}).get("mitochondrial")
        mito = contig["n"][contig["cats"].index(mito_name)] if mito_name and mito_name in (contig.get("cats") or []) else None
        has_cigar = "inserted_bases" in sc
        columns = (sc.get("match_or_mismatch_bases") or 0) + (sc.get("inserted_bases") or 0) + (sc.get("deleted_bases") or 0)
        ref_bases = (sc.get("match_or_mismatch_bases") or 0) + (sc.get("deleted_bases") or 0)
        read_bases = (sc.get("match_or_mismatch_bases") or 0) + (sc.get("inserted_bases") or 0) + (sc.get("soft_clipped_bases") or 0)
        aln["summary"].update({
            "aln_strand_frac": _rate(sc.get("forward_reads") or 0, mapped) if has_cigar else None,
            "aln_mismatch_frac": _rate(sc.get("mismatched_bases") or 0, sc.get("alignment_columns")) if has_cigar else None,
            "aln_insertion_frac": _rate(sc.get("inserted_bases") or 0, columns) if has_cigar else None,
            "aln_deletion_frac": _rate(sc.get("deleted_bases") or 0, columns) if has_cigar else None,
            "aln_soft_clipped_frac": _rate(sc.get("soft_clipped_bases") or 0, read_bases) if has_cigar else None,
            "aln_chrm_frac": _rate(mito, mapped) if mito is not None else None,
            "aln_mean_depth": _rate(ref_bases, (ref or {}).get("total_length"), 4) if ref and ref.get("kind") == "genome" and has_cigar else None})
        aln["params"].update({
            "reference": {k: ref.get(k) for k in ("kind", "evidence", "id", "sequences", "total_length", "mitochondrial")} if ref else None,
            "aligned_length_rule": "read bases inside the alignment (matches, mismatches and insertions); clipped bases are outside",
            "mismatch_rule": "NM minus inserted minus deleted bases, on records with an NM tag; the rate is per alignment column of those records",
            "base_breakdown_rule": "of all alignment columns and clipped bases of the primary mapped records: match, mismatch, insertion, deletion, soft clip, hard clip",
            "soft_clip_rule": "soft-clipped bases of a read / its length",
            "depth_rule": "aligned bases = reference positions under a match, mismatch or deletion column of a primary record; skipped regions (N) are not covered",
            "identity_axis": "identity on the Phred scale, -10 log10(1 - identity): the same number is the quality measured from the alignment"})
        if col.other_reference["genome"]:
            aln["notes"].append({"level": "warning", "code": "several_references",
                                 "message": f"The alignment files of this run name more than one reference; contigs and depth are shown for the first ({ref.get('sequences'):,} sequences) and left out for {len(col.other_reference['genome'])} sample(s) aligned to another."})

        # ---- aln.time: aligned length and identity over time (its own section: quartiles cannot be added up, and the
        # cells of aln are added up by the checks)
        have = {tuple(c["at"]) for c in aln["cells"] if c.get("hists")}
        tcells = [{"at": list(key), "scalars": {"reads_with_start_time": int(clock.rows(col.aln_tables[key].t["aligned_length"], 3600).sum())}}
                  for key in by_order(col.aln_tables) if key in have]
        tcells = [c for c in tcells if c["scalars"]["reads_with_start_time"]]
        if tcells:
            grids = over_time(col.aln_tables, ALN_TABLES, tcells)
            sec = section("aln.time", "genome", tcells, aln_pop, {},
                          {"clock": clock.params(), "identity_axis": aln["params"]["identity_axis"],
                           "quartile_rule": "per clock hour, read from the full table of that hour by rank inside a bin"}, [])
            if grids:
                sec["run_level"] = {"hists2d": grids}
            rd.modules["aln.time"] = sec

        # ---- aln.windows
        if ref and ref.get("kind") == "genome" and col.windows:
            wcells, parts = [], []
            for key in by_order(col.windows):
                if key not in have:
                    continue
                v = col.windows[key]
                if standard:
                    wcells.append({"at": list(key), "scalars": {"aligned_bases": int(v.sum())},
                                   "series": {"aln_depth_windows": build.series("genome_window", {"aligned_bases": v}, unit="bases", n_bins=len(v))}})
                if key[1] == "pass":
                    parts.append(v.reshape(-1, 1))
            if wcells or parts:
                params = {"reference": ref, "depth_rule": aln["params"]["depth_rule"],
                          "window_rule": "one window length for the whole reference, chosen so that it is about 3,000 windows; every contig listed starts a new window, so its last window is shorter; the other contigs are added up and have no window"}
                sec = section("aln.windows", "genome", wcells, aln_pop, {"aln_mean_depth": aln["summary"]["aln_mean_depth"]}, params, [])
                if not standard and parts:
                    total = _sum_rows(parts, 1)[:, 0]
                    sec["run_level"] = {"series": {"aln_depth_windows": build.series("genome_window", {"aligned_bases": total}, unit="bases", n_bins=len(total))}}
                    sec["state"] = aln["state"]
                    sec.pop("partial_reason", None)
                    if aln["state"] == "partial":
                        sec["partial_reason"] = aln.get("partial_reason")
                rd.modules["aln.windows"] = sec

    # ---- dyn.basecalled: output and the read's own numbers over time
    tcells_raw = {tuple(c["at"]): c for c in col.cubes["time"].listed(order)}
    if tcells_raw:
        cells = []
        ten: dict[str, list] = {}
        for key in by_order(tcells_raw):
            c = tcells_raw[key]
            cell = {"at": list(key), "scalars": dict(c.get("scalars") or {})}
            if c.get("hists"):
                cell["hists"] = dict(c["hists"])
            out = col.output.get(key)
            if out is not None and c.get("hists"):
                rows = clock.rows(out, 600)
                if len(rows):
                    ten.setdefault(key[1], []).append(rows)
                    if standard:
                        hours = np.zeros((-(-(len(rows) + clock.bin0_600 % 6) // 6), 2), np.int64)
                        np.add.at(hours, (np.arange(len(rows)) + clock.bin0_600 % 6) // 6, rows)
                        cell["series"] = {"dyn_output_series": build.series("clock_3600s", {"reads": hours[:, 0], "bases": hours[:, 1]})}
            cells.append(cell)
        # whole-number quality tags (qs:i): a 1 Q step of the table holds one value, so the median of an hour is that
        # value itself; reading it from a rank inside the step would put it up to one Q too high
        whole_q = whole_number_qualities(list(col.cubes["qual"].cells.values()))
        grids = over_time(col.read_tables, READ_TABLES, cells, whole=("quality",) if whole_q else ())
        tot = merge.rollup(list(tcells_raw.values()), ["sample", "filter"], [], where={"filter": ["pass"]} if any(k[1] == "pass" for k in tcells_raw) else None)
        sc = (tot[0].get("scalars") if tot else None) or {}
        hs = ((tot[0].get("hists") or {}).get("dyn_speed_basecalled_hist") if tot else None)
        from flowcellar.aggregates import derive
        trimmed = bool((reads_pop or {}).get("trimmed"))
        summary = {"dyn_speed_basecalled": round(derive.median(hs), 1) if hs and derive.total(hs) else None,
                   "dyn_speed_basecalled_mean": _rate((sc.get("bases_with_duration") or 0) * 1000, sc.get("duration_ms"), 1),
                   "dyn_reads_with_start_time_frac": _rate(sc.get("reads_with_start_time") or 0, sc.get("reads"))}
        params = {"clock": clock.params(),
                  "speed_rule": "bases of the read / seconds it spent in the pore (du tag), on reads with du of 1 s or more",
                  "speed_on_trimmed_reads": trimmed,
                  "quartile_rule": "per clock hour, read from the full table of that hour by rank inside a bin",
                  "qual_whole_numbers": whole_q}
        if whole_q:
            params["qual_quartile_rule"] = ("the files give every read quality as a whole number: a quartile of read quality is the value "
                                            "of the 1 Q step that holds its rank, not a point inside the step")
        notes = conflict_notes("reads")
        if trimmed:
            notes.append({"level": "info", "code": "speed_on_trimmed_reads",
                          "message": "These reads were trimmed after basecalling while their duration covers the whole read, so the measured speed is a little lower than the true one."})
        sec = section("dyn.basecalled", "reads", cells, reads_pop, summary, params, notes)
        run_level: dict = {}
        if ten:
            vals = {}
            for flt, parts in sorted(ten.items()):
                total = _sum_rows(parts, 2)
                vals[f"reads_{flt}"], vals[f"bases_{flt}"] = total[:, 0], total[:, 1]
            run_level["series"] = {"dyn_output_series_600s": build.series("clock_600s", vals)}
        if grids:
            run_level["hists2d"] = grids
        if run_level:
            sec["run_level"] = run_level
        rd.modules["dyn.basecalled"] = sec

    # ---- rna.polya
    pcells = col.cubes["polya"].listed(order)
    if pcells:
        cells = []
        for c in pcells:
            cell = {"at": c["at"], "scalars": dict(c.get("scalars") or {})}
            h = (c.get("hists") or {}).get("rna_polya_hist")
            if h:
                if not standard:
                    f = axes.coarsen_factor(h["axis"], "polya_5nt")
                    h = build.hist("polya_5nt", axes.coarsen(build.hist_full(h, "n"), f), missing=h.get("missing"))
                cell["hists"] = {"rna_polya_hist": h}
            cells.append(cell)
        tot = merge.rollup(pcells, ["sample", "filter"], [], where={"filter": ["pass"]} if any(c["at"][1] == "pass" for c in pcells) else None)
        sc = (tot[0].get("scalars") if tot else None) or {}
        h = ((tot[0].get("hists") or {}).get("rna_polya_hist") if tot else None)
        from flowcellar.aggregates import derive
        summary = {"rna_polya_detect_frac": _rate(sc.get("tail_found_reads") or 0, sc.get("reads")),
                   "rna_polya_median": round(derive.median(h), 1) if h and derive.total(h) else None,
                   "rna_polya_tagged_frac": _rate(sc.get("reads_with_tag") or 0, sc.get("reads"))}
        params = {"source": "the pt tag the basecaller writes on each read when it is asked to estimate poly(A)",
                  "classes": "a tail was found (pt above 0); the tag says zero; the tag says none was found (pt = -1); the read has no tag",
                  "histogram": "reads with a tail found; one bin per nucleotide up to 1,000"}
        rd.modules["rna.polya"] = section("rna.polya", "reads", cells, reads_pop, summary, params, conflict_notes("reads"))

    # ---- rna.tx and rna.coverage: the alignment to the transcripts
    tx_aln = {tuple(c["at"]): c for c in col.cubes["aln@tx"].listed(order)}
    tx_body = {tuple(c["at"]): c for c in col.cubes["tx"].listed(order)}
    tref = col.reference.get("transcriptome")
    if tx_aln or conflicts.get("transcriptome") or waiting("transcriptome"):
        from flowcellar.readdata import record
        from flowcellar.readdata.sections import aln_summary, population
        xs = src.get("transcriptome") or []
        unmapped = sum((c.get("scalars") or {}).get("unmapped_reads", 0) for c in tx_aln.values())
        tx_pop = None
        if xs:
            tx_pop = population("align.transcriptome_bam", [x.partial.get("facts") or {} for x in xs if x.partial], barcoded=barcoded,
                                has_unclassified=False, has_fail=False, unmapped_reads=unmapped,
                                labelled=labelled.get("transcriptome", False) or any(x.li.sample and x.li.sample != record.ALL for x in xs))
        cells = []
        pooled_t = pooled_g = None
        pooled_ok = bool(tx_aln)
        for key in by_order(tx_aln):
            c = tx_aln[key]
            scal = c.get("scalars") or {}
            cell = {"at": list(key), "scalars": {k: scal[k] for k in ALN_TX_SCALARS if k in scal}}
            hs = c.get("hists") or {}
            if hs.get("aln_identity_hist"):
                h = hs["aln_identity_hist"]
                if not standard:
                    f = axes.coarsen_factor(h["axis"], "identity_q_0p5")
                    h = build.hist("identity_q_0p5", axes.coarsen(build.hist_full(h, "n"), f), missing=h.get("missing"))
                cell["hists"] = {"aln_identity_hist": h}
            counts = col.tx_counts.get(key)
            usable = counts is not None and key not in col.tx_untold and counts["transcripts"] is not None
            if usable and hs:
                if standard:
                    curves = {"rna_reads_per_transcript": saturation.reads_per_item(counts["transcripts"]),
                              "rna_saturation_curve": saturation.curve(counts["transcripts"], counts["genes"])}
                    cell["curves"] = {k: v for k, v in curves.items() if v}
                if key[1] == "pass":
                    pooled_t = counts["transcripts"] if pooled_t is None else pooled_t + counts["transcripts"]
                    if counts["genes"] is not None:
                        pooled_g = counts["genes"] if pooled_g is None else pooled_g + counts["genes"]
            elif key[1] == "pass" and hs:
                pooled_ok = False
            cells.append(cell)
        tot = merge.rollup(list(tx_aln.values()), ["sample", "filter"], [], where={"filter": ["pass"]} if any(k[1] == "pass" for k in tx_aln) else None)
        tsum = aln_summary(tot[0] if tot else None)
        tsc = (tot[0].get("scalars") if tot else None) or {}
        all_told = pooled_ok and pooled_t is not None
        summary = {"rna_tx_aligned_frac": tsum["aln_mapped_frac"], "rna_tx_identity_median": tsum["aln_identity_median"],
                   "rna_antisense_frac": _rate(tsc.get("reverse_reads") or 0, tsc.get("mapped_reads")) if "reverse_reads" in tsc else None,
                   "rna_tx_transcripts_detected": int((pooled_t > 0).sum()) if all_told else None,
                   "rna_tx_genes_detected": int((pooled_g > 0).sum()) if all_told and pooled_g is not None else None}
        params = {"reference": {k: tref.get(k) for k in ("kind", "evidence", "id", "sequences", "total_length", "length_classes")} if tref else None,
                  "saturation_rule": "the expected number of transcripts (genes) with at least one read if n of the N primary mapped reads had been sequenced, any n of them: the sum over transcripts of 1 - C(N - c, n) / C(N, n), from the read count c of each; no read is drawn and none is identified",
                  "genes_from": "the gene named in each reference name (transcript|gene|...)" if all_told and pooled_g is not None else None,
                  "reads_per_transcript_rule": "transcripts by their number of primary mapped reads, in classes that double; a read counts for the transcript of its primary alignment only",
                  "strand_rule": "a read on the reverse strand of its transcript is antisense"}
        notes = conflict_notes("transcriptome")
        if tx_aln and not all_told:
            notes.append({"level": "info", "code": "reads_per_transcript_not_kept",
                          "message": "Reads per transcript and the saturation curve are not shown for all samples added up: at least one file was read without the count of each transcript."})
        sec = section("rna.tx", "transcriptome", cells, tx_pop, summary, params, notes)
        if all_told:
            curves = {"rna_reads_per_transcript": saturation.reads_per_item(pooled_t), "rna_saturation_curve": saturation.curve(pooled_t, pooled_g)}
            sec["run_level"] = {"curves": {k: v for k, v in curves.items() if v}}
        rd.modules["rna.tx"] = sec

        ccells = []
        class_reads = {cid: 0 for cid, _ in TX_CLASSES}
        merged: dict[tuple, dict] = {}                       # (sample, class) summed over the filters
        for key in by_order(tx_body):
            c = tx_body[key]
            for cid, _ in TX_CLASSES:
                h = (c.get("hists") or {}).get("rna_body_coverage_" + cid)
                n = (c.get("scalars") or {}).get("reads_" + cid, 0)
                class_reads[cid] += n
                if h is None:
                    continue
                cur = merged.get((key[0], cid))
                body = {"scalars": {"reads": int(n)}, "hists": {"rna_body_coverage": h}}
                merged[(key[0], cid)] = body if cur is None else merge.merge_cell(cur, body)
        for (skey, cid), body in merged.items():
            ccells.append({"at": [skey, cid], **body})
        if ccells:
            sec = section("rna.coverage", "transcriptome", ccells, tx_pop,
                          {"rna_body_coverage_reads": int(sum(class_reads.values()))},
                          {"reference": params["reference"],
                           "coverage_rule": "per transcript length class and for each of 100 equal steps from the 5' end to the 3' end of the transcript: the reads whose primary alignment covers the middle of that step",
                           "length_classes": [{"id": cid, "below_bases": hi} for cid, hi in TX_CLASSES]}, [])
            sec["dims"] = ["sample", "read_class"]
            have = {c["at"][0] for c in ccells}
            if sec["state"] == "partial" and sec.get("partial_reason") == "samples_missing" and all(k in have for k in samples.expected_keys()):
                sec["state"] = "ok"
                sec.pop("partial_reason", None)
            rd.modules["rna.coverage"] = sec


def contig_rows(table: dict | None, ref: dict | None) -> list[dict]:
    """The contig table of a cell with the length of every contig: for the index card and for tests."""
    if not table:
        return []
    length = {c["name"]: c["length"] for c in (ref or {}).get("contigs") or []}
    length[OTHER_CONTIGS] = ((ref or {}).get("other") or {}).get("length")
    return [{"contig": name, "reads": n, "aligned_bases": w, "length": length.get(name)} for name, n, w in zip(table["cats"], table["n"], table["w"])]
