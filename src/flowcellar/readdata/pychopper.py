"""Every category of a pychopper statistics table (``Category / Name / Value``), per sample and summed per run.

``tables.read_pychopper`` reduces a table to the three classes read fate needs. This module keeps the whole table:
the section ``cdna.pychopper`` holds one cell per sample and, in ``run_level``, their exact sum.

What each number is a count of was worked out on 55 real tables written over three years (one layout: eleven
categories, the ``ReadStats`` rows in two places because the UMI rows come after ``Hits``). The identities in the
right-hand column held on every one of them; they are tested again on every table that is read, and a table on
which one fails is stored as written, with a note that names the identity.

==================  ==========================================================================  =========================
category            a count of                                                                  adds up to
==================  ==========================================================================  =========================
``ReadStats``       ``PassReads``: reads pychopper classified. ``QcFail``: reads it dropped     PassReads + QcFail =
                    before classifying (mean quality under its limit). ``LenFail``: PIECES it   reads given to pychopper
                    dropped after classifying (shorter than its limit). ``Umi_detected`` and
                    ``Umi_detected_final`` (only when UMI detection was on): pieces in which a
                    UMI was found, and those still counted after pychopper's own checks
``Classification``  ``Primers_found``: reads with one full-length piece. ``Rescue``: PIECES     Primers_found + rescued
                    cut out of fused reads. ``Unusable``: reads without a usable piece          reads + Unusable = PassReads
``Strand``          full-length reads by orientation                                            Primers_found
``RescueStrand``    rescued pieces by orientation                                               Rescue
``RescueSegmentNr`` fused reads by the number of pieces cut out of them (2, 3, ...)             sum of k x n = Rescue
``RescueHitNr``     rescued PIECES by the number of primer hits on the read they came from      Rescue
``UnclassHitNr``    unusable reads by the number of primer hits on them                         Unusable
``Unusable`` rows   usable reads by the number of bases that were NOT kept (primers, adapters   Primers_found + rescued
                    and whatever lies outside the pieces)                                       reads
``Hits``            reads with at least one primer hit, by primer configuration                 PassReads - UnclassHitNr[0]
``AutotuneSample``  for each cut-off q tried on a sample of the reads: reads classified         -
``Parameter``       ``Cutoff(q)``: the cut-off pychopper chose, the q with the most reads       -
                    classified in its tuning sample
==================  ==========================================================================  =========================

Reads written by pychopper (what every later step reads) = Primers_found + Rescue - LenFail.
"""
from __future__ import annotations

import os

import numpy as np

from flowcellar.aggregates import build, merge
from flowcellar.model import axes
from flowcellar.util.errors import InputError, ResourceLimitError
from flowcellar.util.text import clean_text

VERSION = "0.2.1"                    # 0.2.1: a category the table should hold and does not is named (note pychopper_categories_absent)
SECTION = "cdna.pychopper"
MAX_BYTES = 16 << 20                 # a statistics table is 10 to 80 kB; anything far larger is not one
PRIMER_CONFIGS_KEPT = 12             # the most common primer configurations of the run are stored by name
CONFIG_NAME_CHARS = 120
KNOWN_CATEGORIES = ("ReadStats", "Classification", "Strand", "RescueStrand", "UnclassHitNr", "RescueHitNr", "RescueSegmentNr",
                    "Unusable", "Hits", "AutotuneSample", "Parameter")
INT_NAMED = ("UnclassHitNr", "RescueHitNr", "RescueSegmentNr", "Unusable")
#: the one scalar of a cell that is not a sum: two samples share most of their primer configurations, so the number
#: for several samples together is given as the largest of theirs. (The cut-off is one value per table and is stored
#: with the tuning curve, which is never added up.)
SCALAR_MERGE = {"primer_configs": "max"}
CLASS_CATS = ("full_length", "rescued", "unusable")

IDENTITIES = {
    "classes_add_up": "Primers_found + rescued reads + Unusable = PassReads",
    "rescue_pieces": "the pieces of RescueSegmentNr add up to Rescue",
    "strand_sum": "the Strand rows add up to Primers_found",
    "rescue_strand_sum": "the RescueStrand rows add up to Rescue",
    "unusable_hits_sum": "the UnclassHitNr rows add up to Unusable",
    "rescue_hits_sum": "the RescueHitNr rows add up to Rescue",
    "trimmed_sum": "the Unusable rows add up to Primers_found + rescued reads",
    "hits_sum": "the Hits rows add up to PassReads - UnclassHitNr[0]",
    "cutoff_in_curve": "Cutoff(q) is one of the AutotuneSample values",
}


def _int(v) -> int:
    return int(round(float(v or 0)))


def _is(q: float, cutoff: float) -> bool:
    """Whether a tuning value is the cut-off (both are written with 16 digits; allow for the last one)."""
    return abs(q - cutoff) <= 1e-9 * max(1.0, abs(cutoff))


def parse(path: str) -> dict:
    """The table as ``{"cats": {category: {name: value}}, "order": [categories], "rows", "rows_skipped",
    "rows_repeated"}``. A category may be split over the file; a name repeated inside a category is added up."""
    size = os.path.getsize(path)
    if size > MAX_BYTES:
        raise ResourceLimitError(f"{path}: {size >> 20} MB is too large for a pychopper statistics table")
    with open(path, "rb") as fh:
        text = fh.read().decode("utf-8-sig", errors="replace")
    lines = text.splitlines()
    if not lines or [c.strip() for c in lines[0].split("\t")][:3] != ["Category", "Name", "Value"]:
        raise InputError(f"{path}: not a pychopper statistics table (no Category / Name / Value header)")
    cats: dict[str, dict[str, float]] = {}
    rows = skipped = repeated = 0
    for line in lines[1:]:
        p = line.rstrip("\r\n").split("\t")
        if len(p) != 3:
            skipped += bool(line.strip())
            continue
        try:
            v = float(p[2])
        except ValueError:
            skipped += 1
            continue
        if v != v or v in (float("inf"), float("-inf")):
            skipped += 1
            continue
        rows += 1
        d = cats.setdefault(p[0].strip(), {})
        name = p[1].strip()
        if name in d:
            repeated += 1
            d[name] += v
        else:
            d[name] = v
    return {"cats": cats, "order": list(cats), "rows": rows, "rows_skipped": skipped, "rows_repeated": repeated}


def _by_int(d: dict[str, float]) -> tuple[dict[int, int], int]:
    """Rows whose name is a whole number -> {number: count}; the rows whose name is not one are counted."""
    out: dict[int, int] = {}
    bad = 0
    for k, v in d.items():
        try:
            f = float(k)
            if f != int(f) or f < 0:
                raise ValueError
        except ValueError:
            bad += 1
            continue
        out[int(f)] = out.get(int(f), 0) + _int(v)
    return dict(sorted(out.items())), bad


def read_stats(path: str) -> dict:
    """One statistics table, every category.

    ``scalars`` (whole numbers; a number the table does not hold is None): ``input`` (reads given to pychopper),
    ``pass_reads``, ``qc_fail``, ``len_fail`` (pieces), ``primers_found``, ``rescued_reads``, ``rescue_segments``
    (pieces), ``unusable``, ``oriented`` (pieces kept before the length filter = primers_found + rescue_segments),
    ``written``, ``strand_plus``, ``strand_minus``, ``rescue_strand_plus``, ``rescue_strand_minus``,
    ``unusable_no_hit``, ``unusable_multi_hit`` (three hits or more), ``hits_reads`` (reads with a primer hit),
    ``trimmed_reads``, ``trimmed_bases``, ``umi_detected``, ``umi_detected_final``, ``primer_configs`` (distinct).

    Besides: ``cutoff_q`` (the cut-off chosen, a float, or None), ``segments``, ``unusable_hits``, ``rescue_hits``, ``trimmed`` ({whole number: count}), ``configs``
    ([(name, reads)], most common first), ``autotune`` ([(q, reads)] by q), ``parameters``, ``identities``
    ({id: True | False | None}), ``failed`` (ids that are False), ``categories``, ``unknown_categories``, ``rows``,
    ``rows_skipped``.
    """
    t = parse(path)
    cats = t["cats"]
    rs, cl = cats.get("ReadStats", {}), cats.get("Classification", {})
    if "PassReads" not in rs or "Primers_found" not in cl:
        raise InputError(f"{path}: pychopper statistics without ReadStats PassReads or Classification Primers_found")
    bad_names = 0
    ints: dict[str, dict[int, int]] = {}
    for c in INT_NAMED:
        ints[c], bad = _by_int(cats.get(c, {}))
        bad_names += bad
    seg, uh, rh, trimmed = ints["RescueSegmentNr"], ints["UnclassHitNr"], ints["RescueHitNr"], ints["Unusable"]
    strand, rstrand = cats.get("Strand", {}), cats.get("RescueStrand", {})
    hits: dict[str, int] = {}
    for k, v in cats.get("Hits", {}).items():
        name = clean_text(k, CONFIG_NAME_CHARS) or "(no name)"
        hits[name] = hits.get(name, 0) + _int(v)
    configs = sorted(hits.items(), key=lambda kv: (-kv[1], kv[0]))
    autotune = []
    for k, v in cats.get("AutotuneSample", {}).items():
        try:
            autotune.append((float(k), _int(v)))
        except ValueError:
            bad_names += 1
    autotune.sort()
    params = {clean_text(k, 60): float(v) for k, v in cats.get("Parameter", {}).items()}
    cutoff = params.get("Cutoff(q)")

    pf, rescue, unusable = _int(cl.get("Primers_found")), _int(cl.get("Rescue")), _int(cl.get("Unusable"))
    rescued_reads = sum(seg.values())
    sc: dict = {
        "pass_reads": _int(rs.get("PassReads")), "qc_fail": _int(rs.get("QcFail")), "len_fail": _int(rs.get("LenFail")),
        "primers_found": pf, "rescued_reads": rescued_reads, "rescue_segments": rescue, "unusable": unusable,
        "oriented": pf + rescue,
        "strand_plus": _int(strand.get("+")), "strand_minus": _int(strand.get("-")),
        "rescue_strand_plus": _int(rstrand.get("+")), "rescue_strand_minus": _int(rstrand.get("-")),
        "unusable_no_hit": uh.get(0, 0), "unusable_multi_hit": sum(v for k, v in uh.items() if k >= 3),
        "hits_reads": sum(hits.values()), "trimmed_reads": sum(trimmed.values()),
        "trimmed_bases": sum(k * v for k, v in trimmed.items()),
        "umi_detected": _int(rs["Umi_detected"]) if "Umi_detected" in rs else None,
        "umi_detected_final": _int(rs["Umi_detected_final"]) if "Umi_detected_final" in rs else None,
        "primer_configs": len(configs),
    }
    sc["input"] = sc["pass_reads"] + sc["qc_fail"]
    sc["written"] = max(pf + rescue - sc["len_fail"], 0)

    def same(a, b, *needs) -> bool | None:
        return None if any(n not in cats for n in needs) else a == b
    ident = {
        "classes_add_up": same(pf + rescued_reads + unusable, sc["pass_reads"], "RescueSegmentNr") if rescue else pf + unusable == sc["pass_reads"],
        "rescue_pieces": same(sum(k * v for k, v in seg.items()), rescue, "RescueSegmentNr"),
        "strand_sum": same(sc["strand_plus"] + sc["strand_minus"], pf, "Strand"),
        "rescue_strand_sum": same(sc["rescue_strand_plus"] + sc["rescue_strand_minus"], rescue, "RescueStrand"),
        "unusable_hits_sum": same(sum(uh.values()), unusable, "UnclassHitNr"),
        "rescue_hits_sum": same(sum(rh.values()), rescue, "RescueHitNr"),
        "trimmed_sum": same(sc["trimmed_reads"], pf + rescued_reads, "Unusable", "RescueSegmentNr") if rescue else same(sc["trimmed_reads"], pf, "Unusable"),
        "hits_sum": same(sc["hits_reads"], sc["pass_reads"] - uh.get(0, 0), "Hits", "UnclassHitNr"),
        "cutoff_in_curve": None if cutoff is None or not autotune else any(_is(q, cutoff) for q, _ in autotune),
    }
    # A category that must hold rows (its total is not zero) and is not in the file: a table written by a version of
    # pychopper that did not record it. The numbers read from it cannot be known; they are zero in the cell, and the
    # section says so in words (``absent``: category -> what is not known because of it).
    absent: dict[str, str] = {}
    for cat, total, what in (("RescueSegmentNr", rescue, "the number of rescued reads (so the full-length share leaves them out)"),
                             ("UnclassHitNr", unusable, "the primers found on the unusable reads"),
                             ("Strand", pf, "the strand of the full-length reads"),
                             ("RescueStrand", rescue, "the strand of the rescued pieces"),
                             ("RescueHitNr", rescue, "the primers found on the rescued reads"),
                             ("Hits", sc["pass_reads"], "the primer arrangements"),
                             ("Unusable", pf, "the bases cut off each usable read")):
        if total and cat not in cats:
            absent[cat] = what
    return {"scalars": sc, "absent": absent, "cutoff_q": cutoff, "segments": seg, "unusable_hits": uh, "rescue_hits": rh, "trimmed": trimmed, "configs": configs,
            "autotune": autotune, "parameters": params, "identities": ident, "failed": [k for k, v in ident.items() if v is False],
            "categories": [clean_text(c, 60) for c in t["order"]],
            "unknown_categories": {clean_text(c, 60): len(cats[c]) for c in t["order"] if c not in KNOWN_CATEGORIES},
            "rows": t["rows"], "rows_skipped": t["rows_skipped"] + bad_names, "rows_repeated": t["rows_repeated"]}


# ------------------------------------------------------------------------------------------ aggregates
def _count_hist(axis: str, by_int: dict[int, int], *, weighted: bool = False, w_name: str | None = None) -> dict:
    """{whole number: count} -> a stored histogram; with ``weighted`` the weight of a cell is number x count."""
    full = np.zeros(axes.n_cells(axis), np.int64)
    w = np.zeros(axes.n_cells(axis), np.int64) if weighted else None
    if by_int:
        ks = np.fromiter(by_int.keys(), np.int64, len(by_int))
        vs = np.fromiter(by_int.values(), np.int64, len(by_int))
        idx = axes.bin_index(ks, axis)
        np.add.at(full, idx, vs)
        if w is not None:
            np.add.at(w, idx, ks * vs)
    keys = [k for k, v in by_int.items() if v]
    return build.hist(axis, full, w, w_name=w_name, vmin=min(keys) if keys else None, vmax=max(keys) if keys else None)


def cell_of(stats: dict, config_cats: list[str]) -> dict:
    """The aggregates of one table (a cell without ``at``). ``config_cats`` are the primer configurations stored
    by name for the whole run, so that the cells of a run add up category by category."""
    sc = stats["scalars"]
    by_name = dict(stats["configs"])
    named = [int(by_name.get(c, 0)) for c in config_cats]
    other = sc["hits_reads"] - sum(named)
    scalars = dict(sc)
    scalars["primer_config_other"] = other
    table = build.counts(list(config_cats), named, role="categories")
    table["other"] = other
    cell: dict = {
        "scalars": scalars,
        "counts": {
            "cdna_class_table": build.counts(list(CLASS_CATS), [sc["primers_found"], sc["rescued_reads"], sc["unusable"]], role="categories"),
            "cdna_strand_table": build.counts(["+", "-"], [sc["strand_plus"], sc["strand_minus"]], role="categories"),
            "cdna_rescue_strand_table": build.counts(["+", "-"], [sc["rescue_strand_plus"], sc["rescue_strand_minus"]], role="categories"),
            "cdna_primer_config_table": table,
        },
        "hists": {
            "cdna_rescue_segment_hist": _count_hist("count_1", stats["segments"], weighted=True, w_name="pieces"),
            "cdna_unusable_hits_hist": _count_hist("count_1", stats["unusable_hits"]),
            "cdna_rescue_hits_hist": _count_hist("count_1", stats["rescue_hits"]),
            "cdna_trimmed_bases_hist": _count_hist("trimmed_bases_2", stats["trimmed"], weighted=True, w_name="bases"),
        },
    }
    cut = stats["cutoff_q"]
    points = list(stats["autotune"][:256])
    if cut is not None and not any(_is(q, cut) for q, _ in points):
        points = sorted(points + [(cut, None)], key=lambda qn: qn[0])[:256]       # a cut-off set by hand: no reads were counted at it
    if points:
        cell["curves"] = {"cdna_autotune_curve": {"x": [q for q, _ in points], "y": {"reads": [n for _, n in points],
                                                                                    "chosen": [int(cut is not None and _is(q, cut)) for q, _ in points]},
                                                  "x_unit": "cut-off q", "y_unit": "reads classified in the tuning sample",
                                                  "x_scale": "linear"}}
    return cell


def config_cats(all_stats: list[dict], keep: int = PRIMER_CONFIGS_KEPT) -> list[str]:
    """The ``keep`` most common primer configurations of the tables together, most common first."""
    total: dict[str, int] = {}
    for st in all_stats:
        for name, n in st["configs"]:
            total[name] = total.get(name, 0) + n
    return [k for k, _ in sorted(total.items(), key=lambda kv: (-kv[1], kv[0]))[:keep]]


def _exact_median(by_int: dict[int, int]) -> float | None:
    """The median of whole numbers given as {number: count}: the lower of the two middle values for an even count."""
    n = sum(by_int.values())
    if not n:
        return None
    rank, acc = (n + 1) // 2, 0
    for k in sorted(by_int):
        acc += by_int[k]
        if acc >= rank:
            return k
    return None


def _frac(num, den, digits: int = 6):
    return None if num is None or not den else round(num / den, digits)


def summary_of(pooled: dict, all_stats: list[dict]) -> dict:
    """Derived numbers of the tables together (``pooled`` = the summed scalars)."""
    trimmed: dict[int, int] = {}
    configs: set[str] = set()
    for st in all_stats:
        for k, v in st["trimmed"].items():
            trimmed[k] = trimmed.get(k, 0) + v
        configs.update(name for name, _ in st["configs"])
    cut = {st["cutoff_q"] for st in all_stats}
    strands = (pooled.get("strand_plus") or 0) + (pooled.get("strand_minus") or 0)
    return {
        "cdna_pychopper_samples": len(all_stats),
        "cdna_full_length_frac": _frac(pooled["primers_found"] + pooled["rescued_reads"], pooled["pass_reads"]),
        "cdna_rescued_frac": _frac(pooled["rescued_reads"], pooled["pass_reads"]),
        "cdna_unusable_frac": _frac(pooled["unusable"], pooled["pass_reads"]),
        "cdna_no_primer_frac": _frac(pooled["unusable_no_hit"], pooled["pass_reads"]),
        "cdna_multi_primer_frac": _frac(pooled["unusable_multi_hit"], pooled["pass_reads"]),
        "cdna_strand_plus_frac": _frac(pooled.get("strand_plus"), strands),
        "cdna_qc_fail_frac": _frac(pooled["qc_fail"], pooled["input"]),
        "cdna_len_fail_frac": _frac(pooled["len_fail"], pooled["oriented"]),
        "cdna_rescue_pieces_per_read": _frac(pooled["rescue_segments"], pooled["rescued_reads"], 4),
        "cdna_trimmed_bases_median": _exact_median(trimmed),
        "cdna_trimmed_bases_mean": _frac(pooled["trimmed_bases"], pooled["trimmed_reads"], 2),
        "cdna_umi_detected_frac": _frac(pooled.get("umi_detected"), pooled["oriented"]),
        "cdna_umi_final_frac": _frac(pooled.get("umi_detected_final"), pooled["oriented"]),
        "cdna_primer_configs": len(configs),
        "cdna_pychopper_cutoff": next(iter(cut)) if len(cut) == 1 else None,
    }


PARAMS = {
    "counts_of": {
        "input": "reads given to pychopper = PassReads + QcFail",
        "pass_reads": "reads pychopper classified (PassReads)",
        "qc_fail": "reads dropped before classifying, mean quality under pychopper's limit (QcFail)",
        "len_fail": "pieces dropped after classifying, shorter than pychopper's limit (LenFail)",
        "primers_found": "reads with one full-length piece (Primers_found)",
        "rescued_reads": "fused reads that were cut apart = sum of the RescueSegmentNr rows",
        "rescue_segments": "pieces cut out of fused reads (Rescue)",
        "unusable": "reads without a usable piece (Unusable)",
        "oriented": "pieces kept before the length filter = Primers_found + Rescue",
        "written": "reads pychopper wrote = Primers_found + Rescue - LenFail",
        "strand_plus, strand_minus": "full-length reads by orientation (Strand)",
        "rescue_strand_plus, rescue_strand_minus": "rescued pieces by orientation (RescueStrand)",
        "unusable_no_hit": "unusable reads without any primer hit (UnclassHitNr 0)",
        "unusable_multi_hit": "unusable reads with three primer hits or more (UnclassHitNr 3 and above)",
        "hits_reads": "reads with at least one primer hit = sum of the Hits rows",
        "primer_config_other": "reads with a primer hit whose configuration is not one of those stored by name",
        "primer_configs": "distinct primer configurations in the table; for several samples together: the largest of theirs, not a sum",
        "trimmed_reads": "usable reads = sum of the rows of the category Unusable",
        "trimmed_bases": "bases of usable reads that were not kept = sum of name x value over those rows",
        "umi_detected": "pieces in which pychopper found a UMI (Umi_detected); not known when UMI detection was off",
        "umi_detected_final": "of those, the pieces still counted after pychopper's own checks (Umi_detected_final)",
        "cdna_class_table": "a partition of PassReads: full_length = Primers_found, rescued = rescued reads, unusable = Unusable",
        "cdna_primer_config_table": "reads with a primer hit, by configuration; the same named configurations in every cell",
        "cdna_rescue_segment_hist": "fused reads by the number of pieces cut out of them; the weight is pieces",
        "cdna_unusable_hits_hist": "unusable reads by the number of primer hits on them",
        "cdna_rescue_hits_hist": "rescued PIECES by the number of primer hits on the read they came from",
        "cdna_trimmed_bases_hist": "usable reads by the number of bases not kept; the weight is bases",
        "cdna_autotune_curve": "for each cut-off q that pychopper tried on a sample of the reads: reads classified; chosen = 1 at the cut-off it kept (Parameter Cutoff(q))",
    },
    "full_length_rule": "(Primers_found + rescued reads) / PassReads",
    "written_rule": "reads written = Primers_found + Rescue - LenFail; a fused read that was cut apart counts once per piece",
    "umi_rule": "UMI shares are of the pieces kept before the length filter (Primers_found + Rescue)",
    "median_rule": "the median of the bases not kept is exact (read from the rows of the table, not from the histogram)",
}


def build_section(entries: list[dict], *, order: list[str], expected: list[str], conflicts: list[str] = ()) -> dict | None:
    """The section ``cdna.pychopper``.

    ``entries``: one per sample, ``{"key": sample key, "name": display name, "input": input id, "stats": read_stats(),
    "several": runs the table covers when more than one, "weak": linked by a rule only}``. ``order`` is the order of
    the sample keys in the bundle; ``expected`` the keys of the samples that should have a table."""
    if not entries:
        return None
    pos = {k: i for i, k in enumerate(order)}
    entries = sorted(entries, key=lambda e: (pos.get(e["key"], 1 << 30), e["key"]))
    all_stats = [e["stats"] for e in entries]
    cats = config_cats(all_stats)
    cells = [{"at": [e["key"]], **cell_of(e["stats"], cats)} for e in entries]
    # the tuning curves stay with their samples: a curve is not a sum
    pooled = merge.rollup([{k: v for k, v in c.items() if k != "curves"} for c in cells], ["sample"], [], SCALAR_MERGE)[0]
    run_level = {k: v for k, v in pooled.items() if k in ("scalars", "counts", "hists")}
    run_level["counts"]["cdna_primer_config_table"]["other"] = run_level["scalars"]["primer_config_other"]
    notes = [{"level": "info", "code": "source_is_report",
              "message": "Collected from each sample's pychopper statistics file; every number is the one the file holds, or a sum of them."}]
    for e in entries:
        st = e["stats"]
        for ident in st["failed"]:
            notes.append({"level": "warning", "code": "pychopper_layout", "input": e["input"],
                          "message": clean_text(f"{e['name']}: {IDENTITIES[ident]} does not hold in this table; its numbers are shown as written.", 600)})
        if st.get("absent"):
            notes.append({"level": "warning", "code": "pychopper_categories_absent", "input": e["input"], "value": len(st["absent"]),
                          "message": clean_text(f"{e['name']}: the table holds no rows of {', '.join(sorted(st['absent']))} (written by a version of pychopper "
                                                f"that did not record them). Not known for this sample, and shown as zero here: "
                                                f"{'; '.join(st['absent'][k] for k in sorted(st['absent']))}.", 900)})
        if st["unknown_categories"]:
            names = ", ".join(f"{k} ({v} rows)" for k, v in sorted(st["unknown_categories"].items())[:6])
            notes.append({"level": "info", "code": "pychopper_other_categories", "input": e["input"],
                          "message": clean_text(f"{e['name']}: the table holds categories this version does not read: {names}.", 600)})
        if st["rows_skipped"]:
            notes.append({"level": "warning", "code": "rows_skipped", "input": e["input"], "value": st["rows_skipped"],
                          "message": clean_text(f"{e['name']}: {st['rows_skipped']} rows of the table could not be read and were left out.", 600)})
        if e.get("several"):
            notes.append({"level": "info", "code": "covers_several_runs", "input": e["input"], "value": e["several"],
                          "message": clean_text(f"{e['name']}: one pychopper analysis for {e['several']} flow cells; their reads were merged before pychopper ran.", 600)})
        if e.get("weak"):
            notes.append({"level": "warning", "code": "weak_link", "input": e["input"],
                          "message": clean_text(f"{e['name']}: the table is linked by a rule, not by what it contains.", 600)})
    for m in conflicts:
        notes.append({"level": "warning", "code": "two_candidate_sources", "message": clean_text("Not assessed: " + m, 1200)})
    have = {e["key"] for e in entries}
    missing = [k for k in expected if k not in have]
    state = "ok" if not missing and not conflicts else "partial"
    sec = {"version": VERSION, "tier": "assay", "state": state, "inputs": [e["input"] for e in entries], "dims": ["sample"],
           "cells": cells, "run_level": run_level, "scalar_merge": dict(SCALAR_MERGE),
           "summary": summary_of(run_level["scalars"], all_stats),
           "params": {**PARAMS, "primer_configs_named": len(cats)}, "notes": notes,
           "coverage": {"kind": "full", "records": int(run_level["scalars"]["input"]), "records_total": None, "fraction": 1.0}}
    if state == "partial":
        sec["partial_reason"] = "samples_missing"
    return sec
