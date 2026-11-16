"""Read fate: from the reads the flow cell sequenced to the reads counted towards a transcript.

Stages, in the order they are stored: ``sequenced``, ``basecalled``, ``classified`` (given to a sample),
``full_length`` (kept by the cDNA step), ``genome_mapped``, ``quantified`` (counted towards a transcript).

Two rules hold for every bundle:

* **A stage without a source is not known, never zero.** It is listed in ``missing`` with the reason in
  ``sources``; its ``n`` is 0 only because the array needs a number.
* **Stages never grow.** Every stage is a subset of the stages it descends from: ``basecalled`` of ``sequenced``,
  ``classified`` of ``basecalled``, ``full_length`` of ``classified``, and both ``genome_mapped`` and ``quantified``
  of ``full_length`` (those two are not subsets of each other: a read can be counted towards a transcript without
  a genome alignment). A count that is larger than a stage it descends from means the two sources do not describe
  the same reads; it is stored as not known, with both numbers in the reason, and a warning is raised.

Units. ``sequenced`` and ``basecalled`` count reads, a read the basecaller split counting once. From ``full_length``
on the unit is the read as the cDNA step wrote it: a fused read that was cut apart counts once per piece, because
that is what every later file holds.
"""
from __future__ import annotations

import re

from flowcellar.aggregates import build as agg
from flowcellar.readdata import tables

VERSION = "0.1.1"                    # 0.1.1: a table of reads merged from several flow cells is said as that (wording of sources and notes)
STAGES = ("sequenced", "basecalled", "classified", "full_length", "genome_mapped", "quantified")
SAMPLE_STAGES = ("classified", "full_length", "genome_mapped", "quantified")
#: the stages a stage descends from, nearest first
ANCESTORS = {"sequenced": (), "basecalled": ("sequenced",), "classified": ("basecalled", "sequenced"),
             "full_length": ("classified", "basecalled", "sequenced"),
             "genome_mapped": ("full_length", "classified", "basecalled", "sequenced"),
             "quantified": ("full_length", "classified", "basecalled", "sequenced")}
WORDS = {"sequenced": "sequenced", "basecalled": "basecalled", "classified": "given to a sample", "full_length": "full-length",
         "genome_mapped": "mapped to the genome", "quantified": "counted towards a transcript"}
UNCLASSIFIED = "unclassified"


def sequenced_reads(tier1: dict | None) -> int | None:
    """Reads the flow cell produced, from the run-health part of the bundle."""
    y = ((tier1 or {}).get("modules") or {}).get("yield") or {}
    v = (y.get("summary") or {}).get("yield_reads")
    if v is None:
        for c in y.get("cells") or []:
            v = (c.get("scalars") or {}).get("reads")
            if v is not None:
                break
    return int(v) if v is not None else None


def never_grow(counts: dict[str, int | None], stages=STAGES) -> list[tuple[str, int, str, int]]:
    """The violations of "stages never grow" among known counts: (stage, its count, ancestor, the ancestor's count)."""
    bad = []
    for s in stages:
        v = counts.get(s)
        if v is None:
            continue
        for a in ANCESTORS[s]:
            if counts.get(a) is not None:
                if v > counts[a]:
                    bad.append((s, v, a, counts[a]))
                break
    return bad


def _enforce(counts: dict, sources: dict, stages, who: str, warnings: list[dict]) -> None:
    """Store a stage that grew as not known, with both numbers in its reason. Repeats until nothing grows."""
    for _ in range(len(stages)):
        bad = never_grow(counts, stages)
        if not bad:
            return
        s, v, a, av = bad[0]
        was = sources.get(s, "")
        counts[s] = None
        sources[s] = (f"not known: the count from its source ({v:,}; {was}) is larger than the reads {WORDS[a]} ({av:,}), "
                      "so the two files do not describe the same reads")
        warnings.append({"level": "warning", "code": "funnel_stage_grew",
                         "message": f"Read fate{who}: {v:,} reads {WORDS[s]} against {av:,} {WORDS[a]}. A later stage cannot hold more reads "
                                    "than an earlier one, so the later count is shown as not known. Two files that describe the same reads may both be linked."})


def _once_per_parent(sc: dict) -> tuple[int, bool]:
    """Reads of a cell with a split read counted once: (count, exact)."""
    reads = int(sc.get("reads", 0))
    children, parents = sc.get("split_children"), sc.get("split_parents")
    if children and parents is not None:
        return reads - (int(children) - int(parents)), True
    return reads, not children


def build_stage_table(stages, counts: dict, sources: dict) -> dict:
    missing = [s for s in stages if counts.get(s) is None]
    return agg.counts(list(stages), [int(counts[s]) if counts.get(s) is not None else 0 for s in stages], role="funnel",
                        missing=missing, sources={s: sources.get(s, "not known: no source file for this stage") for s in stages})


def build(*, run, tier1, samples, len_cells, aln_cells, reads_pop, aln_pop, collected, assay_group, barcoded, reads_inputs,
          genome_inputs, reads_complete, genome_complete, conflicts, reads_waiting=False, genome_waiting=False,
          reads_roles=(), sample_state=None) -> dict:
    """The sections ``funnel`` and ``funnel.samples``, the warnings they raise and which input fed which.

    ``sample_state``: ``{need: {"read": set of sample keys whose files are read whole, "conflict": set of sample
    keys with two candidate sources (None in it: the conflict concerns every sample)}}``. With it "not known" is
    decided sample by sample: a conflict or an unread file of one sample does not blank the others."""
    warnings: list[dict] = []
    sample_state = sample_state or {}

    def whole(need: str, k: str, run_wide: bool) -> bool:
        st = sample_state.get(need)
        return run_wide if st is None else (k in st["read"] and not in_conflict(need, k))

    def in_conflict(need: str, k: str) -> bool:
        st = sample_state.get(need)
        if st is None:
            return bool(conflicts.get(need))
        return k in st["conflict"] or None in st["conflict"]
    why: dict[str, dict[str, str]] = {}                 # every sample's reasons, also of a sample with no number at all
    used: dict[str, set] = {}
    pych, quant = collected["pychopper"], collected["quant"]
    stage_pop = (reads_pop or {}).get("stage")
    keys = [k for k in samples.order() if k != UNCLASSIFIED]
    expected = samples.expected_keys()

    def cell_reads(cells, key, filters=("pass",)):
        sel = [c for c in cells if c["at"][0] == key and (filters is None or c["at"][1] in filters)]
        return sel

    # ---- per sample
    per: dict[str, dict] = {}
    for k in keys:
        counts: dict[str, int | None] = {s: None for s in SAMPLE_STAGES}
        srcs: dict[str, str] = {}
        ids: set[str] = set()
        py = pych.get(k)
        py_ok = bool(py) and not py["several"] and not py["weak"] and py["consistent"]
        lc = cell_reads(len_cells, k)
        ac = cell_reads(aln_cells, k, None)
        # classified
        if lc and stage_pop in ("classified", "basecalled_pass", "basecalled_all") and whole("reads", k, reads_complete):
            n, exact = 0, True
            for c in lc:
                v, ex = _once_per_parent(c.get("scalars") or {})
                n += v
                exact = exact and ex
            counts["classified"] = n
            srcs["classified"] = "reads of this sample in the basecalled files, every record read" + ("" if exact else "; split reads could not be counted once")
            ids.update(reads_inputs)
        elif py_ok:
            counts["classified"] = py["scalars"]["input"]
            srcs["classified"] = "pychopper table: reads given to pychopper (PassReads + QcFail)"
            ids.add(py["input"])
        elif in_conflict("reads", k):
            srcs["classified"] = "not assessed: two candidate sources"
        elif reads_waiting and not all(r.startswith("align.") for r in reads_roles):
            srcs["classified"] = "not known yet: the basecalled files have not been read"
        # full-length
        fl_bam = None
        for cells, pop in ((aln_cells, aln_pop), (len_cells, reads_pop)):
            sel = cell_reads(cells, k, None)
            if sel and (pop or {}).get("stage") == "full_length_trimmed":
                unm = sum((c.get("scalars") or {}).get("unmapped_reads", 0) for c in cell_reads(aln_cells, k, None))
                if unm > 0 and (whole("genome", k, genome_complete) if cells is aln_cells else whole("reads", k, reads_complete)):
                    fl_bam = sum((c.get("scalars") or {}).get("reads", 0) for c in sel)
                    break
        if py_ok:
            counts["full_length"] = py["scalars"]["written"]
            srcs["full_length"] = "pychopper table: reads written (Primers_found + Rescue - LenFail)"
            ids.add(py["input"])
            if fl_bam is not None and fl_bam != py["scalars"]["written"]:
                warnings.append({"level": "info", "code": "full_length_counts_differ",
                                 "message": f"{samples.name(k)}: the pychopper table says {py['scalars']['written']:,} reads were written, the alignment file holds {fl_bam:,} reads."})
        elif fl_bam is not None:
            counts["full_length"] = fl_bam
            srcs["full_length"] = "reads of this run in the alignment file of full-length reads (it kept its unmapped reads), every record read"
            ids.update(genome_inputs or reads_inputs)
        elif py and py["several"]:
            srcs["full_length"] = ("not known for this flow cell alone: the reads of this sample from several flow cells were "
                                   "merged before pychopper ran, so there is one pychopper table for all of them")
        # mapped to the genome
        if ac and whole("genome", k, genome_complete):
            counts["genome_mapped"] = sum((c.get("scalars") or {}).get("mapped_reads", 0) for c in ac)
            srcs["genome_mapped"] = "primary mapped records of the genome alignment file, every record read"
            ids.update(genome_inputs)
        elif in_conflict("genome", k):
            srcs["genome_mapped"] = "not assessed: two candidate sources"
        elif ac or genome_waiting:
            srcs["genome_mapped"] = "not known yet: the genome alignment file is not read whole"
        # counted towards a transcript
        q = quant.get(k) or {}
        for role in tables.QUANT_PREFERENCE:
            d = q.get(role)
            if not d or d["weak"]:
                continue
            if d["several"]:
                srcs["quantified"] = (f"not known for this flow cell alone: the reads of this sample from several flow cells were "
                                      f"merged before {d['tool']} ran, so there is one {d['tool']} table for all of them")
                break
            counts["quantified"] = d["counted"]
            ver = f" {d['version']}" if d.get("version") else ""
            proc = f" ({d['processed']:,} processed)" if d.get("processed") is not None else ""
            srcs["quantified"] = f"{d['tool']}{ver}: {d['counted_field']}{proc}"
            ids.add(d["input"])
            break
        why[k] = dict(srcs)
        if all(v is None for v in counts.values()):
            continue
        _enforce(counts, srcs, SAMPLE_STAGES, f" of {samples.name(k)}", warnings)
        for s in SAMPLE_STAGES:
            if counts[s] is None:
                srcs.setdefault(s, "not known: no source file for this stage on this machine")
        per[k] = {"counts": counts, "sources": srcs, "inputs": ids}

    # ---- the flow cell
    seq = sequenced_reads(tier1)
    counts = {s: None for s in STAGES}
    srcs: dict[str, str] = {}
    counts["sequenced"] = seq
    srcs["sequenced"] = "MinKNOW: reads counted during the run" if seq is not None else "not known: no run folder of this run is in the library"
    all_len = [c for c in len_cells if c["at"][1] == "pass"]
    if all_len and stage_pop in ("basecalled_all", "basecalled_pass") and reads_complete:
        n, exact = 0, True
        for c in all_len:
            v, ex = _once_per_parent(c.get("scalars") or {})
            n += v
            exact = exact and ex
        counts["basecalled"] = n
        minq = (reads_pop or {}).get("min_qscore")
        has_fail = any(c["at"][1] == "fail" for c in len_cells)
        srcs["basecalled"] = ("reads the basecaller kept" + (f" (Q {minq:g} or better)" if minq is not None else " (passed its filter)" if has_fail else "")
                              + ", every record of the basecalled files read; a split read counts once" + ("" if exact else " where its parent is known"))
    elif conflicts.get("reads"):
        srcs["basecalled"] = "not assessed: two candidate sources"
    elif stage_pop == "classified":
        srcs["basecalled"] = "not known: the linked files hold only the reads given to samples; the reads without a barcode are not on this machine"
    elif stage_pop in ("full_length_trimmed", "aligned_primary"):
        srcs["basecalled"] = "not known: no basecalled file is linked; the reads on this machine were already filtered and trimmed"
    elif reads_roles and all(r.startswith("align.") for r in reads_roles):
        srcs["basecalled"] = "not known: no basecalled file is linked; the reads on this machine come from alignment files"
    elif (len_cells or reads_waiting) and not reads_complete:
        srcs["basecalled"] = "not known yet: the basecalled files are not read whole"
    else:
        srcs["basecalled"] = "not known: no basecalled file is linked to this run"
    want = [k for k in expected if k != UNCLASSIFIED] or list(per)
    for s in SAMPLE_STAGES:
        have = [k for k in want if k in per and per[k]["counts"].get(s) is not None]
        if want and len(have) == len(want):
            counts[s] = sum(per[k]["counts"][s] for k in want)
            srcs[s] = (f"sum over the {len(want)} samples" if len(want) > 1 else "the run's one sample") + ": " + _common(per, want, s)
        elif have:
            srcs[s] = f"not known for the whole flow cell: {len(have)} of {len(want)} samples have a source"
        else:
            # the samples' own reasons, also of samples that have no number at all (a table shared by two runs)
            reasons = {(per[k]["sources"] if k in per else why.get(k, {})).get(s) for k in want} - {None}
            needs_of = {"classified": ("reads",), "full_length": ("pychopper.stats",), "genome_mapped": ("genome",),
                        "quantified": tables.QUANT_PREFERENCE}[s]
            if any(conflicts.get(n) for n in needs_of):
                srcs[s] = "not assessed: two candidate sources"
            elif (s == "genome_mapped" and genome_waiting) or (s == "classified" and reads_waiting and not reasons):
                srcs[s] = "not known yet: the files are not read whole"
            else:
                srcs[s] = reasons.pop() if len(reasons) == 1 else "not known: no source file for this stage on this machine"
    stages = [s for s in STAGES if _applies(s, assay_group, barcoded, counts)]
    _enforce(counts, srcs, stages, "", warnings)
    usable_stage = "quantified" if (assay_group in ("cdna", "drna") or (assay_group is None and counts.get("quantified") is not None)) else "genome_mapped"
    usable = None
    reason = None
    if seq is None:
        reason = "the reads sequenced are not known (no run folder of this run is in the library)"
    elif counts.get(usable_stage) is None:
        reason = f"the reads {WORDS[usable_stage]} are not known for every sample of this flow cell"
    else:
        usable = round(counts[usable_stage] / seq, 6) if seq else None
    loss = None
    if seq and counts.get("basecalled") is not None:
        loss = round(1.0 - counts["basecalled"] / seq, 6)
    used_ids: set[str] = set()
    for k in per:
        used_ids |= per[k]["inputs"]
    if counts.get("basecalled") is not None:
        used_ids |= set(reads_inputs)
    table = build_stage_table(stages, counts, srcs)
    missing = table.get("missing") or []
    sec = {"version": VERSION, "tier": "read_level", "state": "partial" if missing else "ok", "inputs": sorted(used_ids, key=_idnum), "dims": [],
           "cells": [{"at": [], "counts": {"funnel_stage_table": table}}],
           "summary": {"funnel_reads_sequenced": seq, "funnel_usable_frac": usable, "funnel_usable_stage": usable_stage,
                       "funnel_basecall_loss_frac": loss, "funnel_reads_basecalled": counts.get("basecalled")},
           "params": {"usable_not_assessed_reason": reason,
                      "branches": "mapped to the genome and counted towards a transcript both descend from the full-length reads; "
                                  "neither is a part of the other, so either may be the larger",
                      "units": "sequenced and basecalled count reads, a split read once; from full-length on a fused read that was cut apart counts once per piece",
                      "never_grow": "a stage larger than a stage it descends from is stored as not known"},
           "notes": []}
    if missing:
        sec["partial_reason"] = "inputs_missing"
    sections = {"funnel": sec}
    for i in used_ids:
        used.setdefault(i, set()).add("funnel")
    if per:
        pos = {k: i for i, k in enumerate(samples.order())}
        cells = []
        ids: set[str] = set()
        for k in sorted(per, key=lambda k: pos.get(k, 1 << 30)):
            st = [s for s in SAMPLE_STAGES if _applies(s, assay_group, barcoded, per[k]["counts"], sample=True)]
            cells.append({"at": [k], "counts": {"funnel_stage_table": build_stage_table(st, per[k]["counts"], per[k]["sources"])}})
            ids |= per[k]["inputs"]
        any_missing = any(c["counts"]["funnel_stage_table"].get("missing") for c in cells) or any(k not in per for k in want)
        s2 = {"version": VERSION, "tier": "read_level", "state": "partial" if any_missing else "ok", "inputs": sorted(ids, key=_idnum), "dims": ["sample"],
              "cells": cells, "params": {"units": sec["params"]["units"]},
              "notes": [{"level": "info", "code": "records_and_reads",
                         "message": "Stages count reads: primary records only. Secondary and supplementary alignment records are never counted as reads."}]}
        if any_missing:
            s2["partial_reason"] = "inputs_missing" if all(k in per for k in want) else "samples_missing"
        sections["funnel.samples"] = s2
        for i in ids:
            used.setdefault(i, set()).add("funnel.samples")
    return {"sections": sections, "warnings": warnings, "used": used}


def _idnum(i: str) -> tuple:
    d = "".join(ch for ch in i if ch.isdigit())
    return (int(d) if d else 0, i)


def _common(per: dict, want: list[str], stage: str) -> str:
    """The source of a stage in words when every sample has the same one (a per-sample count in brackets left out)."""
    texts = {re.sub(r"\s*\([\d,]+ processed\)", "", per[k]["sources"].get(stage, "")) for k in want}
    return texts.pop() if len(texts) == 1 else "sources differ between samples, see the per-sample table"


def _applies(stage: str, assay_group: str | None, barcoded: bool, counts: dict, sample: bool = False) -> bool:
    """Whether a stage belongs to this kind of run. A stage that has a number is always shown."""
    if counts.get(stage) is not None or stage in ("sequenced", "basecalled", "genome_mapped"):
        return True
    if stage == "classified":
        return barcoded or sample
    if stage == "full_length":
        return assay_group == "cdna"
    if stage == "quantified":
        return assay_group in ("cdna", "drna")
    return True
