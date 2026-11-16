"""The read-data checks of a bundle, in the stored rule form.

Limits come from the catalogue (``thresholds`` with ``when`` clauses on assay group, flow cell family, chemistry and
model tier). Rules:

* a check is evaluated only on a read population it applies to; on another one it is stored as not assessed with
  the reason, never judged on reads it was not written for;
* per-sample checks are evaluated on every sample and stored as one rule: the worst sample decides, the sentence
  names the samples;
* a check is not assessed while the files behind it are only partly counted, when they are linked by an
  unconfirmed rule, or when two candidate sources describe the same reads;
* read fate is assessed only when every stage it needs has a source;
* a library may change the limits of a check for one assay group, as it can for run health:
  ``[profiles.<assay group>.rules.<check id>]`` with ``warn``, ``fail`` or ``enabled = false`` in ``flowcellar.toml``.

Nothing here reads a file or the library: the bundle and the ``[profiles]`` table handed over are the only inputs.
"""
from __future__ import annotations

from flowcellar import catalogue
from flowcellar.aggregates import merge
from flowcellar.readdata import copy_pending as copy
from flowcellar.readdata import sections as S
from flowcellar.util.text import fmt_bases, fmt_pct

UNCLASSIFIED = "unclassified"
SWITCHED_OFF = "this check is switched off in the library's flowcellar.toml"


# ------------------------------------------------------------------------------------------ wording
def _fmt(metric: str, v) -> str:
    kind = copy.CHECK_COPY[metric]["fmt"]
    if v is None:
        return "not known"
    if kind == "pct":
        return fmt_pct(v, 1)
    if kind == "pct2":
        return fmt_pct(v, 2)
    if kind == "q":
        return f"Q {v:.1f}"
    if kind == "length":
        return f"{int(round(v)):,} bases"
    if kind == "bases":
        return fmt_bases(v)
    return str(v)


def _limit_sentence(metric: str, kind: str, warn) -> str:
    lt = catalogue.limit_text()
    return lt[kind].replace("{limit}", _fmt(metric, warn))


def _context(bundle: dict) -> dict:
    idn = bundle.get("identity") or {}
    assay = idn.get("assay") or {}
    bc = [b for b in idn.get("basecalling") or [] if b.get("primary")] or list(idn.get("basecalling") or [])
    model = (bc[0].get("model") if bc else None) or ""
    chem = str(assay.get("chemistry") or "").lower()
    if not chem:
        chem = "rna004" if model.startswith("rna004") else "r10" if "r10" in model else "r9" if "r9" in model else ""
    chem = "rna004" if "rna004" in chem else "r10" if chem.startswith("r10") else "r9" if chem.startswith("r9") else (chem or None)
    return {"assay": assay.get("group"), "family": (idn.get("flow_cell") or {}).get("family"), "chemistry": chem,
            "model_tier": bc[0].get("model_tier") if bc else None,
            "both_ends": any(b.get("barcode_both_ends") for b in idn.get("basecalling") or [])}


def _override(metric_id: str, ctx: dict) -> dict:
    """What the library's configuration says about one check for the assay group of this run."""
    ov = (ctx.get("overrides") or {}).get(metric_id)
    return ov if isinstance(ov, dict) else {}


def _limits(metric_id: str, ctx: dict) -> tuple[str, float | None, float | None, str] | None:
    """(kind, warn, fail, basis): the first ``when`` clause of the catalogue that matches, with the levels the
    library's configuration sets instead; None when neither gives a limit."""
    m = catalogue.metric(metric_id) or {}
    th = m.get("thresholds") or {}
    if th.get("kind") not in ("min", "max"):
        return None
    ov = _override(metric_id, ctx)
    numbers = {k: ov[k] for k in ("warn", "fail") if isinstance(ov.get(k), (int, float)) and not isinstance(ov.get(k), bool)}
    found = None
    for r in th.get("rules") or []:
        when = r.get("when") or {}
        if all(ctx.get(k) is not None and ctx.get(k) == v for k, v in when.items()):
            words = [copy.ASSAY_WORDS.get(when["assay"], when["assay"]) + " runs" if "assay" in when else None,
                     when.get("family"), {"r10": "R10", "r9": "R9", "rna004": "RNA004"}.get(when.get("chemistry")),
                     f"{when['model_tier']} basecalling" if "model_tier" in when else None]
            basis = "built-in limit" + (" for " + ", ".join(w for w in words if w) if any(words) else "")
            found = (th["kind"], r.get("warn"), r.get("fail"), basis)
            break
    if not numbers:
        return found
    kind, warn, fail, _ = found or (th["kind"], None, None, "")
    return kind, numbers.get("warn", warn), numbers.get("fail", fail), "limit set in this library's flowcellar.toml"


def _level(kind: str, v: float, warn, fail) -> str:
    if kind == "min":
        return "fail" if fail is not None and v < fail else "warn" if warn is not None and v < warn else "pass"
    return "fail" if fail is not None and v > fail else "warn" if warn is not None and v > warn else "pass"


def _rule(check_id: str, metric: str, module: str, *, outcome: str, level: str = "pass", value=None, threshold=None, message: str,
          reason: str | None = None, population=None, sample=None) -> dict:
    m = catalogue.metric(check_id) or {}
    role = m.get("verdict_role") if m.get("verdict_role") in ("key", "supporting") else "info"
    sev = "none"
    if outcome == "triggered":
        sev = "fail" if (level == "fail" and role == "key") else "warn"
    r = {"id": check_id, "role": role, "tier": "read_level", "module": module, "severity": sev, "outcome": outcome,
         "scope": {"sample": sample}, "metric": metric, "value": value, "threshold": threshold or {"kind": "rule"},
         "message": message, "not_evaluated_reason": reason, "panel": None, "population": population}
    return r


def _not_assessed(check_id: str, metric: str, module: str, reason: str, population=None) -> dict:
    msg = catalogue.limit_text()["not_assessed"].replace("{reason}", reason)
    return _rule(check_id, metric, module, outcome="not_evaluated", message=msg, reason=reason, population=population,
                 threshold={"kind": "rule", "basis": reason})


def _blocked(sec: dict) -> str | None:
    """Why no check may be assessed on a section, or None."""
    codes = {n.get("code") for n in sec.get("notes") or []}
    if "weak_link" in codes:
        return copy.REASONS["weak_link"]
    if "mixed_sources" in codes:
        return copy.REASONS["mixed_sources"]
    if sec.get("state") == "partial" and sec.get("partial_reason") in ("partly_read", "waiting"):
        return copy.REASONS["not_counted_yet"]
    if not sec.get("cells"):
        return copy.REASONS["two_sources"] if "two_candidate_sources" in codes else copy.REASONS["no_value"]
    return None


def _population_ok(check_id: str, pop: dict | None) -> str | None:
    wanted = (catalogue.metric(check_id) or {}).get("populations")
    if not wanted or wanted == "any" or pop is None:
        return None
    if pop.get("stage") in wanted:
        return None
    return copy.REASONS["population"].format(wanted=" or ".join(copy.POPULATION_WORDS.get(w, w) for w in wanted),
                                             have=copy.POPULATION_WORDS.get(pop.get("stage"), "another set of reads"))


def _per_sample(bundle: dict, sec: dict, value_of) -> tuple[dict[str, float], float | None]:
    """({sample key: value}, value over all samples) from the cells of a section with dims sample x filter."""
    cells = sec.get("cells") or []
    dims = sec.get("dims") or []
    out: dict[str, float] = {}
    if dims == ["sample", "filter"]:
        for key in dict.fromkeys(c["at"][0] for c in cells):
            if key == UNCLASSIFIED:
                continue
            mine = [c for c in cells if c["at"][0] == key]
            tot = S._total(mine, pass_only=True)
            v = value_of(tot)
            if v is not None:
                out[key] = v
        named = [c for c in cells if c["at"][0] != UNCLASSIFIED]
        allv = value_of(S._total(named, pass_only=True)) if named else None
    else:
        for c in cells:
            if c["at"][0] == UNCLASSIFIED:
                continue
            v = value_of(c)
            if v is not None:
                out[c["at"][0]] = v
        rolled = merge.rollup([c for c in cells if c["at"][0] != UNCLASSIFIED], dims, []) if cells else []
        allv = value_of(rolled[0]) if rolled else None
    return out, allv


def _names(bundle: dict, keys: list[str]) -> str:
    by = {s.get("key"): (s.get("name") or s.get("alias") or s.get("key")) for s in bundle.get("samples") or []}
    names = [str(by.get(k, k)) for k in keys]
    if len(names) <= 3:
        return ", ".join(names[:-1]) + (" and " if len(names) > 1 else "") + names[-1]
    return f"{names[0]}, {names[1]} and {len(names) - 2} more"


def _sample_check(bundle: dict, ctx: dict, check_id: str, metric: str, module: str, value_of, *, extra_block=None) -> dict | None:
    sec = (bundle.get("modules") or {}).get(module)
    if not sec:
        return None
    pop = sec.get("population")
    if _override(check_id, ctx).get("enabled") is False:
        return _not_assessed(check_id, metric, module, SWITCHED_OFF, pop)
    why = _blocked(sec) or _population_ok(check_id, pop) or (extra_block(sec) if extra_block else None)
    if why:
        return _not_assessed(check_id, metric, module, why, pop)
    lim = _limits(check_id, ctx)
    if lim is None:
        return _not_assessed(check_id, metric, module, copy.REASONS["no_limit"], pop)
    kind, warn, fail, basis = lim
    values, overall = _per_sample(bundle, sec, value_of)
    if not values:
        return _not_assessed(check_id, metric, module, copy.REASONS["no_value"], pop)
    levels = {k: _level(kind, v, warn, fail) for k, v in values.items()}
    worst_key = (min if kind == "min" else max)(values, key=lambda k: values[k])
    worst = values[worst_key]
    off = [k for k in values if levels[k] != "pass"]
    level = "fail" if "fail" in levels.values() else "warn" if "warn" in levels.values() else "pass"
    label = copy.CHECK_COPY[metric]["label"]
    lim_s = _limit_sentence(metric, kind, warn)
    extreme = "lowest" if kind == "min" else "highest"
    if not off:
        if len(values) == 1:
            msg = copy.SENTENCES["one"].format(label=label, value=_fmt(metric, worst))
        else:
            msg = copy.SENTENCES["all_samples"].format(label=label, value=_fmt(metric, overall if overall is not None else worst),
                                                       extreme=extreme, worst=_fmt(metric, worst))
    elif len(off) == 1 and len(values) == 1:
        msg = copy.SENTENCES["one"].format(label=label, value=_fmt(metric, worst))
    elif len(off) == 1:
        msg = copy.SENTENCES["one_sample_off"].format(label=label, worst=_fmt(metric, worst), names=_names(bundle, off))
    else:
        order = sorted(off, key=lambda k: values[k], reverse=(kind == "max"))
        msg = copy.SENTENCES["samples_off"].format(label=label, side="below" if kind == "min" else "above", k=len(off), n=len(values),
                                                   names=_names(bundle, order), extreme=extreme, worst=_fmt(metric, worst))
    return _rule(check_id, metric, module, outcome="triggered" if off else "ok", level=level, value=round(float(worst), 6),
                 threshold={"kind": kind, "warn": warn, "fail": fail, "basis": basis + f"; evaluated on each of {len(values)} sample{'s' if len(values) != 1 else ''}", "applied": "fixed"},
                 message=msg + " " + lim_s, population=pop, sample=worst_key if len(off) == 1 else None)


def _run_check(check_id: str, metric: str, module: str, ctx: dict, value, pop, *, sentence: str | None = None, block: str | None = None) -> dict:
    if _override(check_id, ctx).get("enabled") is False:
        return _not_assessed(check_id, metric, module, SWITCHED_OFF, pop)
    if block:
        return _not_assessed(check_id, metric, module, block, pop)
    lim = _limits(check_id, ctx)
    if lim is None:
        return _not_assessed(check_id, metric, module, copy.REASONS["no_limit"], pop)
    if value is None:
        return _not_assessed(check_id, metric, module, copy.REASONS["no_value"], pop)
    kind, warn, fail, basis = lim
    level = _level(kind, value, warn, fail)
    msg = (sentence or copy.SENTENCES["one"].format(label=copy.CHECK_COPY[metric]["label"], value=_fmt(metric, value)))
    return _rule(check_id, metric, module, outcome="ok" if level == "pass" else "triggered", level=level, value=round(float(value), 6),
                 threshold={"kind": kind, "warn": warn, "fail": fail, "basis": basis, "applied": "fixed"},
                 message=msg + " " + _limit_sentence(metric, kind, warn), population=pop)


def read_checks(bundle: dict, profiles: dict | None = None) -> list[dict]:
    """``profiles`` is the ``[profiles]`` table of the library's configuration (assay group -> ``rules`` -> check id)."""
    mods = bundle.get("modules") or {}
    ctx = _context(bundle)
    group = (profiles or {}).get(ctx.get("assay") or "")
    ctx["overrides"] = (group.get("rules") if isinstance(group, dict) else None) or {}
    rules: list[dict] = []

    def add(r):
        if r is not None:
            rules.append(r)

    # ---- basecalled output
    yb = mods.get("yield.basecalled")
    if yb:
        pop = yb.get("population")
        blk = _blocked(yb)
        add(_run_check("yield_pass_frac", "yield_pass_frac", "yield.basecalled", ctx, (yb.get("summary") or {}).get("yield_pass_frac"), pop,
                       block=blk or _population_ok("yield_pass_frac", pop)))
        cells = yb.get("cells") or []
        passed = sum((c.get("scalars") or {}).get("bases", 0) for c in cells if c["at"][1] == "pass") if cells else None
        run = (mods.get("run") or {}).get("summary") or {}
        early = None
        done = run.get("run_completion_frac")
        hours = run.get("run_duration_h")
        if not ((done is not None and done >= 0.9) or (hours is not None and hours >= 48)):
            early = copy.REASONS["stopped_early"]
        add(_run_check("yield_basecalled_bases", "yield_basecalled_bases", "yield.basecalled", ctx, passed, pop,
                       block=blk or _population_ok("yield_basecalled_bases", pop) or early))
    # ---- length and quality, per sample
    add(_sample_check(bundle, ctx, "len_n50", "len_n50", "len", lambda c: S.len_summary(c)["len_n50"]))
    add(_sample_check(bundle, ctx, "len_median", "len_median", "len", lambda c: S.len_summary(c)["len_median"]))
    add(_sample_check(bundle, ctx, "len_short_read_frac", "len_short_read_frac", "len", lambda c: S.len_summary(c)["len_short_read_frac"]))
    add(_sample_check(bundle, ctx, "qual_median", "qual_median", "qual", S.qual_median))
    # ---- alignment, per sample
    mapped_only = lambda sec: copy.REASONS["mapped_only"] if (sec.get("population") or {}).get("stage") == "aligned_primary" else None  # noqa: E731
    add(_sample_check(bundle, ctx, "aln_mapped_frac", "aln_mapped_frac", "aln", lambda c: S.aln_summary(c)["aln_mapped_frac"], extra_block=mapped_only))
    add(_sample_check(bundle, ctx, "aln_identity_hist", "aln_identity_median", "aln", lambda c: S.aln_summary(c)["aln_identity_median"]))
    add(_sample_check(bundle, ctx, "aln_supplementary_frac", "aln_supplementary_frac", "aln", lambda c: S.aln_summary(c)["aln_supplementary_frac"]))
    # ---- cDNA library, per sample
    cd = mods.get("cdna")
    if cd:
        def full(c):
            sc = (c or {}).get("scalars") or {}
            return (sc.get("primers_found", 0) + sc.get("rescued_reads", 0)) / sc["pass_reads"] if sc.get("pass_reads") else None

        def rescued(c):
            sc = (c or {}).get("scalars") or {}
            return sc.get("rescued_reads", 0) / sc["pass_reads"] if sc.get("pass_reads") else None
        both = None
        if ctx["both_ends"] or ((mods.get("len") or {}).get("population") or {}).get("both_end_barcodes"):
            v = (cd.get("summary") or {}).get("cdna_full_length_frac")
            both = lambda sec: copy.REASONS["both_ends"].format(value=fmt_pct(v, 1) if v is not None else "not known")  # noqa: E731
        add(_sample_check(bundle, ctx, "cdna_full_length_frac", "cdna_full_length_frac", "cdna", full, extra_block=both))
        add(_sample_check(bundle, ctx, "cdna_rescued_frac", "cdna_rescued_frac", "cdna", rescued))
    # ---- read fate
    fu = mods.get("funnel")
    if fu:
        sm = fu.get("summary") or {}
        why = (fu.get("params") or {}).get("usable_not_assessed_reason")
        stage = copy.STAGE_WORDS.get(sm.get("funnel_usable_stage"), "usable")
        v = sm.get("funnel_usable_frac")
        add(_run_check("funnel_usable_frac", "funnel_usable_frac", "funnel", ctx, v, None, block=why if v is None else None,
                       sentence=copy.SENTENCES["usable"].format(stage=stage, value=fmt_pct(v, 1)) if v is not None else None))
        loss = sm.get("funnel_basecall_loss_frac")
        if loss is not None:
            add(_run_check("funnel_basecall_loss_frac", "funnel_basecall_loss_frac", "funnel", ctx, loss, None))
    return rules
