"""The verdict of an analysis: its checks, their limits and the sentences stored with the result.

An analysis has its own verdict. It is about the tissue, the library and the workflow, and it **never changes
the verdict of a flow cell run**: a good flow cell with badly permeabilised tissue is still a good flow cell.

Three kinds of limit, in this order of trust:

1. the pipeline's own alarms, shown word for word (WARN -> warning, ERROR -> fail);
2. vendor conventions, used where an alarm cannot be read;
3. provisional defaults for the numbers no vendor publishes a limit for.

The numbers are the catalogue's (``flowcellar.sc.catalogue``); a library overrides them per assay in
``flowcellar.toml`` under ``[profiles.<assay>.rules.<check id>]`` with ``warn``, ``fail`` or ``enabled``.
A check whose number has no source is "not assessed", never a failure. Pure: no file is read or written.
"""
from __future__ import annotations

from flowcellar.model import status as status_model
from flowcellar.sc import MODULE_ID, PIPELINE_WORD
from flowcellar.sc import catalogue as cat
from flowcellar.util.timeutil import now_z

NOT_ABOUT_RUNS = "This verdict is about the analysis. It does not change the verdict of any flow cell."
PROFILE_ID = "sc-default"
_COUNT_KEYS = frozenset(("spatial_units_under_tissue", "barcodes_cells", "workflow_failed_tasks"))


def _pct(x) -> str:
    v = float(x) * 100.0
    return (f"{v:.0f}" if abs(v - round(v)) < 0.05 and abs(v) >= 10 else f"{v:.1f}") + " %"


def _fmt(key: str | None, x) -> str:
    if x is None:
        return "–"
    if key in _COUNT_KEYS or (isinstance(x, int) and not isinstance(x, bool)):
        return f"{int(x):,}"
    if key == "spatial_wf_sr_umi_ratio":
        return f"{float(x):.2f}"
    return _pct(x)


def _limits(check: dict, override: dict | None) -> dict:
    """The limits of one check: the metric's numbers, then the check's own, then the library's."""
    lim = {k: v for k, v in cat.limits(check["id"]).items() if isinstance(v, (int, float)) and not isinstance(v, bool)}
    # a limit the catalogue words in a sentence has its number on the check (warn_at, fail_at)
    lim.update({k: check[k + "_at"] for k in ("warn", "fail") if isinstance(check.get(k + "_at"), (int, float))})
    for k in ("warn", "fail", "warn_low", "warn_high"):
        if override and isinstance(override.get(k), (int, float)) and not isinstance(override.get(k), bool):
            lim[k] = override[k]
    return lim


def _rule(check: dict, *, outcome: str, severity: str = "none", value=None, message: str = "", threshold: dict | None = None,
          reason: str | None = None, role: str | None = None) -> dict:
    r: dict = {"id": check["id"], "module": MODULE_ID, "role": role or check["role"], "severity": severity, "outcome": outcome}
    if check.get("value"):
        r["metric"] = check["value"]
    if value is not None:
        r["value"] = value if isinstance(value, (int, float, str)) and not isinstance(value, bool) else str(value)
    if threshold:
        r["threshold"] = threshold
    r["message"] = message
    if reason:
        r["not_evaluated_reason"] = reason
    return r


def judge(bundle: dict, *, profiles: dict | None = None, in_library: bool = False, sibling: dict | None = None) -> dict:
    """The ``status`` of an analysis bundle. ``profiles`` is ``config["profiles"]`` of the library; ``sibling``
    is the summary of the analysis of the same reads by the other pipeline, when the library knows the pair."""
    ana = (bundle.get("identity") or {}).get("analysis") or {}
    sec = ((bundle.get("modules") or {}).get(MODULE_ID)) or {}
    s = sec.get("summary") or {}
    notes = sec.get("notes") or []
    assay, pipeline = ana.get("assay"), ana.get("pipeline")
    words = cat.unit_words(assay)
    over_all = (((profiles or {}).get(assay or "") or {}).get("rules")) or {}
    custom_any = False
    rules: list[dict] = []
    for check in cat.checks():
        cid = check["id"]
        if check.get("pipeline") and check["pipeline"] != pipeline:
            continue
        if check.get("assays") and assay not in check["assays"]:
            continue
        over = over_all.get(cid) if isinstance(over_all.get(cid), dict) else None
        if over and over.get("enabled") is False:
            continue
        custom = bool(over and any(k in over for k in ("warn", "fail", "warn_low", "warn_high")))
        custom_any = custom_any or custom
        lim = _limits(check, over)
        basis = "set in the library configuration" if custom else (cat.metric(cid) or {}).get("threshold_basis")
        kind = check["kind"]

        def thr(kind=kind, lim=lim, basis=basis) -> dict:
            t: dict = {"kind": kind, "applied": "rule" if kind == "rule" else "fixed"}
            if kind in ("min", "max"):
                t["warn"], t["fail"] = lim.get("warn"), lim.get("fail")
            elif kind == "range" and lim.get("warn_low") is not None and lim.get("warn_high") is not None:
                t["pass_band"] = [float(lim["warn_low"]), float(lim["warn_high"])]
            if basis:
                t["basis"] = basis
            return t

        def text(template: str, value=None, check=check, lim=lim, **extra) -> str:
            key = check.get("value")
            return template.format(value=_fmt(key, value), warn=_fmt(key, lim.get("warn")), fail=_fmt(key, lim.get("fail")),
                                   lo=_fmt(key, lim.get("warn_low")), hi=_fmt(key, lim.get("warn_high")),
                                   unit=words["name"], units=words["plural"],
                                   a_unit=("an " if words["name"][:1] in "aeiou8" else "a ") + words["name"], **extra)

        # ---- the checks that are not one number against a limit -------------------------------------------
        if cid == "sp_pipeline_alarms":
            if s.get("spatial_alarms_known") is not True:
                rules.append(_rule(check, outcome="not_evaluated", message="The pipeline's alarms could not be read from its report.",
                                   reason="report_unreadable", threshold=thr()))
                continue
            alarms = [n for n in notes if n.get("code") == "pipeline_alarm" and n.get("level") in ("warning", "error")]
            if not alarms:
                rules.append(_rule(check, outcome="ok", message=check["ok"], threshold=thr()))
            for a in alarms:
                sev = "fail" if a["level"] == "error" else "warn"
                value = a.get("value")
                part = f" ({value})" if value not in (None, "", "True", "true") else ""
                msg = check[sev].format(title=a.get("title") or "", value_part=part, message=a.get("message") or "").strip()
                rules.append(_rule(check, outcome="triggered", severity=sev, value=a.get("title"), message=msg, threshold=thr()))
            continue
        if cid == "sp_registration":
            flag = s.get("spatial_registration_flag")
            if flag is None:
                continue                                    # this slide type has no such check
            if not flag:
                rules.append(_rule(check, outcome="ok", message=check["ok"], threshold=thr()))
            else:
                # the pipeline raises its own alarm for it; then this is the same finding, listed once as a reason
                covered = any(n.get("code") == "pipeline_alarm" and "registration" in str(n.get("title", "")).lower() for n in notes)
                rules.append(_rule(check, outcome="triggered", severity="info" if covered else "warn", message=check["warn"],
                                   threshold=thr(), role="info" if covered else None))
            continue
        if cid == "sc_link_runs":
            links = ana.get("links") or []
            if not in_library:
                rules.append(_rule(check, outcome="not_evaluated", role="info", reason="no_library", threshold=thr(),
                                   message="Made without a library: the flow cell runs are named, not looked up."))
            elif not links:
                rules.append(_rule(check, outcome="triggered", severity="info", role="info", message=check["warn_none"], threshold=thr()))
            else:
                missing = [ln for ln in links if not ln.get("resolved")]
                if missing:
                    names = ", ".join(str(ln.get("flow_cell_id")) for ln in missing)
                    msg = check["warn_unresolved"].format(missing=len(missing), n=len(links), names=names,
                                                          are="is" if len(missing) == 1 else "are")
                    rules.append(_rule(check, outcome="triggered", severity="info", role="info", value=len(missing), message=msg, threshold=thr()))
                else:
                    rules.append(_rule(check, outcome="ok", role="info", value=len(links), threshold=thr(),
                                       message=check["ok"].format(n=len(links), s="" if len(links) == 1 else "s")))
            continue
        if cid == "sc_wf_execution":
            failed, tasks = s.get("workflow_failed_tasks"), s.get("workflow_tasks")
            if failed is None:
                rules.append(_rule(check, outcome="not_evaluated", reason="no_source", threshold=thr(),
                                   message="The workflow's trace file is not in the folder."))
            elif failed > 0:
                rules.append(_rule(check, outcome="triggered", severity="warn", value=int(failed), threshold=thr(),
                                   message=check["warn"].format(value=f"{int(failed):,}", tasks=f"{int(tasks or 0):,}")))
            else:
                rules.append(_rule(check, outcome="ok", value=0, threshold=thr(), message=check["ok"].format(tasks=f"{int(tasks or 0):,}")))
            continue
        if cid == "sc_cells":
            cells, expected = s.get("barcodes_cells"), s.get("expected_cells")
            if cells is None or not isinstance(expected, (int, float)) or expected <= 0:
                continue
            frac = float(check.get("warn_fraction_of_expected") or 0.5)
            low = cells < frac * expected
            rules.append(_rule(check, outcome="triggered" if low else "ok", severity="warn" if low else "none", value=int(cells), threshold=thr(),
                               message=check["warn" if low else "ok"].format(value=f"{int(cells):,}", expected=f"{int(expected):,}")))
            continue
        if cid == "sp_units_under_tissue":
            v = s.get("spatial_units_under_tissue")
            if v is None:
                rules.append(_rule(check, outcome="not_evaluated", reason="no_source", message="The number of units under tissue is not known.", threshold=thr()))
                continue
            t = {"kind": "min", "applied": "fixed", "fail": lim.get("fail"),
                 "warn": lim.get("warn") if assay in (check.get("warn_assays") or [assay]) else None}
            if basis:
                t["basis"] = basis
            if lim.get("fail") is not None and v < lim["fail"]:
                rules.append(_rule(check, outcome="triggered", severity="fail", value=int(v), message=text(check["fail"], v), threshold=t))
            elif t["warn"] is not None and v < t["warn"]:
                rules.append(_rule(check, outcome="triggered", severity="warn", value=int(v), message=text(check["warn"], v), threshold=t))
            else:
                rules.append(_rule(check, outcome="ok", value=int(v), message=text(check["ok"], v), threshold=t))
            continue

        # ---- one number against a limit ---------------------------------------------------------------------
        key = check.get("value")
        v = s.get(key) if key else None
        if cid == "sp_wf_sr_umi_ratio":
            mine, other = s.get("barcodes_unique_molecules"), (sibling or {}).get("barcodes_unique_molecules")
            if not mine or not other or pipeline == "spaceranger":
                continue                                    # stated once, on the workflow's side of a pair
            v = mine / other
        if not isinstance(v, (int, float)) or isinstance(v, bool):
            rules.append(_rule(check, outcome="not_evaluated", reason="no_source", threshold=thr(),
                               message=f"{cat.label(cid)}: not available for this analysis."))
            continue
        sev = "none"
        if kind == "min":
            sev = "fail" if (lim.get("fail") is not None and v < lim["fail"]) else "warn" if (lim.get("warn") is not None and v < lim["warn"]) else "none"
        elif kind == "max":
            sev = "fail" if (lim.get("fail") is not None and v > lim["fail"]) else "warn" if (lim.get("warn") is not None and v > lim["warn"]) else "none"
        elif kind == "range":
            lo, hi = lim.get("warn_low"), lim.get("warn_high")
            sev = "warn" if ((lo is not None and v < lo) or (hi is not None and v > hi)) else "none"
        if sev == "fail" and check["role"] != "key":
            sev = "warn"
        if sev == "fail" and not check.get("fail"):
            sev = "warn"
        template = check[sev] if sev != "none" else check["ok"]
        rules.append(_rule(check, outcome="triggered" if sev != "none" else "ok", severity=sev, value=float(v), message=text(template, v),
                           threshold=thr()))

    judged = [r for r in rules if r["outcome"] != "not_evaluated" and r["role"] != "info"]
    verdict = status_model.worst_of_rules(rules) if judged else "unknown"
    trig = [r for r in rules if r["outcome"] == "triggered" and r["severity"] in ("warn", "fail") and r["role"] != "info"]
    alarms = [r for r in trig if r["id"] == "sp_pipeline_alarms"]
    others = [r for r in trig if r["id"] != "sp_pipeline_alarms"]
    pipe = PIPELINE_WORD.get(pipeline, "The pipeline")
    n_ok = sum(1 for r in judged if r["outcome"] == "ok")

    def count(n: int, word: str) -> str:
        return f"{n} {word}{'' if n == 1 else 's'}"
    if verdict == "unknown":
        headline = "No check could be assessed: the files that carry the numbers were not found."
    elif not trig:
        headline = f"All {count(len(judged), 'check')} passed."
    else:
        parts = []
        if alarms:
            parts.append(f"{pipe} raised {'one alarm' if len(alarms) == 1 else count(len(alarms), 'alarm')}.")
        if others:
            failed = [r for r in others if r["severity"] == "fail"]
            if failed:
                parts.append(f"{count(len(failed), 'check')} failed" + (f" and {len(others) - len(failed)} more need a look." if len(others) > len(failed) else "."))
            else:
                more = " more" if alarms else ""
                parts.append(f"One{more} check needs a look." if len(others) == 1 else f"{len(others)}{more} checks need a look.")
        else:
            parts.append("Every other check is within its limit." if n_ok else "")
        headline = " ".join(p for p in parts if p)
    status = {
        "verdict": verdict, "headline": headline, "detail": NOT_ABOUT_RUNS,
        "profile": {"id": PROFILE_ID, "version": cat.version(), "source": "library_config" if custom_any else "builtin"},
        "evaluated": now_z(), "rules": rules, "cause": None, "modifiers": [], "provisional": False, "pending": [],
        "judged": {"at": now_z(), "as_of": None, "history_hash": None, "library": bool(in_library), "frozen": False},
    }
    part = sec.get("state") == "partial"
    if part:
        status["modifiers"].append({"code": "partly_read", "message": "Part of the analysis could not be read or was left out; "
                                                                      "the notes of the analysis say what."})
    return status


def reasons(status: dict, limit: int = 3) -> list[dict]:
    """The triggered checks behind a verdict, the pipeline's alarms first, as the index stores them."""
    trig = [r for r in status.get("rules", []) if r.get("outcome") == "triggered" and r.get("severity") in ("warn", "fail") and r.get("role") != "info"]
    trig.sort(key=lambda r: (r["id"] != "sp_pipeline_alarms", r.get("severity") != "fail", r.get("role") != "key"))
    return [{"check": r["id"], "severity": status_model.rule_level(r), "message": r.get("message", "")} for r in trig[:limit]]
