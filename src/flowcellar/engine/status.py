"""The status engine: every tier 1 key and supporting check, the two-part verdict, history as of
the run's start, flow cell use, and the stored sentences.

One implementation: rules are evaluated here, once, and the viewer displays what is stored. Sentences of single
checks are the catalogue's display copy (``display.message`` with ``limit_text``); the few sentences the catalogue
does not hold yet are in ``copy_pending``.
"""
from __future__ import annotations

from typing import Sequence

from flowcellar import catalogue
from flowcellar.contracts import Card, Declared
from flowcellar.engine import baselines as BL
from flowcellar.engine import copy_pending as CP
from flowcellar.engine import flowcell_use as FU
from flowcellar.engine.card import facts_of
from flowcellar.engine.limits import BAD_STOPS, FLOOR, HORIZONS, PENDING, SPEED_FAIL_BELOW, SPEED_TOLERANCE, Limits
from flowcellar.model import bundle as B
from flowcellar.model import status as S
from flowcellar.util.fingerprint import hash_obj
from flowcellar.util.text import fmt_bases, fmt_count, fmt_pct
from flowcellar.util.timeutil import now_z

ORDER = ("pass", "warn", "fail")


# ------------------------------------------------------------------------------------------------ wording
def _display(mid: str) -> dict:
    return catalogue.display(mid)


def label(mid: str) -> str:
    return _display(mid).get("label") or mid


def panel(mid: str) -> str | None:
    return _display(mid).get("panel")


def message(mid: str, value_text: str, tail: str = "", *, key: str = "message", **fields) -> str:
    """The catalogue's sentence of a check with its value, and a second sentence about the limit."""
    tpl = _display(mid).get(key) or _display(mid).get("message") or (label(mid) + ": {value}.")
    s = tpl.replace("{value}", value_text)
    for k, v in fields.items():
        s = s.replace("{" + k + "}", str(v))
    return (s + " " + tail).strip() if tail else s


def limit_sentence(kind: str, **fields) -> str:
    s = catalogue.limit_text()[kind]
    for k, v in fields.items():
        s = s.replace("{" + k + "}", str(v))
    return s


def pct_vs(v: float, limit: float) -> str:
    """A share worded with as few decimals as keep it apart from its limit (69.8 %, not 70 %, against 70 %)."""
    for d in (0, 1, 2):
        if f"{v * 100:.{d}f}" != f"{limit * 100:.{d}f}":
            return fmt_pct(v, d)
    return fmt_pct(v, 2)


def _pct_limit(x: float) -> str:
    return fmt_pct(x, 0 if abs(x * 100 - round(x * 100)) < 1e-9 else 1)


def duration_text(hours: float) -> str:
    return f"{round(hours * 60)} minutes" if hours < 1 else f"{hours:.1f} hours" if hours < 10 else f"{hours:.0f} h"


def _num(x):
    """Whole numbers as integers in stored thresholds."""
    if isinstance(x, float) and x.is_integer() and abs(x) < 1e15:
        return int(x)
    return x


# ----------------------------------------------------------------------------------------------- history
def _prepare(r: dict, history: Sequence[Card], baseline_groups, flow_cell: dict | None, declared: Declared | None,
             *, stored_groups: bool = False) -> list[dict]:
    """The facts of the other runs, each with its comparison group and the use of its flow cell. The group comes
    from the configuration's dated rules; with ``stored_groups`` (no configuration at hand) it is the one each
    run was last judged with, as kept in its card."""
    H = []
    for c in history or ():
        f = (c or {}).get("facts") if isinstance(c, dict) else None
        if not isinstance(f, dict) or not f.get("run_key") or f.get("run_key") == r.get("run_key"):
            continue
        H.append(dict(f))
    if stored_groups:
        r["group"] = r.get("group") or BL.group_key(r, None)
        for p in H:
            p["group"] = p.get("group") or BL.group_key(p, None)
    else:
        lab = BL.group_label(r, baseline_groups, (declared.baseline_group if declared else None))
        r["group_label"] = lab
        r["group"] = BL.group_key(r, lab)
        for p in H:
            p["group"] = BL.group_key(p, BL.group_label(p, baseline_groups))
    uses = FU.annotate(H + [r], flow_cell)
    for p in H:
        p["use_index"] = uses.get(p["run_key"], {}).get("index", 1)
    return H


def _horizon(hours: float | None) -> int | None:
    return next((h for h in HORIZONS if hours is not None and hours >= h), None)


def _explanation_table(r: dict, hz: int | None, with_speed: bool) -> list[tuple]:
    """(metric id, facts column, floor, aspect, formatter, extra peer filter, absolute floor) of every number that is
    placed against history without judging a first-use run."""
    hzc = hz if hz in (24, 48, 72) else None
    t = [("fc_occupancy_median_24h", "occ_median_24h", FLOOR["occ"], "loading", lambda v: fmt_pct(v), None, 0.0),
         ("fc_pores_first_scan", "mux0_pores", FLOOR["pores"], "flow_cell", lambda v: f"{v:,.0f}", None, 0.0),
         ("fc_pore_survival_24h_frac", "pore_survival_24h", FLOOR["other"], "other", lambda v: fmt_pct(v), None, 0.0),
         ("fc_sequencing_duty", "seq_duty", FLOOR["other"], "other", lambda v: fmt_pct(v), None, 0.0)]
    if hzc:
        t.append((f"yield_per_pore_mb_{hzc}h", f"mbpp_{hzc}h", FLOOR["per_pore"], "other", lambda v: f"{v:.1f} Mb", None, 0.0))
        k = 48 if hzc >= 48 else 24
        t.append((f"yield_kb_per_active_channel_hour_{k}h", f"kb_per_active_ch_hour_{k}h", FLOOR["kb"], "library", lambda v: f"{v:,.0f}", None, 0.0))
    if with_speed:
        t.append(("dyn_speed_est_median", "speed_med", FLOOR["speed"], "other", lambda v: f"{v:.0f} bases per second", None, 0.0))
    mk = r.get("mk")
    t.append(("fc_adapter_frac", "adapter_frac", 0.25, "other", lambda v: fmt_pct(v, 1), (lambda p: p.get("mk") == mk), 0.03))
    return t


def _history_hash(r: dict, hist: BL.History) -> str | None:
    """A hash of the keys and values of the earlier runs this run can be compared with (output at its horizon and
    the numbers of quick detail). It does not change when a later run arrives, nor when an earlier run gets its
    read count."""
    hz = _horizon(r.get("run_hours"))
    used: set = set()
    if hz:
        for unit in ("reads", "bases"):
            used.update((p["run_key"], f"{unit}_{hz}h", p[f"{unit}_{hz}h"]) for p in hist.peers(r, f"{unit}_{hz}h"))
    for _mid, col, _fl, _a, _fmt, extra, _fa in _explanation_table(r, hz, with_speed=False):
        used.update((p["run_key"], col, p[col]) for p in hist.peers(r, col, extra))
    if not used:
        return None
    return "h1:" + hash_obj(sorted([list(u) for u in used], key=lambda u: (u[0], u[1])), 16)


# ------------------------------------------------------------------------------------------------- judge
class _Rules:
    def __init__(self, limits: Limits):
        self.rules: list[dict] = []
        self.levels: dict[int, str] = {}
        self.limits = limits

    def add(self, mid, role, level, value, threshold, msg, module, *, tier="run_health", metric=None, **extra) -> dict:
        if role != "key" and level == "fail":
            level = "warn"
        if role == "info":
            level = "pass"
        rule = {"id": mid, "module": module, "role": role, "tier": tier, "severity": level if level != "pass" else "none",
                "outcome": "ok" if level == "pass" else "triggered", "metric": metric or mid, "value": value,
                "threshold": threshold, "message": msg, "panel": panel(mid)}
        rule.update(extra)
        self.rules.append(rule)
        self.levels[id(rule)] = level
        return rule

    def not_assessed(self, mid, role, module, reason, *, tier="run_health", metric=None, **extra) -> None:
        rule = {"id": mid, "module": module, "role": role, "tier": tier, "severity": "none", "outcome": "not_evaluated",
                "metric": metric or mid, "value": None, "threshold": {"kind": "rule", "basis": reason, "applied": "rule"},
                "message": limit_sentence("not_assessed", reason=reason), "not_evaluated_reason": reason, "panel": panel(mid)}
        rule.update(extra)
        self.rules.append(rule)
        self.levels[id(rule)] = "pass"

    def level(self, rule) -> str:
        return self.levels.get(id(rule), "pass")

    def worst(self, tier: str) -> str:
        v = "pass"
        for x in self.rules:
            if x["tier"] == tier and x.get("role") != "info" and ORDER.index(self.level(x)) > ORDER.index(v):
                v = self.level(x)
        return v


def _fixed_threshold(lim: dict, role: str, basis_key: str) -> dict:
    basis = CP.BASIS["configured"] if lim.get("configured") else CP.BASIS.get(basis_key, CP.BASIS["default"])
    return {"kind": lim["kind"] if lim["kind"] in ("min", "max") else "min", "warn": _num(lim.get("warn")),
            "fail": _num(lim.get("fail")) if role == "key" else None, "basis": basis, "applied": "fixed"}


def _level_fixed(v: float, lim: dict, role: str) -> str:
    warn, fail = lim.get("warn"), lim.get("fail") if role == "key" else None
    if lim["kind"] == "max":
        return "fail" if fail is not None and v > fail else "warn" if warn is not None and v > warn else "pass"
    return "fail" if fail is not None and v < fail else "warn" if warn is not None and v < warn else "pass"


def _empty_status(profile: dict, when: str, **kw) -> dict:
    st = {"profile": profile, "evaluated": when, "rules": [], "modifiers": [], "estimates": [], "provisional": False, "pending": [],
          "explanations": [], "detail": ""}
    st.update(kw)
    return st


def judge_facts(r: dict, hist: BL.History, use: dict, limits: Limits, *, declared: Declared, in_library: bool,
                read_rules: Sequence[dict] = (), read_level: dict | None = None, read_files_waiting: int = 0) -> dict:
    """The ``status`` member for one run. ``r``: its facts (``card.facts_of``) with ``group`` set; ``hist``: the other
    runs; ``use``: the use of its flow cell."""
    now = now_z()
    profile = limits.profile()
    if declared.flow_cell_use is not None or declared.flow_cell_check_pores is not None or declared.baseline_group:
        profile = dict(profile, source="command_line") if not limits.custom else profile
    as_of = (r["started"] or "")[:19] + "Z" if r.get("started") else None
    judged = {"at": now, "as_of": as_of, "history_hash": None, "library": bool(in_library), "frozen": True}
    read_rules = [dict(x) for x in read_rules or ()]
    has_read = bool(read_level and read_level.get("state") in ("partial", "complete")) or bool(read_rules)
    t2_assessed = [x for x in read_rules if x.get("outcome") != "not_evaluated" and x.get("role") != "info"]
    v2 = S.worst_of_rules(t2_assessed, None) if t2_assessed else ("unknown" if has_read else "not_available")
    fresh = use["index"] == 1
    fam, ag = r.get("family"), r.get("agroup") or "unknown"
    hours = r.get("run_hours")

    def reused_mod() -> dict:
        return {"code": "reused_flow_cell", "message": CP.MODIFIER["reused_flow_cell"].format(index=use["index"])}

    def read_detail() -> list[str]:
        if not has_read:
            return [CP.DETAIL["read_data_absent"]]
        rl = read_level or {}
        out = [CP.DETAIL["read_data_samples"].format(with_data=rl.get("samples_with_data"), expected=rl.get("samples_expected"))
               if rl.get("samples_expected") else CP.DETAIL["read_data_present"]]
        if v2 == "unknown":
            out.append(CP.DETAIL["read_data_not_assessed"])
        return out

    def finish(st: dict, v1: str) -> dict:
        st["by_tier"] = {"run_health": v1, "read_level": v2}
        st["verdict"] = S.combine_verdicts(v1, v2)
        st["rules"] = st["rules"] + read_rules
        st["judged"] = judged
        if read_files_waiting and has_read:
            # the read-data part is a verdict on what has been read so far: say so wherever the status is shown
            total = [x for x in read_rules if x.get("role") != "info"]
            key = "read_data_so_far_one" if read_files_waiting == 1 else "read_data_so_far"
            st["modifiers"].append({"code": "read_data_so_far", "message": CP.MODIFIER[key].format(
                files=read_files_waiting, assessed=len(t2_assessed), total=len(total))})
        if limits.custom:
            st["modifiers"].append({"code": "custom_profile", "message": CP.MODIFIER["custom_profile"]})
        return st

    # ---- runs that cannot be judged on their flow cell at all ----------------------------------------
    if r.get("run_state") == "raw_only" or (not r.get("has_tier1") and r.get("run_state") not in ("running", "interrupted")):
        cause = "raw_only" if r.get("run_state") == "raw_only" else "inputs_missing"
        st = _empty_status(profile, now, cause=cause, headline=CP.HEADLINE[cause], detail=" ".join(read_detail()))
        return finish(st, "incomplete")
    if r.get("run_state") == "running":
        judged["frozen"] = False
        st = _empty_status(profile, now, cause=None, headline=CP.HEADLINE["running"], detail=" ".join(read_detail()))
        return finish(st, "running")
    interrupted = r.get("run_state") == "interrupted"
    if hours is None or not r.get("has_tier1"):
        cause = "no_final_summary" if interrupted else "inputs_missing"
        st = _empty_status(profile, now, cause=cause, headline=CP.HEADLINE[cause], detail=" ".join(read_detail()))
        return finish(st, "incomplete")
    dur_txt = duration_text(hours)
    if hours < 1:
        st = _empty_status(profile, now, cause="stopped_early", headline=CP.HEADLINE["stopped_early"].format(duration=dur_txt),
                           detail=" ".join(read_detail()) if has_read else "")
        # too short to be a loading or to be judged as one: the use of the flow cell stays in the relations
        return finish(st, "incomplete")

    R = _Rules(limits)
    notes: set[str] = set()

    # ---- how the run ended, messages, files (written-out rules) -----------------------------------------
    if r.get("has_run_section"):
        state, stop = str(r.get("state") or ""), str(r.get("stop_reason") or "")
        if limits.enabled("run_stop_cause"):
            bad = stop in BAD_STOPS or "ERROR" in state
            by_person = "USER" in state
            ok_msg = CP.RULE["stop_by_person" if by_person else "stop_planned" if "COMPLETED" in state else "stop_other"].format(duration=dur_txt)
            R.add("run_stop_cause", limits.role("run_stop_cause"), "fail" if bad else "pass", r.get("state"),
                  {"kind": "rule", "basis": CP.BASIS["run_stop_cause"], "applied": "rule"},
                  ok_msg if not bad else message("run_stop_cause", (stop if stop in BAD_STOPS else state).replace("STOPPED_", "").replace("PROTOCOL_", "").replace("_", " ").lower()),
                  "run", metric="run_stop_cause")
        if limits.enabled("run_messages"):
            n_err, n_warn, n_disk = r.get("msg_error") or 0, r.get("msg_warn") or 0, r.get("msg_warn_disk") or 0
            nw = n_warn - n_disk
            lvl = "fail" if n_err > 0 else "warn" if nw > 0 else "pass"
            txt = CP.RULE["messages_none"] if lvl == "pass" else CP.RULE["messages_error"].format(n=n_err) if n_err else \
                (CP.RULE["messages_warning_one"] if nw == 1 else CP.RULE["messages_warning"].format(n=nw))
            if lvl == "pass" and n_disk:
                txt += CP.RULE["messages_computer"].format(n=n_disk)
            # the stored value is the count the sentence speaks of: errors when one is marked, else the warnings that
            # are not notices about the computer
            R.add("run_messages", limits.role("run_messages"), lvl, n_err if n_err > 0 else max(nw, 0),
                  {"kind": "rule", "basis": CP.BASIS["run_messages"], "applied": "rule"}, message("run_messages", txt), "run",
                  metric="run_messages_error" if n_err > 0 or nw <= 0 else "run_messages_warning")
        if limits.enabled("run_output_files"):
            fb, wf, wb = r.get("fallback_files") or 0, r.get("write_failed") or 0, r.get("write_backlog") or 0
            files_bad = fb > 0 or wf > 0 or wb > 0
            R.add("run_output_files", limits.role("run_output_files"), "warn" if files_bad else "pass", fb,
                  {"kind": "rule", "basis": CP.BASIS["run_output_files"], "applied": "rule"},
                  message("run_output_files", CP.RULE["files_value"].format(n=fb) if fb else CP.RULE["files_bytes"]) if files_bad else CP.RULE["files_ok"],
                  "run", metric="files_in_fallback")
    else:
        for mid, metric in (("run_stop_cause", "run_stop_cause"), ("run_messages", "run_messages_error"), ("run_output_files", "files_in_fallback")):
            if limits.enabled(mid):
                R.not_assessed(mid, limits.role(mid), "run", CP.REASON["report_unreadable"], metric=metric)

    # ---- fixed limits -------------------------------------------------------------------------------------
    w = limits.warranty()
    pores = r.get("mux0_pores")
    if fresh and pores is not None and w is not None and limits.enabled("fc_pores_first_scan"):
        lim = limits.fixed("fc_pores_first_scan")
        role = limits.role("fc_pores_first_scan")
        if declared.flow_cell_check_pores is not None:
            v = int(declared.flow_cell_check_pores)
            th = _fixed_threshold(lim, role, "fc_pores_flow_cell_check")
            R.add("fc_pores_first_scan", role, _level_fixed(v, lim, role), pores, th,
                  CP.RULE["flow_cell_check"].format(value=f"{v:,}", limit=f"{w:,.0f}"), "fc")
        else:
            lvl = _level_fixed(pores, lim, role)
            key = "message_use_unknown" if (lvl != "pass" and use["source"] != "declared") else "message"
            R.add("fc_pores_first_scan", role, lvl, pores, _fixed_threshold(lim, role, "fc_pores_first_scan"),
                  message("fc_pores_first_scan", f"{pores:,}", key=key, limit=f"{w:,.0f}"), "fc")
    if fresh and r.get("chan_pore_frac") is not None and limits.enabled("fc_channels_with_pore_first_scan"):
        lim = limits.fixed("fc_channels_with_pore_first_scan")
        if lim:
            role = limits.role("fc_channels_with_pore_first_scan")
            v = round(r["chan_pore_frac"], 4)
            val = CP.RULE["count_of"].format(count=f"{r['msg0_num_pores']:,}", total=f"{r['channel_count']:,}", share=pct_vs(v, lim["warn"]))
            R.add("fc_channels_with_pore_first_scan", role, _level_fixed(v, lim, role), v,
                  _fixed_threshold(lim, role, "fc_channels_with_pore_first_scan"),
                  message("fc_channels_with_pore_first_scan", val, limit_sentence("min", limit=_pct_limit(lim["warn"]))), "fc",
                  metric="fc_channels_with_pore_first_scan_frac")

    def share_check(mid, col, need_24h=False, digits=None, basis=None, tail_key=None):
        v = r.get(col)
        if v is None or (need_24h and hours < 24) or not limits.enabled(mid):
            return
        lim = limits.fixed(mid)
        if not lim or lim.get("warn") is None:
            return
        role = limits.role(mid)
        v = round(v, 4)
        val = pct_vs(v, lim["warn"]) if digits is None else fmt_pct(v, digits)
        own = _display(mid).get("limit_text")
        tail = own.replace("{limit}", _pct_limit(lim["warn"])) if own else limit_sentence(lim["kind"] if lim["kind"] in ("min", "max") else "min", limit=_pct_limit(lim["warn"]))
        R.add(mid, role, _level_fixed(v, lim, role), v, _fixed_threshold(lim, role, basis or mid), message(mid, val, tail), "fc")
    share_check("fc_occupancy_first_hour", "occ_30_90")
    share_check("fc_occupancy_median", "occ_median_run", need_24h=True)
    share_check("fc_recovering_frac", "rec_of_working", digits=1)
    share_check("fc_recovering_last_third_frac", "rec_last_third_of_working", need_24h=True, digits=1)
    if r.get("temp_max_dev_c") is not None and limits.enabled("fc_temperature_max_dev_c"):
        lim = limits.fixed("fc_temperature_max_dev_c")
        if lim:
            role = limits.role("fc_temperature_max_dev_c")
            v = round(r["temp_max_dev_c"], 2)
            R.add("fc_temperature_max_dev_c", role, _level_fixed(v, lim, role), v, _fixed_threshold(lim, role, "fc_temperature_max_dev_c"),
                  message("fc_temperature_max_dev_c", f"{v:.1f}", limit_sentence("max", limit=f"{lim['warn']:g} °C")) if v > 0 else CP.RULE["temperature_ok"],
                  "fc")

    # ---- checks that need the read count -------------------------------------------------------------------
    full = r.get("full")
    pending: list[str] = []
    applies = [x for x in PENDING if limits.enabled(x) and not (x == "dyn_speed_drop_frac" and hours < 24)
               and not (x in ("fc_active_channels", "len_est_n50") and not fresh)]
    module_of = {"dyn_speed_est_median": "dyn", "dyn_speed_drop_frac": "dyn", "dyn_signal_positive_frac": "dyn", "dyn_unblock_frac": "dyn",
                 "fc_active_channels": "fc", "len_est_n50": "len.est"}
    metric_of = {"fc_active_channels": "fc_active_channels_frac"}
    summary_missing = False
    if full and not full.get("error"):
        part = "" if full["complete"] else " " + CP.RULE["from_part"].format(share=fmt_pct(full.get("coverage_frac") or 0.0))
        if "fc_active_channels" in applies:
            lim = limits.fixed("fc_active_channels")
            role = limits.role("fc_active_channels")
            if full["complete"] and full.get("active_ch_frac") is not None and lim:
                v = round(full["active_ch_frac"], 4)
                val = CP.RULE["count_of"].format(count=f"{full['active_channels']:,}", total=f"{r['channel_count']:,}", share=pct_vs(v, lim["warn"]))
                R.add("fc_active_channels", role, _level_fixed(v, lim, role), v, _fixed_threshold(lim, role, "fc_active_channels"),
                      message("fc_active_channels", val, limit_sentence("min", limit=_pct_limit(lim["warn"]))), "fc", metric="fc_active_channels_frac")
            elif not full["complete"]:
                R.not_assessed("fc_active_channels", role, "fc", CP.REASON["channels_need_every_read"], metric="fc_active_channels_frac")
        if "dyn_speed_est_median" in applies:
            lo, hi, s = r.get("speed_min"), r.get("speed_max"), full.get("speed_med")
            role = limits.role("dyn_speed_est_median")
            if lo and hi and s is not None:
                lvl = "pass" if lo * (1 - SPEED_TOLERANCE) <= s <= hi * (1 + SPEED_TOLERANCE) else "fail" if s < lo * (1 - SPEED_FAIL_BELOW) else "warn"
                R.add("dyn_speed_est_median", role, lvl, s,
                      {"kind": "run_config", "pass_band": [round(lo * (1 - SPEED_TOLERANCE), 1), round(hi * (1 + SPEED_TOLERANCE), 1)],
                       "fail": round(lo * (1 - SPEED_FAIL_BELOW), 1) if role == "key" else None, "basis": CP.BASIS["dyn_speed_est_median"],
                       "applied": "run_config"},
                      message("dyn_speed_est_median", CP.RULE["speed_value"].format(value=f"{s:.0f}"),
                              limit_sentence("run_config", low=f"{lo:.0f}", high=f"{hi:.0f}") + part), "dyn")
            elif s is not None:
                R.not_assessed("dyn_speed_est_median", role, "dyn", CP.REASON["no_target_speed"])

        def dyn_share(mid, col, module="dyn"):
            v = full.get(col)
            lim = limits.fixed(mid)
            if mid not in applies or v is None or not lim or lim.get("warn") is None:
                return
            role = limits.role(mid)
            v = round(v, 4)
            R.add(mid, role, _level_fixed(v, lim, role), v, _fixed_threshold(lim, role, "default"),
                  message(mid, fmt_pct(v, 1), limit_sentence(lim["kind"], limit=_pct_limit(lim["warn"])) + part), module)
        dyn_share("dyn_speed_drop_frac", "speed_drop")
        dyn_share("dyn_signal_positive_frac", "sp_frac")
        dyn_share("dyn_unblock_frac", "unblock_frac")
        if "len_est_n50" in applies and full.get("est_n50") is not None:
            lim = limits.fixed("len_est_n50")
            if lim and lim.get("warn") is not None:
                role = limits.role("len_est_n50")
                v = round(full["est_n50"])
                R.add("len_est_n50", role, _level_fixed(v, lim, role), v, _fixed_threshold(lim, role, "default"),
                      message("len_est_n50", CP.RULE["bases_value"].format(value=f"{v:,}"),
                              limit_sentence("min", limit=CP.RULE["bases_value"].format(value=f"{lim['warn']:,.0f}")) + part), "len.est")
    elif full and full.get("error"):
        for x in applies:
            R.not_assessed(x, limits.role(x), module_of[x], CP.REASON["summary_unreadable"], metric=metric_of.get(x))
    elif not r.get("summary_present", True):
        summary_missing = True
        for x in applies:
            R.not_assessed(x, limits.role(x), module_of[x], CP.REASON["summary_missing"], metric=metric_of.get(x))
    else:
        pending = list(applies)

    # ---- outcome: output at the longest fixed horizon reached ---------------------------------------------------
    hz = _horizon(hours)
    outcome = None
    capped = False
    can_fail_output = use["source"] == "declared" or w is None or (pores or 0) >= w
    if hz and fresh:
        for unit, role in limits.output_units():
            col = f"{unit}_{hz}h"
            x = r.get(col)
            mid = f"yield_{'reads' if unit == 'reads' else 'estimated_bases'}_at_{hz}h"
            if x is None or not limits.enabled(mid):
                continue
            role = limits.overrides.get(mid, {}).get("role", role)
            fmt = fmt_count if unit == "reads" else fmt_bases
            b, trust, n = hist.base(r, col, FLOOR["yield"])
            fl = limits.output_levels(unit, hz)
            level, tail, th = "pass", "", {}
            applied = "fixed"
            if b:
                z = b["z"]
                if role == "key":
                    level = ("fail" if trust == "trusted" else "warn") if z <= -2.5 else "warn" if z <= -1.5 else "pass"
                    blk = BL.block(r, b, -1.5, -2.5, hz, r.get("group_label"))
                else:
                    level = "warn" if z <= -2.5 else "pass"
                    blk = BL.block(r, b, -2.5, None, hz, r.get("group_label"))
                if level == "fail" and not can_fail_output:
                    level, capped = "warn", True
                applied = "baseline"
                th = {"kind": "baseline", "warn": None, "fail": None, "peers": b["n"],
                      "basis": CP.BASIS["history"].format(n=b["n"], horizon=hz), "baseline": blk}
                lo_, hi_ = _low_of_range(fmt(max(blk["band"][0], 0)), fmt(blk["band"][1])), fmt(blk["band"][1])
                if level != "pass":
                    tail = limit_sentence("history_warn", horizon=f"{hz} h", level=fmt(blk["warn_level"]), n=b["n"])
                else:
                    tail = limit_sentence("history_ok", low=lo_, high=hi_, n=b["n"])
                if b["guard"]:
                    tail += " " + limit_sentence("history_spread", n=b["n"])
                elif trust == "provisional":
                    tail += " " + limit_sentence("history_early", n=b["n"])
            if (not b) or trust != "trusted":
                if fl:
                    lf = "fail" if (fl[1] is not None and x < fl[1]) else "warn" if x < fl[0] else "pass"
                    if lf == "fail" and not can_fail_output:
                        lf, capped = "warn", True
                    if role != "key" and lf == "fail":
                        lf = "warn"
                    if b:
                        th["warn"], th["fail"] = _num(fl[0]), (_num(fl[1]) if role == "key" else None)
                        th["basis"] += CP.BASIS["history_fixed_also"].format(horizon=hz, level=fmt(fl[0]))
                        if ORDER.index(lf) > ORDER.index(level):
                            level, applied = lf, "fixed"
                            tail = CP.RULE["fixed_level_also"].format(horizon=hz, level=fmt(fl[0])) + " " + \
                                limit_sentence("history_spread" if b["guard"] else "history_early", n=b["n"])
                        elif level == "pass":
                            applied = "fixed"
                    else:
                        level, applied = lf, "fixed"
                        th = {"kind": "min", "warn": _num(fl[0]), "fail": _num(fl[1]) if role == "key" else None, "peers": n,
                              "basis": CP.BASIS["fixed_only"].format(assay=ag, family=fam, horizon=hz)}
                        tail = limit_sentence("history_none", horizon=f"{hz} h", level=fmt(fl[0]))
                elif not b:
                    th = {"kind": "baseline", "warn": None, "fail": None, "peers": n, "basis": CP.BASIS["no_level"]}
                    tail = CP.RULE["output_no_level"]
                    applied = "baseline"
            th["applied"] = applied
            rule = R.add(mid, role, level, x, th, message(mid, fmt(x), tail), "yield")
            if role == "key" and outcome is None:
                bl = th.get("baseline") or {}
                outcome = {"metric": mid, "horizon_h": hz, "value": x, "n": (b or {}).get("n", 0), "trust": trust, "unit": unit,
                           "level": bl.get("warn_level") if applied == "baseline" else th.get("warn"),
                           "guard": bool(bl.get("spread_guard")), "has_level": bool(b or fl), "rule": rule}

    # ---- explanations: numbers placed against history that never change a first-use verdict -------------------
    expl = []
    reused_outcome = False
    for mid, col, fl_, aspect, fmtv, extra, fl_abs in _explanation_table(r, hz, with_speed=bool(full and full.get("complete"))):
        if not fresh and mid in ("fc_pores_first_scan", "fc_pore_survival_24h_frac", "fc_sequencing_duty"):
            continue
        x = r.get(col)
        if x is None:
            continue
        b, trust, n = hist.base(r, col, fl_, extra, fl_abs)
        if not b:
            continue
        blk = BL.block(r, b, label=r.get("group_label"))
        z = b["z"]
        direction = "low" if z <= -1.5 else "high" if z >= 1.5 else "typical"
        lo_b = max(blk["band"][0], 0.0)                  # none of these numbers can be negative
        strip = fmtv(lo_b).replace(" %", "").replace(" Mb", "").replace(" bases per second", "")
        m_ = f"{label(mid)}: {fmtv(x)}. " + limit_sentence("history_ok", low=strip, high=fmtv(blk["band"][1]), n=b["n"])
        if mid == "fc_pores_first_scan" and z <= -4 and trust == "trusted":
            m_ += " " + CP.MODIFIER["flow_cell_use_possibly_reused"]
            notes.add("possibly_reused")
        expl.append({"id": mid, "aspect": aspect, "value": BL.sig(x, 6), "direction": direction, "baseline": blk, "message": m_})
        if (not fresh) and mid.startswith("yield_per_pore_mb_") and limits.enabled(mid):
            level = ("fail" if trust == "trusted" else "warn") if z <= -2.5 else "warn" if z <= -1.5 else "pass"
            b2 = BL.block(r, b, -1.5, -2.5, hz, r.get("group_label"))
            rule = R.add(mid, "key", level, round(x, 2),
                         {"kind": "baseline", "warn": None, "fail": None, "peers": b["n"], "basis": CP.BASIS["reused"].format(horizon=hz, n=b["n"]),
                          "baseline": b2, "applied": "baseline"},
                         message(mid, f"{x:.1f} Mb", limit_sentence("history_ok", low=f"{max(b2['band'][0], 0):.1f}", high=f"{b2['band'][1]:.1f} Mb", n=b["n"])),
                         "yield")
            reused_outcome = True
            outcome = {"metric": mid, "horizon_h": hz, "value": round(x, 2), "n": b["n"], "trust": trust, "unit": "per_pore",
                       "level": b2["warn_level"], "guard": bool(b2["spread_guard"]), "has_level": True, "rule": rule}

    # ---- roll-up ------------------------------------------------------------------------------------------------
    v1 = R.worst("run_health")
    cause = None
    judged_output = (outcome is not None and outcome.get("has_level")) or reused_outcome
    if v1 == "pass" and interrupted:
        v1, cause = "incomplete", "no_final_summary"
    elif v1 == "pass" and not judged_output and hours < 6:
        v1, cause = "incomplete", "output_not_judged"
    mods = []
    if not fresh:
        mods.append(reused_mod())
    if "possibly_reused" in notes:
        mods.append({"code": "flow_cell_use_unknown", "message": CP.MODIFIER["flow_cell_use_possibly_reused"]})
    if capped:
        mods.append({"code": "flow_cell_use_unknown", "message": CP.MODIFIER["flow_cell_use_unknown"]})
    first_outcome = outcome if outcome and outcome["unit"] != "per_pore" else None
    if first_outcome and first_outcome["trust"] == "provisional":
        mods.append({"code": "baseline_spread_guard", "message": limit_sentence("history_spread", n=first_outcome["n"])} if first_outcome["guard"] else
                    {"code": "baseline_provisional", "message": limit_sentence("history_early", n=first_outcome["n"])})
    if first_outcome and first_outcome["trust"] == "none" and fresh:
        mods.append({"code": "no_history", "message": CP.MODIFIER["no_history" if first_outcome["has_level"] else "no_history_no_level"]})
    if pending:
        mods.append({"code": "quick_detail", "message": CP.MODIFIER["quick_detail"]})
    if full and not full.get("error") and not full["complete"]:
        mods.append({"code": "reads_partly_counted", "message": CP.MODIFIER["reads_partly_counted"].format(share=fmt_pct(full.get("coverage_frac") or 0.0))})
    if full and full.get("error"):
        mods.append({"code": "reads_not_counted", "message": CP.MODIFIER["reads_not_counted"]})
    trig = [x for x in R.rules if x["outcome"] == "triggered"]
    trig.sort(key=lambda x: (x["tier"] != "run_health", x["role"] != "key", R.level(x) != "fail"))
    so_far = " so far" if pending else ""
    if cause == "no_final_summary":
        head = CP.HEADLINE["no_final_summary"]
    elif v1 == "incomplete":
        head = CP.HEADLINE["output_not_judged"].format(duration=dur_txt)
    elif not trig:
        facts = []
        if first_outcome:
            facts.append(CP.HEADLINE_FACT["reads" if first_outcome["unit"] == "reads" else "bases"].format(
                value=fmt_count(first_outcome["value"]) if first_outcome["unit"] == "reads" else fmt_bases(first_outcome["value"]),
                horizon=first_outcome["horizon_h"]))
        if r.get("occ_30_90") is not None:
            facts.append(CP.HEADLINE_FACT["occupancy"].format(value=fmt_pct(r["occ_30_90"])))
        head = CP.HEADLINE["ran_well"].format(so_far=so_far, facts=", ".join(facts)) if facts else CP.HEADLINE["ran_well_plain"].format(so_far=so_far)
    else:
        def names(rules, skip=()) -> list[str]:
            out: list[str] = []
            for x in rules:
                nme = label(x["id"])
                nme = nme[0].lower() + nme[1:] if not nme.startswith(("ONT", "MinKNOW")) else nme
                if nme not in out and nme not in skip:
                    out.append(nme)
            return out

        def listed(items: list[str], n: int) -> str:
            return ", ".join(items[:n]) + (f" and {len(items) - n} more" if len(items) > n else "")
        # "Failed on" names only the checks that failed; a check that only warned (a supporting check can never
        # fail) follows as "also needs a look"
        failed = names([x for x in trig if R.level(x) == "fail"]) if v1 == "fail" else []
        if failed:
            also = names([x for x in trig if R.level(x) != "fail"], skip=failed)
            head = ""
            for n_fail, n_also in ((3, 2), (2, 1), (1, 1), (1, 0)):
                head = CP.HEADLINE["failed_on"].format(list=listed(failed, n_fail))
                if also and n_also:
                    head += CP.HEADLINE["also_needs_a_look"].format(list=listed(also, n_also))
                elif also:
                    head += CP.HEADLINE["also_needs_a_look"].format(list=f"{len(also)} more check{'s' if len(also) != 1 else ''}")
                if len(head) <= 140:
                    break
        else:
            seen = names(trig)
            head = CP.HEADLINE["needs_a_look"].format(list=listed(seen, 3))
            if len(head) > 140:
                head = CP.HEADLINE["needs_a_look"].format(list=", ".join(seen[:2]) + f" and {len(seen) - 2} more")
            if len(head) > 140:
                head = CP.HEADLINE["many"].format(n=len(seen), first=seen[0])
    det = []
    if pending:
        det.append(CP.DETAIL["pending"].format(count=CP.COUNT_WORDS.get(len(pending), f"{len(pending)} checks wait")))
    if summary_missing:
        det.append(CP.DETAIL["summary_missing"])
    det.extend(read_detail())
    judged["history_hash"] = _history_hash(r, hist)
    judged["frozen"] = not pending
    st = {"profile": profile, "evaluated": now, "headline": head, "rules": R.rules, "detail": " ".join(det), "cause": cause,
          "modifiers": mods, "estimates": [], "provisional": bool(pending), "pending": pending, "explanations": expl}
    return finish(st, v1)


def _low_of_range(low: str, high: str) -> str:
    """The lower end of "low to high" as it is written before the higher end: without its unit when both ends carry
    the same one ("63.2 to 140 M", "52 to 142 Gb"), whole when they differ ("727 Mb to 1.1 Gb")."""
    lo, hi = low.rsplit(" ", 1), high.rsplit(" ", 1)
    return lo[0] if len(lo) == 2 and len(hi) == 2 and lo[1] == hi[1] else low


def judge(bundle: dict, history: Sequence[Card], *, declared: Declared | None = None, profile: dict | None = None,
          baseline_groups: Sequence[dict] = (), flow_cell: dict | None = None, in_library: bool = True) -> dict:
    """ENGINE ``judge``: the status file of one bundle. Pure: reads and writes no file."""
    declared = declared or Declared()
    r = facts_of(bundle)
    H = _prepare(r, history if in_library else (), baseline_groups, flow_cell, declared)
    use = FU.use_of(r, H, flow_cell, declared.flow_cell_use)
    limits = Limits(r.get("family"), r.get("agroup"), {k: v for k, v in (profile or {}).items() if not k.startswith("_")})
    hist = BL.History(H, (profile or {}).get("_baselines"))
    try:
        from flowcellar.readdata import api as readdata_api
        own = {k: v for k, v in (profile or {}).items() if not k.startswith("_")}       # the library's [profiles], when it has any
        read_rules = list((readdata_api.read_checks(bundle, own) if own else readdata_api.read_checks(bundle)) or [])
    except Exception:  # noqa: BLE001  read-data checks never take the run-health verdict with them
        read_rules = []
    waiting = sum(1 for i in ((bundle.get("provenance") or {}).get("inputs") or [])
                  if any(n.get("code") == "waiting" for n in i.get("notes") or []))
    st = judge_facts(r, hist, use, limits, declared=declared, in_library=in_library, read_rules=read_rules,
                     read_level=(bundle.get("tiers") or {}).get("read_level"), read_files_waiting=waiting)
    relations = B.default_relations()
    relations["flow_cell_use"] = FU.relation(use)
    sha = bundle.get("content_sha256") or B.content_sha256(bundle)
    sf = S.new_status_file(bundle.get("bundle_id") or r.get("run_key"), sha, st, relations)
    sf["engine"] = {"group": r.get("group"), "group_label": r.get("group_label"), "restart": bool(use.get("restart"))}
    return sf


def history_hash(card_: Card, history: Sequence[Card]) -> str | None:
    """What ``judge`` stores in ``status.judged.history_hash`` for this run against ``history``."""
    r = dict((card_ or {}).get("facts") or {})
    if not r.get("run_key") or not r.get("has_tier1") or r.get("run_hours") is None or r["run_hours"] < 1 \
            or r.get("run_state") in ("running", "raw_only"):
        return None
    H = _prepare(r, history, (), None, None, stored_groups=True)
    return _history_hash(r, BL.History(H))


def baselines_today(cards: Sequence[Card]) -> dict:
    """The ``baselines`` table of the index: per comparison group and metric as of today. It drives the overlay
    "against the lab today" in the viewer and is never used for a verdict."""
    import numpy as np
    metrics = [("yield_reads_at_72h", "reads_72h", FLOOR["yield"]), ("yield_reads_at_48h", "reads_48h", FLOOR["yield"]),
               ("yield_reads_at_24h", "reads_24h", FLOOR["yield"]), ("yield_estimated_bases_at_72h", "bases_72h", FLOOR["yield"]),
               ("yield_estimated_bases_at_48h", "bases_48h", FLOOR["yield"]), ("yield_estimated_bases_at_24h", "bases_24h", FLOOR["yield"]),
               ("fc_pores_first_scan", "mux0_pores", FLOOR["pores"]), ("fc_occupancy_median_24h", "occ_median_24h", FLOOR["occ"]),
               ("yield_kb_per_active_channel_hour_48h", "kb_per_active_ch_hour_48h", FLOOR["kb"]),
               ("yield_per_pore_mb_48h", "mbpp_48h", FLOOR["per_pore"]), ("fc_pore_survival_24h_frac", "pore_survival_24h", FLOOR["other"]),
               ("fc_sequencing_duty", "seq_duty", FLOOR["other"])]
    facts = [dict(c["facts"]) for c in cards or () if isinstance(c, dict) and isinstance(c.get("facts"), dict) and c["facts"].get("run_key")]
    uses = FU.annotate(facts, None)
    cfg = BL.settings(None)
    groups: dict[str, list[dict]] = {}
    for f in facts:
        f["group"] = f.get("group") or BL.group_key(f, None)
        if f.get("has_tier1") and f.get("t1_complete", True) and f["group"] and uses.get(f["run_key"], {}).get("index", 1) == 1 \
                and (f.get("run_hours") or 0) >= 24 and f.get("started") and not f.get("exclude"):
            groups.setdefault(f["group"], []).append(f)
    out = []
    now = now_z()
    for gname, rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: r["started"])[-int(cfg["window"]):]
        entry = {"group": gname, "n_runs": len(rs), "as_of": now, "window_first": rs[0]["started"][:10], "window_last": rs[-1]["started"][:10],
                 "metrics": {}}
        for mid, col, floor in metrics:
            vals = [r[col] for r in rs if r.get(col) is not None]
            if not vals:
                continue
            b = BL.robust(vals[0], vals, floor)
            if b is None:
                continue
            n = len(vals)
            trust = "trusted" if n >= cfg["trusted_runs"] else "provisional" if n >= cfg["min_runs"] else "none"
            rel = b["rel"] if b["rel"] == b["rel"] else None
            guard = bool(rel is not None and rel > cfg["spread_guard"])
            if guard and trust == "trusted":
                trust = "provisional"
            entry["metrics"][mid] = {"n": n, "median": BL.sig(b["median"], 6), "mad_sd": BL.sig(b["rel"] * abs(b["median"]), 6) if rel is not None else None,
                                     "floor": floor, "scale": BL.sig(b["scale"], 6), "trust": trust,
                                     "spread": round(rel, 3) if rel is not None else None, "spread_guard": guard}
        for name, col in (("output_gb_per_3h", "spark_output_gb_per_3h"), ("working_channels", "spark_working_channels")):
            rows = [r for r in rs if isinstance(r.get(col), list) and len(r[col]) == 32]
            if not rows:
                continue
            mat = np.array([r[col] for r in rows], float)
            live = np.array([[1.0 if 3 * i < (r.get("run_hours") or 0) else np.nan for i in range(32)] for r in rows])
            mat = mat * live
            import warnings as _w
            with np.errstate(all="ignore"), _w.catch_warnings():
                _w.simplefilter("ignore")
                med, lo, hi = np.nanmedian(mat, 0), np.nanquantile(mat, 0.25, 0), np.nanquantile(mat, 0.75, 0)

            def r3(a):
                return [None if not np.isfinite(v) else round(float(v), 3) for v in a]
            entry[f"curve_{name}"] = {"median": r3(med), "q25": r3(lo), "q75": r3(hi), "points_h": 3}
        out.append(entry)
    return {"groups": out}
