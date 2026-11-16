"""Shared by the ENGINE tests: the invariants every tier 1 bundle must satisfy beyond the schema (totals that must
agree, what a status may say), and helpers that make variations of a measured bundle for the status tests."""
from __future__ import annotations

import copy
import json
import re

from flowcellar.aggregates import derive
from flowcellar.engine import api
from flowcellar.engine.copy_pending import BANNED_WORDS
from flowcellar.model import schema


def banned(text: str) -> list[str]:
    return [w for w in BANNED_WORDS if re.search(r"(?<![a-z])" + re.escape(w) + r"(?![a-z])", text.lower())]


def status_sentences(status: dict) -> list[str]:
    return ([status.get("headline") or "", status.get("detail") or ""] + [r["message"] for r in status.get("rules", [])]
            + [m["message"] for m in status.get("modifiers", [])] + [e["message"] for e in status.get("explanations", [])])


def problems(b: dict) -> list[str]:
    """Every invariant that does not hold, in words; empty for a good bundle. ``b`` may be measured or composed."""
    out: list[str] = []

    def check(ok, msg):
        if not ok:
            out.append(msg)
    out.extend("schema: " + e for e in schema.schema_errors(b))
    M, axes = b["modules"], b["axes"]
    inputs = {i["id"] for i in b["provenance"]["inputs"]}
    check(list(b)[:4] == ["format", "schema_version", "kind", "bundle_id"], "format, schema_version, kind, bundle_id come first")
    check(b["bundle_id"] == b["identity"]["run_key"], "bundle_id equals identity.run_key")
    ch_count = b["identity"]["flow_cell"].get("channel_count")
    for mid, sec in M.items():
        check(set(sec.get("inputs", [])) <= inputs, f"{mid}: every input id exists")
        check((sec["state"] == "partial") == ("partial_reason" in sec), f"{mid}: partial_reason exactly when the state is partial")
        seen = set()
        for c in sec["cells"] + ([dict(sec["run_level"], at=["<run level>"] * len(sec["dims"]))] if sec.get("run_level") else []):
            at = tuple(c["at"])
            check(len(at) == len(sec["dims"]) and at not in seen, f"{mid}: cell {at} matches dims and is unique")
            seen.add(at)
            for name, h in (c.get("hists") or {}).items():
                check(h["axis"] in axes, f"{mid}/{name}: axis {h['axis']} is defined")
                if "w" in h:
                    check(len(h["w"]) == len(h["n"]) and "w_under" in h and "w_over" in h, f"{mid}/{name}: w is parallel to n and has both ends")
                check(h["start"] + len(h["n"]) <= axes[h["axis"]].get("n", 10 ** 9), f"{mid}/{name}: bins lie inside the axis")
            for name, h in (c.get("hists2d") or {}).items():
                check(h["x"] in axes and h["y"] in axes, f"{mid}/{name}: both axes are defined")
                check(len(h["i"]) == len(h["j"]) == len(h["n"]) and "outside" in h, f"{mid}/{name}: i, j, n parallel; outside stored")
                pairs = list(zip(h["i"], h["j"]))
                check(sorted(pairs) == pairs and len(set(pairs)) == len(pairs), f"{mid}/{name}: cells sorted and unique")
            for name, s in (c.get("series") or {}).items():
                check(s["axis"] in axes, f"{mid}/{name}: axis {s['axis']} is defined")
                check(len({len(x) for x in s["values"].values()}) == 1, f"{mid}/{name}: value arrays have equal length")
                check(set(s.get("levels", [])) <= set(s["values"]) and set(s.get("levels", [])) <= set(s.get("merge", {})),
                      f"{mid}/{name}: every level is a value with a merge rule")
            for name, k in (c.get("counts") or {}).items():
                check(len(k["cats"]) == len(k["n"]), f"{mid}/{name}: cats and n have equal length")
            for name, gr in (c.get("grids") or {}).items():
                for cn, ch in gr["channels"].items():
                    check(len(ch["data"]) == gr["rows"] * gr["cols"], f"{mid}/{name}/{cn}: one character per cell")
                check(gr.get("rows_axis") in axes and gr.get("cols_axis") in axes, f"{mid}/{name}: grid axes are defined")
            for name, cm in (c.get("channels") or {}).items():
                check(all(len(x) == ch_count * cm.get("wells", 1) for x in cm["values"].values()), f"{mid}/{name}: one value per channel (and well)")
    # ---- yield
    if "yield" in M and M["yield"]["cells"][0].get("series"):
        y, ys = M["yield"]["cells"][0], M["yield"]["summary"]
        vv = y["series"]["yield_series"]["values"]
        ratio = M["yield"]["params"]["events_to_base_ratio"]
        check(sum(vv["reads"]) == y["scalars"]["reads"] == ys["yield_reads"], "yield: series reads sum to the scalar")
        check(sum(vv["selected_events"]) == y["scalars"]["selected_events"], "yield: series events sum to the scalar")
        check(all(x == round(e * ratio) for e, x in zip(vv["selected_events"], vv["estimated_bases"])), "yield: estimated bases = round(events x ratio) per bin")
        hs = [ys.get(f"yield_estimated_bases_at_{h}h") for h in (6, 12, 24, 48, 72)]
        hs = [x for x in hs if x is not None]
        check(hs == sorted(hs) and (not hs or hs[-1] <= sum(vv["estimated_bases"])), "yield: fixed-horizon values are cumulative")
        for h in (6, 12, 24, 48, 72):
            if ys.get(f"yield_reads_at_{h}h") is not None:
                check(ys[f"yield_reads_at_{h}h"] == sum(vv["reads"][: 6 * h]), f"yield: reads at {h} h equal the first {6 * h} bins")
    # ---- fc
    if "fc" in M and "fc_channel_state_series" in M["fc"]["cells"][0].get("series", {}):
        fcell = M["fc"]["cells"][0]
        st = fcell["series"]["fc_channel_state_series"]
        check(st["unit"] == "channel_seconds" and len(st["values"]) >= 18, "fc: 18 channel states as channel-seconds")
        nb = len(next(iter(st["values"].values())))
        per_bin = [sum(s[i] for s in st["values"].values()) for i in range(nb)]
        if nb > 2 and ch_count:
            check(all(x == ch_count * 600 for x in per_bin[:-1]), "fc: every full 600 s bin sums to channels x 600 s")
            check(per_bin[-1] <= ch_count * 600, "fc: the last bin holds no more than channels x 600 s")
        check(sum(st["values"]["strand"]) == fcell["scalars"]["strand_channel_seconds"], "fc: strand series sums to its scalar")
        if "fc_pore_scan_series" in fcell["series"]:
            ps = fcell["series"]["fc_pore_scan_series"]
            v = ps["values"]
            check(v["single_pore"][0] + v["reserved_pore"][0] == M["fc"]["summary"]["fc_pores_first_scan"], "fc: first scan equals fc_pores_first_scan")
            check(v["single_pore"][-1] + v["reserved_pore"][-1] == M["fc"]["summary"]["fc_pores_last_scan"], "fc: last scan equals fc_pores_last_scan")
            check(ps.get("levels") == ["t_s"] and v["t_s"] == sorted(v["t_s"]), "fc: scan times are a level, in order")
        if "fc_temperature_series" in fcell["series"]:
            ts = fcell["series"]["fc_temperature_series"]
            check(sorted(ts["levels"]) == sorted(ts["values"]), "fc: all four temperature values are levels")
        tot = sum(per_bin)
        if tot and M["fc"]["summary"].get("fc_sequencing_duty") is not None:
            duty = sum(sum(st["values"][s]) for s in ("strand", "adapter")) / tot
            check(abs(duty - M["fc"]["summary"]["fc_sequencing_duty"]) < 5e-4, "fc: fc_sequencing_duty follows from the stored series")
    # ---- full detail
    detail = b["tiers"]["run_health"].get("detail")
    if "dyn" in M and M["dyn"]["cells"] and "yield" in M:
        reads = M["yield"]["cells"][0]["scalars"]["reads"]
        complete = M["dyn"]["state"] == "ok"
        dr = sum(c["scalars"]["reads"] for c in M["dyn"]["cells"])
        if complete:
            check(dr == reads, f"dyn: reads over all end reasons ({dr}) equal MinKNOW's read count ({reads})")
            check(M["dyn"]["coverage"]["kind"] == "full" and M["dyn"]["params"]["completeness"]["equal"] is True, "dyn: a complete pass says so")
        else:
            check(M["dyn"]["coverage"]["kind"] in ("byte_ranges", "truncated") and M["dyn"]["partial_reason"] == "partly_read",
                  "dyn: a partly read summary says how much was covered")
        for c in M["dyn"]["cells"]:
            n = c["scalars"]["reads"]
            check(derive.total(c["hists"]["dyn_duration_hist"]) == n, f"dyn {c['at']}: duration histogram holds every read")
            s = c["hists"]["dyn_speed_est_hist"]
            check(derive.total(s) + s["missing"] == n and s["missing"] == c["scalars"]["reads_under_1s"], f"dyn {c['at']}: speed histogram plus reads under 1 s holds every read")
            sv = c["series"]["dyn_end_reason_series"]["values"]
            check(sum(sv["reads"]) == n and sum(sv["duration_ms"]) == c["scalars"]["duration_ms"], f"dyn {c['at']}: hourly series sum to the scalars")
            if "hists2d" in c:
                h2 = c["hists2d"]["dyn_speed_est_series"]
                check(sum(h2["n"]) + h2["outside"] == derive.total(s), f"dyn {c['at']}: speed over time holds the reads of 1 s or longer")
        for c in M["len.est"]["cells"]:
            h = c["hists"]["len_est_hist"]
            check(derive.total(h) == c["scalars"]["reads"], f"len.est {c['at']}: histogram holds every read")
            check(abs(derive.total(h, "w") - c["scalars"]["estimated_bases"]) <= len(h["w"]) + 2, f"len.est {c['at']}: bases per bin sum to the cell's bases")
        check({tuple(c["at"]) for c in M["dyn"]["cells"]} == {tuple(c["at"]) for c in M["len.est"]["cells"]}, "dyn and len.est hold the same cells")
        fcell = M["fc"]["cells"][0]
        if complete and "fc_channel_map" in fcell.get("channels", {}):
            cm = fcell["channels"]["fc_channel_map"]["values"]
            check(sum(cm["reads"]) == fcell["scalars"]["reads_in_channel_map"] <= reads, "fc: channel map reads sum to the reads with a valid channel")
            check(sum(1 for x in cm["reads"] if x > 0) == M["fc"]["summary"]["fc_active_channels"], "fc: fc_active_channels equals the channels with a read")
            wr = fcell["channels"]["fc_channel_well_reads"]["values"]["reads"]
            check(sum(wr) + fcell["scalars"]["reads_mux_unknown"] == sum(cm["reads"]), "fc: reads per well and reads without a well add up to the map")
            mu = fcell["counts"]["fc_mux_usage"]
            check(sum(mu["n"]) + mu["other"] == dr, "fc: reads per well number add up to all reads")
        if not complete:
            check("fc_channel_map" not in fcell.get("channels", {}) and "grids" not in fcell, "fc: no channel map from a part of the reads")
    check((detail == "full") == ("dyn" in M and M["dyn"]["state"] == "ok"), "run_health.detail is full exactly with a complete read pass")
    # ---- status
    S = b.get("status")
    if S:
        R = S["rules"]

        def worst(rules):
            t = [r for r in rules if r["outcome"] == "triggered" and r.get("role") != "info"]
            return "fail" if any(r["severity"] == "fail" for r in t) else "warn" if t else "pass"
        h = S["by_tier"]["run_health"]
        t1 = [r for r in R if r.get("tier") == "run_health"]
        if h == "incomplete":
            check(S.get("cause") in ("stopped_early", "output_not_judged", "no_final_summary", "raw_only", "inputs_missing") and worst(t1) == "pass",
                  "an incomplete run names its cause and has no triggered check")
        elif h != "running":
            check(h == worst(t1), f"by_tier.run_health '{h}' equals the worst triggered check")
        check(all(r["severity"] == "none" for r in R if r["outcome"] != "triggered"), "severity is none unless a check triggered")
        check(all(not (r["outcome"] == "triggered" and r["severity"] == "fail" and r.get("role") != "key") for r in R), "only key checks can fail")
        check(len({(r["id"], json.dumps(r.get("scope"))) for r in R}) == len(R), "check ids are unique per scope")
        check(all((r.get("threshold") or {}).get("applied") in ("fixed", "baseline", "run_config", "rule") for r in t1), "every rule names the limits that applied")
        summ = {k: x for sec in M.values() for k, x in (sec.get("summary") or {}).items()}
        scal = {k: x for sec in M.values() if not sec["dims"] for c in sec["cells"] for k, x in (c.get("scalars") or {}).items()}
        for r in t1:
            bl = (r.get("threshold") or {}).get("baseline")
            if bl and isinstance(r["value"], (int, float)):
                check(abs((r["value"] - bl["median"]) / bl["scale"] - bl["z"]) < 0.02, f"{r['id']}: z follows from value, median and scale")
                if bl["trust"] != "trusted":
                    check(not (r["outcome"] == "triggered" and r["severity"] == "fail" and r["threshold"]["applied"] == "baseline"),
                          f"{r['id']}: an untrusted history cannot fail")
            if r["outcome"] != "not_evaluated":
                m = r.get("metric")
                if m in summ and r["id"] != "fc_pores_first_scan":
                    v = summ[m]
                    check(r["value"] == v or (isinstance(v, float) and abs(r["value"] - v) < 0.006), f"{r['id']}: value equals the summary {m}")
                elif m in scal:
                    check(r["value"] == scal[m], f"{r['id']}: value equals the scalar {m}")
                elif m not in summ:
                    check(False, f"{r['id']}: metric {m} is in no section summary or run-level scalar")
        check(bool(S.get("provisional")) == bool(S.get("pending")), "provisional exactly while checks are pending")
        check(not (set(S.get("pending", [])) & {r["id"] for r in R}), "a pending check is not also evaluated")
        check(not (S["judged"]["frozen"] and S.get("provisional")), "a provisional status is not frozen")
        started = b["identity"]["protocol_run"].get("started") or ""
        for r in R + S.get("explanations", []):
            bl = (r.get("threshold") or {}).get("baseline") or r.get("baseline")
            if bl:
                check(bl["as_of"][:19] == started[:19], f"{r['id']}: history is taken as of the run's start")
                check(bl["trust"] in ("trusted", "provisional") and not (bl.get("spread_guard") and bl["trust"] == "trusted"),
                      f"{r['id']}: trust is consistent with the spread guard")
                check(bl["n"] >= 5 and (bl["trust"] != "trusted" or bl["n"] >= 8), f"{r['id']}: trust is consistent with the number of runs")
        for e in S.get("explanations", []):
            z = e["baseline"]["z"]
            check(e["direction"] == ("low" if z <= -1.5 else "high" if z >= 1.5 else "typical"), f"explanation {e['id']}: direction follows from z")
        outcome_rules = [r for r in R if re.match(r"yield_(reads|estimated_bases)_at_\d+h$", r["id"]) and r["role"] == "key"]
        reused = any(m["code"] == "reused_flow_cell" for m in S.get("modifiers", []))
        check(h != "pass" or bool(outcome_rules) or reused, "Run health is Pass only with an output check")
        for t in status_sentences(S):
            check(not banned(t), f"a stored sentence uses a word the interface never shows {banned(t)}: {t[:90]!r}")
            check(not re.search(r"\d%", t), f"percent sign without a space before it: {t[:60]!r}")
            check("{" not in t and "}" not in t, f"an unfilled placeholder: {t[:90]!r}")
        check(len(S.get("headline") or "") <= 140, "the headline is at most 140 characters")
        lib = copy.copy(b)
        for k in ("status", "user", "relations"):
            lib.pop(k, None)
        check(not schema.schema_errors(lib), "the schema accepts the bundle without status, user and relations")
    return out


def clone(bundle: dict, *, key: str | None = None, flow_cell: str | None = None, started: str | None = None, **summaries) -> dict:
    """A copy of a measured bundle with another identity and other summary numbers (``section__key=value``)."""
    b = copy.deepcopy(bundle)
    b.pop("content_sha256", None)
    if flow_cell:
        b["identity"]["flow_cell"]["id"] = flow_cell
    if key:
        b["identity"]["run_key"] = b["bundle_id"] = key
        b["identity"]["protocol_run"]["short_id"] = key.split("_")[1]
    if started:
        b["identity"]["protocol_run"]["started"] = started
    for name, value in summaries.items():
        section, _, k = name.partition("__")
        b["modules"][section.replace("len_est", "len.est")]["summary"][k] = value
    return b


def history_of(base: dict, n: int, *, first_day: int = 1, month: str = "2026-01", **per_run) -> list[dict]:
    """``n`` cards of earlier first-use runs of the same group as ``base``, on ``n`` different days. ``per_run`` maps
    a summary (``section__key``) to a list of ``n`` values."""
    cards = []
    for i in range(n):
        vals = {k: v[i] for k, v in per_run.items()}
        b = clone(base, key=f"TST{20000 + i:05d}_{i + 1:08x}", flow_cell=f"TST{20000 + i:05d}",
                  started=f"{month}-{first_day + i:02d}T09:00:00Z", **vals)
        cards.append(api.card(b, None))
    return cards


def judged(bundle: dict, history=(), **kw) -> dict:
    """The status of a bundle; ``rules`` by id are added under ``_rules`` for the tests."""
    sf = api.judge(bundle, history, in_library=kw.pop("in_library", bool(history)), **kw)
    st = sf["status"]
    st["_rules"] = {r["id"]: r for r in st["rules"]}
    st["_file"] = sf
    return st
