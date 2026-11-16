"""A synthetic library for the viewer's tests: composed bundles, status files and an index in the layout
the library writes, written straight from a few numbers per run.

Nothing here reads a run folder and nothing here is the engine: the numbers are made up so that every state the
viewer must draw exists once (full and quick detail, raw data only, a run that stopped after minutes, a re-used
flow cell with a restart, complete and partial read data, a root that is offline, a project in two places, a run
without a project, and names that try to be markup). Every name is invented; flow cells are TST00001 onwards.

    python tests/viewer/viewer_synth.py OUT_DIR            write the library (bundles/, status/, index/, ...)
    python tests/viewer/viewer_synth.py OUT_DIR --empty    a library without runs
"""
from __future__ import annotations

import datetime as dt
import json
import math
import os
import re
import sys
import warnings
from dataclasses import dataclass, field

import numpy as np

from flowcellar import SCHEMA_VERSION, __version__
from flowcellar.aggregates import build as agg
from flowcellar.model import axes as ax
from flowcellar.model import bundle as bundle_model
from flowcellar.model import status as status_model
from flowcellar.util import atomic

HOSTILE_RUN_NAME = '<script>window.__pwned=1</script><img src=x onerror="window.__pwned=2">'
HOSTILE_SAMPLE = '"><svg onload=window.__pwned=3>&amp;%%JS%%</script><script>window.__pwned=4</script>'
HOSTILE_PROJECT = "</title><script>window.__pwned=5</script>"
HOSTILE_MESSAGE = "Disk <b>full</b> </script><script>window.__pwned=6</script> %%PAYLOAD%%"

END_REASONS = [("signal_positive", 0.955), ("unblock_mux_change", 0.028), ("signal_negative", 0.012), ("mux_change", 0.005)]
STATES = ["strand", "adapter", "pore", "unavailable", "unblocking", "no_pore", "zero", "saturated", "multiple",
          "unknown_positive", "unknown_negative", "disabled", "locked", "membrane", "unclassified",
          "unclassified_following_reset", "pending_manual_reset", "pending_mux_change"]
HORIZONS = (6, 12, 24, 48, 72)


@dataclass
class Spec:
    fc: str
    short: str
    started: str                                  # RFC 3339 UTC
    hours: float = 72.0
    planned_h: float = 72.0
    family: str = "promethion"                    # promethion | minion
    assay: str = "cdna"                           # cdna | gdna | drna
    detail: str | None = "full"                   # full | quick | None (raw data only)
    read_data: str = "absent"                     # absent | complete | partial
    partial_reason: str | None = None             # samples_missing | partly_read | waiting
    samples: tuple[str, ...] = ()
    pores: int = 7200
    occ: float = 0.9
    half_life_h: float = 70.0
    experiment: str = "demo_experiment"
    run_name: str = "demo_run"
    project: str | None = "Demo Alpha"
    sample_project: str | None = None             # a second project, through the samples
    roots: tuple[str, ...] = ("sequencer",)
    position: str = "TSTP2-00001-A"
    use: tuple[int, int] = (1, 1)                 # loading index, loadings known
    restart_of: str | None = None
    peers: int = 12                               # comparable runs before this one
    run_state: str = "finished"
    stop: str = "PROTOCOL_COMPLETED"
    source_state: str = "ok"
    error_message: str | None = None
    seed: int = 1
    extra: dict = field(default_factory=dict)

    @property
    def key(self) -> str:
        return f"{self.fc}_{self.short}"

    @property
    def channels(self) -> int:
        return 3000 if self.family == "promethion" else 512


def _iso(t: dt.datetime) -> str:
    return t.strftime("%Y-%m-%dT%H:%M:%SZ")


def _t(s: str) -> dt.datetime:
    return dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ")


def _lognormal_hist(axis: str, rng, n_reads: int, median: float, sigma: float, weight: bool = True):
    """An internal histogram (under, bins, over) of n_reads log-normal values on a registry axis, and its weights."""
    edges = ax.EDGES[axis]
    if n_reads <= 0:
        z = np.zeros(len(edges) + 1, np.int64)
        return z, z.copy()
    mids = np.sqrt(edges[:-1] * edges[1:]) if edges[0] > 0 else (edges[:-1] + edges[1:]) / 2
    with np.errstate(divide="ignore"):
        lx = np.log(np.maximum(mids, 1e-9))
    p = np.exp(-0.5 * ((lx - math.log(median)) / sigma) ** 2)
    p /= p.sum()
    n = np.array(agg.largest_remainder(list(p * n_reads), n_reads), np.int64)
    full = np.concatenate([[0], n, [0]])
    w = np.concatenate([[0], np.rint(n * mids).astype(np.int64), [0]])
    return full, (w if weight else None)


def _normal_hist(axis: str, n_reads: int, mean: float, sd: float):
    edges = ax.EDGES[axis]
    mids = (edges[:-1] + edges[1:]) / 2
    p = np.exp(-0.5 * ((mids - mean) / sd) ** 2)
    p /= p.sum()
    n = np.array(agg.largest_remainder(list(p * n_reads), n_reads), np.int64) if n_reads > 0 else np.zeros(len(mids), np.int64)
    return np.concatenate([[0], n, [0]])


def _quantile(edges, counts, q, log=False):
    total = counts.sum()
    if total <= 0:
        return None
    acc, want = 0.0, q * total
    for i, c in enumerate(counts):
        if acc + c >= want and c > 0:
            f = (want - acc) / c
            if log:
                return float(10 ** (math.log10(edges[i]) + f * (math.log10(edges[i + 1]) - math.log10(edges[i]))))
            return float(edges[i] + f * (edges[i + 1] - edges[i]))
        acc += c
    return float(edges[-1])


def _n50(edges, w):
    total = w.sum()
    if total <= 0:
        return None
    acc = 0.0
    for i in range(len(w) - 1, -1, -1):
        if acc + w[i] >= total / 2 and w[i] > 0:
            f = (total / 2 - acc) / w[i]
            return float(10 ** (math.log10(edges[i + 1]) - f * (math.log10(edges[i + 1]) - math.log10(edges[i]))))
        acc += w[i]
    return float(edges[0])


def _fmt_count(x):
    return f"{x / 1e6:.1f} M" if x >= 1e6 else f"{int(round(x)):,}"


def _fmt_bases(x):
    return f"{x / 1e9:.0f} Gb" if x >= 1e11 else f"{x / 1e9:.1f} Gb" if x >= 1e9 else f"{x / 1e6:.0f} Mb"


# ---------------------------------------------------------------------------------------------------------
# one bundle
# ---------------------------------------------------------------------------------------------------------

def make_bundle(spec: Spec, library: dict | None = None) -> dict:
    """A composed bundle (measured content, status, user, relations) for one synthetic run."""
    rng = np.random.default_rng(spec.seed)
    lib = library or {"id": "lib_0123456789ab", "name": "Demo library"}
    n_ch = spec.channels
    product = "FLO-PRO114M" if spec.family == "promethion" else "FLO-MIN114"
    kit = {"cdna": "SQK-PCB114-24", "gdna": "SQK-NBD114-24", "drna": "SQK-RNA004"}[spec.assay]
    device = "P2_SOLO" if spec.family == "promethion" else "MINION_MK1D"
    t0 = _t(spec.started)
    raw_only = spec.detail is None
    pid = f"{spec.short}-1111-4222-8333-444455556666"
    identity = {
        "run_key": spec.key,
        "flow_cell": {"id": spec.fc, "user_specified_id": spec.fc, "product_code": None if raw_only else product,
                      "channel_count": n_ch, "wells_per_channel": None if raw_only else (4 if spec.family == "promethion" else 4),
                      "channel_layout": "promethion_3000" if spec.family == "promethion" else "minion_512", "family": spec.family},
        "protocol_run": {"id": None if raw_only else pid, "short_id": spec.short, "started": spec.started,
                         "ended": None if raw_only or spec.run_state == "running" else _iso(t0 + dt.timedelta(hours=spec.hours, seconds=120)),
                         "time_origin": None if raw_only else _iso(t0 + dt.timedelta(seconds=110)),
                         "time_origin_source": None if raw_only else "minknow_report",
                         "run_state": "raw_only" if raw_only else spec.run_state},
        "acquisitions": [], "device": {"type": None if raw_only else device, "id": spec.position.rsplit("-", 1)[0], "position": spec.position},
        "kit": {}, "assay": {"type": "unknown", "group": "unknown", "molecule": "unknown", "source": "unknown", "confidence": "low", "evidence": []},
        "software": {}, "basecalling": [],
    }
    rel_dir = "/".join(x for x in ((spec.project or "").replace(" ", "_") or None, spec.experiment, spec.run_name[:20].replace("/", "_").replace("<", "_").replace(">", "_"),
                                    f"20260101_1200_{spec.position}_{spec.key}") if x)
    root0 = spec.roots[0]
    run_dir = f"/data/{root0}/{rel_dir}"
    b = bundle_model.new_bundle(identity)
    b["provenance"].update({"created": "2026-11-13T12:00:00Z", "updated": "2026-11-13T12:00:00Z", "revision": 1,
                            "library": {"id": lib["id"], "name": lib["name"]},
                            "run_dir": {"path": run_dir, "root": root0, "relpath": rel_dir, "layout": "raw_only" if raw_only else "standard", "alternates": []},
                            "profile": "standard", "partials": [], "warnings": []})
    relations = bundle_model.default_relations()
    relations["flow_cell_use"].update({"index": spec.use[0], "of": spec.use[1], "source": "library" if spec.use[1] > 1 else "unknown",
                                       "restart_of": spec.restart_of})
    projects = []
    for name in (spec.project, spec.sample_project):
        if name:
            projects.append({"id": project_id(name), "name": name, "source": "path", "roots": list(spec.roots)})
    relations["projects"] = projects
    relations["paths"] = [{"root": r, "relpath": rel_dir} for r in spec.roots]

    if raw_only:
        b["provenance"]["warnings"].append({"level": "warning", "code": "raw_only",
                                            "message": "Only raw data was found in this folder (sub-folders: pod5); identity comes from the folder name."})
        b["tiers"] = {"run_health": {"state": "absent", "detail": None, "missing_roles": ["minknow.final_summary", "minknow.report_json"]},
                      "read_level": {"state": "absent"}}
        status = {"profile": {"id": "unknown", "version": "0.1.1", "source": "builtin"}, "evaluated": "2026-11-13T12:00:00Z",
                  "verdict": "incomplete", "by_tier": {"run_health": "incomplete", "read_level": "not_available"}, "cause": "raw_only",
                  "headline": "Only raw data was found for this run: no report, so the flow cell cannot be judged.", "detail": "",
                  "rules": [], "modifiers": [], "estimates": [], "provisional": False, "pending": [], "explanations": [],
                  "judged": {"at": "2026-11-13T12:00:00Z", "as_of": spec.started, "history_hash": None, "library": True, "frozen": True}}
        return bundle_model.compose(b, status=status, relations=relations)

    # ---- identity of a run with a MinKNOW report ----
    full = spec.detail == "full"
    speed = 400 if spec.assay != "drna" else 130
    identity["flow_cell"]["asic_id"] = "00A1B2C3"
    identity["protocol_run"].update({"protocol": f"sequencing/sequencing_{product}:{kit}", "experiment_id": spec.experiment, "sample_id": spec.run_name,
                                     "utc_offset_min": 0, "state": spec.stop, "planned_duration_s": spec.planned_h * 3600.0})
    identity["device"].update({"host": "fakehost-01"})
    identity["kit"] = {"sequencing_kit": kit, "barcoding_kits": [], "speed_bps": speed, "kit_category": ["DNA"]}
    identity["assay"] = {"type": {"cdna": "cdna_pcr", "gdna": "gdna", "drna": "direct_rna"}[spec.assay], "group": spec.assay,
                         "chemistry": "r10.4.1" if spec.assay != "drna" else "rna004", "molecule": "rna" if spec.assay == "drna" else "dna",
                         "barcoded": spec.assay != "drna", "source": "run_files", "confidence": "high", "evidence": [f"kit {kit}", f"flow cell product {product}"]}
    identity["software"] = {"minknow": "6.10.12", "bream": "8.10.13", "distribution": "26.01.15"}
    bins = max(1, int(round(spec.hours * 6)))
    t_mid = (np.arange(bins) + 0.5) / 6.0

    # ---- pore scans ----
    period_s = 5400 + 290
    scan_t = [276 + i * period_s for i in range(int((spec.hours * 3600 - 276) // period_s) + 1)] if spec.hours * 3600 > 300 else [min(200, int(spec.hours * 3600 * 0.8))]
    wells = n_ch * 4
    pores_t = [max(0, int(round(spec.pores * math.exp(-math.log(2) * (s / 3600) / spec.half_life_h)))) for s in scan_t]
    single = [int(round(p * 0.78)) for p in pores_t]
    reserved = [p - s for p, s in zip(pores_t, single)]
    scan = {"single_pore": single, "reserved_pore": reserved,
            "unavailable": [int(wells * 0.03)] * len(scan_t), "saturated": [int(wells * 0.02 + i * 3) for i in range(len(scan_t))],
            "zero": [0] * len(scan_t), "multiple": [int(wells * 0.01)] * len(scan_t), "other": [0] * len(scan_t)}
    for i in range(len(scan_t)):
        scan["zero"][i] = wells - sum(scan[k][i] for k in scan if k != "zero")
    scan_series = agg.series("mux_scan_index", {**scan, "t_s": scan_t}, unit="wells", levels={"t_s": "same"})

    # ---- channel states: channel-seconds per 10 minutes ----
    work0 = min(0.92, spec.pores / wells * 1.35) * n_ch
    working = work0 * np.exp(-math.log(2) * t_mid / (spec.half_life_h * 1.05))
    occ = np.clip(spec.occ - 0.02 * t_mid / max(spec.hours, 1) + rng.normal(0, 0.004, bins), 0.02, 0.995)
    occ[:1] *= 0.8
    scan_bins = {int(s // 600) for s in scan_t if s // 600 < bins}
    states = {k: np.zeros(bins) for k in STATES}
    for i in range(bins):
        w, o, f = working[i], occ[i], (0.6 if i in scan_bins else 1.0)
        states["strand"][i] = w * o * 0.92 * 600 * f
        states["adapter"][i] = w * o * 0.08 * 600 * f
        states["pore"][i] = w * (1 - o) * 600 * f
        states["unavailable"][i] = w * (0.015 + 0.04 * t_mid[i] / max(spec.hours, 1)) * 600
        states["unblocking"][i] = n_ch * 0.005 * 600
        states["unclassified"][i] = n_ch * 0.006 * 600
        if i in scan_bins:
            lost = w * 600 * 0.4
            for k, share in (("zero", 0.4), ("saturated", 0.1), ("unknown_negative", 0.15), ("locked", 0.2), ("unclassified", 0.15)):
                states[k][i] += lost * share
    total = n_ch * 600
    st_int = {k: np.floor(v).astype(np.int64) for k, v in states.items()}
    used = sum(st_int.values())
    st_int["no_pore"] = np.maximum(0, total - (used - st_int["no_pore"]))
    over = sum(st_int.values()) - total
    st_int["pore"] = np.maximum(0, st_int["pore"] - np.maximum(0, over))
    state_series = agg.series("run_time_600s", st_int, unit="channel_seconds", n_bins=bins)

    # ---- output ----
    mean_len = {"cdna": 1150.0, "gdna": 6200.0, "drna": 900.0}[spec.assay]
    est_speed = speed * 0.92
    est_bases = np.rint(st_int["strand"] * est_speed).astype(np.int64)
    reads_bin = np.rint(est_bases / mean_len).astype(np.int64)
    events = np.rint(est_bases / 0.6).astype(np.int64)
    n_reads, n_bases = int(reads_bin.sum()), int(est_bases.sum())
    cum_r, cum_b = np.cumsum(reads_bin), np.cumsum(est_bases)
    y_sum = {"yield_reads": n_reads, "yield_estimated_bases": n_bases}
    for h in HORIZONS:
        if spec.hours >= h - 1e-6:
            y_sum[f"yield_reads_at_{h}h"] = int(cum_r[h * 6 - 1])
            y_sum[f"yield_estimated_bases_at_{h}h"] = int(cum_b[h * 6 - 1])
            y_sum[f"yield_per_pore_mb_{h}h"] = round(float(cum_b[h * 6 - 1]) / max(spec.pores, 1) / 1e6, 2)
    half = int(np.searchsorted(cum_b, n_bases / 2)) if n_bases else 0
    y_sum.update({"yield_half_time_h": round((half + 1) / 6, 2), "yield_per_pore_mb": round(n_bases / max(spec.pores, 1) / 1e6, 2),
                  "yield_rate_peak_gb_per_h": round(float(est_bases.max()) * 6 / 1e9, 3) if bins else 0.0,
                  "yield_first_24h_frac": round(float(cum_b[min(bins, 144) - 1]) / max(n_bases, 1), 4)})
    b["modules"]["yield"] = {
        "version": "0.1.0", "tier": "run_health", "state": "ok", "inputs": ["in3"], "dims": [],
        "cells": [{"at": [], "scalars": {"reads": n_reads, "selected_events": int(events.sum()), "selected_raw_samples": int(events.sum()) * 8},
                   "series": {"yield_series": agg.series("run_time_600s", {"reads": reads_bin, "selected_events": events, "estimated_bases": est_bases}, n_bins=bins)}}],
        "summary": y_sum, "params": {"events_to_base_ratio": 0.6},
    }

    # ---- run ----
    msgs = [[_iso(t0 + dt.timedelta(seconds=20)), "info", "disk_space_info_user", "Disk /data has 4000 GB space remaining"],
            [_iso(t0 + dt.timedelta(seconds=130)), "info", "sequencing.start", "Starting sequencing procedure"]]
    if spec.hours > 30:
        msgs.append([_iso(t0 + dt.timedelta(hours=26)), "warning", "disk_space_warning", "Disk space is getting low on the output drive"])
    if spec.error_message:
        msgs.append([_iso(t0 + dt.timedelta(hours=min(spec.hours * 0.9, 20))), "error", "device_error", spec.error_message])
    n_warn = sum(1 for m in msgs if m[1] == "warning")
    n_err = sum(1 for m in msgs if m[1] == "error")
    b["modules"]["run"] = {
        "version": "0.1.0", "tier": "run_health", "state": "ok", "inputs": ["in1", "in2"], "dims": [],
        "cells": [{"at": [], "scalars": {"messages_info": len(msgs) - n_warn - n_err, "messages_warning": n_warn, "messages_error": n_err, "messages_warning_computer": 0,
                                           "files_in_fallback": 0, "raw_bytes_completed": int(n_reads * 31000), "raw_bytes_failed": 0, "raw_bytes_backlog": 0, "pod5_files": max(1, n_reads // 200000)},
                   "tables": {"run_messages": {"columns": [{"id": "time", "type": "datetime"}, {"id": "severity", "type": "str"}, {"id": "identifier", "type": "str"}, {"id": "user_message", "type": "str"}],
                                               "key": ["time"], "merge": "concat", "rows": msgs, "total_rows": len(msgs), "truncated": False}}}],
        "summary": {"run_duration_h": round(spec.hours, 3), "run_completion_frac": round(min(1.0, spec.hours / spec.planned_h), 4), "run_stop_cause": spec.stop,
                    "run_stop_reason": "STOPPED_PROTOCOL_ENDED" if spec.stop == "PROTOCOL_COMPLETED" else "STOPPED_USER_REQUESTED",
                    "run_messages_error": n_err, "run_messages_warning": n_warn},
        "notes": [{"level": "info", "code": "mux_scan_messages_moved", "message": f"{2 * len(scan_t)} messages are pore-scan start and result messages; they feed the pore-scan series and are not repeated here."}],
    }

    # ---- flow cell ----
    busy = st_int["strand"] + st_int["adapter"]
    workt = busy + st_int["pore"]
    with np.errstate(divide="ignore", invalid="ignore"):
        occ_bin = np.where(workt > 0, busy / np.maximum(workt, 1), np.nan)
    first_hour = occ_bin[3:9] if bins >= 9 else occ_bin
    recov = st_int["unavailable"] + st_int["unblocking"]
    third = bins - bins // 3
    n_temp = max(1, int(spec.hours // 2))
    temp = {"flowcell_c": [round(34 + float(rng.normal(0, 0.03)), 2) for _ in range(n_temp)], "chamber_c": [round(33.2 + float(rng.normal(0, 0.05)), 2) for _ in range(n_temp)],
            "target_min_c": [33.5] * n_temp, "target_max_c": [34.5] * n_temp}
    volt = [float(-170 - 5 * (i // 4)) for i in range(len(scan_t))]
    fc_sum = {
        "fc_pores_first_scan": pores_t[0], "fc_pores_last_scan": pores_t[-1], "fc_pore_half_life_h": round(spec.half_life_h, 1),
        "fc_channels_with_pore_first_scan": int(min(n_ch, round(work0 * 1.03))), "fc_channels_with_pore_first_scan_frac": round(min(1.0, work0 * 1.03 / n_ch), 4),
        "fc_occupancy_first_hour": round(float(np.nanmean(first_hour)), 4), "fc_occupancy_median": round(float(np.nanmedian(occ_bin)), 4),
        "fc_occupancy_median_24h": round(float(np.nanmedian(occ_bin[:144])), 4), "fc_sequencing_duty": round(float(busy.sum()) / (total * bins), 4),
        "fc_adapter_frac": round(float(st_int["adapter"].sum()) / max(1, int(busy.sum())), 4),
        "fc_recovering_frac": round(float(recov.sum()) / max(1, int((workt + recov).sum())), 4),
        "fc_recovering_last_third_frac": round(float(recov[third:].sum()) / max(1, int((workt + recov)[third:].sum())), 4),
        "fc_temperature_max_dev_c": 0.0, "fc_active_channels_peak": int(round(float(workt.max()) / 600)),
        "fc_bias_voltage_start_mv": volt[0], "fc_bias_voltage_end_mv": volt[-1],
    }
    if spec.hours >= 24:
        fc_sum["fc_pore_survival_24h_frac"] = round(math.exp(-math.log(2) * 24 / spec.half_life_h), 4)
    fc_cell = {"at": [], "scalars": {"strand_channel_seconds": int(st_int["strand"].sum()), "adapter_channel_seconds": int(st_int["adapter"].sum()),
                                      "pore_channel_seconds": int(st_int["pore"].sum()), "total_channel_seconds": total * bins, "pore_scans": len(scan_t), "temperature_samples": n_temp},
               "series": {"fc_channel_state_series": state_series, "fc_pore_scan_series": scan_series,
                          "fc_temperature_series": agg.series("run_time_7200s", temp, unit="degC", levels={k: "mean" for k in temp})}}
    fc_params = {"sample_rate_hz": 5000, "mux_scan_period_h": 1.5, "occupancy_first_hour_window_min": [30, 90], "temperature_sample_period_min": 120,
                 "target_speed_min_bps": speed * 0.875, "target_speed_max_bps": speed * 1.125, "temperature_stop_sample_dropped": True}
    fc_notes = [{"level": "info", "code": "temperature_stop_sample_dropped", "message": "The temperature sample taken as the run stopped is left out."}]
    ch_reads = None
    if full:
        fc_cell["series"]["fc_bias_voltage_series"] = agg.series("mux_scan_index", {"voltage_mv": volt, "t_s": scan_t}, unit="mV", levels={"voltage_mv": "same", "t_s": "same"})
        alive = rng.random(n_ch) < min(0.97, work0 / n_ch + 0.04)
        alive[int(n_ch * 0.30): int(n_ch * 0.30) + n_ch // 40] = False               # a patch that never worked
        weights = rng.gamma(4.0, 1.0, n_ch) * alive
        weights /= weights.sum()
        ch_reads = np.array(agg.largest_remainder(list(weights * n_reads), n_reads), np.int64)
        ch_bases = np.rint(ch_reads * (n_bases / max(n_reads, 1))).astype(np.int64)
        ch_seq = np.rint(weights * float(st_int["strand"].sum())).astype(np.int64)
        good = np.where(alive, rng.integers(1, 5, n_ch), 0)
        last_scan = np.where(alive, rng.integers(max(1, len(scan_t) // 2), len(scan_t) + 1, n_ch), 0)
        layout = identity["flow_cell"]["channel_layout"]
        fc_cell["channels"] = {
            "fc_channel_map": {"layout": layout, "values": {"reads": ch_reads.tolist(), "estimated_bases": ch_bases.tolist(), "sequencing_seconds": ch_seq.tolist()}},
            "fc_well_map": {"layout": layout, "values": {"good_wells_first_scan": good.tolist(), "last_good_scan": last_scan.tolist()},
                            "merge": {"good_wells_first_scan": "first", "last_good_scan": "last"}},
        }
        cols = max(1, int(math.ceil(spec.hours / 2)))
        decay = np.exp(-math.log(2) * (np.arange(cols) * 2 + 1) / (spec.half_life_h * 1.05))
        sec = np.outer(weights * n_ch, decay) * 7200 * spec.occ * 0.8
        sec[:, :] = np.minimum(sec, 7200)
        lv = np.rint(np.sqrt(sec) / (math.sqrt(7200) / 63)).astype(np.int64)
        fc_cell["grids"] = {"fc_channel_time_grid": {"rows": n_ch, "cols": cols, "rows_axis": "channel", "cols_axis": "run_time_7200s", "layout": {"kind": "matrix"},
                                                    "channels": {"strand_seconds": {"enc": "u6", "data": agg.u6_encode(lv), "scale": round(math.sqrt(7200) / 63, 5),
                                                                                    "transform": "sqrt", "levels": 64, "unit": "s", "nodata": False}}}}
        mux_n = agg.largest_remainder([n_reads * s for s in (0.27, 0.25, 0.24, 0.24)], n_reads)
        fc_cell["counts"] = {"fc_mux_usage": {**agg.counts(["1", "2", "3", "4"], mux_n, w=[int(v * n_bases / max(n_reads, 1)) for v in mux_n], role="categories"), "w_name": "estimated_bases"}}
        mat = np.zeros((len(scan_t), len(ax.EDGES["current_2pa"]) - 1), np.int64)
        for i in range(len(scan_t)):
            full_h = _normal_hist("current_2pa", pores_t[i], 190 - i * 0.6, 9)
            mat[i] = full_h[1:-1]
        fc_cell["hists2d"] = {"fc_open_pore_current": agg.hist2d("mux_scan_index", "current_2pa", mat, outside=3)}
        fc_sum.update({"fc_active_channels": int(alive.sum()), "fc_active_channels_frac": round(float(alive.sum()) / n_ch, 4)})
    b["modules"]["fc"] = {"version": "0.1.0", "tier": "run_health", "state": "ok" if full else "partial", "inputs": ["in2", "in4"], "dims": [],
                          "cells": [fc_cell], "summary": fc_sum, "params": fc_params, "notes": fc_notes}
    if not full:
        b["modules"]["fc"]["partial_reason"] = "inputs_missing"

    # ---- reads through the pore ----
    dyn_sum: dict = {}
    le_sum: dict = {}
    if full:
        dyn_cells, le_cells = [], []
        hours_n = int(math.ceil(spec.hours))
        per_h = np.add.reduceat(reads_bin, np.arange(0, bins, 6))
        per_h_b = np.add.reduceat(est_bases, np.arange(0, bins, 6))
        len_axis, tot_n, tot_w = "read_length_log40", None, None
        for reason, share in END_REASONS:
            r_n = int(round(n_reads * share)) if reason != "signal_positive" else n_reads - sum(int(round(n_reads * s)) for _, s in END_REASONS[1:])
            r_bases = int(round(n_bases * share))
            med = mean_len * (0.85 if reason == "signal_positive" else 0.35)
            hn, hw = _lognormal_hist(len_axis, rng, r_n, med, 0.75)
            dur = _lognormal_hist("duration_log40", rng, r_n, med / est_speed, 0.8, weight=False)[0]
            cell = {"at": [reason], "scalars": {"reads": r_n, "minknow_events": int(r_bases / 0.6), "estimated_bases": int(hw.sum()), "duration_ms": int(r_n * med / est_speed * 1000),
                                                "reads_under_1s": int(r_n * 0.004)},
                    "hists": {"dyn_duration_hist": agg.hist("duration_log40", dur)},
                    "series": {"dyn_end_reason_series": agg.series("run_time_3600s", {"reads": np.rint(per_h * share).astype(np.int64), "estimated_bases": np.rint(per_h_b * share).astype(np.int64)}, n_bins=hours_n)}}
            if reason == "signal_positive":
                sp = _normal_hist("speed_5bps", r_n - int(r_n * 0.004), est_speed, 22)
                cell["hists"]["dyn_speed_est_hist"] = agg.hist("speed_5bps", sp, missing=int(r_n * 0.004))
                mat = np.zeros((hours_n, len(ax.EDGES["speed_5bps"]) - 1), np.int64)
                for h in range(hours_n):
                    mat[h] = _normal_hist("speed_5bps", int(per_h[h] * share) if h < len(per_h) else 0, est_speed * (1 - 0.03 * h / max(hours_n, 1)), 22)[1:-1]
                cell["hists2d"] = {"dyn_speed_est_series": agg.hist2d("run_time_3600s", "speed_5bps", mat, outside=0)}
            dyn_cells.append(cell)
            le_cells.append({"at": [reason], "scalars": {"reads": r_n, "estimated_bases": int(hw.sum())},
                             "hists": {"len_est_hist": agg.hist(len_axis, hn, hw, w_name="estimated_bases")}})
            tot_n = hn if tot_n is None else tot_n + hn
            tot_w = hw if tot_w is None else tot_w + hw
        edges = ax.EDGES[len_axis]
        dyn_sum = {"dyn_speed_est_median": round(est_speed, 1), "dyn_speed_drop_frac": 0.03, "dyn_signal_positive_frac": END_REASONS[0][1], "dyn_unblock_frac": END_REASONS[1][1],
                   "dyn_signal_negative_frac": END_REASONS[2][1], "dyn_duration_median_s": round(mean_len * 0.8 / est_speed, 2), "dyn_reads": n_reads}
        le_sum = {"len_est_n50": round(_n50(edges, tot_w[1:-1])), "len_est_median": round(_quantile(edges, tot_n[1:-1], 0.5, log=True)), "len_est_mean": round(n_bases / max(n_reads, 1), 1)}
        cov = {"kind": "full", "records": n_reads, "records_total": n_reads, "fraction": 1.0}
        b["modules"]["dyn"] = {"version": "0.1.0", "tier": "run_health", "state": "ok", "inputs": ["in5"], "dims": ["end_reason"], "cells": dyn_cells, "summary": dyn_sum,
                               "run_level": {"counts": {"dyn_pore_type": agg.counts(["not_set"], [n_reads], role="categories")}}, "coverage": cov,
                               "params": {"events_to_base_ratio": 0.6}}
        b["modules"]["len.est"] = {"version": "0.1.0", "tier": "run_health", "state": "ok", "inputs": ["in5"], "dims": ["end_reason"], "cells": le_cells, "summary": le_sum, "coverage": cov}
    b["tiers"]["run_health"] = {"state": "complete" if full else "partial", "detail": spec.detail, "missing_roles": [] if full else ["minknow.sequencing_summary"]}
    b["axes"] = {k: {a: v for a, v in ax.AXES[k].items() if a != "note"} for k in ("run_time_600s", "run_time_3600s", "run_time_7200s", "mux_scan_index", "channel", "read_length_log40", "read_length_log20",
                                           "qscore_0p25", "qscore_1", "duration_log40", "speed_5bps", "identity_q_0p1", "mapq_1", "current_2pa")}
    inputs = [("in1", "minknow.final_summary", 840), ("in2", "minknow.report_json", 2_400_000), ("in3", "minknow.throughput", 320_000), ("in4", "minknow.pore_activity", 2_300_000)]
    if full:
        inputs += [("in5", "minknow.sequencing_summary", n_reads * 290), ("in6", "minknow.pore_scan_data", 61_000_000)]
    b["provenance"]["inputs"] = [{"id": i, "role": role, "path": f"{run_dir}/{role.split('.')[1]}_{spec.key}.txt", "root": root0, "relpath": f"{rel_dir}/{role.split('.')[1]}_{spec.key}.txt",
                                  "size": size, "mtime": "2026-11-13T11:00:00Z", "state": "ok", "scope": "run", "link": {"method": "run_dir", "confidence": "certain"},
                                  "records": n_reads if role == "minknow.sequencing_summary" else None} for i, role, size in inputs]

    # ---- read data ----
    rl: dict = {"state": "absent"}
    rd_summary: dict = {}
    sample_rows: list[dict] = []
    if spec.read_data != "absent" and spec.samples:
        rd_summary, rl, sample_rows = _read_data(b, spec, rng, n_reads, run_dir, rel_dir)
    b["tiers"]["read_level"] = rl

    status, relations = _judge(spec, b, relations, y_sum, fc_sum, dyn_sum, le_sum, rd_summary, n_err)
    out = bundle_model.compose(b, status=status, relations=relations)
    out["content_sha256"] = bundle_model.content_sha256(out)
    return bundle_model.order_members(out)


def _read_data(b, spec, rng, n_reads, run_dir, rel_dir):
    """The read-data sections of one run; exact unless a section is marked partial."""
    partial = spec.read_data == "partial"
    reason = spec.partial_reason if partial else None
    names = list(spec.samples)
    with_data = names[: max(1, len(names) // 2)] if reason == "samples_missing" else names
    frac = 0.46 if reason == "partly_read" else 1.0
    basecalled = int(n_reads * 0.985)
    classified_share = 0.9
    pop = {"stage": "basecalled_all", "min_qscore": 9, "trimmed": False, "both_end_barcodes": False, "source_role": "dorado.bam",
           "evidence": ["@PG dorado basecaller with --min-qscore 9", "the file holds reads below Q 9 as well"]}
    aln_pop = {"stage": "aligned_primary", "min_qscore": 9, "trimmed": True, "both_end_barcodes": False, "source_role": "align.genome_bam", "evidence": ["@PG minimap2"]}
    samples, rows = [], []
    for i, name in enumerate(names):
        proj = spec.sample_project if (spec.sample_project and i >= len(names) // 2) else spec.project
        samples.append({"key": f"s{i + 1}", "barcode": f"barcode{i + 1:02d}", "alias": name, "name": name,
                        "sample_uid": f"{project_id(proj or 'none')}/{name}"[:256], "type": "test_sample", "expected": True, "source": "basecaller_sample_sheet", "project": proj})
    samples.append({"key": "unclassified", "barcode": None, "alias": "unclassified", "name": "No barcode found", "expected": False, "source": "bam_read_group"})
    b["samples"] = samples
    coarse = len(names) >= 13
    len_axis = "read_length_log20" if coarse else "read_length_log40"
    mean_len = {"cdna": 1100.0, "gdna": 6000.0, "drna": 880.0}[spec.assay]
    shares = rng.dirichlet(np.full(len(names), 8.0)) * classified_share
    yb_cells, len_cells, qual_cells, aln_cells, cdna_cells, fs_cells = [], [], [], [], [], []
    tot_n = np.zeros(len(ax.EDGES[len_axis]) + 1, np.int64)
    tot_w = tot_n.copy()
    tot_q = np.zeros(len(ax.EDGES["qscore_0p25"]) + 1, np.int64)
    pass_reads = pass_bases = all_reads = all_bases = mapped = aln_reads = 0
    covered = 0
    keys = [(s["key"], s["name"]) for s in samples if s["key"] != "unclassified" and s["name"] in with_data] + [("unclassified", None)]
    per_sample: dict = {}
    for key, name in keys:
        share = shares[names.index(name)] if name else (1 - classified_share)
        s_reads = int(basecalled * share * frac)
        for flt, part, qmean, lfac in (("pass", 0.91, 17.5, 1.0), ("fail", 0.09, 6.5, 0.7)):
            n = int(s_reads * part)
            hn, hw = _lognormal_hist(len_axis, rng, n, mean_len * lfac * (0.9 + 0.2 * rng.random()), 0.7)
            hq = _normal_hist("qscore_0p25", n, qmean + rng.normal(0, 0.3), 3.2 if flt == "pass" else 1.4)
            bases = int(hw.sum())
            yb_cells.append({"at": [key, flt], "scalars": {"reads": n, "bases": bases}})
            len_cells.append({"at": [key, flt], "scalars": {"reads": n, "bases": bases}, "hists": {"len_hist": agg.hist(len_axis, hn, hw, w_name="bases")}})
            qcell = {"at": [key, flt], "hists": {"qual_hist": agg.hist("qscore_0p25", hq, missing=0)}}
            if not coarse:
                lx, qx = ax.EDGES["read_length_log20"], ax.EDGES["qscore_1"]
                ln = (hn[1:-1].reshape(-1, 2).sum(axis=1)).astype(np.float64)
                lw = (hw[1:-1].reshape(-1, 2).sum(axis=1)).astype(np.float64)
                qn = hq[1:-1].reshape(-1, 4).sum(axis=1).astype(np.float64)
                m = np.outer(ln, qn / max(qn.sum(), 1))
                mi = np.floor(m).astype(np.int64)
                mw = np.floor(np.outer(lw, qn / max(qn.sum(), 1))).astype(np.int64)
                assert mi.shape == (len(lx) - 1, len(qx) - 1)
                qcell["hists2d"] = {"qual_vs_length": agg.hist2d("read_length_log20", "qscore_1", mi, mw, outside=0, missing=0, w_name="bases")}
            qual_cells.append(qcell)
            all_reads += n
            all_bases += bases
            if flt == "pass":
                pass_reads += n
                pass_bases += bases
                tot_n += hn
                tot_w += hw
                tot_q += hq
                if name:
                    covered += n
                    per_sample[key] = {"reads": n, "bases": bases, "n50": round(_n50(ax.EDGES[len_axis], hw[1:-1]) or 0),
                                       "q_med": round(_quantile(ax.EDGES["qscore_0p25"], hq[1:-1], 0.5) or 0, 2)}
        if name:
            n = int(s_reads * 0.91)
            m_reads = int(n * 0.972)
            ident = _normal_hist("identity_q_0p1", m_reads - int(m_reads * 0.03), 19.5, 3.5)
            ident[-1] = int(m_reads * 0.03)
            mq = np.zeros(len(ax.EDGES["mapq_1"]) + 1, np.int64)
            mq[1 + 60] = int(m_reads * 0.9)
            mq[1 + 1] = int(m_reads * 0.04)
            mq[1 + 30] = m_reads - int(mq.sum())
            aln_cells.append({"at": [key, "pass"], "scalars": {"reads": n, "mapped_reads": m_reads, "unmapped_reads": n - m_reads, "perfect_reads": int(m_reads * 0.03),
                                                              "alignment_columns": int(m_reads * mean_len), "edit_distance": int(m_reads * mean_len * 0.012),
                                                              "secondary_records": int(m_reads * 0.05), "supplementary_records": int(m_reads * 0.02)},
                              "hists": {"aln_identity_hist": agg.hist("identity_q_0p1", ident), "aln_mapq_hist": agg.hist("mapq_1", mq)}})
            mapped += m_reads
            aln_reads += n
            if spec.assay == "cdna":
                inp = n
                found, rescued, unusable = int(inp * 0.86), int(inp * 0.01), int(inp * 0.10)
                cdna_cells.append({"at": [key], "scalars": {"input": inp, "pass_reads": inp, "primers_found": found, "rescued_reads": rescued, "rescue_segments": rescued * 2,
                                                            "unusable": unusable, "len_fail": int(inp * 0.02), "qc_fail": inp - found - rescued - unusable - int(inp * 0.02), "umi_detected": int(inp * 0.8)}})
                fl = found + rescued
                fs_cells.append({"at": [key], "counts": {"funnel_stage_table": agg.counts(["classified", "full_length", "genome_mapped", "quantified"], [n, fl, 0, int(fl * 0.95)], role="funnel", missing=["genome_mapped"],
                                                                                    sources={"classified": "Dorado BAM: reads with this barcode", "full_length": "pychopper report", "quantified": "quantification result"})}})
            else:
                fs_cells.append({"at": [key], "counts": {"funnel_stage_table": agg.counts(["classified", "genome_mapped"], [n, m_reads], role="funnel",
                                                                                    sources={"classified": "Dorado BAM: reads with this barcode", "genome_mapped": "alignment file: primary mapped reads"})}})
    state = "partial" if partial else "ok"
    records_total = int(basecalled)
    cov = {"kind": "full" if not reason or reason != "partly_read" else "byte_ranges", "records": all_reads, "records_total": records_total if reason == "partly_read" else all_reads,
           "fraction": round(frac, 4)}
    edges = ax.EDGES[len_axis]
    summary = {"yield_basecalled_reads": all_reads, "yield_basecalled_bases": all_bases, "yield_pass_frac": round(pass_reads / max(all_reads, 1), 4),
               "len_n50": round(_n50(edges, tot_w[1:-1])), "len_median": round(_quantile(edges, tot_n[1:-1], 0.5, log=True)), "len_mean": round(pass_bases / max(pass_reads, 1), 1),
               "qual_median": round(_quantile(ax.EDGES["qscore_0p25"], tot_q[1:-1], 0.5), 2), "aln_mapped_frac": round(mapped / max(aln_reads, 1), 4), "aln_identity_median": 0.989}

    def section(mod_id, dims, cells, summ, population=None, tier="read_level"):
        s = {"version": "0.1.0", "tier": tier, "state": state, "inputs": ["rd1"], "dims": dims, "cells": cells, "summary": summ, "coverage": cov, "detail": "coarse" if coarse else "standard"}
        if population:
            s["population"] = population
        if partial:
            s["partial_reason"] = reason
        b["modules"][mod_id] = s

    section("yield.basecalled", ["sample", "filter"], yb_cells, {k: summary[k] for k in ("yield_basecalled_reads", "yield_basecalled_bases", "yield_pass_frac")}, pop)
    section("len", ["sample", "filter"], len_cells, {k: summary[k] for k in ("len_n50", "len_median", "len_mean")}, pop)
    section("qual", ["sample", "filter"], qual_cells, {"qual_median": summary["qual_median"]}, pop)
    section("aln", ["sample", "filter"], aln_cells, {"aln_mapped_frac": summary["aln_mapped_frac"], "aln_identity_median": summary["aln_identity_median"]}, aln_pop)
    classified = covered
    stages = ["sequenced", "basecalled", "classified"]
    counts = [n_reads, int(all_reads / frac) if frac < 1 else all_reads, classified]
    missing = ["basecalled"] if reason == "partly_read" else []
    sources = {"sequenced": "MinKNOW report: read count", "basecalled": "Dorado BAM: every record", "classified": "Dorado BAM: reads with a barcode of the sample sheet"}
    if spec.assay == "cdna":
        fl = sum(c["scalars"]["primers_found"] + c["scalars"]["rescued_reads"] for c in cdna_cells)
        stages += ["full_length", "genome_mapped", "quantified"]
        counts += [fl, 0, int(fl * 0.95)]
        missing.append("genome_mapped")
        sources.update({"full_length": "pychopper reports", "quantified": "quantification results"})
        usable, usable_stage = int(fl * 0.95), "quantified"
        b["modules"]["cdna"] = {"version": "0.1.0", "tier": "assay", "state": "ok", "inputs": ["rd2"], "dims": ["sample"], "cells": cdna_cells,
                                "summary": {"cdna_full_length_frac": round(fl / max(1, sum(c["scalars"]["input"] for c in cdna_cells)), 4), "cdna_samples": len(cdna_cells)},
                                "notes": [{"level": "info", "code": "source_is_report", "message": "Collected from each sample's pychopper statistics file; exact."}]}
    else:
        stages += ["genome_mapped"]
        counts += [mapped]
        sources["genome_mapped"] = "alignment files: primary mapped reads"
        usable, usable_stage = mapped, "genome_mapped"
    usable_frac = round(usable / max(n_reads, 1), 4) if not partial else None
    f_sum = {"funnel_reads_sequenced": n_reads, "funnel_usable_stage": usable_stage}
    if usable_frac is not None:
        f_sum["funnel_usable_frac"] = usable_frac
    summary["funnel_usable_frac"] = usable_frac
    b["modules"]["funnel"] = {"version": "0.1.0", "tier": "read_level", "state": state, "inputs": ["rd1"], "dims": [],
                              "cells": [{"at": [], "counts": {"funnel_stage_table": agg.counts(stages, counts, role="funnel", missing=missing, sources=sources)}}], "summary": f_sum, "coverage": cov}
    b["modules"]["funnel.samples"] = {"version": "0.1.0", "tier": "read_level", "state": state, "inputs": ["rd1"], "dims": ["sample"], "cells": fs_cells, "coverage": cov}
    if partial:
        b["modules"]["funnel"]["partial_reason"] = reason
        b["modules"]["funnel.samples"]["partial_reason"] = reason
    b["identity"]["basecalling"] = [{"basecaller": "dorado", "version": "1.1.0", "model": "dna_r10.4.1_e8.2_400bps_sup@v5.2.0", "modbase_models": [], "live": False,
                                     "min_qscore": 9, "inputs": ["rd1"], "model_tier": "sup", "variant": None, "primary": True}]
    b["provenance"]["inputs"].append({"id": "rd1", "role": "dorado.bam", "path": f"{run_dir}/basecalled/{spec.key}.bam", "root": spec.roots[0], "relpath": f"{rel_dir}/basecalled/{spec.key}.bam",
                                      "size": all_bases * 2, "mtime": "2026-11-13T11:30:00Z", "state": "partially_read" if reason == "partly_read" else "ok", "scope": "run",
                                      "link": {"method": "bam_header_rg", "confidence": "certain", "evidence": "@RG PU names the flow cell and DS the run"}, "records": all_reads, "coverage": cov})
    rl = {"state": "partial" if partial else "complete", "origin": "offline_basecalling", "samples_expected": len(names), "samples_with_data": len(with_data),
          "reads_expected": n_reads, "reads_covered": covered, "missing_roles": [], "population": pop, "coverage": cov}
    rows = [{"key": s["key"], "name": s["name"], "alias": s["alias"], "barcode": s["barcode"], "project": s.get("project"), "sample_uid": s["sample_uid"],
             **(per_sample.get(s["key"]) or {"reads": None, "bases": None, "n50": None, "q_med": None})} for s in samples if s["key"] != "unclassified"]
    summary.update({"_pass_bases": pass_bases, "_covered": covered, "_waiting": 3 if reason == "waiting" else 0})
    return summary, rl, rows


def _judge(spec, b, relations, y, fc, dyn, le, rd, n_err):
    """A status in the stored form, with sentences in the catalogue's style. The numbers decide; nothing is forced."""
    full, quick = spec.detail == "full", spec.detail == "quick"
    rna = spec.assay != "gdna"
    group = " | ".join((b["identity"]["device"]["type"], b["identity"]["flow_cell"]["product_code"], b["identity"]["kit"]["sequencing_kit"]))
    rules: list[dict] = []

    def rule(rid, module, role, panel, value, outcome, message, threshold, severity=None, tier="run_health", metric=None, reason=None):
        r = {"id": rid, "module": module, "role": role, "tier": tier, "severity": severity or "none", "outcome": outcome, "metric": metric or rid, "value": value,
             "threshold": threshold, "message": message, "panel": panel}
        if reason:
            r["not_evaluated_reason"] = reason
        rules.append(r)

    def minimum(rid, module, role, panel, value, warn, fail, text, limit_word, metric=None):
        if fail is not None and value < fail and role == "key":
            sev, out = "fail", "triggered"
        elif value < warn:
            sev, out = "warn", "triggered"
        else:
            sev, out = None, "ok"
        rule(rid, module, role, panel, value, out, f"{text} {limit_word}", {"kind": "min", "warn": warn, "fail": fail, "basis": "built-in profile", "applied": "fixed"}, sev, metric=metric)

    stop_ok = spec.stop in ("PROTOCOL_COMPLETED", "PROTOCOL_STOPPED_BY_USER")
    rule("run_stop_cause", "run", "key", "run/events", spec.stop, "ok" if stop_ok else "triggered",
         f"The run ended as planned after {spec.hours:.0f} h." if spec.stop == "PROTOCOL_COMPLETED" else f"The run was stopped by a person after {spec.hours:.1f} h." if stop_ok else "The run ended with an error.",
         {"kind": "rule", "basis": "pass when the protocol completed or a person stopped it", "applied": "rule"}, None if stop_ok else "fail")
    rule("run_messages", "run", "key", "run/events", n_err, "triggered" if n_err else "ok",
         "MinKNOW messages: 1 error message needs a look." if n_err else "MinKNOW messages: none that need a look.",
         {"kind": "rule", "basis": "warn on an error message of the device", "applied": "rule"}, "warn" if n_err else None)
    too_short = spec.hours < 1.0
    if not too_short:
        if spec.use[0] > 1:
            rule("fc_pores_first_scan", "fc", "key", "fc/pores_per_scan", fc["fc_pores_first_scan"], "not_evaluated",
                 f"Not assessed: a re-used flow cell is not judged on its pore count ({fc['fc_pores_first_scan']:,} working pores at the first scan).",
                 {"kind": "rule", "basis": "re-used flow cell", "applied": "rule"}, reason="re-used flow cell")
        else:
            minimum("fc_pores_first_scan", "fc", "key", "fc/pores_per_scan", fc["fc_pores_first_scan"], 5000 if spec.family == "promethion" else 800, 2500 if spec.family == "promethion" else 400,
                    f"Working pores at the first scan after loading: {fc['fc_pores_first_scan']:,}.", f"The limit is {5000 if spec.family == 'promethion' else 800:,} or more.")
    if not too_short and spec.use[0] == 1:
        minimum("fc_channels_with_pore_first_scan", "fc", "supporting", "fc/first_scan", fc["fc_channels_with_pore_first_scan_frac"], 0.65, None,
                f"Channels with a working pore at the first scan: {fc['fc_channels_with_pore_first_scan']:,} of {spec.channels:,} ({fc['fc_channels_with_pore_first_scan_frac'] * 100:.0f} %).",
                "The limit is 65 % or more.", metric="fc_channels_with_pore_first_scan_frac")
    if not too_short:
        if spec.hours >= 1.5:
            minimum("fc_occupancy_first_hour", "fc", "key", "fc/occupancy", fc["fc_occupancy_first_hour"], 0.7, 0.4,
                    f"Occupancy in the first hour: {fc['fc_occupancy_first_hour'] * 100:.0f} %.", "ONT's guide is 70 % or more.")
            minimum("fc_occupancy_median", "fc", "supporting", "fc/occupancy", fc["fc_occupancy_median"], 0.7, None,
                    f"Occupancy over the run: {fc['fc_occupancy_median'] * 100:.0f} %.", "The limit is 70 % or more.")
        v = fc["fc_recovering_frac"]
        rule("fc_recovering_frac", "fc", "key", "fc/channel_activity", v, "triggered" if v > 0.12 else "ok", f"Blocked pores: {v * 100:.1f} % of working-pore time. The limit is 12 % or less.",
             {"kind": "max", "warn": 0.12, "fail": 0.2, "basis": "built-in profile", "applied": "fixed"}, "warn" if v > 0.12 else None)
        rule("fc_temperature_max_dev_c", "fc", "supporting", "fc/temperature", 0.0, "ok", "Temperature: inside the target band for the whole run.",
             {"kind": "max", "warn": 0.5, "fail": None, "basis": "built-in profile", "applied": "fixed"})
    pending: list[str] = []
    if full and not too_short:
        sp = dyn["dyn_speed_est_median"]
        lo = b["modules"]["fc"]["params"]["target_speed_min_bps"]
        rule("dyn_speed_est_median", "dyn", "key", "dyn/speed_over_time", sp, "ok", f"Speed through the pore: {sp:.0f} bases per second, estimated. The target set by this run's own kit is {lo:.0f} to {lo / 0.875 * 1.125:.0f}.",
             {"kind": "run_config", "warn": None, "fail": lo * 0.85, "pass_band": [lo, lo / 0.875 * 1.125], "basis": "the run's own kit", "applied": "run_config"})
        minimum("dyn_signal_positive_frac", "dyn", "key", "dyn/why_reads_ended", dyn["dyn_signal_positive_frac"], 0.9, 0.8,
                f"Reads that ended normally: {dyn['dyn_signal_positive_frac'] * 100:.1f} %.", "The limit is 90 % or more.")
        minimum("fc_active_channels", "fc", "supporting", "fc/channel_map", fc["fc_active_channels_frac"], 0.65, None,
                f"Channels that produced reads: {fc['fc_active_channels']:,} of {spec.channels:,} ({fc['fc_active_channels_frac'] * 100:.0f} %).", "The limit is 65 % or more.", metric="fc_active_channels_frac")
        minimum("len_est_n50", "len.est", "supporting", "len/estimated_length", le["len_est_n50"], 500, None, f"Estimated read N50: {le['len_est_n50']:,} bases.", "The limit is 500 bases or more.")
    elif quick and not too_short:
        pending = ["dyn_speed_est_median", "dyn_speed_drop_frac", "dyn_signal_positive_frac", "dyn_unblock_frac", "fc_active_channels", "len_est_n50"]

    # the output check at the judged horizon, against history when there is enough of it
    horizon = next((h for h in reversed(HORIZONS) if spec.hours >= h - 1e-6), None)
    explanations, modifiers = [], []
    out_rule = None
    if horizon:
        unit = "reads" if rna else "estimated_bases"
        rid = f"yield_{unit}_at_{horizon}h"
        value = y[rid]
        fixed = (9e6 if rna else 9e9) * horizon / 24 * (1 if spec.family == "promethion" else 0.12)
        fmt = _fmt_count if rna else _fmt_bases
        label = f"{'Reads' if rna else 'Estimated bases'} at {horizon} h: {fmt(value)}."
        if spec.peers >= 5 and spec.use[0] == 1:
            median = fixed * 2.4
            scale = median * 0.2
            band = [median - 1.5 * scale, median + 1.5 * scale]
            warn_level, fail_level = median - 2 * scale, (median - 3.2 * scale if spec.peers >= 8 else None)
            z = (value - median) / scale
            base = {"group": group, "trust": "trusted" if spec.peers >= 8 else "provisional", "n": spec.peers, "median": median, "scale": scale, "z": round(z, 2), "rank_below": 3,
                    "band": band, "warn_z": -2.0, "fail_z": -3.2 if spec.peers >= 8 else None, "horizon_h": horizon, "as_of": spec.started, "spread": 0.2, "spread_guard": False,
                    "warn_level": warn_level, "fail_level": fail_level}
            if fail_level is not None and value < fail_level:
                sev, outc, msg = "fail", "triggered", f"{label} This lab's fail level at {horizon} h is {fmt(fail_level)} ({spec.peers} comparable runs)."
            elif value < warn_level:
                sev, outc, msg = "warn", "triggered", f"{label} This lab's warning level at {horizon} h is {fmt(warn_level)} ({spec.peers} comparable runs)."
            else:
                sev, outc, msg = None, "ok", f"{label} Typical here is {fmt(band[0])} to {fmt(band[1])} ({spec.peers} comparable runs)."
            rule(rid, "yield", "key", "yield/output_over_time", value, outc, msg, {"kind": "baseline", "warn": fixed, "fail": fixed * 0.5, "basis": f"{spec.peers} comparable runs before this one",
                                                                                   "peers": spec.peers, "baseline": base, "applied": "baseline"}, sev)
            if spec.peers < 8:
                modifiers.append({"code": "baseline_provisional", "message": f"Lab typical rests on {spec.peers} comparable runs: it can warn but not fail a run until there are 8."})
            for eid, aspect, val, med, sc in (("fc_pores_first_scan", "flow_cell", fc["fc_pores_first_scan"], 6800 if spec.family == "promethion" else 1200, 700),
                                              ("fc_occupancy_median_24h", "loading", fc["fc_occupancy_median_24h"], 0.86, 0.05),
                                              ("yield_kb_per_active_channel_hour_48h", "library", 1200.0 * (spec.occ / 0.9), 1150.0, 110.0)):
                zz = (val - med) / sc
                explanations.append({"id": eid, "aspect": aspect, "value": val, "direction": "low" if zz <= -1.5 else "high" if zz >= 1.5 else "typical",
                                     "baseline": {"group": group, "trust": "trusted" if spec.peers >= 8 else "provisional", "n": spec.peers, "median": med, "scale": sc, "z": round(zz, 2),
                                                  "band": [med - 1.5 * sc, med + 1.5 * sc], "as_of": spec.started}})
        else:
            sev, outc = (("warn", "triggered") if value < fixed else (None, "ok"))
            why = "No lab history for this combination yet" if spec.use[0] == 1 else "A re-used flow cell is not compared with the lab's first-use runs"
            rule(rid, "yield", "key", "yield/output_over_time", value, outc, f"{label} {why}; the fixed level at {horizon} h is {fmt(fixed)}.",
                 {"kind": "min", "warn": fixed, "fail": fixed * 0.5, "basis": "built-in profile", "applied": "fixed"}, sev)
            if spec.use[0] == 1:
                modifiers.append({"code": "no_history", "message": "No lab history for this combination yet: the output was judged on fixed levels."})
        out_rule = rules[-1]
    if spec.use[0] > 1:
        modifiers.append({"code": "reused_flow_cell", "message": f"Re-used flow cell, use {spec.use[0]} of {spec.use[1]}: pore count and total output are shown but not judged."})
    if quick:
        modifiers.append({"code": "quick_detail", "message": "Made from MinKNOW's report files. Speed, why reads ended, read length and the channel map appear when the read count finishes."})

    # read data
    has_rd = b["tiers"]["read_level"]["state"] != "absent"
    if has_rd:
        partial = b["tiers"]["read_level"]["state"] == "partial"
        if rd.get("funnel_usable_frac") is not None:
            v = rd["funnel_usable_frac"]
            sev = "fail" if v < 0.35 else "warn" if v < 0.55 else None
            rule("funnel_usable_frac", "funnel", "key", None, v, "triggered" if sev else "ok", f"Usable reads: {v * 100:.1f} % of the reads sequenced. The limit is 55 % or more.",
                 {"kind": "min", "warn": 0.55, "fail": 0.35, "basis": "built-in profile", "applied": "fixed"}, sev, tier="read_level")
        else:
            rule("funnel_usable_frac", "funnel", "key", None, None, "not_evaluated", "Not assessed: read data is not complete for this run yet.",
                 {"kind": "rule", "basis": "needs complete read data", "applied": "rule"}, tier="read_level", reason="read data is not complete for this run yet")
        v = rd["aln_mapped_frac"]
        if partial:
            rule("aln_mapped_frac", "aln", "key", None, None, "not_evaluated", "Not assessed: read data is not complete for this run yet.",
                 {"kind": "rule", "basis": "needs complete read data", "applied": "rule"}, tier="read_level", reason="read data is not complete for this run yet")
        else:
            rule("aln_mapped_frac", "aln", "key", None, v, "triggered" if v < 0.9 else "ok", f"Mapped reads: {v * 100:.1f} %. The limit is 90 % or more.",
                 {"kind": "min", "warn": 0.9, "fail": 0.7, "basis": "built-in profile", "applied": "fixed"}, "warn" if v < 0.9 else None, tier="read_level")
            q = rd["qual_median"]
            rule("qual_median", "qual", "supporting", None, q, "triggered" if q < 12 else "ok", f"Median read quality: Q {q:.1f}. The limit is Q 12 or more.",
                 {"kind": "min", "warn": 12, "fail": None, "basis": "built-in profile", "applied": "fixed"}, "warn" if q < 12 else None, tier="read_level")

    rh = status_model.worst_of_rules(rules, "run_health")
    cause = None
    if too_short:
        rh, cause = "incomplete", "stopped_early"
    elif spec.hours < 6 and rh == "pass":
        rh, cause = "incomplete", "output_not_judged"
    if spec.run_state == "running":
        rh = "running"
    rd_rules = [r for r in rules if r["tier"] == "read_level"]
    data = "not_available" if not has_rd else ("unknown" if all(r["outcome"] == "not_evaluated" for r in rd_rules) else status_model.worst_of_rules(rules, "read_level"))
    verdict = status_model.combine_verdicts(rh, data)
    trig = [r for r in rules if r["outcome"] == "triggered" and r["tier"] == "run_health"]
    if cause == "stopped_early":
        headline = f"Stopped after {max(1, round(spec.hours * 60))} minutes. Too short to judge."
    elif cause == "output_not_judged":
        headline = f"Stopped after {spec.hours:.1f} hours. No check failed, but the run is too short to judge its output."
    elif spec.run_state == "running":
        headline = f"Sequencing since {spec.hours:.0f} h. Checks are made as the numbers arrive."
    elif trig:
        names = {"fc_pores_first_scan": "pores at the first scan", "fc_occupancy_first_hour": "occupancy in the first hour", "fc_occupancy_median": "occupancy over the run",
                 "fc_channels_with_pore_first_scan": "channels with a pore at the first scan", "fc_recovering_frac": "blocked pores", "run_messages": "a MinKNOW error message",
                 "run_stop_cause": "how the run ended"}
        words = [names.get(r["id"], "output") for r in trig]
        headline = ("Failed on: " if rh == "fail" else "Needs a look: ") + ", ".join(dict.fromkeys(words)) + "."
    elif out_rule:
        headline = f"The flow cell ran well ({_fmt_count(y[out_rule['id']]) + ' reads' if rna else _fmt_bases(y[out_rule['id']])} by {horizon} h, {fc['fc_occupancy_first_hour'] * 100:.0f} % of pores busy in the first hour)."
    else:
        headline = "No check triggered."
    rl = b["tiers"]["read_level"]
    detail = ""
    if has_rd:
        detail = f"Read data for {rl['samples_with_data']} of {rl['samples_expected']} samples" + ("." if rl["state"] == "complete" else {
            "samples_missing": "; the others have no files yet.", "partly_read": f"; {rl['coverage']['fraction'] * 100:.0f} % of the file has been read so far.",
            "waiting": "; 3 files wait to be read."}.get(spec.partial_reason or "", "."))
    status = {"profile": {"id": spec.assay, "version": "0.1.1", "source": "builtin"}, "evaluated": "2026-11-13T12:00:00Z", "verdict": verdict, "headline": headline,
              "by_tier": {"run_health": rh, "read_level": data}, "rules": rules, "detail": detail, "cause": cause, "modifiers": modifiers, "estimates": [],
              "provisional": quick, "pending": pending, "explanations": explanations,
              "judged": {"at": "2026-11-13T12:00:00Z", "as_of": spec.started, "history_hash": "h1:0011223344556677", "library": True, "frozen": True}}
    return status, relations


def project_id(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_") or "project"


# ---------------------------------------------------------------------------------------------------------
# the index
# ---------------------------------------------------------------------------------------------------------

def _columns(spec: Spec, b: dict) -> dict:
    st, ident, tiers = b["status"], b["identity"], b["tiers"]
    s: dict = {}
    for m in b["modules"].values():
        s.update(m.get("summary") or {})
    rs = status_model.reasons(st)
    expl = {e["aspect"]: {"low": -1, "typical": 0, "high": 1}.get(e["direction"]) for e in st.get("explanations") or []}
    out_rule = next((r for r in st["rules"] if r["id"].startswith("yield_") and "_at_" in r["id"]), None)
    raw = spec.detail is None
    hours = None if raw else spec.hours
    spark, spark_max = "", 0.0
    if not raw:
        v = b["modules"]["yield"]["cells"][0]["series"]["yield_series"]["values"]["estimated_bases"]
        per3 = [sum(v[i:i + 18]) / 1e9 for i in range(0, len(v), 18)]
        spark, spark_max = agg.sparkline(per3)
        spark = spark[: len(per3)].ljust(32, "A") if len(per3) < 32 else spark
        spark = spark[: len(per3)]
    rl = tiers["read_level"]
    rd = spec.extra.get("rd") or {}
    projects = b["relations"]["projects"]
    cols = {
        "id": spec.key, "fc": spec.fc, "use": spec.use[0], "use_of": spec.use[1], "use_source": "library" if spec.use[1] > 1 else "unknown",
        "restart": bool(spec.restart_of), "restart_of": spec.restart_of, "started": spec.started, "hours": hours, "planned_h": None if raw else spec.planned_h,
        "run_state": "raw_only" if raw else spec.run_state, "device_type": ident["device"].get("type"), "position": spec.position,
        "product": ident["flow_cell"].get("product_code"), "kit": ident["kit"].get("sequencing_kit"), "assay": ident["assay"]["type"], "assay_group": ident["assay"]["group"],
        "chemistry": ident["assay"].get("chemistry"), "minknow": ident["software"].get("minknow"), "experiment": ident["protocol_run"].get("experiment_id"),
        "run_name": ident["protocol_run"].get("sample_id"), "channels": spec.channels,
        "group": None if raw else " | ".join((ident["device"]["type"], ident["flow_cell"]["product_code"], ident["kit"]["sequencing_kit"])),
        "folder": b["provenance"]["run_dir"]["relpath"].rsplit("/", 1)[-1], "roots": list(spec.roots),
        "project": projects[0]["name"] if projects else None, "project_ids": [p["id"] for p in projects], "source_state": spec.source_state,
        "verdict": st["verdict"], "health": st["by_tier"]["run_health"], "data": st["by_tier"]["read_level"], "provisional": st["provisional"],
        "pending": len(st["pending"]), "checks_done": len(st["rules"]), "cause": st["cause"], "modifiers": [m["code"] for m in st["modifiers"]],
        "headline": st["headline"], "detail": st["detail"], "reasons": [r["message"] for r in rs], "reason_levels": [r["severity"] for r in rs],
        "reason_checks": [r["check"] for r in rs], "reason_tiers": [r["tier"] for r in rs],
        "outcome_metric": out_rule["id"] if out_rule else None, "outcome_h": int(out_rule["id"].rsplit("_at_", 1)[1][:-1]) if out_rule else None,
        "outcome_value": out_rule["value"] if out_rule else None, "outcome_level": None, "outcome_n": spec.peers if out_rule else None, "outcome_trust": None,
        "mark_flow_cell": expl.get("flow_cell"), "mark_loading": expl.get("loading"), "mark_library": expl.get("library"),
        "reads": s.get("yield_reads"), "gb": s["yield_estimated_bases"] / 1e9 if "yield_estimated_bases" in s else None,
        "pores": s.get("fc_pores_first_scan"), "pores_last": s.get("fc_pores_last_scan"), "active_peak": s.get("fc_active_channels_peak"),
        "pore_half_h": s.get("fc_pore_half_life_h"), "duty_half_h": None, "occ_first": s.get("fc_occupancy_first_hour"), "occ_24h": s.get("fc_occupancy_median_24h"),
        "occ_run": s.get("fc_occupancy_median"), "recovering": s.get("fc_recovering_frac"), "kb_per_ch_h_48h": None, "mb_per_pore_48h": s.get("yield_per_pore_mb_48h"),
        "speed": s.get("dyn_speed_est_median"), "n50_est": s.get("len_est_n50"), "normal_end_frac": s.get("dyn_signal_positive_frac"), "active_channels": s.get("fc_active_channels"),
        "detail_t1": "none" if raw else spec.detail, "spark_output": spark, "spark_output_max_gb_per_3h": spark_max, "spark_channels": "", "spark_channels_max": 0,
        "history_hash": st["judged"]["history_hash"], "limits": "builtin-0.1.1",
        "read_data": rl["state"], "t2_population": (rl.get("population") or {}).get("stage"), "t2_samples": rl.get("samples_expected"),
        "t2_samples_with_data": rl.get("samples_with_data"), "t2_reads": rl.get("reads_covered"),
        "t2_passed_gb": rd["_pass_bases"] / 1e9 if rd else None, "t2_n50": rd.get("len_n50"), "t2_q_med": rd.get("qual_median"),
        "t2_ident_med": rd.get("aln_identity_median"), "t2_mapped_frac": rd.get("aln_mapped_frac"), "t2_usable": rd.get("funnel_usable_frac"),
        "coverage": (rl.get("coverage") or {}).get("fraction"), "t2_files_waiting": rd.get("_waiting", 0),
        "bundle": bundle_model.bundle_filename(spec.key), "bundle_bytes": 0, "bundle_gzip_bytes": 0,
    }
    for h in HORIZONS:
        cols[f"reads_{h}h"] = s.get(f"yield_reads_at_{h}h")
        cols[f"gb_{h}h"] = s[f"yield_estimated_bases_at_{h}h"] / 1e9 if f"yield_estimated_bases_at_{h}h" in s else None
    return cols


def _encode(values: list):
    """A string column with at most 64 distinct values is dictionary-encoded, as LIBRARY does."""
    if all(v is None or isinstance(v, str) for v in values):
        distinct = list(dict.fromkeys(values))
        if len(distinct) <= 64 and len(values) > len(distinct):
            return {"dict": distinct, "codes": [distinct.index(v) for v in values]}
    return values


def make_index(specs: list[Spec], bundles: list[dict], library: dict, roots: list[dict], *, written: str = "2026-11-13T12:05:00Z",
               progress: dict | None = None, notices: list[dict] | None = None, sample_rows: dict | None = None) -> dict:
    order = sorted(range(len(specs)), key=lambda i: specs[i].started, reverse=True)
    rows = [_columns(specs[i], bundles[i]) for i in order]
    names = list(rows[0].keys()) if rows else []
    columns = {k: [r[k] for r in rows] for k in names}
    for k in ("verdict", "health", "data", "assay", "assay_group", "kit", "product", "device_type", "position", "chemistry", "minknow", "run_state", "group",
              "use_source", "detail_t1", "read_data", "t2_population", "source_state", "limits", "project", "experiment"):
        if k in columns:
            columns[k] = _encode(columns[k])
    verdict_keys = ("pass", "warn", "fail", "incomplete", "running", "unknown")

    def tally(rs):
        t = {k: 0 for k in verdict_keys}
        for r in rs:
            t[r["verdict"] if r["verdict"] in t else "unknown"] += 1
        return t

    samples = []
    for i in order:
        for s in (sample_rows or {}).get(specs[i].key, []):
            samples.append({"sample_uid": s["sample_uid"], "name": s["name"], "alias": s["alias"], "barcode": s["barcode"], "project": s["project"], "run": specs[i].key,
                            "key": s["key"], "reads": s["reads"], "bases": s["bases"], "n50": s["n50"], "q_med": s["q_med"]})
    proj: dict = {}
    for r in rows:
        ids = r["project_ids"] or ["_none"]
        b = bundles[[s.key for s in specs].index(r["id"])]
        for p in b["relations"]["projects"] or [{"id": "_none", "name": None}]:
            e = proj.setdefault(p["id"], {"id": p["id"], "name": p["name"], "rows": [], "roots": set(), "samples": set()})
            e["rows"].append(r)
            e["roots"].update(r["roots"])
        del ids
    for s in samples:
        if s["project"]:
            pe = proj.get(project_id(s["project"]))
            if pe:
                pe["samples"].add(s["sample_uid"])
    projects = []
    for e in proj.values():
        rs = e["rows"]
        projects.append({"id": e["id"], "name": e["name"], "runs": len(rs), "flow_cells": len({r["fc"] for r in rs}), "samples": len(e["samples"]),
                         "first": min(r["started"] for r in rs), "last": max(r["started"] for r in rs), "verdicts": tally(rs),
                         "estimated_bases": int(sum((r["gb"] or 0) * 1e9 for r in rs)), "reads": int(sum(r["reads"] or 0 for r in rs)), "roots": sorted(e["roots"])})
    projects.sort(key=lambda p: (p["id"] == "_none", p["last"]), reverse=False)
    projects = sorted([p for p in projects if p["id"] != "_none"], key=lambda p: p["last"], reverse=True) + [p for p in projects if p["id"] == "_none"]
    root_rows = []
    for root in roots:
        rs = [r for r in rows if root["id"] in r["roots"]]
        root_rows.append({"id": root["id"], "path": root["path"], "live": root.get("live", False), "project_depth": 1, "state": root.get("state", "reachable"),
                          "offline_since": root.get("offline_since"), "last_scanned": root.get("last_scanned", written), "runs": len(rs),
                          "detail": {"full": sum(1 for r in rs if r["detail_t1"] == "full"), "quick": sum(1 for r in rs if r["detail_t1"] == "quick"),
                                     "raw_only": sum(1 for r in rs if r["detail_t1"] == "none")},
                          "read_data_files": {"found": sum(1 for r in rs if r["read_data"] != "absent") * 2, "read": sum(1 for r in rs if r["read_data"] == "complete") * 2}})
    months: list[str] = []
    if rows:
        first, last = min(r["started"] for r in rows)[:7], max(r["started"] for r in rows)[:7]
        y, m = int(first[:4]), int(first[5:7])
        while f"{y:04d}-{m:02d}" <= last:
            months.append(f"{y:04d}-{m:02d}")
            y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    activity = {"months": months, "runs": [sum(1 for r in rows if r["started"][:7] == mo) for mo in months],
                "estimated_bases": [int(sum((r["gb"] or 0) * 1e9 for r in rows if r["started"][:7] == mo)) for mo in months],
                "reads": [int(sum(r["reads"] or 0 for r in rows if r["started"][:7] == mo)) for mo in months]}
    groups = []
    for g in sorted({r["group"] for r in rows if r["group"]}):
        rs = [r for r in rows if r["group"] == g and r["detail_t1"] in ("full", "quick") and r["use"] == 1]
        if not rs:
            continue
        curves = []
        for r in rs:
            b = bundles[[s.key for s in specs].index(r["id"])]
            v = b["modules"]["yield"]["cells"][0]["series"]["yield_series"]["values"]["estimated_bases"]
            curves.append([sum(v[i:i + 18]) / 1e9 if i < len(v) else None for i in range(0, 32 * 18, 18)])
        arr = np.array([[np.nan if x is None else x for x in c] for c in curves], dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med, q25, q75 = np.nanmedian(arr, axis=0), np.nanpercentile(arr, 25, axis=0), np.nanpercentile(arr, 75, axis=0)
        clean = lambda a: [None if not np.isfinite(x) else round(float(x), 4) for x in a]   # noqa: E731
        n_peers = max(len(rs), max(s.peers for s in specs if s.key in {r["id"] for r in rs}))
        groups.append({"group": g, "n_runs": n_peers, "as_of": written, "window_first": min(r["started"] for r in rs)[:10], "window_last": max(r["started"] for r in rs)[:10],
                       "metrics": {"yield_estimated_bases_at_48h": {"n": n_peers, "median": 9e10, "spread": 0.2, "spread_guard": False}},
                       "curve_output_gb_per_3h": {"median": clean(med), "q25": clean(q25), "q75": clean(q75)}})
    by_detail = {"full": 0, "quick": 0, "none": 0}
    for r in rows:
        by_detail[r["detail_t1"]] += 1

    def facet(k):
        return sorted({r[k] for r in rows if r[k] is not None})

    return {
        "format": "flowcellar.index", "schema": SCHEMA_VERSION, "written": written,
        "meta": {"library": library["name"], "runs": len(rows), "generator": {"name": "flowcellar", "version": __version__}, "catalogue": "0.1.1",
                 "limits": {"id": "builtin", "version": "0.1.1", "hash": "0" * 16},
                 "layout": "columnar: one array per field; low-cardinality strings as {dict, codes}",
                 "sparklines": "32 points on a fixed 0 to 96 h axis (3 h per point), 6 bits per point, base64url alphabet, scaled to the row's own maximum",
                 "marks": "-1 low, 0 typical, 1 high, null unknown"},
        "library": {"id": library["id"], "name": library["name"], "path": library["path"], "created": "2026-10-01T09:00:00Z", "scanned": written,
                    "counts": {"runs": len(rows), "flow_cells": len({r["fc"] for r in rows}), "projects": len([p for p in projects if p["id"] != "_none"]),
                               "samples": len({s["sample_uid"] for s in samples}), "by_verdict": tally(rows), "by_detail": by_detail,
                               "with_read_data": sum(1 for r in rows if r["read_data"] != "absent")},
                    "estimated_bases": int(sum((r["gb"] or 0) * 1e9 for r in rows)), "reads": int(sum(r["reads"] or 0 for r in rows)),
                    "first_run": min((r["started"] for r in rows), default=None), "last_run": max((r["started"] for r in rows), default=None)},
        "runs": {"n": len(rows), "columns": columns}, "samples": samples, "projects": projects, "roots": root_rows, "activity": activity,
        "baselines": {"groups": groups},
        "facets": {k: facet(k) for k in ("verdict", "health", "data", "assay_group", "kit", "product", "device_type", "project", "detail_t1")} | {"root": [r["id"] for r in roots]},
        "notices": notices or [],
        "progress": progress or {"state": "idle", "started": None, "host": None, "time_box_minutes": 240,
                                 "tier1": {"runs_waiting": 0, "bytes_waiting": 0, "bytes_done": 0, "mb_per_s": None, "eta_minutes": None},
                                 "read_data": {"files_total": 0, "files_read": 0, "files_waiting": 0, "bytes_waiting": 0, "eta_minutes": None}, "message": None},
    }


# ---------------------------------------------------------------------------------------------------------
# the default library
# ---------------------------------------------------------------------------------------------------------

ROOTS = [
    {"id": "sequencer", "path": "/data/sequencer", "live": True},
    {"id": "projects", "path": "/data/projects"},
    {"id": "archive", "path": "/mnt/archive/nanopore", "state": "offline", "offline_since": "2026-11-12T02:10:00Z", "last_scanned": "2026-11-11T22:00:00Z"},
]


def default_specs() -> list[Spec]:
    four = ("sample_A_control_24h", "sample_B_treated_24h", "sample_C_control_48h", "sample_D_treated_48h")
    return [
        Spec("TST00001", "0a1b2c3d", "2026-11-04T15:42:00Z", samples=four, read_data="complete", pores=8100, occ=0.95, run_name="demo_cdna_plate1", experiment="demo_experiment_one", seed=1),
        Spec("TST00002", "1b2c3d4e", "2026-11-02T10:05:00Z", samples=four, read_data="partial", partial_reason="samples_missing", pores=6900, occ=0.64, run_name="demo_cdna_plate2",
             experiment="demo_experiment_one", position="TSTP2-00001-B", seed=2),
        Spec("TST00003", "2c3d4e5f", "2026-10-29T09:00:00Z", family="minion", assay="gdna", pores=1420, occ=0.93, run_name="demo_genome_isolates", experiment="demo_isolates",
             project="Demo Beta", roots=("projects",), position="TSTMN-00001", peers=1, seed=3),
        Spec("TST00004", "3d4e5f6a", "2026-09-03T14:00:00Z", hours=88, planned_h=96, detail="quick", pores=7300, occ=0.91, run_name="demo_cdna_quick", experiment="demo_experiment_two",
             project="Demo Beta", roots=("projects",), peers=6, seed=4),
        Spec("TST00005", "4e5f6a7b", "2026-07-10T08:00:00Z", hours=48, planned_h=72, pores=7000, occ=0.9, run_name="demo_first_loading", experiment="demo_reuse", use=(1, 2),
             stop="PROTOCOL_STOPPED_BY_USER", seed=5),
        Spec("TST00005", "5f6a7b8c", "2026-07-12T09:00:00Z", hours=0.08, planned_h=72, pores=2400, occ=0.7, run_name="demo_second_loading", experiment="demo_reuse", use=(2, 2),
             restart_of="TST00005_6a7b8c9d", stop="PROTOCOL_STOPPED_BY_USER", seed=6),
        Spec("TST00005", "6a7b8c9d", "2026-07-12T09:20:00Z", hours=40, planned_h=72, pores=2500, occ=0.84, half_life_h=30, run_name="demo_second_loading", experiment="demo_reuse",
             use=(2, 2), seed=7),
        Spec("TST00006", "7b8c9d0e", "2026-03-13T14:45:00Z", detail=None, run_name="", experiment="demo_old", project="Demo Gamma", roots=("archive",), source_state="offline", seed=8),
        Spec("TST00007", "8c9d0e1f", "2026-06-01T11:00:00Z", hours=3.5, planned_h=72, pores=6100, occ=0.88, run_name="demo_stopped_early", experiment="demo_experiment_two",
             stop="PROTOCOL_STOPPED_BY_USER", project=None, seed=9),
        Spec("TST00008", "9d0e1f2a", "2026-05-20T10:00:00Z", samples=four[:2], read_data="partial", partial_reason="partly_read", pores=2100, occ=0.35, half_life_h=40,
             run_name="demo_low_pores", experiment="demo_experiment_two", error_message=HOSTILE_MESSAGE, seed=10),
        Spec("TST00009", "0e1f2a3b", "2026-04-11T10:00:00Z", samples=(HOSTILE_SAMPLE, "sample_plain"), read_data="complete", pores=7600, occ=0.92, run_name=HOSTILE_RUN_NAME,
             experiment="demo <b>experiment</b>", project=HOSTILE_PROJECT, roots=("sequencer", "projects"), seed=11),
        Spec("TST00010", "1f2a3b4c", "2026-11-12T08:00:00Z", hours=14, planned_h=72, detail="quick", run_state="running", pores=7900, occ=0.93, run_name="demo_running_now",
             experiment="demo_experiment_three", stop="RUNNING", seed=12),
        Spec("TST00011", "2a3b4c5d", "2026-08-15T10:00:00Z", samples=tuple(f"sample_{i:02d}_long_descriptive_name" for i in range(1, 15)), read_data="complete", assay="gdna",
             pores=7700, occ=0.9, run_name="demo_many_samples_with_a_long_run_name", experiment="demo_experiment_three", project="Demo Alpha", sample_project="Demo Beta",
             roots=("sequencer",), peers=3, seed=13),
        Spec("TST00012", "3b4c5d6e", "2026-08-20T10:00:00Z", samples=four[:3], read_data="partial", partial_reason="waiting", pores=7100, occ=0.89, run_name="demo_waiting_files",
             experiment="demo_experiment_three", seed=14),
    ]


def make_library(out_dir: str, specs: list[Spec] | None = None, *, name: str = "Demo library", empty: bool = False,
                 scanning: bool = False) -> dict:
    """Write a library folder in the layout of a library (bundles, status files, index, marker, configuration).
    Returns ``{"path", "specs", "bundles", "index"}``."""
    out_dir = os.path.abspath(out_dir)
    library = {"id": "lib_0123456789ab", "name": name, "path": out_dir}
    specs = [] if empty else (default_specs() if specs is None else specs)
    for d in ("bundles", "status", "index", "cards", "state"):
        os.makedirs(os.path.join(out_dir, d), exist_ok=True)
    atomic.atomic_write_json(os.path.join(out_dir, "LIBRARY"), {"format": "flowcellar.library", "layout_version": 1, "library_id": library["id"],
                                                               "created": "2026-10-01T09:00:00Z", "generator": __version__})
    toml = [f'[library]\nname = "{name}"\n'] + [f'[[roots]]\nid = "{r["id"]}"\npath = "{r["path"]}"\n' for r in ROOTS]
    atomic.atomic_write_text(os.path.join(out_dir, "flowcellar.toml"), "\n".join(toml))
    bundles, sample_rows = [], {}
    for spec in specs:
        b = make_bundle(spec, library)
        rd = {}
        for mid in ("yield.basecalled", "len", "qual", "aln", "funnel"):
            rd.update((b["modules"].get(mid) or {}).get("summary") or {})
        if b["tiers"]["read_level"]["state"] != "absent":
            pass_bases = sum(c["scalars"]["bases"] for c in b["modules"]["yield.basecalled"]["cells"] if c["at"][1] == "pass")
            rd["_pass_bases"] = pass_bases
            rd["_waiting"] = 3 if spec.partial_reason == "waiting" else 0
            spec.extra["rd"] = rd
            by_key = {c["at"][0]: c for c in b["modules"]["len"]["cells"] if c["at"][1] == "pass"}
            rows = []
            for s in b["samples"]:
                if s["key"] == "unclassified":
                    continue
                c = by_key.get(s["key"])
                rows.append({"key": s["key"], "name": s["name"], "alias": s["alias"], "barcode": s["barcode"], "project": s.get("project"), "sample_uid": s["sample_uid"],
                             "reads": c["scalars"]["reads"] if c else None, "bases": c["scalars"]["bases"] if c else None, "n50": None, "q_med": None})
            sample_rows[spec.key] = rows
        measured, composed = bundle_model.split(b)
        bundle_model.write_bundle(os.path.join(out_dir, "bundles", bundle_model.bundle_filename(spec.key)), measured)
        status_model.write_status(os.path.join(out_dir, "status", status_model.status_filename(spec.key)),
                                  status_model.new_status_file(spec.key, b.get("content_sha256"), composed["status"], composed["relations"]))
        bundles.append(b)
    notices = [] if empty else [{"code": "root_offline", "level": "warn", "message": "/mnt/archive/nanopore has not answered since 12 Nov 2026, 02:10. Runs stored there are shown as last seen.",
                                 "root": "archive", "since": "2026-11-12T02:10:00Z"}]
    progress = None
    if scanning:
        progress = {"state": "scanning", "started": "2026-11-13T11:00:00Z", "host": "fakehost-01", "time_box_minutes": 240,
                    "tier1": {"runs_waiting": 2, "bytes_waiting": 40_000_000_000, "bytes_done": 10_000_000_000, "mb_per_s": 160, "eta_minutes": 5},
                    "read_data": {"files_total": 12, "files_read": 7, "files_waiting": 5, "bytes_waiting": 90_000_000_000, "eta_minutes": 40},
                    "message": "Counting reads: 1 of 3 runs done, about 5 min left."}
    index = make_index(specs, bundles, library, [] if empty else ROOTS, notices=notices, progress=progress, sample_rows=sample_rows)
    atomic.atomic_write_json(os.path.join(out_dir, "index", "index.json.gz"), index)
    return {"path": out_dir, "specs": specs, "bundles": bundles, "index": index}


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    made = make_library(sys.argv[1], empty="--empty" in sys.argv, scanning="--scanning" in sys.argv)
    print(json.dumps({"path": made["path"], "runs": [s.key for s in made["specs"]]}))
