"""Synthetic MinKNOW run folders for tests, written from scratch with invented names only.

    from flowcellar.testing.fixtures import RunSpec, make_run
    run = make_run(tmp_path, RunSpec(reads=4000, hours=8))
    run.run_dir, run.run_key, run.expected["reads"]

    python -m flowcellar.testing.fixtures OUT_DIR [--reads 4000] [--hours 8] [--family minion] [--seed 1]

This first version writes the newest file layout (MinKNOW 6.10-like): final summary, sample sheet, report JSON,
throughput, pore activity, the pore-scan table and a 16-column sequencing summary with basecalling off.
What it promises (and ``FixtureRun.expected`` records, so tests do not recompute it):

* the sequencing summary, the last row of the throughput file and the report's ``yield_summary`` agree exactly on
  the number of reads and on the sum of ``minknow_events``;
* per minute, the pore-activity file gives every channel exactly ``60 x sample rate`` samples, split over 18 states;
* the pore-scan table holds one row per channel, well and scan, and its counts per scan equal ``mux_scan_results``
  in the report;
* every identifier is invented and derived from the seed: equal specs give byte-identical folders.

It does not promise that channel-state time and read durations describe the same molecules (they come from
separate models). Builders extend this module: older layouts (15 columns, ``run_id`` = acquisition id,
``instrument`` = device), running, interrupted and raw-only folders, look-alikes, basecalling on, and the
read-data files of ``flowcellar.testing`` siblings.
"""
from __future__ import annotations

import csv
import io
import json
import os
import sys
from dataclasses import asdict, dataclass, field

import numpy as np

CHANNEL_STATES = ("strand", "adapter", "pore", "unavailable", "unblocking", "no_pore", "zero", "saturated", "multiple",
                  "unknown_positive", "unknown_negative", "disabled", "locked", "membrane", "unclassified",
                  "unclassified_following_reset", "pending_manual_reset", "pending_mux_change")
SCAN_CATEGORIES = ("single_pore", "reserved_pore", "unavailable", "saturated", "zero", "multiple", "other")
END_REASONS = ("signal_positive", "signal_negative", "unblock_mux_change", "mux_change")
_END_REASON_P = (0.955, 0.02, 0.015, 0.01)
_FAMILIES = {
    # family: (channels, device_type, product, protocol script, position, temperature key, target temperature)
    "minion": (512, "MINION_MK1D", "FLO-MIN114", "sequencing_MIN114_DNA_e8_2_400K", "MD-000001", "minion", 35.0),
    "promethion": (3000, "P2_SOLO", "FLO-PRO114M", "sequencing_PRO114_DNA_e8_2_400K", "P2S-00000-A", "promethion", 34.0),
}
_KIT_CATEGORY = {"SQK-PCB": ["RNA", "PCR", "Multiplexing"], "SQK-PCS": ["RNA", "PCR"], "SQK-NBD": ["DNA", "PCR-Free", "Multiplexing"],
                 "SQK-LSK": ["DNA", "PCR-Free"], "SQK-RBK": ["DNA", "Multiplexing"]}


@dataclass
class RunSpec:
    """What to fake. Every default is an invented value."""
    flow_cell_id: str = "TST00001"
    family: str = "minion"                  # minion (512 channels) or promethion (3,000)
    kit: str = "SQK-PCB114-24"              # a cDNA kit by default: the judged outcome is then reads
    experiment: str = "demo_experiment"
    sample: str = "demo_sample"
    start: str = "2026-03-02T09:00:00Z"     # protocol start, UTC
    hours: float = 8.0                      # sequencing time; the planned duration is the same unless planned_hours is set
    planned_hours: float | None = None
    reads: int = 4000
    seed: int = 1
    minknow: str = "6.10.12"
    sample_rate: int = 5000
    events_to_base_ratio: float = 0.6
    first_scan_pore_frac: float = 0.72      # wells with a pore at the first scan
    occupancy: float = 0.9                  # strand + adapter over strand + adapter + pore
    scan_period_h: float = 1.5
    state: str = "PROTOCOL_COMPLETED"       # how the protocol ended
    disk_warning: bool = False              # add a MinKNOW warning about disk space (a notice about the computer)
    host: str = "fakehost-01"
    # what the family table fixes, for a run of another kind (a direct RNA flow cell, a second position); None = the family's own
    product: str | None = None              # flow cell product code, for example FLO-MIN004RA
    protocol_script: str | None = None      # the name of the sequencing script in the protocol identifier
    position: str | None = None             # the device position, which is also part of the folder name
    device_type: str | None = None
    speed_bps: int = 400                    # the speed the kit aims at; the target band is this value plus and minus an eighth


@dataclass
class FixtureRun:
    """Where the fake run is and what a correct reader must find in it."""
    run_dir: str
    run_key: str
    flow_cell_id: str
    protocol_run_id: str
    acquisition_id: str
    files: dict[str, str]
    expected: dict
    spec: RunSpec = field(repr=False, default_factory=RunSpec)


def _hex(rng: np.random.Generator, n: int) -> str:
    return "".join(rng.choice(list("0123456789abcdef"), n))


def _uuid(rng: np.random.Generator) -> str:
    h = _hex(rng, 32)
    return f"{h[:8]}-{h[8:12]}-4{h[13:16]}-a{h[17:20]}-{h[20:32]}"


def _iso(epoch: float, frac: bool = True) -> str:
    import datetime as dt
    t = dt.datetime.fromtimestamp(epoch, dt.UTC)
    return t.strftime("%Y-%m-%dT%H:%M:%S") + (".%06dZ" % t.microsecond if frac else "Z")


def _epoch(iso: str) -> float:
    import datetime as dt
    return dt.datetime.strptime(iso, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.UTC).timestamp()


def _write(path: str, text: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write(text)


def _largest_remainder(shares: np.ndarray, total: int) -> np.ndarray:
    exact = shares / shares.sum() * total
    fl = np.floor(exact).astype(np.int64)
    rem = int(total - fl.sum())
    order = np.argsort(-(exact - fl), kind="stable")
    fl[order[:rem]] += 1
    return fl


def make_run(out_dir: str | os.PathLike, spec: RunSpec | None = None, **overrides) -> FixtureRun:
    """Write one fake run folder under ``out_dir`` (``<experiment>/<sample>/<date>_<time>_<position>_<fc>_<run8>/``)."""
    spec = spec or RunSpec()
    if overrides:
        spec = RunSpec(**{**asdict(spec), **overrides})
    if spec.family not in _FAMILIES:
        raise ValueError(f"family must be one of {sorted(_FAMILIES)}")
    channels, device_type, product, script, position, temp_key, temp_target = _FAMILIES[spec.family]
    device_type, product = spec.device_type or device_type, spec.product or product
    script, position = spec.protocol_script or script, spec.position or position
    rng = np.random.default_rng(spec.seed)
    prot_id = _uuid(rng)
    acq_id, cal_id = _hex(rng, 40), _hex(rng, 40)
    run8, acq8 = prot_id[:8], acq_id[:8]
    fc = spec.flow_cell_id
    t_start = _epoch(spec.start)
    t0 = t_start + 252.0                      # data reading starts after the first pore scan
    minutes = int(round(spec.hours * 60))
    run_s = minutes * 60
    planned_s = int(round((spec.planned_hours or spec.hours) * 3600))
    stamp = spec.start[:10].replace("-", "") + "_" + spec.start[11:16].replace(":", "")
    run_dir = os.path.join(os.fspath(out_dir), spec.experiment, spec.sample, f"{stamp}_{position}_{fc}_{run8}")
    os.makedirs(os.path.join(run_dir, "other_reports"), exist_ok=True)
    os.makedirs(os.path.join(run_dir, "pod5"), exist_ok=True)
    protocol = f"sequencing/{script}:{product}:{spec.kit}:{spec.speed_bps}"
    files: dict[str, str] = {}

    # ---------------------------------------------------------------- reads -> sequencing summary
    n = spec.reads
    decay = np.exp(-np.linspace(0, 1.2, 4096))            # fewer reads late in the run
    start = np.sort(rng.choice(np.linspace(0, run_s - 30, 4096), size=n, p=decay / decay.sum()) + rng.uniform(0, 5, n))
    reason = rng.choice(len(END_REASONS), size=n, p=_END_REASON_P)
    length = np.clip(rng.lognormal(mean=np.log(900), sigma=0.75, size=n), 60, 60000)       # estimated bases
    length[reason == 2] = np.clip(length[reason == 2] * 0.15, 20, None)                    # ejected reads are short
    events = np.maximum(np.rint(length / spec.events_to_base_ratio), 1).astype(np.int64)
    speed = np.clip(rng.normal(370, 18, n), 250, 480)                                      # estimated bases per second
    duration = np.round(events * spec.events_to_base_ratio / speed + rng.uniform(0.05, 0.4, n), 4)
    duration = np.minimum(duration, run_s - start - 0.5).clip(0.2)
    channel = rng.integers(1, channels + 1, size=n)
    mux = rng.choice([1, 2, 3, 4], size=n, p=[0.4, 0.3, 0.2, 0.1])
    ss_name = f"sequencing_summary_{fc}_{run8}_{acq8}.txt"
    buf = io.StringIO()
    buf.write("\t".join(["filename_fastq", "filename_fast5", "filename_pod5", "filename_bam", "parent_read_id", "read_id",
                         "run_id", "channel", "mux", "minknow_events", "start_time", "duration", "pore_type",
                         "experiment_id", "sample_id", "end_reason"]) + "\n")
    read_ids = [_uuid(rng) for _ in range(n)]
    for i in range(n):
        buf.write(f"-\t-\t{fc}_{run8}_{acq8}_{i // 4000}.pod5\t-\t{read_ids[i]}\t{read_ids[i]}\t{prot_id}\t{channel[i]}\t{mux[i]}\t"
                  f"{events[i]}\t{start[i]:.6f}\t{duration[i]:.6f}\tnot_set\t{spec.experiment}\t{spec.sample}\t{END_REASONS[reason[i]]}\n")
    files["minknow.sequencing_summary"] = os.path.join(run_dir, ss_name)
    _write(files["minknow.sequencing_summary"], buf.getvalue())
    _write(os.path.join(run_dir, "pod5", f"{fc}_{run8}_{acq8}_0.pod5"), "")

    # ---------------------------------------------------------------- throughput (cumulative, one row per minute)
    end_min = np.minimum(((start + duration) // 60).astype(np.int64), minutes - 1)
    reads_min = np.bincount(end_min, minlength=minutes)
    events_min = np.bincount(end_min, weights=events, minlength=minutes).astype(np.int64)
    samples_min = np.bincount(end_min, weights=np.rint(duration * spec.sample_rate), minlength=minutes).astype(np.int64)
    cum_r, cum_e, cum_s = np.cumsum(reads_min), np.cumsum(events_min), np.cumsum(samples_min)
    cum_b = np.rint(cum_e * spec.events_to_base_ratio).astype(np.int64)
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["Experiment Time (minutes)", "Reads", "Basecalled Reads Passed", "Basecalled Reads Failed", "Basecalled Reads Skipped",
                "Selected Raw Samples", "Selected Events", "Estimated Bases", "Basecalled Bases", "Basecalled Samples"])
    for m in range(minutes):
        w.writerow([m, cum_r[m], 0, 0, 0, cum_s[m], cum_e[m], cum_b[m], 0, 0])
    files["minknow.throughput"] = os.path.join(run_dir, f"throughput_{fc}_{run8}_{acq8}.csv")
    _write(files["minknow.throughput"], buf.getvalue())

    # ---------------------------------------------------------------- pore scans (report + table)
    wells = 4
    n_scans = int(spec.hours / spec.scan_period_h) + 1
    scan_t = [int(k * spec.scan_period_h * 3600) for k in range(n_scans)]
    total_wells = channels * wells
    well_state = rng.choice(len(SCAN_CATEGORIES), size=total_wells,
                            p=_scan_p(spec.first_scan_pore_frac))        # state at the first scan
    scan_counts, scan_rows = [], []
    ps_cols = ["channel", "well", "flow_cell_id", "mux_scan_assessment", "pore_median", "pore_sd", "scan_number",
               "scan_offset_voltage_mv", "scan_running_voltage_mv", "seconds_since_start_of_run", "total_good_samples"]
    ch_idx = np.repeat(np.arange(1, channels + 1), wells)
    well_idx = np.tile(np.arange(1, wells + 1), channels)
    for k in range(n_scans):
        if k:
            alive = well_state <= 1
            die = alive & (rng.random(total_wells) < 0.06)
            well_state = np.where(die, rng.choice([2, 4, 6], size=total_wells), well_state)
        counts = np.bincount(well_state, minlength=len(SCAN_CATEGORIES))
        scan_counts.append({c: int(counts[i]) for i, c in enumerate(SCAN_CATEGORIES)})
        volt = -180 - 5 * (k // 3)
        good = well_state <= 1
        med = np.where(good, np.round(rng.normal(195 - k, 9, total_wells), 3), np.nan)
        sd = np.where(good, np.round(np.abs(rng.normal(2.2, 0.4, total_wells)), 3), np.nan)
        for i in range(total_wells):
            scan_rows.append((ch_idx[i], well_idx[i], fc, SCAN_CATEGORIES[well_state[i]],
                              "" if np.isnan(med[i]) else med[i], "" if np.isnan(sd[i]) else sd[i], k + 1, 10 + k, volt,
                              scan_t[k] + 3, 50000 if good[i] else 0))
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(ps_cols)
    w.writerows(scan_rows)
    files["minknow.pore_scan"] = os.path.join(run_dir, "other_reports", f"pore_scan_data_{fc}_{run8}_{acq8}.csv")
    _write(files["minknow.pore_scan"], buf.getvalue())
    pores = [c["single_pore"] + c["reserved_pore"] for c in scan_counts]

    # ---------------------------------------------------------------- pore activity (samples per state per minute)
    per_min_total = channels * 60 * spec.sample_rate
    t = np.arange(minutes) / max(minutes - 1, 1)
    working = min(pores[0] / total_wells * 1.25, 0.92) * (1 - 0.35 * t)      # share of channels with a working pore
    seq = working * spec.occupancy
    shares = np.zeros((len(CHANNEL_STATES), minutes))
    shares[0] = seq * 0.93                       # strand
    shares[1] = seq * 0.07                       # adapter
    shares[2] = working * (1 - spec.occupancy)   # pore (open, waiting)
    shares[3] = 0.03 + 0.02 * t                  # unavailable
    shares[4] = 0.004                            # unblocking
    shares[6] = 0.02 + 0.10 * t                  # zero
    shares[7] = 0.01 + 0.05 * t                  # saturated
    shares[8] = 0.01                             # multiple
    shares[14] = 0.002                           # unclassified
    shares[5] = np.clip(1 - shares.sum(0), 0.0, None)   # no_pore takes the rest
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["Channel State", "Experiment Time (minutes)", "State Time (samples)"])
    state_samples = np.zeros((len(CHANNEL_STATES), minutes), np.int64)
    for m in range(minutes):
        state_samples[:, m] = _largest_remainder(shares[:, m], per_min_total)
    for si in np.argsort(CHANNEL_STATES, kind="stable"):                 # MinKNOW writes the states in name order
        for m in range(minutes):
            w.writerow([CHANNEL_STATES[si], m, int(state_samples[si, m])])
    files["minknow.pore_activity"] = os.path.join(run_dir, f"pore_activity_{fc}_{run8}_{acq8}.csv")
    _write(files["minknow.pore_activity"], buf.getvalue())
    seq_s = state_samples[0] + state_samples[1]
    pore_s = state_samples[2]
    lo, hi = min(30, minutes), min(90, minutes)
    occ_first = float(seq_s[lo:hi].sum() / max((seq_s[lo:hi] + pore_s[lo:hi]).sum(), 1)) if hi > lo else None

    # ---------------------------------------------------------------- small text files
    ended = t0 + run_s
    final = {
        "instrument": spec.host, "position": position, "flow_cell_id": fc, "sample_id": spec.sample,
        "protocol_group_id": spec.experiment, "protocol": protocol, "protocol_run_id": prot_id, "acquisition_run_id": acq_id,
        "started": _iso(t_start).replace("Z", "+00:00"), "acquisition_stopped": _iso(ended).replace("Z", "+00:00"),
        "processing_stopped": _iso(ended + 4).replace("Z", "+00:00"), "basecalling_enabled": "0",
        "sequencing_summary_file": ss_name, "fast5_files_in_final_dest": "0", "fast5_files_in_fallback": "0",
        "pod5_files_in_final_dest": str(n // 4000 + 1), "pod5_files_in_fallback": "0", "fastq_files_in_final_dest": "0",
        "fastq_files_in_fallback": "0", "bam_files_in_final_dest": "0", "bam_files_in_fallback": "0",
        "bai_files_in_final_dest": "0", "bai_files_in_fallback": "0",
    }
    files["minknow.final_summary"] = os.path.join(run_dir, f"final_summary_{fc}_{run8}_{acq8}.txt")
    _write(files["minknow.final_summary"], "".join(f"{k}={v}\n" for k, v in final.items()))
    files["minknow.sample_sheet"] = os.path.join(run_dir, f"sample_sheet_{fc}_{stamp}_{run8}.csv")
    _write(files["minknow.sample_sheet"],
           "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,flow_cell_product_code,kit\n"
           f"{prot_id},{position},{fc},{spec.sample},{spec.experiment},{product},{spec.kit}\n")
    files["minknow.barcode_alignment"] = os.path.join(run_dir, f"barcode_alignment_{fc}_{run8}_{acq8}.tsv")
    _write(files["minknow.barcode_alignment"],
           "barcode\talias\ttype\ttarget_unclassified\tacquisition_run_id\tprotocol_group_id\tsample_id\tflow_cell_id\tstarted\n"
           f"unclassified\tunclassified\tna\t{n}\t{acq_id}\t{spec.experiment}\t{spec.sample}\t{fc}\t{_iso(t_start).replace('Z', '+00:00')}\n")

    # ---------------------------------------------------------------- report JSON
    major, minor, patch = (int(x) for x in spec.minknow.split("."))
    kit_cat = next((v for k, v in _KIT_CATEGORY.items() if spec.kit.startswith(k)), ["DNA"])
    messages = [
        {"time": _iso(t_start - 2), "severity": "MESSAGE_SEVERITY_INFO", "identifier": "disk_space_info_user",
         "user_message": "Disk /data has 4000 GB space remaining",
         "extra_data": {"space_remaining": "4000", "filesystem_id": "/data", "unit": "GB"}},
        {"time": _iso(t_start), "severity": "MESSAGE_SEVERITY_INFO", "identifier": "base_device.wait_for_temperature",
         "user_message": f"Waiting up to 300 seconds for temperature to stabilise at {temp_target}°C",
         "extra_data": {"timeout": "300", "target": str(temp_target)}},
        {"time": _iso(t_start + 40), "severity": "MESSAGE_SEVERITY_INFO", "identifier": "sequencing.start",
         "user_message": "Starting sequencing procedure"},
    ]
    if spec.disk_warning:
        messages.append({"time": _iso(t0 + 3600), "severity": "MESSAGE_SEVERITY_WARNING", "identifier": "disk_usage_alert",
                         "user_message": "Disk /data is running low on space",
                         "extra_data": {"space_remaining": "90", "filesystem_id": "/data", "unit": "GB"}})
    ch_with_pore = []
    ws = np.asarray([SCAN_CATEGORIES.index(r[3]) for r in scan_rows], dtype=np.int64).reshape(n_scans, channels, wells)
    for k in range(n_scans):
        ch_with_pore.append(int(((ws[k] <= 1).any(axis=1)).sum()))
        messages.append({"time": _iso(t0 + scan_t[k] - 200), "severity": "MESSAGE_SEVERITY_INFO",
                         "identifier": "sequencing.mux_scan", "user_message": "Performing Pore Scan"})
        messages.append({"time": _iso(t0 + scan_t[k]), "severity": "MESSAGE_SEVERITY_INFO", "identifier": "mux_scan_result",
                         "user_message": f"Pore scan for flow cell {fc} has found a total of {pores[k]} pores. "
                                         f"{ch_with_pore[k]} pores available for sequencing immediately",
                         "extra_data": {"flow_cell_id": fc, "num_pores": str(ch_with_pore[k]), "total_pores": str(pores[k])}})
    bucket_min = 120
    n_temp = int(run_s // (bucket_min * 60)) + 1
    temps = []
    for k in range(n_temp):
        v = temp_target + float(rng.normal(0, 0.004))
        if k == n_temp - 1:
            v = temp_target - 2.1            # the sample taken as the run stops reads low: a reader must leave it out
        second = "chamber_temperature" if temp_key == "promethion" else "asic_temperature"
        temps.append({temp_key: {("flowcell_temperature" if temp_key == "promethion" else "heatsink_temperature"): v,
                                 second: 37.7 if temp_key == "promethion" else 455.0},
                      "target_temperature": {"minimum": temp_target, "maximum": temp_target}})
    report = {
        "host": {"serial": spec.host},
        "protocol_run_info": {
            "run_id": prot_id, "protocol_id": protocol, "args": ["--pod5=on", "--base_calling=off", f"--mux_scan_period={spec.scan_period_h}"],
            "state": spec.state, "start_time": _iso(t_start), "end_time": _iso(ended + 6),
            "acquisition_run_ids": [cal_id, acq_id],
            "user_info": {"protocol_group_id": spec.experiment, "sample_id": spec.sample, "user_specified_flow_cell_id": "",
                          "user_specified_product_code": "", "kit_info": {"sequencing_kit": spec.kit}},
            "output_path": run_dir,
            "device": {"device_id": position, "device_type": device_type, "max_channel_count": channels,
                       "max_wells_per_channel": wells, "can_set_temperature": True, "digitisation": 2048},
            "flow_cell": {"has_flow_cell": True, "channel_count": channels, "wells_per_channel": wells, "flow_cell_id": fc,
                          "product_code": product, "user_specified_flow_cell_id": "", "user_specified_product_code": "",
                          "asic_id_str": "00" + _hex(rng, 14)},
            "meta_info": {"identifier": protocol, "name": script, "tags": {
                "barcoding": {"bool_value": spec.kit.endswith(("-24", "-96"))}, "experiment type": {"string_value": "sequencing"},
                "barcoding kits": {"array_value": "[]"}, "flow cell": {"string_value": product},
                "kit_category": {"array_value": json.dumps(kit_cat, separators=(",", ":"))},
                "sample rate": {"int_value": str(spec.sample_rate)}, "kit": {"string_value": spec.kit}, "speed": {"int_value": str(spec.speed_bps)}}},
            "software_versions": {"minknow": {"major": major, "minor": minor, "patch": patch, "full": spec.minknow},
                                  "bream": "8.10.13", "distribution_version": "26.01.15", "distribution_status": "STABLE",
                                  "protocol_configuration": spec.minknow, "basecaller_connected_version": "7.13.6+0000000"},
            "stop_origin": {"run_until": {"criteria": "runtime"}} if spec.state == "PROTOCOL_COMPLETED" else {"user": {}},
            "phase_history": [{"phase": "PHASE_INITIALISING", "timestamp": _iso(t_start), "origin": {"protocol_phase_management": {}}},
                              {"phase": "PHASE_MUX_SCAN", "timestamp": _iso(t_start + 40), "origin": {"protocol_phase_management": {}}},
                              {"phase": "PHASE_SEQUENCING", "timestamp": _iso(t0), "origin": {"protocol_phase_management": {}}}],
        },
        "user_messages": sorted(messages, key=lambda m: m["time"]),
        "acquisitions": [
            {"acquisition_run_info": {"run_id": cal_id, "state": "ACQUISITION_COMPLETED", "stop_reason": "STOPPED_USER_REQUESTED",
                                      "start_time": _iso(t_start + 5), "data_read_start_time": _iso(t_start + 6),
                                      "data_read_end_time": _iso(t_start + 12), "end_time": _iso(t_start + 13), "yield_summary": {},
                                      "config_summary": {"sample_rate": spec.sample_rate, "channel_count": channels, "purpose": "CALIBRATION"}},
             "bucket_size": 1},
            {"acquisition_run_info": {
                "run_id": acq_id, "state": "ACQUISITION_COMPLETED",
                "stop_reason": "STOPPED_PROTOCOL_ENDED" if spec.state == "PROTOCOL_COMPLETED" else "STOPPED_USER_REQUESTED",
                "start_time": _iso(t0 - 1), "data_read_start_time": _iso(t0), "data_read_end_time": _iso(ended), "end_time": _iso(ended + 2),
                "yield_summary": {"read_count": str(n), "selected_raw_samples": str(int(cum_s[-1])),
                                  "selected_events": str(int(cum_e[-1])), "estimated_selected_bases": str(int(cum_b[-1]))},
                "config_summary": {"reads_directory": run_dir, "events_to_base_ratio": spec.events_to_base_ratio,
                                   "sample_rate": spec.sample_rate, "channel_count": channels, "purpose": "SEQUENCING",
                                   "pod5_reads_enabled": True},
                "writer_summary": {"bytes_to_write_produced": str(n * 9000), "bytes_to_write_completed": str(n * 9000)},
                "bream_info": {
                    "mux_scan_metadata": {"auto_mux_scan_period_hours": spec.scan_period_h, "enable_reserve_pore": True},
                    "mux_scan_results": [{"counts": scan_counts[k], "mux_scan_timestamp": str(scan_t[k])} for k in range(n_scans)],
                    "target_translocation_speed": {"minimum": round(spec.speed_bps * 0.875), "maximum": round(spec.speed_bps * 1.125)},
                    "target_temperature": {"minimum": temp_target - 0.1, "maximum": temp_target + 0.1}},
                "target_run_until_criteria": {"stop_criteria": {"criteria": {"runtime": {
                    "@type": "type.googleapis.com/google.protobuf.UInt64Value", "value": str(planned_s)}}}}},
             "bucket_size": bucket_min,
             "temperature": [{"temperatures": temps}],
             "bias_voltage": [{"bias_voltages": [{}] + [{"acquisition_index": str(scan_t[k] * spec.sample_rate),
                                                         "bias_voltage": -180 - 5 * (k // 3), "time_seconds": str(scan_t[k])}
                                                        for k in range(n_scans)]}]},
        ],
        "report_data_generation_time": _iso(ended + 8),
    }
    files["minknow.report_json"] = os.path.join(run_dir, f"report_{fc}_{stamp}_{run8}.json")
    _write(files["minknow.report_json"], json.dumps(report, indent=1))

    def at_h(cum, h):
        return int(cum[h * 60 - 1]) if minutes >= h * 60 else None
    expected = {
        "run_key": f"{fc}_{run8}", "reads": n, "minknow_events": int(events.sum()),
        "estimated_bases": int(cum_b[-1]), "selected_raw_samples": int(cum_s[-1]),
        "channel_count": channels, "family": spec.family, "device_type": device_type, "product_code": product, "kit": spec.kit,
        "assay_group": "cdna" if spec.kit.startswith(("SQK-PCB", "SQK-PCS", "SQK-DCS")) else "gdna",
        "run_hours": minutes / 60, "planned_s": planned_s, "time_origin": _iso(t0), "started": _iso(t_start),
        "pores_per_scan": pores, "pores_first_scan": pores[0], "channels_with_pore_first_scan": ch_with_pore[0],
        "pore_scans": n_scans, "pore_scan_rows": len(scan_rows), "occupancy_first_hour": occ_first,
        "reads_by_end_reason": {END_REASONS[i]: int((reason == i).sum()) for i in range(len(END_REASONS))},
        "reads_at_h": {h: at_h(cum_r, h) for h in (6, 12, 24, 48, 72)},
        "channels_with_reads": int(len(np.unique(channel))),
        "channel_samples_per_minute": per_min_total, "minutes": minutes,
        "temperature_samples": n_temp, "temperature_last_sample_low": True,
        "messages_warning": 1 if spec.disk_warning else 0, "state": spec.state,
        "summary_columns": 16, "basecalling": False,
    }
    return FixtureRun(run_dir=run_dir, run_key=f"{fc}_{run8}", flow_cell_id=fc, protocol_run_id=prot_id, acquisition_id=acq_id,
                      files=files, expected=expected, spec=spec)


def _scan_p(pore_frac: float) -> list[float]:
    """Probabilities of the seven scan results at the first scan, given the share of wells with a pore."""
    rest = 1.0 - pore_frac
    return [pore_frac * 0.97, pore_frac * 0.03, rest * 0.30, rest * 0.25, rest * 0.30, rest * 0.10, rest * 0.05]


def make_library_tree(out_dir: str | os.PathLike, *, seed: int = 1) -> list[FixtureRun]:
    """A small tree with two roots and three runs on two flow cells, for discovery and scan tests:

    ``root_a/project_one/…`` holds two runs of one flow cell (a first loading and a later one);
    ``root_b/project_two/…`` holds one run of another flow cell.
    """
    out = os.fspath(out_dir)
    a = make_run(os.path.join(out, "root_a"), RunSpec(flow_cell_id="TST00001", experiment="project_one", sample="pool_1",
                                                      start="2026-03-02T09:00:00Z", hours=8, reads=3000, seed=seed))
    b = make_run(os.path.join(out, "root_a"), RunSpec(flow_cell_id="TST00001", experiment="project_one", sample="pool_1_reload",
                                                      start="2026-03-03T10:00:00Z", hours=7, reads=2000, seed=seed + 1,
                                                      first_scan_pore_frac=0.35))
    c = make_run(os.path.join(out, "root_b"), RunSpec(flow_cell_id="TST00002", experiment="project_two", sample="pool_2",
                                                      start="2026-03-05T08:30:00Z", hours=8, reads=3000, seed=seed + 2,
                                                      kit="SQK-NBD114-24"))
    return [a, b, c]


def main(argv: list[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(prog="python -m flowcellar.testing.fixtures",
                                description="Write a fake MinKNOW run folder (invented names only) for tests and demonstrations.")
    p.add_argument("out_dir")
    p.add_argument("--reads", type=int, default=4000)
    p.add_argument("--hours", type=float, default=8.0)
    p.add_argument("--family", choices=sorted(_FAMILIES), default="minion")
    p.add_argument("--kit", default="SQK-PCB114-24")
    p.add_argument("--flow-cell-id", default="TST00001")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--tree", action="store_true", help="write the two-root tree of three runs instead of one run")
    a = p.parse_args(argv)
    if a.tree:
        runs = make_library_tree(a.out_dir, seed=a.seed)
    else:
        runs = [make_run(a.out_dir, RunSpec(reads=a.reads, hours=a.hours, family=a.family, kit=a.kit,
                                            flow_cell_id=a.flow_cell_id, seed=a.seed))]
    for r in runs:
        print(json.dumps({"run_key": r.run_key, "run_dir": r.run_dir, "reads": r.expected["reads"],
                          "pores_first_scan": r.expected["pores_first_scan"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
