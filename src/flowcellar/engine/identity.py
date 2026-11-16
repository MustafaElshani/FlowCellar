"""Which run a folder holds: recognition by its files, identity fields with their precedence, assay detection,
run state.

Precedence inside a run folder: report JSON, final summary, sample sheet, file names, folder name. The date in a
folder name is never used for a run that has identity files: it is local time.
"""
from __future__ import annotations

import json
import os
import re
import time

from flowcellar.contracts import RunIdentity
from flowcellar.model.runkey import FLOW_CELL_RE
from flowcellar.sources import minknow as mk
from flowcellar.sources.minknow import g, to_float, to_int
from flowcellar.util.errors import InputError
from flowcellar.util.text import clean_text
from flowcellar.util.timeutil import from_epoch, iso_z, parse_ts

HEX_ID = re.compile(r"^[0-9a-f]{8}[0-9a-f-]*$")
_FOLDER_TIME = re.compile(r"^(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})_")

# kit prefix -> (assay type, group, subtype, molecule); first match wins, the direct RNA kits first
_KIT_RULES = (
    (r"^SQK-(RNA|DRB)\d{3}", "direct_rna", "drna", "native", "rna"),
    (r"^SQK-(PCS|PCB)\d{3}", "cdna_pcr", "cdna", "pcr_cdna", "dna"),
    (r"^SQK-DCS\d{3}", "cdna_direct", "cdna", "direct_cdna", "dna"),
    (r"^SQK-(PBK|RPB|16S|RAB)", "amplicon", "gdna", "amplicon", "dna"),
    (r"^SQK-(LSK|NBD|RBK|RAD|ULK|MLK|CAS)", "gdna", "gdna", "native", "dna"),
)
_OVERRIDE = {"gdna": ("gdna", "native", "dna"), "amplicon": ("gdna", "amplicon", "dna"), "cdna_pcr": ("cdna", "pcr_cdna", "dna"),
             "cdna_direct": ("cdna", "direct_cdna", "dna"), "direct_rna": ("drna", "native", "rna")}
_LAYOUTS = {"promethion": "promethion_3000", "minion": "minion_512", "flongle": "flongle_126"}
_PURPOSES = {"SEQUENCING": "sequencing", "CALIBRATION": "calibration", "FLOW_CELL_CHECK": "flow_cell_check"}


def family_of(channels: int | None) -> str | None:
    if not channels:
        return None
    return "promethion" if channels >= 2000 else "minion" if channels >= 400 else "flongle"


def chemistry_of(product: str | None) -> str | None:
    p = product or ""
    if re.search(r"004RA", p):
        return "rna004"
    if "114" in p:
        return "r10.4.1"
    if re.search(r"MIN106|MIN107|PRO002|FLG001", p):
        return "r9.4.1"
    return None


def detect_assay(kit: str | None, product: str | None, protocol: str | None, kit_category: list[str],
                 barcoding: bool | None, override: str | None = None) -> dict:
    """A stated assay (``--assay``) always wins."""
    ev: list[str] = []
    a = {"type": "unknown", "group": "unknown", "chemistry": chemistry_of(product), "molecule": "unknown", "barcoded": None,
         "source": "unknown", "confidence": "low", "evidence": ev, "subtype": None}
    script = (protocol or "").split(":")[0].split("/")[-1]
    rna_evidence = bool(re.search(r"_RNA(_|$)", script)) or bool(product and product.endswith("RA"))
    cat0 = kit_category[0] if kit_category else None
    signals = 0
    for pattern, typ, grp, sub, mol in _KIT_RULES:
        if kit and re.search(pattern, kit):
            a.update(type=typ, group=grp, subtype=sub, molecule=mol, source="kit_table")
            ev.append(f"kit {kit}")
            signals = 1
            break
    if a["type"] != "unknown":
        if cat0 in ("RNA", "DNA"):
            if (cat0 == "RNA") == (a["group"] in ("cdna", "drna")):
                signals += 1
                a["source"] = "run_files"
                ev.append("kit category " + "/".join(kit_category))
            else:
                ev.append("kit category " + "/".join(kit_category) + " does not match the kit table")
        if a["group"] == "drna" and rna_evidence:
            signals += 1
    elif rna_evidence:
        a.update(type="direct_rna", group="drna", subtype="native", molecule="rna", source="run_files")
        signals = 1 + (cat0 == "RNA")
    elif cat0 == "RNA":
        a.update(type="cdna_pcr", group="cdna", subtype="pcr_cdna", molecule="dna", source="run_files")
        ev.append("kit category " + "/".join(kit_category))
        signals = 1
    elif cat0 == "DNA":
        a.update(type="gdna", group="gdna", subtype="native", molecule="dna", source="run_files")
        ev.append("kit category " + "/".join(kit_category))
        signals = 1
    if script:
        ev.append(f"protocol {script}")
    if product:
        ev.append(f"flow cell product {product}")
    if a["type"] != "unknown":
        a["confidence"] = "high" if signals >= 2 else "medium"
    if kit:
        a["barcoded"] = bool(re.search(r"-(24|96)$", kit)) or bool(barcoding)
    elif barcoding is not None:
        a["barcoded"] = bool(barcoding)
    if override:
        grp, sub, mol = _OVERRIDE[override]
        a.update(type=override, group=grp, subtype=sub, molecule=mol, source="user", confidence="certain")
        ev.append("stated on the command line")
    return a


def _stop_origin_text(so) -> str | None:
    if not so:
        return None
    if isinstance(so, dict):
        k = next(iter(so), None)
        v = so.get(k)
        if isinstance(v, dict) and v:
            inner = next(iter(v.values()))
            return f"{k}:{inner}" if isinstance(inner, str) and len(inner) < 64 and "[" not in inner else str(k)
        return str(k) if k else None
    return clean_text(so, 200)


def _tag_array(tags: dict, key: str) -> list[str]:
    try:
        v = json.loads(g(tags, key, "array_value", default="[]") or "[]")
    except (ValueError, TypeError):
        return []
    return [clean_text(x, 200) for x in v] if isinstance(v, list) else []


def _protocol_fields(protocol: str | None) -> tuple[str | None, str | None]:
    """(flow cell product, kit) from ``sequencing/<script>:<product>:<kit>:<speed>``."""
    parts = (protocol or "").split(":")
    product = parts[1] if len(parts) > 1 and parts[1].startswith("FLO-") else None
    kit = parts[2] if len(parts) > 2 and parts[2].startswith("SQK-") else None
    return product, kit


def _valid_fc(x) -> str | None:
    x = (x or "").strip() if isinstance(x, str) else None
    return x if x and FLOW_CELL_RE.match(x) else None


def _clean_id(x) -> str | None:
    x = (x or "").strip().lower() if isinstance(x, str) else None
    return x if x and HEX_ID.match(x) else None


def build_identity(listing: mk.Listing, rep: dict, fs: dict, sheet: list[dict], *, assay_override: str | None = None,
                   stale_after_hours: float = 6.0, now: float | None = None) -> tuple[dict, list[dict]]:
    """(the ``identity`` member of a bundle, warnings). ``rep`` may be a whole report, only its head, or empty."""
    warnings: list[dict] = []
    pri = rep.get("protocol_run_info") if isinstance(rep.get("protocol_run_info"), dict) else {}
    sw = mk.software_versions(rep)
    fc = pri.get("flow_cell") if isinstance(pri.get("flow_cell"), dict) else {}
    dev = pri.get("device") if isinstance(pri.get("device"), dict) else {}
    row = sheet[0] if sheet else {}
    names_key = mk.key_from_names(listing)

    # ---- flow cell and run ids, by precedence; a disagreement is kept as a warning and the files win
    fc_sources = [("report", _valid_fc(fc.get("flow_cell_id")) or _valid_fc(fc.get("user_specified_flow_cell_id"))),
                  ("final summary", _valid_fc(fs.get("flow_cell_id"))), ("sample sheet", _valid_fc(row.get("flow_cell_id"))),
                  ("file names", names_key[0] if names_key else None)]
    fc_id = next((v for _, v in fc_sources if v), None)
    prot_sources = [("report", _clean_id(pri.get("run_id"))), ("final summary", _clean_id(fs.get("protocol_run_id"))),
                    ("sample sheet", _clean_id(row.get("protocol_run_id")))]
    prot_id = next((v for _, v in prot_sources if v), None)
    short = prot_id.replace("-", "")[:8] if prot_id else (names_key[1] if names_key else None)
    seen_fc = {v for _, v in fc_sources if v}
    if len(seen_fc) > 1:
        warnings.append({"level": "warning", "code": "identity_sources_differ",
                         "message": "The files of this folder name more than one flow cell (" + ", ".join(sorted(seen_fc)) +
                                    f"); {fc_id} is used, in the order report, final summary, sample sheet, file names."})
    if names_key and short and names_key[1] != short:
        warnings.append({"level": "warning", "code": "identity_sources_differ",
                         "message": f"File names carry the run id {names_key[1]}, the files themselves {short}; the files are used."})
    if not fc_id or not short:
        raise InputError(f"{listing.run_dir}: the flow cell and the run could not be told from the files of this folder",
                         hint="a run folder holds a report, a final summary or a sample sheet written by MinKNOW")

    channels = to_int(fc.get("channel_count"), 0) or to_int(dev.get("max_channel_count"), 0) or None
    wells = to_int(fc.get("wells_per_channel"), 0) or to_int(dev.get("max_wells_per_channel"), 0) or None
    protocol = fs.get("protocol") or g(pri, "meta_info", "identifier") or pri.get("protocol_id") or None
    p_product, p_kit = _protocol_fields(protocol)
    product = fc.get("product_code") or fc.get("user_specified_product_code") or p_product or row.get("flow_cell_product_code") or None
    tags = g(pri, "meta_info", "tags", default={}) or {}
    kit = g(pri, "user_info", "kit_info", "sequencing_kit") or g(tags, "kit", "string_value") or p_kit or row.get("kit") or None
    kit_category = _tag_array(tags, "kit_category")
    family = family_of(channels)
    seq, n_seq = mk.sequencing_acquisition(rep)
    ari = (seq or {}).get("acquisition_run_info") or {}
    if n_seq > 1:
        warnings.append({"level": "warning", "code": "several_sequencing_acquisitions",
                         "message": f"The report lists {n_seq} sequencing acquisitions; the one with the most reads is read."})
    if channels is None:
        channels = to_int(g(ari, "config_summary", "channel_count"), 0) or None
        family = family_of(channels)
    t0 = ari.get("data_read_start_time") or ari.get("start_time")
    t0s = parse_ts(t0)
    started = iso_z(pri.get("start_time")) or iso_z(fs.get("started"))
    if not started:
        ref = listing.mtimes.get("minknow.sample_sheet") or listing.oldest_mtime
        if ref:
            started = from_epoch(ref)
            warnings.append({"level": "info", "code": "start_time_from_file_time",
                             "message": "No report and no final summary give the start of this run; the time the sample sheet was written is used."})
    utc_off = None
    m = re.search(r"([+-])(\d{2}):?(\d{2})$", fs.get("started", "") or "")
    if m:
        utc_off = (1 if m.group(1) == "+" else -1) * (int(m.group(2)) * 60 + int(m.group(3)))
    planned = to_float(g(ari, "target_run_until_criteria", "stop_criteria", "criteria", "runtime", "value"))
    if planned is None:
        for a in pri.get("args") or []:
            mm = re.match(r"--experiment_time=(\d+(?:\.\d+)?)", str(a))
            if mm:
                planned = float(mm.group(1)) * 3600
    acq_list = []
    for a in rep.get("acquisitions") or []:
        r = a.get("acquisition_run_info") if isinstance(a, dict) else None
        if not isinstance(r, dict) or not r.get("run_id"):
            continue
        ds = r.get("data_read_start_time")
        purpose = g(r, "config_summary", "purpose")
        acq_list.append({
            "id": str(r["run_id"]), "short_id": str(r["run_id"])[:8],
            "purpose": _PURPOSES.get(purpose, "other" if purpose else "unknown"),
            "started": iso_z(r.get("start_time")), "data_read_started": iso_z(ds), "ended": iso_z(r.get("end_time")),
            "offset_s": round(parse_ts(ds) - t0s, 3) if parse_ts(ds) is not None and t0s is not None else None,
            "stop_reason": r.get("stop_reason"), "sample_rate_hz": to_float(g(r, "config_summary", "sample_rate")),
            "read_count": to_int(g(r, "yield_summary", "read_count"), None)})
    dev_id = dev.get("device_id") or None
    fs_position, fs_instrument = (fs.get("position") or "").strip(), (fs.get("instrument") or "").strip()
    # MinKNOW 6.10 writes instrument = host and position = device; before that instrument = device, position empty
    position = fs_position or dev_id or fs_instrument or row.get("position_id") or None
    serial = g(rep, "host", "serial")
    host = (fs_instrument if fs_position else None) or (serial if isinstance(serial, str) and serial not in ("", "SKU", "unknown") else None)
    device_type = dev.get("device_type") or None
    device_id = dev_id or position
    if device_id and (str(device_type or "").startswith("P2") or re.match(r"^P2[S-]", device_id)):
        device_id = re.sub(r"-[A-Z0-9]$", "", device_id)
    barcoding = g(tags, "barcoding", "bool_value")
    assay = detect_assay(kit, product, protocol, kit_category, barcoding if isinstance(barcoding, bool) else None, assay_override)
    has_final = "minknow.final_summary" in listing.files
    now = time.time() if now is None else now
    if has_final:
        run_state = "finished"
    else:
        run_state = "running" if (now - listing.newest_mtime) < stale_after_hours * 3600 else "interrupted"
    speed = g(tags, "speed", "int_value")
    identity = {
        "run_key": f"{fc_id}_{short}",
        "flow_cell": {"id": fc_id, "user_specified_id": fc.get("user_specified_flow_cell_id") or None, "product_code": product,
                      "asic_id": fc.get("asic_id_str") or None, "channel_count": channels, "wells_per_channel": wells,
                      "channel_layout": _LAYOUTS.get(family) if family else None, "family": family},
        "protocol_run": {"id": prot_id, "short_id": short, "protocol": protocol,
                         "experiment_id": fs.get("protocol_group_id") or g(pri, "user_info", "protocol_group_id") or row.get("experiment_id") or None,
                         "sample_id": fs.get("sample_id") or g(pri, "user_info", "sample_id") or row.get("sample_id") or None,
                         "started": started, "ended": iso_z(pri.get("end_time")) or iso_z(fs.get("processing_stopped")),
                         "time_origin": iso_z(t0), "time_origin_source": "minknow_report" if t0 else None, "utc_offset_min": utc_off,
                         "state": pri.get("state") or None, "stop_origin": _stop_origin_text(pri.get("stop_origin")),
                         "run_state": run_state, "planned_duration_s": planned},
        "acquisitions": acq_list,
        "device": {"type": device_type, "id": device_id, "position": position, "host": host},
        "kit": {"sequencing_kit": kit, "barcoding_kits": _tag_array(tags, "barcoding kits"),
                "speed_bps": to_float(speed), "kit_category": kit_category},
        "assay": assay,
        "software": {"minknow": mk.minknow_version(sw), "bream": sw.get("bream") or None,
                     "distribution": sw.get("distribution_version") or None,
                     "basecaller_live": sw.get("basecaller_connected_version") or sw.get("guppy_connected_version") or None},
        "basecalling": [],
    }
    for sect in ("flow_cell", "protocol_run", "device", "software"):
        for k, v in identity[sect].items():
            if isinstance(v, str):
                identity[sect][k] = clean_text(v, 256)
    return identity, warnings


def raw_only_identity(listing: mk.Listing) -> tuple[dict, list[dict]]:
    """The stub identity of a folder that holds raw data only: flow cell and run from names, nothing else."""
    key = mk.key_from_names(listing)
    if not key:
        raise InputError(f"{listing.run_dir}: holds raw data, but no name in it gives a flow cell and a run id")
    fc_id, short = key
    m = _FOLDER_TIME.match(listing.name)
    if m:
        started = f"{m.group(1)}-{m.group(2)}-{m.group(3)}T{m.group(4)}:{m.group(5)}:00Z"
        how = "the folder name, whose time is the sequencer's local time"
    elif listing.oldest_mtime:
        started = from_epoch(listing.oldest_mtime)
        how = "the oldest file time in the folder"
    else:
        started, how = None, "nothing"
    family = "promethion" if fc_id.startswith("P") else "minion" if fc_id.startswith("F") else None
    channels = {"promethion": 3000, "minion": 512}.get(family)
    parts = listing.name.split("_")
    position = "_".join(parts[2:-2]) if m and len(parts) >= 5 else None
    identity = {
        "run_key": f"{fc_id}_{short}",
        "flow_cell": {"id": fc_id, "product_code": None, "channel_count": channels, "wells_per_channel": None,
                      "channel_layout": _LAYOUTS.get(family) if family else None, "family": family},
        "protocol_run": {"id": None, "short_id": short, "started": started, "ended": None, "time_origin": None,
                         "time_origin_source": None, "run_state": "raw_only"},
        "acquisitions": [], "device": {"type": None, "id": None, "position": position},
        "kit": {}, "assay": {"type": "unknown", "group": "unknown", "molecule": "unknown", "source": "unknown",
                             "confidence": "low", "evidence": []},
        "software": {}, "basecalling": [],
    }
    subs = ", ".join(listing.subdirs) or "none"
    warnings = [{"level": "warning", "code": "raw_only",
                 "message": f"Only raw data was found in this folder (sub-folders: {subs}). The flow cell and the run come from "
                            f"file and folder names; the start time comes from {how}; the device family is taken from the flow cell id."}]
    return identity, warnings


def identify(run_dir: str) -> RunIdentity:
    """ENGINE ``identify``: small files only. Raises InputError when the folder is not a run folder."""
    listing = mk.list_run_folder(run_dir)
    if listing.kind == "lookalike":
        raise InputError(f"{listing.run_dir}: not a run folder: it is named like one but holds only basecaller output",
                         hint="give the folder MinKNOW wrote (with its report and final summary), or pass this one as read data")
    if listing.kind is None:
        raise InputError(f"{listing.run_dir}: not a MinKNOW run folder (no report, final summary, sample sheet, throughput, "
                         "pore activity or sequencing summary, and no raw data)")
    if listing.kind == "run_raw_only":
        ident, _ = raw_only_identity(listing)
        fcell = ident["flow_cell"]
        return RunIdentity(run_key=ident["run_key"], flow_cell_id=fcell["id"], run_short=ident["protocol_run"]["short_id"],
                           run_dir=listing.run_dir, kind="run_raw_only", run_state="raw_only", protocol_run_id=None,
                           started=ident["protocol_run"]["started"], files=dict(listing.files), sizes=dict(listing.sizes),
                           newest_mtime=listing.newest_mtime,
                           facts={"device_type": None, "product_code": None, "kit": None, "experiment_id": None, "sample_id": None,
                                  "channel_count": fcell["channel_count"], "family": fcell["family"], "folder": listing.name})
    fs = mk.read_final_summary(listing.files["minknow.final_summary"]) if "minknow.final_summary" in listing.files else {}
    sheet = mk.read_sample_sheet(listing.files["minknow.sample_sheet"]) if "minknow.sample_sheet" in listing.files else []
    head = mk.read_report_head(listing.files["minknow.report_json"]) if "minknow.report_json" in listing.files else {}
    ident, _ = build_identity(listing, head, fs, sheet)
    pr, fcell = ident["protocol_run"], ident["flow_cell"]
    return RunIdentity(run_key=ident["run_key"], flow_cell_id=fcell["id"], run_short=pr["short_id"], run_dir=listing.run_dir,
                       kind="run", run_state=pr["run_state"], protocol_run_id=pr["id"], started=pr["started"],
                       files=dict(listing.files), sizes=dict(listing.sizes), newest_mtime=listing.newest_mtime,
                       facts={"device_type": ident["device"]["type"], "product_code": fcell["product_code"],
                              "kit": ident["kit"]["sequencing_kit"], "experiment_id": pr["experiment_id"], "sample_id": pr["sample_id"],
                              "channel_count": fcell["channel_count"], "family": fcell["family"], "folder": listing.name,
                              "assay_group": ident["assay"]["group"]})


def run_dir_layout(run_dir: str, experiment: str | None, sample: str | None) -> str:
    """``standard`` when the folder still sits under ``<experiment>/<sample>/`` as MinKNOW wrote it, else ``moved``."""
    parts = os.path.abspath(run_dir).rstrip("/").split(os.sep)
    return "standard" if len(parts) >= 3 and parts[-2] == sample and parts[-3] == experiment else "moved"
