"""Section ``run``: how the run went as MinKNOW tells it (messages, files written, how it ended)."""
from __future__ import annotations

from flowcellar.modules.common import VERSIONS, note
from flowcellar.sources.minknow import to_int
from flowcellar.util.text import clean_text
from flowcellar.util.timeutil import iso_z

MAX_MESSAGES = 500
#: MinKNOW messages that only feed the pore-scan series
_SCAN_MESSAGES = ("sequencing.mux_scan", "mux_scan_result")


def is_computer_notice(identifier: str) -> bool:
    """A notice about the computer, not about the flow cell (disk space): listed in the run log, never judged."""
    return "disk" in identifier.lower()


def build_section(rep: dict, ari: dict, fs: dict, *, run_hours: float | None, planned_s: float | None, report_problem: str | None,
                  inputs: list[str], in_rj: str | None, run_minutes: int | None = None, duration_source: str | None = None) -> tuple[dict | None, dict]:
    """(section, numbers for the status engine). Without a readable report the section says so and holds what the
    final summary gives."""
    pri = rep.get("protocol_run_info") if isinstance(rep.get("protocol_run_info"), dict) else {}
    msgs = [m for m in (rep.get("user_messages") or []) if isinstance(m, dict)]
    kept: list[list] = []
    n_info = n_warn = n_err = n_disk = 0
    ids: dict[str, int] = {}
    for m in msgs:
        idn, sev = str(m.get("identifier", "")), str(m.get("severity", ""))
        if idn in _SCAN_MESSAGES:
            continue
        short = sev.replace("MESSAGE_SEVERITY_", "").lower() or "unknown"
        kept.append([iso_z(m.get("time")), short, clean_text(idn, 128), clean_text(m.get("user_message", ""), 500)])
        if sev.endswith("WARNING"):
            n_warn += 1
            ids[idn] = ids.get(idn, 0) + 1
            if is_computer_notice(idn):
                n_disk += 1
        elif sev.endswith("ERROR"):
            n_err += 1
            ids["ERROR:" + idn] = ids.get("ERROR:" + idn, 0) + 1
        elif sev.endswith("INFO"):
            n_info += 1
    ws = ari.get("writer_summary") if isinstance(ari.get("writer_summary"), dict) else {}
    write_failed = to_int(ws.get("bytes_to_write_failed"), 0)
    produced, completed = to_int(ws.get("bytes_to_write_produced"), 0), to_int(ws.get("bytes_to_write_completed"), 0)
    fallback = sum(int(v) for k, v in fs.items() if k.endswith("_in_fallback") and v.isdigit())
    file_counts = {k[: -len("_files_in_final_dest")] + "_files": int(v) for k, v in fs.items()
                   if k.endswith("_files_in_final_dest") and v.isdigit()}
    backlog = max(produced - completed, 0)
    derived = {"msg_warn": n_warn, "msg_error": n_err, "msg_warn_disk": n_disk, "fallback_files": fallback,
               "write_failed": write_failed, "write_backlog": backlog, "state": pri.get("state"), "stop_reason": ari.get("stop_reason")}
    if not rep and not fs and run_hours is None:
        return None, derived
    notes: list[dict] = []
    scalars = {"files_in_fallback": fallback, **file_counts}
    tables = {}
    summary: dict = {"run_duration_h": round(run_hours, 2) if run_hours is not None else None,
                     "run_completion_frac": round(min(run_hours * 3600 / planned_s, 1.0), 4) if planned_s and run_hours is not None else None}
    state, reason = "ok", None
    if rep:
        scalars.update(messages_info=n_info, messages_warning=n_warn, messages_error=n_err, messages_warning_computer=n_disk,
                       raw_bytes_completed=completed, raw_bytes_failed=write_failed, raw_bytes_backlog=backlog)
        tables["run_messages"] = {"columns": [{"id": "time", "type": "datetime"}, {"id": "severity", "type": "str"},
                                              {"id": "identifier", "type": "str"}, {"id": "user_message", "type": "str"}],
                                  "key": ["time"], "merge": "concat", "rows": kept[:MAX_MESSAGES], "total_rows": len(kept),
                                  "truncated": len(kept) > MAX_MESSAGES}
        summary.update(run_stop_cause=pri.get("state"), run_stop_reason=ari.get("stop_reason"),
                       run_messages_error=n_err, run_messages_warning=n_warn - n_disk)
        notes.append(note("info", "mux_scan_messages_moved",
                          f"{len(msgs) - len(kept)} of the {len(msgs)} messages are pore-scan start and result messages; "
                          "they feed the pore-scan series of the flow cell section and are not repeated here."))
    else:
        state = "error" if report_problem else "partial"
        reason = "inputs_missing"
        notes.append(note("error" if report_problem else "warning", "report_unreadable" if report_problem else "report_missing",
                          ("MinKNOW's report could not be used: " + report_problem + ".") if report_problem else
                          "This folder holds no report from MinKNOW: messages and how the run ended are not known.", input=in_rj))
    cell: dict = {"at": [], "scalars": scalars}
    if tables:
        cell["tables"] = tables
    sec = {"version": VERSIONS["run"], "tier": "run_health", "state": state, "inputs": list(inputs), "dims": [], "cells": [cell],
           "summary": summary, "notes": notes,
           "params": {"duration_minutes": run_minutes, "duration_source": duration_source,
                      "duration_rule": "the last minute of MinKNOW's throughput file (of the pore-activity file when there is none)"}}
    if state == "partial":
        sec["partial_reason"] = reason
    return sec, derived
