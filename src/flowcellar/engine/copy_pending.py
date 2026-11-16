"""Sentences of the status engine that the catalogue does not hold yet.

A sentence that is missing in the catalogue is kept in this
table, so that every one of them can be found and moved into ``catalogue/metrics_nanopore.yaml`` in one go.
Placeholders are filled with ``str.format``. Sentences of single checks come from the catalogue's display copy
(``display.message`` and ``limit_text``) and are not repeated here.
"""
from __future__ import annotations

HEADLINE = {
    "raw_only": "Only raw data was found for this run: no report, so the flow cell cannot be judged.",
    "running": "This run is still going: MinKNOW writes its report when the run ends.",
    "no_final_summary": "MinKNOW wrote no final summary for this run: it was interrupted, or its files are incomplete.",
    "inputs_missing": "The files MinKNOW writes at the end of a run are missing or cannot be read, so this run cannot be judged.",
    "stopped_early": "Stopped after {duration}. Too short to judge.",
    "output_not_judged": "Stopped after {duration}. No check failed, but the run is too short to judge its output.",
    "ran_well": "The flow cell ran well{so_far} ({facts}).",
    "ran_well_plain": "The flow cell ran well{so_far}.",
    "failed_on": "Failed on {list}.",
    "also_needs_a_look": " Also needs a look: {list}.",
    "needs_a_look": "Needs a look: {list}.",
    "many": "{n} checks need a look, first {first}.",
}
HEADLINE_FACT = {
    "reads": "{value} reads by {horizon} h",
    "bases": "{value} estimated by {horizon} h",
    "occupancy": "{value} of pores busy in the first hour",
}
MODIFIER = {
    "reused_flow_cell": "Re-used flow cell (use {index}): judged on output per starting pore, occupancy, speed and read quality only.",
    "flow_cell_use_possibly_reused": "Possibly re-used: the first scan is far below this lab's fresh flow cells.",
    "flow_cell_use_unknown": "The use of this flow cell is not declared and the first scan is below the warranty figure, "
                             "so the output check can warn but not fail.",
    "no_history": "No lab history for this combination yet: the output was judged on fixed levels.",
    "no_history_no_level": "No lab history for this combination yet, and no fixed level for this assay: the output was not judged.",
    "quick_detail": "Made from MinKNOW's report files. Speed, why reads ended, read length and the channel map appear when the read count finishes.",
    "reads_partly_counted": "The read count is not complete: speed, why reads ended and read length rest on {share} of the reads.",
    "reads_not_counted": "The reads were not counted: the sequencing summary could not be read, so speed, why reads ended and "
                         "read length are not known.",
    "custom_profile": "Limits changed in this library's configuration were used.",
    "read_data_so_far_one": "Read data so far: 1 file of this run has not been read yet; {assessed} of {total} read-data checks "
                            "could be assessed.",
    "read_data_so_far": "Read data so far: {files} files of this run have not been read yet; {assessed} of {total} read-data "
                        "checks could be assessed.",
}
RULE = {
    "stop_planned": "The run ended as planned after {duration}.",
    "stop_by_person": "The run was stopped by a person after {duration}.",
    "stop_other": "The run ended after {duration}.",
    "messages_none": "none that need a look",
    "messages_error": "{n} marked as an error",
    "messages_warning": "{n} warnings",
    "messages_warning_one": "1 warning",
    "messages_computer": " ({n} notices about the computer's disk are in the run log)",
    "files_ok": "Raw data written: every file reached the output folder.",
    "files_value": "{n} files",
    "files_bytes": "a part of the raw data",
    "temperature_ok": "Temperature: inside the target band for the whole run.",
    "flow_cell_check": "Pores at the flow cell check before loading: {value}. ONT's warranty minimum is {limit}.",
    "pores_ok": "Working pores at the first scan after loading: {value}.",
    "speed_value": "{value} bases per second",
    "count_of": "{count} of {total} ({share})",
    "bases_value": "{value} bases",
    "fixed_level_also": "The fixed level at {horizon} h is {level}.",
    "from_part": "From {share} of the reads.",
    "output_no_level": "No lab history for this combination yet, and no fixed level for this assay.",
    "ont_guide": "ONT's guide is {limit} or more.",
}
REASON = {
    "report_unreadable": "MinKNOW's report is missing or could not be read",
    "summary_missing": "the sequencing summary is not in the run folder",
    "summary_unreadable": "the sequencing summary could not be read",
    "channels_need_every_read": "channels can only be counted when every read has been read",
    "no_target_speed": "the report gives no target speed for this run",
}
BASIS = {
    "run_stop_cause": "pass when the protocol completed or a person stopped it; fail on a device, disk, temperature or internal error",
    "run_messages": "fail on any MinKNOW message marked as an error; warn on a warning that is not a notice about the computer (disk space)",
    "run_output_files": "warn when any file was written to a fall-back location or raw data could not be written",
    "fc_pores_first_scan": "built-in profile. ONT's warranty minimum belongs to the flow cell check before loading; here it is a warning level for the first scan after loading",
    "fc_pores_flow_cell_check": "ONT's warranty minimum, compared with the declared result of the flow cell check before loading",
    "fc_channels_with_pore_first_scan": "built-in profile: share of channels with a working pore at the first scan",
    "fc_occupancy_first_hour": "built-in profile: ONT's guide is 70 % in the first hour (minutes 30 to 90); the fail level is this tool's own",
    "fc_active_channels": "built-in profile: share of channels that produced one read or more",
    "fc_temperature_max_dev_c": "built-in profile: distance from the run's own target band",
    "dyn_speed_est_median": "the run's own target speed band, widened by 5 % for the estimate's own bias; fail more than 15 % below it",
    "default": "built-in profile",
    "configured": "this library's configuration",
    "history": "this lab's history as of the run's start: {n} first-use runs of the same device type, flow cell product and kit that "
               "started before this run and reached {horizon} h; a run is never compared with itself or with a run started on the same day",
    "history_fixed_also": "; the fixed level at {horizon} h ({level}) also applies until the history is trusted",
    "fixed_only": "fixed level for {assay} on {family} at {horizon} h: fewer than 5 comparable first-use runs started before this one",
    "no_level": "no fixed level exists for this assay and device, and fewer than 5 comparable first-use runs started before this one",
    "reused": "re-used flow cell: output per starting pore at {horizon} h against {n} first-use runs of the same device type, "
              "flow cell product and kit that started before this run",
}
DETAIL = {
    "pending": "{count} for the read count.",
    "read_data_absent": "Read data: not available yet.",
    "read_data_samples": "Read data for {with_data} of {expected} samples.",
    "read_data_present": "Read data is available.",
    "read_data_not_assessed": "None of the read-data checks could be assessed.",
    "summary_missing": "The sequencing summary is not in the run folder, so speed, why reads ended and read length are not known.",
}
COUNT_WORDS = {1: "One check waits", 2: "Two checks wait", 3: "Three checks wait", 4: "Four checks wait", 5: "Five checks wait",
               6: "Six checks wait"}
#: words the interface never shows; no stored sentence may hold one (tests/engine checks every sentence)
BANNED_WORDS = ("tier", "rule", "threshold", "cut-off", "baseline", "cohort", "z-score", "bundle", "card", "partial", "sidecar",
                "aborted", "duty time", "funnel", "sankey", "unclassified", "merged", "aggregated", "stale", "override", "n/a",
                "skipped", "bad", "poor", "critical", "abnormal", "read-level", "!")
