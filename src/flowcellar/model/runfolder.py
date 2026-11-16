"""Recognising a MinKNOW run folder by the files in it, never by its name.

One rule, used by the discovery walk (LIBRARY) and by ``flowcellar run`` (ENGINE).
"""
from __future__ import annotations

import fnmatch
import re
from typing import Iterable

from flowcellar.model.runkey import parse_minknow_name, run_key_from_folder_name

# role -> file name pattern, relative to the run folder
TIER1_ROLES: dict[str, str] = {
    "minknow.report_json": "report_*.json",
    "minknow.final_summary": "final_summary_*.txt",
    "minknow.sample_sheet": "sample_sheet_*.csv",
    "minknow.throughput": "throughput_*.csv",
    "minknow.pore_activity": "pore_activity_*.csv",
    "minknow.sequencing_summary": "sequencing_summary_*.txt",
    "minknow.barcode_alignment": "barcode_alignment_*.tsv",
    "minknow.output_hash": "output_hash_*.csv",
    "minknow.pore_scan": "other_reports/pore_scan_data_*.csv",
}
#: any one of these makes a folder a run folder
IDENTITY_PATTERNS = ("final_summary_*.txt", "report_*.json", "sequencing_summary_*.txt", "sample_sheet_*.csv",
                     "throughput_*.csv", "pore_activity_*.csv")
RAW_PATTERNS = ("pod5*", "fast5*", "other_reports")
BASECALLER_OUTPUT = ("bam_pass", "bam_fail", "fastq_pass", "fastq_fail", "fastq_skip", "bam_skip")
#: folders the generic walk never lists; the configured ``prune`` list is added to these
DEFAULT_PRUNE = ("pod5*", "fast5*", "fastq_pass", "fastq_fail", "fastq_skip", "bam_pass", "bam_fail", "bam_skip", "work",
                 "slurm_logs", "logs", "node_modules", "__pycache__", "env", "envs", "conda", "venv", "renv",
                 "site-packages", ".*")
_LOOKALIKE = re.compile(r"_[A-Z]{2,3}[0-9]{4,6}_[0-9a-f]{8}$")


def _any(names: Iterable[str], patterns: Iterable[str]) -> bool:
    pats = tuple(patterns)
    return any(fnmatch.fnmatchcase(n, p) for n in names for p in pats)


def classify_listing(dirname: str, names: Iterable[str], dirs: Iterable[str] | None = None) -> str | None:
    """What a directory is, from its own listing (file and folder names, not paths).

    ``"run"``: holds at least one identity file (and not the final summaries of several runs). ``"run_raw_only"``: none of those, but raw data and a name or a
    raw file name that gives ``{flow cell}_{8 hex}``. ``"lookalike"``: named like a run but holds only basecaller
    output. ``None``: anything else.

    ``dirs`` names the entries that are folders. Give it whenever the caller knows: a raw-data marker must then
    be a FOLDER (``pod5*``, ``fast5*``, ``other_reports``) or a raw file (``*.pod5``, ``*.fast5``), not a script
    that happens to be called ``fast5_to_pod5.sh``, and the run key must come from the folder's own name or from
    a raw file in it, never from the name of a sub-folder. Without it the folder ABOVE a run folder could be
    taken for a raw-only run, and the run inside it would never be reached.
    """
    names = [n for n in names if not n.startswith(".")]
    if _any(names, IDENTITY_PATTERNS):
        # MinKNOW writes one run per folder. A folder with the final summaries of several runs is a collection of
        # copied reports (somebody's working folder), not a run: taking it for one would give that run a second
        # place and a project that is none.
        finals = {p.run_key for p in (parse_minknow_name(n) for n in names if fnmatch.fnmatchcase(n, "final_summary_*.txt")) if p}
        if len(finals) > 1:
            return None
        return "run"
    if _any(names, RAW_PATTERNS):
        if dirs is None:
            if run_key_from_folder_name(dirname) or any(parse_minknow_name(n) for n in names):
                return "run_raw_only"
            return None
        folders = set(dirs)
        raw_files = [n for n in names if n not in folders and n.lower().endswith((".pod5", ".fast5"))]
        raw_dirs = [n for n in names if n in folders and _any([n], RAW_PATTERNS)]
        if not (raw_dirs or raw_files):
            return None
        if run_key_from_folder_name(dirname) or any(parse_minknow_name(n) for n in raw_files):
            return "run_raw_only"
        return None
    if names and _LOOKALIKE.search(dirname) and all(n in BASECALLER_OUTPUT for n in names):
        return "lookalike"
    return None


def is_pruned(name: str, extra: Iterable[str] = ()) -> bool:
    return _any([name], tuple(DEFAULT_PRUNE) + tuple(extra))
