"""Read-data files for ``flowcellar run RUN_DIR [PATH...] --input ROLE=PATH --sample-sheet CSV`` (no library).

Three sources of candidates:

* the run folder itself and every folder given: the built-in patterns (``bam_pass/``, ``bam_fail/``, the names
  dorado demux writes), resolved by directed descent;
* every file and folder given: recognition by content (a BAM is unaligned or aligned by its header; pychopper
  statistics, quantification metadata, read summaries and sample sheets by their first bytes);
* ``--input ROLE=PATH``: the role is stated, the path is a file or a folder of such files.

A file named on the command line that carries no identifier is linked to the run as stated by the person
(method ``manual``); identifiers inside a file still win.
"""
from __future__ import annotations

import os
import re
import stat as _stat
from typing import Mapping, Sequence

from flowcellar.contracts import Candidate, RunIdentity
from flowcellar.model import runfolder
from flowcellar.readdata import patterns
from flowcellar.readdata.sample_sheet import sniff

COMMAND_LINE = "named on the command line"
MAX_DEPTH = 6
MAX_FILES = 5000
_HINTS = (".bam", ".csv", ".tsv", ".json", ".txt")
_STRIP = (("pychopper.stats", re.compile(r"(.+?)_pychopper_stats\.tsv$")), ("quant.oarfish", re.compile(r"(.+?)\.meta_info\.json$")),
          ("polya.read_table", re.compile(r"(.+?)_readid_annotated\.csv$")),
          # the endings pipelines give an alignment, in the spellings met so far; without them the genome and the
          # transcriptome alignment of one sample would be taken for two samples
          ("align.genome_bam", re.compile(r"(.+?)(?:[._](?:genome|gnome|genomic))?(?:[._]sorted)?\.bam$")),
          ("align.transcriptome_bam", re.compile(r"(.+?)(?:[._](?:transcriptome|txome|transcripts))?(?:[._]sorted)?\.bam$")))


def guess_sample(role: str, path: str, top: str | None = None) -> str | None:
    """The sample a file belongs to, from the naming habits of the tools: the folder under ``bam_pass``, the name
    before ``_pychopper_stats.tsv``, the salmon output folder, the kallisto output folder."""
    name = os.path.basename(path)
    parent = os.path.basename(os.path.dirname(path))
    if role == "dorado.bam":
        grand = os.path.basename(os.path.dirname(os.path.dirname(path)))
        if grand in ("bam_pass", "bam_fail") or (top and os.path.dirname(os.path.dirname(path)) == os.path.normpath(top)):
            return parent
        return None
    if role == "quant.salmon":
        d = os.path.basename(os.path.dirname(os.path.dirname(path))) if parent == "aux_info" else parent
        return re.sub(r"_salmon$", "", d) or None
    if role == "quant.kallisto":
        return parent or None
    for r, rx in _STRIP:
        if r == role:
            m = rx.match(name)
            return m.group(1) if m else None
    return None


def _candidate(role: str, path: str, *, sample: str | None, rule: str, exp: str | None, top: str | None = None) -> Candidate | None:
    try:
        st = os.lstat(path)
    except OSError:
        return None
    if not _stat.S_ISREG(st.st_mode):
        return None
    return Candidate(role=role, path=path, size=st.st_size, mtime_ns=st.st_mtime_ns, root_id=None,
                     relpath=os.path.relpath(path, top) if top else os.path.basename(path), sample=sample, rule=rule, experiment_dir=exp)


def _files_below(top: str) -> list[str]:
    """Plain files below a folder that may be read-data files, by name (pruned, at most six levels, links skipped)."""
    out: list[str] = []
    stack = [(top, 0)]
    while stack and len(out) < MAX_FILES:
        d, depth = stack.pop()
        try:
            entries = sorted(os.scandir(d), key=lambda e: e.name)
        except OSError:
            continue
        for e in entries:
            if e.name.startswith(".") or e.is_symlink():
                continue
            if e.is_dir(follow_symlinks=False):
                keep = e.name in ("bam_pass", "bam_fail") or not runfolder.is_pruned(e.name)
                if depth < MAX_DEPTH and keep:
                    stack.append((e.path, depth + 1))
            elif e.name.lower().endswith(_HINTS):
                out.append(e.path)
    return sorted(out)


def candidates_for_run(run: RunIdentity, paths: Sequence[str] = (), explicit: Mapping[str, Sequence[str]] | None = None,
                       sample_sheet: str | None = None) -> list[Candidate]:
    from flowcellar.readdata.api import default_input_rules
    rules = [r for r in default_input_rules() if r.glob]
    exp = os.path.dirname(os.path.normpath(run.run_dir)) if run.run_dir else None
    found: dict[str, Candidate] = {}

    def add(c: Candidate | None) -> None:
        if c is not None and c.path not in found:
            found[c.path] = c

    tier1 = {os.path.normpath(p) for p in (run.files or {}).values()}
    # ---- stated on the command line: the role is given
    for role, plist in (explicit or {}).items():
        for p in plist:
            p = os.path.abspath(p)
            if os.path.isdir(p):
                for f in _files_below(p):
                    if role.endswith("_bam") or role == "dorado.bam":
                        ok = f.lower().endswith(".bam")
                    else:
                        ok = sniff(f) == role
                    if ok:
                        add(_candidate(role, f, sample=guess_sample(role, f, p), rule=f"{COMMAND_LINE} as {role}", exp=exp, top=p))
            else:
                add(_candidate(role, p, sample=guess_sample(role, p), rule=f"{COMMAND_LINE} as {role}", exp=exp))
    if sample_sheet:
        add(_candidate("basecaller.sample_sheet", os.path.abspath(sample_sheet), sample=None, rule=f"{COMMAND_LINE} as the sample sheet", exp=exp))
    # ---- the run folder: the built-in patterns only (its own small files belong to run health)
    if run.run_dir and os.path.isdir(run.run_dir):
        got, _ = patterns.find(run.run_dir, rules, experiment_dir=exp)
        for c in got:
            add(c)
        # with live basecalling on, MinKNOW's own summary holds length, quality and barcode of every read
        ss = (run.files or {}).get("minknow.sequencing_summary")
        if ss and os.path.isfile(ss) and sniff(ss) == "dorado.summary":
            add(_candidate("dorado.summary", ss, sample=None, rule="the run's own sequencing summary, written with basecalling on", exp=exp, top=run.run_dir))
    # ---- the paths given: the built-in patterns, then recognition by content
    for p in paths:
        p = os.path.abspath(p)
        if os.path.isdir(p):
            got, _ = patterns.find(p, rules, experiment_dir=exp)
            for c in got:
                add(Candidate(**{**c.__dict__, "rule": f"{COMMAND_LINE}; {c.rule}"}))
            for f in _files_below(p):
                if f in found or os.path.normpath(f) in tier1:
                    continue
                role = sniff(f)
                if role and role != "flowcellar.partial":
                    add(_candidate(role, f, sample=guess_sample(role, f, p), rule=f"{COMMAND_LINE}; told by its content", exp=exp, top=p))
        else:
            role = sniff(p)
            if role and role != "flowcellar.partial":
                add(_candidate(role, p, sample=guess_sample(role, p), rule=f"{COMMAND_LINE}; told by its content", exp=exp))
    return [found[k] for k in sorted(found)]
