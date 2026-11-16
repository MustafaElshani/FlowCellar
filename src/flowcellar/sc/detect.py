"""Which folders are analyses. A folder is recognised by the files in it, never by its name.

Three layouts are known:

``spaceranger``
    a folder with ``metrics_summary.csv`` and ``web_summary.html`` and either ``spatial/`` or
    ``binned_outputs/`` (the ``outs/`` folder of Space Ranger). It is Visium HD when ``feature_slice.h5`` or
    ``binned_outputs/`` is there.
``wf_v3``
    a workflow output folder with ``wf-single-cell-report.html`` and ``params.json``; one analysis per
    sub-folder that holds ``<alias>.read_summary.tsv`` or an ``<alias>.gene_raw_feature_bc_matrix*`` folder.
``wf_v2``
    a workflow output folder with ``wf-single-cell-report.html`` and no ``params.json``; one analysis per
    sub-folder that holds ``read_tags.tsv`` or ``config_stats.json``.

Only listings and ``stat`` calls happen here: no file is opened. Standard library only.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

from flowcellar.contracts import InputRule
from flowcellar.sc import PIPELINE_SPACERANGER, PIPELINE_WF, ROLE_SPACERANGER, ROLE_WF

REPORT_NAME = "wf-single-cell-report.html"
METRICS_NAME = "metrics_summary.csv"
#: folders that hold copies or scratch of an analysis, never an analysis of their own
SKIP_NAMES = frozenset((".ipynb_checkpoints", "SPATIAL_RNA_COUNTER_CS", "work", "execution"))
_SAMPLE_SKIP = frozenset(("execution", "work", "fusions", "fastqs"))


@dataclass(frozen=True)
class Claim:
    """One analysis folder with the files that were found in it (role -> absolute path)."""
    pipeline: str                       # "spaceranger" or "wf-single-cell"
    layout: str                         # "spaceranger", "wf_v2" or "wf_v3"
    path: str                           # the outs/ folder, or the sample folder of the workflow
    hd: bool = False                    # the file layout of Visium HD
    alias: str | None = None            # wf-single-cell: the name of the sample folder
    top: str | None = None              # wf-single-cell: the workflow output folder (report, params.json, execution/)
    files: dict[str, str] = field(default_factory=dict)

    @property
    def variant(self) -> str:
        """The folder name that tells one output from a re-run of it."""
        if self.pipeline == PIPELINE_SPACERANGER:
            base = os.path.basename(self.path.rstrip("/"))
            return os.path.basename(os.path.dirname(self.path.rstrip("/"))) if base == "outs" else base
        return os.path.basename((self.top or self.path).rstrip("/"))

    def to_json(self) -> dict:
        return {"pipeline": self.pipeline, "layout": self.layout, "path": self.path, "hd": self.hd, "alias": self.alias,
                "top": self.top, "files": dict(self.files)}


def input_rules() -> list[InputRule]:
    """How the discovery walk finds analyses by itself: the two generic anchor files. A library adds patterns of
    its own with ``[[inputs]]`` entries of the roles ``analysis.spaceranger`` (the ``metrics_summary.csv`` of an
    ``outs/`` folder) and ``analysis.wf_single_cell`` (the report of a workflow output folder).

    Both are FILE names looked for in every folder the walk lists, never a folder name: a pattern that starts at
    a folder (``**/outs/...``) makes every folder of that name a place the generic walk does not enter, and runs
    or read-data files below a folder called ``outs`` would no longer be found. Whether a ``metrics_summary.csv``
    is Space Ranger's is decided from the folder it lies in (``is_spaceranger_outs``), not from that folder's name."""
    return [InputRule(role=ROLE_SPACERANGER, glob="**/" + METRICS_NAME),
            InputRule(role=ROLE_WF, glob="**/" + REPORT_NAME)]


def skipped(path: str) -> bool:
    """True when a path lies inside a copy or scratch folder (notebook checkpoints, a pipestance, work/)."""
    return any(p in SKIP_NAMES for p in os.path.normpath(path).split(os.sep))


def _listing(path: str) -> tuple[set[str], set[str]]:
    files: set[str] = set()
    dirs: set[str] = set()
    try:
        with os.scandir(path) as it:
            for e in it:
                try:
                    if e.is_symlink():
                        continue
                    (dirs if e.is_dir(follow_symlinks=False) else files).add(e.name)
                except OSError:
                    continue
    except OSError:
        pass
    return files, dirs


def _first(path: str, *names: str) -> str | None:
    for n in names:
        p = os.path.join(path, n)
        if os.path.isfile(p):
            return p
    return None


# ---------------------------------------------------------------------------------------------------------
# Space Ranger
# ---------------------------------------------------------------------------------------------------------

def is_spaceranger_outs(path: str) -> bool:
    files, dirs = _listing(path)
    return METRICS_NAME in files and "web_summary.html" in files and bool({"spatial", "binned_outputs"} & dirs)


def _pipestance_metrics(parent: str) -> str | None:
    """``metrics_summary_json.json`` of the pipestance beside ``outs/``, when it was kept. A bounded look: the
    file sits under ``SPATIAL_RNA_COUNTER_CS/**/SPATIAL_REPORTER/fork0/join-*/files/`` (one level deeper in
    Space Ranger 4 than in 3)."""
    top = os.path.join(parent, "SPATIAL_RNA_COUNTER_CS")
    if not os.path.isdir(top):
        return None
    stack = [(top, 0)]
    listed = 0
    while stack and listed < 400:
        d, depth = stack.pop()
        files, dirs = _listing(d)
        listed += 1
        if os.path.basename(d) == "SPATIAL_REPORTER":
            fork = os.path.join(d, "fork0")
            _, joins = _listing(fork)
            for j in sorted(joins):
                if j.startswith("join-"):
                    hit = os.path.join(fork, j, "files", "metrics_summary_json.json")
                    if os.path.isfile(hit):
                        return hit
            continue
        if depth >= 4:
            continue
        # the stages are named in capitals; forks and chunks are not on the way
        stack.extend((os.path.join(d, n), depth + 1) for n in sorted(dirs, reverse=True) if n.isupper() or n.startswith("_"))
    return None


def spaceranger_claim(outs: str) -> Claim | None:
    outs = os.path.abspath(outs)
    if not is_spaceranger_outs(outs):
        return None
    files, dirs = _listing(outs)
    f: dict[str, str] = {"metrics_summary": os.path.join(outs, METRICS_NAME), "web_summary": os.path.join(outs, "web_summary.html")}
    for role, name in (("molecule_info", "molecule_info.h5"), ("feature_slice", "feature_slice.h5"),
                       ("raw_matrix", "raw_feature_bc_matrix.h5"), ("filtered_matrix", "filtered_feature_bc_matrix.h5"),
                       ("bam", "possorted_genome_bam.bam"), ("bam_index", "possorted_genome_bam.bam.bai")):
        if name in files:
            f[role] = os.path.join(outs, name)
    hd = "feature_slice" in f or "binned_outputs" in dirs
    if "binned_outputs" in dirs:
        for um in ("002", "008", "016"):
            sq = os.path.join(outs, "binned_outputs", f"square_{um}um")
            for kind in ("raw", "filtered"):
                p = os.path.join(sq, f"{kind}_feature_bc_matrix.h5")
                if os.path.isfile(p):
                    f[f"{kind}_matrix_{um}um"] = p
            p = os.path.join(sq, "spatial", "scalefactors_json.json")
            if os.path.isfile(p):
                f[f"scalefactors_{um}um"] = p
    if "spatial" in dirs:
        sp = os.path.join(outs, "spatial")
        for role, names in (("tissue_positions", ("tissue_positions.csv", "tissue_positions_list.csv")),
                            ("scalefactors", ("scalefactors_json.json",))):
            p = _first(sp, *names)
            if p:
                f[role] = p
    parent = os.path.dirname(outs)
    for role, name in (("versions", "_versions"), ("cmdline", "_cmdline"), ("invocation", "_invocation")):
        p = os.path.join(parent, name)
        if os.path.isfile(p):
            f[role] = p
    pm = _pipestance_metrics(parent)
    if pm:
        f["pipestance_metrics"] = pm
    return Claim(PIPELINE_SPACERANGER, "spaceranger", outs, hd=hd, files=f)


# ---------------------------------------------------------------------------------------------------------
# wf-single-cell
# ---------------------------------------------------------------------------------------------------------

def _matrix_dirs(sample_dir: str, alias: str, dirs: set[str]) -> dict[str, str]:
    """The raw count matrices of one sample folder by role: ``gene_matrix``, ``gene_matrix_8um``, ...
    (the processed matrices are normalised and are never read)."""
    out: dict[str, str] = {}
    for what in ("gene", "transcript"):
        for suffix, role in (("", f"{what}_matrix"), ("_8um", f"{what}_matrix_8um"), ("_2um", f"{what}_matrix_2um")):
            for name in (f"{alias}.{what}_raw_feature_bc_matrix{suffix}", f"{what}_raw_feature_bc_matrix{suffix}"):
                if name in dirs and os.path.isfile(os.path.join(sample_dir, name, "matrix.mtx.gz")):
                    out[role] = os.path.join(sample_dir, name)
                    break
    return out


def _wf_sample_claim(top: str, sample_dir: str, top_files: set[str]) -> Claim | None:
    alias = os.path.basename(sample_dir.rstrip("/"))
    files, dirs = _listing(sample_dir)
    v3_summary = f"{alias}.read_summary.tsv"
    mats = _matrix_dirs(sample_dir, alias, dirs)
    is_v3 = v3_summary in files or any(d.startswith(alias + ".") and "_feature_bc_matrix" in d for d in dirs)
    is_v2 = "read_tags.tsv" in files or "config_stats.json" in files
    if not (is_v3 or is_v2 or mats):
        return None
    layout = "wf_v3" if is_v3 or ("params.json" in top_files and not is_v2) else "wf_v2"
    f: dict[str, str] = {"report": os.path.join(top, REPORT_NAME)} if REPORT_NAME in top_files else {}
    if "params.json" in top_files:
        f["params"] = os.path.join(top, "params.json")
    trace = os.path.join(top, "execution", "trace.txt")
    if os.path.isfile(trace):
        f["trace"] = trace
    for role, names in (
            ("config_stats", (f"{alias}.config_stats.json", "config_stats.json")),
            ("read_summary", (v3_summary, "read_tags.tsv")),
            ("mito", (f"{alias}.gene_expression_mito_per_cell.tsv", "gene.expression.mito-per-cell.tsv")),
            ("whitelist", (f"{alias}.whitelist.tsv", "whitelist.tsv")),
            ("bam", (f"{alias}.tagged.bam", "tagged.bam")),
            ("bc_assignment", (f"{alias}.bc_assignment_summary.tsv",))):
        hit = next((n for n in names if n in files), None)
        if hit:
            f[role] = os.path.join(sample_dir, hit)
    if "bam" in f and os.path.isfile(f["bam"] + ".bai"):
        f["bam_index"] = f["bam"] + ".bai"
    f.update(mats)
    hd = "gene_matrix_8um" in f or "gene_matrix_2um" in f
    return Claim(PIPELINE_WF, layout, os.path.abspath(sample_dir), hd=hd, alias=alias, top=os.path.abspath(top), files=f)


def wf_claims(top: str) -> list[Claim]:
    """The analyses of one workflow output folder (the folder that holds the report): one per sample folder."""
    top = os.path.abspath(top)
    top_files, top_dirs = _listing(top)
    if REPORT_NAME not in top_files:
        return []
    out = []
    for name in sorted(top_dirs):
        if name.startswith(".") or name in _SAMPLE_SKIP or name in SKIP_NAMES:
            continue
        c = _wf_sample_claim(top, os.path.join(top, name), top_files)
        if c is not None:
            out.append(c)
    return out


# ---------------------------------------------------------------------------------------------------------
# any folder
# ---------------------------------------------------------------------------------------------------------

def claims_at(path: str) -> list[Claim]:
    """The analyses a folder stands for: an ``outs/`` folder, a folder that holds one, a workflow output folder, or
    one sample folder of a workflow. An empty list when it is none of these."""
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        return []
    c = spaceranger_claim(path)
    if c:
        return [c]
    if os.path.isdir(os.path.join(path, "outs")):
        c = spaceranger_claim(os.path.join(path, "outs"))
        if c:
            return [c]
    found = wf_claims(path)
    if found:
        return found
    top = os.path.dirname(path)
    top_files, _ = _listing(top)
    if REPORT_NAME in top_files or "params.json" in top_files:
        c = _wf_sample_claim(top, path, top_files)
        if c:
            return [c]
    return []


def claims_from_anchor(role: str, anchor_path: str, relpath: str | None = None) -> list[Claim]:
    """The analyses behind one anchor file the discovery walk reported. ``relpath`` is its path below the scanned
    folder: copies and scratch folders are told by the names below that folder only, so that a library whose data
    lies under ``/work/...`` or in a folder named ``work`` keeps its analyses."""
    if skipped(relpath if relpath else anchor_path):
        return []
    folder = os.path.dirname(anchor_path)
    if role == ROLE_SPACERANGER:
        c = spaceranger_claim(folder)
        return [c] if c else []
    if role == ROLE_WF:
        return wf_claims(folder)
    return []


def spaceranger_outs_below(top: str, *, max_depth: int = 3, max_listings: int = 400) -> list[str]:
    """``outs/`` folders of Space Ranger inside a workflow output folder (the older layout keeps them there)."""
    out, stack, listed = [], [(os.path.abspath(top), 0)], 0
    while stack and listed < max_listings:
        d, depth = stack.pop()
        _, dirs = _listing(d)
        listed += 1
        for n in sorted(dirs):
            if n.startswith(".") or n in SKIP_NAMES or "_feature_bc_matrix" in n:
                continue
            p = os.path.join(d, n)
            if n == "outs" and is_spaceranger_outs(p):
                out.append(p)
            elif depth < max_depth:
                stack.append((p, depth + 1))
    return sorted(out)
