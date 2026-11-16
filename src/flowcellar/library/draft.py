"""``flowcellar init``: draft a full ``flowcellar.toml``, optionally from what a walk of the roots sees.

``init`` is optional (``scan`` creates a library by itself). It exists for the person who wants to shape the
configuration first: every section is written out with its default, and ``--from-scan`` adds, as comments to
check and switch on, the read-data file patterns the folders suggest. Listings only: no file is opened.
"""
from __future__ import annotations

import os
from collections import Counter, defaultdict

from flowcellar.config import loader
from flowcellar.discovery import discover
from flowcellar.model import runfolder

WIDE = 24               # a folder with more sub-folders than this is sampled, not walked
SAMPLED = 3
MAX_LISTINGS = 4000     # per root, for the look at read-data files
QUANT = {"meta_info.json": ("quant.salmon", "aux_info"), "run_info.json": ("quant.kallisto", None)}


def _look(root: str) -> tuple[list[tuple[str, str]], bool]:
    """(folder relative to the root, file name) of files that look like read data, from a bounded, sampled walk."""
    hits: list[tuple[str, str]] = []
    listed = 0
    cut = False
    stack = [(root, 0)]
    while stack:
        d, depth = stack.pop()
        if listed >= MAX_LISTINGS:
            cut = True
            break
        try:
            with os.scandir(d) as it:
                entries = sorted((e.name, e.is_dir(follow_symlinks=False), e.is_file(follow_symlinks=False)) for e in it)
        except OSError:
            continue
        listed += 1
        names = [n for n, _, _ in entries if not n.startswith(".")]
        if runfolder.classify_listing(os.path.basename(d), names, [n for n, is_dir, _ in entries if is_dir]) == "run":
            continue
        rel = os.path.relpath(d, root)
        for n, _, is_file in entries:
            if is_file and not n.startswith(".") and (n.endswith(".bam") or n in QUANT or n.endswith(".meta_info.json")):
                hits.append((rel, n))
        if depth >= 9:
            continue
        subs = [n for n, is_dir, _ in entries if is_dir and not n.startswith(".")
                and not runfolder.is_pruned(n, ("bam_pass", "bam_fail"))]
        subs += [n for n, is_dir, _ in entries if is_dir and n in ("bam_pass",)]
        if len(subs) > WIDE and depth >= 1:
            subs, cut = subs[:SAMPLED], True
        stack.extend((os.path.join(d, n), depth + 1) for n in reversed(subs))
    return hits, cut


def _tail_groups(names: list[str]) -> dict[str, int]:
    """Group file names by what they end in: ``a_x_sorted.bam``, ``b_x_sorted.bam`` -> ``{sample}_x_sorted.bam``."""
    stems = [n.rsplit(".", 1) for n in names]
    for k in (3, 2, 1):
        groups: Counter = Counter()
        for stem, ext in stems:
            toks = stem.split("_")
            if len(toks) <= k:
                groups = Counter()
                break
            groups["_" + "_".join(toks[-k:]) + "." + ext] += 1
        if groups and min(groups.values()) >= 2:
            return {"{sample}" + tail: n for tail, n in groups.items()}
    ext = Counter(e for _, e in stems)
    return {"{sample}." + e: n for e, n in ext.items()}


def _common_tail(names: list[str]) -> str:
    """What all names end in, from a ``_`` or ``.`` on: ``a_salmon``, ``b_salmon`` -> ``_salmon``."""
    if len(names) < 2:
        return ""
    rev = [n[::-1] for n in names]
    n_same = 0
    for chars in zip(*rev):                              # character by character from the end, while all agree
        if len(set(chars)) != 1:
            break
        n_same += 1
    tail = rev[0][:n_same][::-1]
    cut = min((i for i in (tail.find("_"), tail.find(".")) if i >= 0), default=-1)
    return tail[cut:] if cut >= 0 and len(tail) < min(len(n) for n in names) else ""


def propose_inputs(roots: list[str]) -> tuple[list[dict], bool]:
    """Read-data patterns the folders suggest: ``[{"role", "glob", "files", "folders", "check"}]``, most files first."""
    by_dir: dict[str, list[str]] = defaultdict(list)
    sampled = False
    for root in roots:
        hits, cut = _look(root)
        sampled = sampled or cut
        for rel, name in hits:
            by_dir[rel].append(name)
    seen_in: Counter = Counter()                       # folder name -> in how many of the hit folders it occurs
    for d in by_dir:
        for c in set(d.split("/")):
            seen_in[c] += 1
    found: dict[tuple[str, str], dict] = {}
    quant: dict[tuple[str, str, str], list[str]] = defaultdict(list)

    def add(role: str, glob: str, n: int, check: str | None) -> None:
        f = found.setdefault((role, glob), {"role": role, "glob": glob, "files": 0, "folders": 0, "check": check})
        f["files"] += n
        f["folders"] += 1

    def generic(comps: list[str]) -> str:
        """The last folder, with the one above it when that name is not particular to one place."""
        keep = comps[-2:] if len(comps) >= 2 and seen_in[comps[-2]] >= 2 else comps[-1:]
        return "/".join(keep)

    for rel, names in by_dir.items():
        comps = [c for c in rel.split("/") if c and c != "."]
        if not comps:
            continue
        bams = [n for n in names if n.endswith(".bam")]
        if "bam_pass" in comps and bams:
            inner = "/".join(["{sample}"] * (len(comps) - comps.index("bam_pass") - 1))
            add("dorado.bam", "**/bam_pass/" + (inner + "/" if inner else "") + "*.bam", len(bams), None)
        elif bams:
            for pat, n in _tail_groups(bams).items():
                low = pat.lower()
                role = "align.transcriptome_bam" if ("txome" in low or "transcript" in low) else "align.genome_bam"
                add(role, "**/" + generic(comps) + "/" + pat, n,
                    "check the role: align.genome_bam, align.transcriptome_bam or dorado.bam")
        for n in names:
            if n in QUANT and (QUANT[n][1] is None or comps[-1] == QUANT[n][1]):
                up = comps[:-1] if QUANT[n][1] else comps              # ends in the folder named after the sample
                if up:
                    quant[(QUANT[n][0], generic(up[:-1]) if len(up) > 1 else "", ("/" + QUANT[n][1] if QUANT[n][1] else "") + "/" + n)
                          ].append(up[-1])
            elif n.endswith(".meta_info.json") and n != "meta_info.json":
                add("quant.oarfish", "**/" + generic(comps) + "/{sample}.meta_info.json", 1, None)
    for (role, where, tail), sample_dirs in quant.items():
        glob = "**/" + (where + "/" if where else "") + "{sample}" + _common_tail(sample_dirs) + tail
        f = found.setdefault((role, glob), {"role": role, "glob": glob, "files": 0, "folders": 0,
                                            "check": "check where {sample} sits in the path"})
        f["files"] += len(sample_dirs)
        f["folders"] += len(sample_dirs)
    return sorted(found.values(), key=lambda f: (-f["files"], f["glob"])), sampled


def _example(root: dict, found) -> str | None:
    """One run folder of this root with its path components numbered, so that project_depth can be chosen."""
    rels = sorted((f.relpath for rid, f in found.runs if rid == root["id"] and f.relpath), key=lambda r: (-r.count("/"), r))
    if not rels:
        return None
    comps = rels[0].split("/")
    return " / ".join([f"{i}={c}" for i, c in enumerate(comps[:-1], 1)] + ["<run folder>"])


def draft_config(name: str, root_paths: list[str], *, look: bool = True) -> tuple[str, dict]:
    """The text of a full flowcellar.toml and what the look at the roots found (``{"runs", "raw_only", ...}``)."""
    paths = [loader.normalise_root_path(p) for p in root_paths]
    ids = loader.root_ids_for(paths, set())
    roots = [{"id": i, "path": p, "live": False, "project_depth": 1, "project_pattern": None} for i, p in zip(ids, paths)]
    facts: dict = {"roots": len(roots), "runs": 0, "raw_only": 0, "proposals": 0, "sampled": False}
    found = None
    proposals: list[dict] = []
    if roots and look:
        cfg = loader.merge({}).config
        cfg["roots"] = roots
        found = discover(cfg)
        c = found.counts()
        facts.update(runs=c["runs"], raw_only=c["raw_only"], lookalikes=c["lookalikes"], listings=c["listings"],
                     offline=[rid for rid, o in found.roots.items() if o.state != "reachable"])
        proposals, facts["sampled"] = propose_inputs([r["path"] for r in roots if r["id"] not in facts["offline"]])
        facts["proposals"] = len(proposals)
    d = loader.defaults()
    out = [loader._HEAD.format(name=loader.toml_str(name))]
    if not roots:
        out.append('\n# [[roots]]\n# id = "minknow"\n# path = "/data/minknow"\n# live = true\n# project_depth = 1\n')
    for r in roots:
        n = sum(1 for rid, _ in found.runs if rid == r["id"]) if found else None
        out.append("\n")
        if n is not None:
            out.append(f"# {n} run folder{'s' if n != 1 else ''} found here.\n")
            ex = _example(r, found)
            if ex:
                out.append(f"# One of them: {ex}\n# Set project_depth to the number in front of the folder that names the project.\n")
        out.append(loader.root_block(r))
    dv = d["discovery"]
    out.append(f"""
# ---- How folders are walked -----------------------------------------------------------------------------
[discovery]
min_age_seconds = {dv['min_age_seconds']}             # a file changed more recently than this is left for the next scan
stale_after_hours = {dv['stale_after_hours']}             # a run without a final summary counts as running while its files still change
prune = []                        # folder names (shell patterns) never entered, added to the built-in list
subtree_listing_budget = {dv['subtree_listing_budget']}    # folders listed under one first-level folder before it is given up, and named
root_check_seconds = {dv['root_check_seconds']}           # a root that does not answer in this time is skipped; nothing under it changes
stall_timeout_seconds = {dv['stall_timeout_seconds']}

# ---- Read data ------------------------------------------------------------------------------------------
# FlowCellar finds by itself what MinKNOW and Dorado write inside a run folder. Files that a pipeline writes
# elsewhere need one [[inputs]] rule each: a role and a pattern. "{{sample}}" captures the sample name,
# "{{variant:a|b}}" one of several versions of a file, "**" crosses folders. See what the rules find with
#     flowcellar discover -L <this library> --explain
""")
    if proposals:
        out.append("#\n# Patterns the folders suggest" + (" (wide folders were sampled, so the counts are a lower bound)"
                                                         if facts["sampled"] else "") + ". Check each, then remove the '#':\n")
        for p in proposals[:24]:
            out.append(f"\n# {p['files']} file{'s' if p['files'] != 1 else ''} in {p['folders']} folder{'s' if p['folders'] != 1 else ''}"
                       + (f"; {p['check']}" if p["check"] else "") + "\n"
                       f"# [[inputs]]\n# role = {loader.toml_str(p['role'])}\n# glob = {loader.toml_str(p['glob'])}\n")
    else:
        out.append('#\n# [[inputs]]\n# role = "align.genome_bam"\n# glob = "**/analysis/alignments/{sample}_sorted.bam"\n')
    sv, bv, fv = d["scan"], d["baselines"], d["flow_cell"]
    out.append(f"""
# [variants]
# primary = "main"                # the version that counts when a sample has several

[samples]
unclassified = ["unclassified"]   # names of the reads that were given to no sample
# [samples.alias]
# "name as recorded" = "name to show"

# ---- Verdicts -------------------------------------------------------------------------------------------
[baselines]
window = {bv['window']}                       # comparable earlier runs a run is compared with
min_runs = {bv['min_runs']}                      # below this many, fixed limits are used
trusted_runs = {bv['trusted_runs']}
spread_guard = {bv['spread_guard']}

[flow_cell]
restart_max_hours = {fv['restart_max_hours']}             # a shorter earlier start on the same flow cell is a restart, not a loading
restart_max_output_frac = {fv['restart_max_output_frac']}

# ---- Scan, hub, serve -----------------------------------------------------------------------------------
[scan]
jobs = {sv['jobs']}                          # worker processes for large files; 0 chooses by itself
threads = {sv['threads']}
time_box_minutes = {sv['time_box_minutes']}            # stop starting new passes after this long; 0 = no limit
live_workers = {sv['live_workers']}                  # workers used on a root marked live = true

[hub]
sidecars = true                   # one data file per run beside the hub page, so that it opens without a server
reports = false                   # also write one offline report per run

[serve]
host = "127.0.0.1"
port = 8737
""")
    return "".join(out), facts
