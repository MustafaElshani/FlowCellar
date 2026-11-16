"""One source of record per run, sample and need (finding R2).

A fragment is never added to another fragment that describes the same reads. For every run, sample and need exactly
one *input set* is counted:

=============  =================================================  ==================================================
need           feeds                                              role preference, first that exists
=============  =================================================  ==================================================
``reads``      ``yield.basecalled``, ``len``, ``qual``, read fate  ``dorado.bam``; ``dorado.summary``;
               up to "given to a sample"                           ``align.genome_bam``; ``align.transcriptome_bam``
``genome``     ``aln``, ``aln.windows``, "mapped to the genome"    ``align.genome_bam``
``transcriptome``  ``rna.tx``, ``rna.coverage``                    ``align.transcriptome_bam``
the role       ``cdna``, the later stages of read fate             pychopper statistics, quantification metadata
=============  =================================================  ==================================================

Inside one role: a file linked by its content beats one linked by a path rule; a rule marked ``primary = true``
wins; then the configured primary variant; then the file without a variant when variant files exist beside it.
What is left must be **one set**: the files one rule found in one folder (the chunks a basecaller writes; the
``bam_pass`` and ``bam_fail`` folders of one output are one set). Two sets that nothing separates are *two candidate
sources*: neither is counted, both are named, and the sections that need them say "not assessed". Nothing is summed.
"""
from __future__ import annotations

import os
import re

from flowcellar.contracts import LinkedInput

READS_PREFERENCE = ("dorado.bam", "dorado.summary", "align.genome_bam", "align.transcriptome_bam")
ROLE_NEEDS = {"dorado.bam": ("reads",), "dorado.summary": ("reads",), "align.genome_bam": ("reads", "genome"),
              "align.transcriptome_bam": ("reads", "transcriptome")}
#: needs no section of this build reads: an input wanted for nothing else is not passed over
#: (none: the alignment to the transcripts feeds the sections rna.tx and rna.coverage)
INACTIVE_NEEDS: tuple[str, ...] = ()
NEED_WORDS = {"reads": "the reads", "genome": "the genome alignment", "transcriptome": "the transcriptome alignment",
              "pychopper.stats": "the cDNA classification", "quant.salmon": "the salmon result", "quant.oarfish": "the oarfish result",
              "quant.kallisto": "the kallisto result", "polya.read_table": "the annotated read table"}
_PASS_FAIL = ("bam_pass", "bam_fail", "bam_skip", "fastq_pass", "fastq_fail", "fastq_skip", "pass", "fail")
ALL = "*"                   # the "sample" of a file that holds several samples, or whose reads name none


def needs_of(role: str) -> tuple[str, ...]:
    return ROLE_NEEDS.get(role, (role,))


def bundle_need(need: str) -> str:
    """The word the bundle schema uses for a need."""
    return need if need in ("reads", "genome", "transcriptome", "mods") else "collected"


def set_key(li: LinkedInput) -> tuple:
    """Files with an equal key are parts of one output and are added up; different keys are different sources."""
    c = li.candidate
    parts = [("*" if p in _PASS_FAIL else p) for p in os.path.dirname(c.path).split(os.sep)]
    rule = re.sub(r"\b(bam|fastq)_(pass|fail|skip)\b", r"\1_*", c.rule or "")      # the pass and the fail rule find one output
    return (c.role, li.variant, rule, os.sep.join(parts))


def _sample(li: LinkedInput) -> str:
    return li.sample or ALL


def _strong(li: LinkedInput) -> bool:
    return li.link.get("confidence") in ("certain", "high") or bool(li.link.get("confirmed"))


def _pick(group: list[LinkedInput], need: str, primary_variant: str | None) -> tuple[list[LinkedInput], dict[int, str], list[list[LinkedInput]] | None]:
    """(inputs counted, reason by id() for the inputs set aside, the candidate sets when nothing separates them)."""
    aside: dict[int, str] = {}
    pool = group
    if need == "reads":
        role = next(r for r in READS_PREFERENCE if any(li.candidate.role == r for li in group))
        pool = [li for li in group if li.candidate.role == role]
    if any(_strong(li) for li in pool) and not all(_strong(li) for li in pool):
        for li in pool:
            if not _strong(li):
                aside[id(li)] = "a file linked by what it contains is counted instead of this one, which is linked by a path rule"
        pool = [li for li in pool if _strong(li)]
    flagged = [li for li in pool if li.candidate.primary]
    if flagged and len(flagged) < len(pool):
        for li in pool:
            if not li.candidate.primary:
                aside[id(li)] = "another file of this sample comes from a rule marked primary and is the one counted"
        pool = flagged
    variants = {li.variant for li in pool}
    if len(variants) > 1:
        keep = primary_variant if primary_variant in variants else None if None in variants else False
        if keep is not False:
            for li in pool:
                if li.variant != keep:
                    aside[id(li)] = (f"another version of this file is the one counted (the primary variant is "
                                     f"{keep or 'the file without a variant name'}); this one is version '{li.variant}'")
            pool = [li for li in pool if li.variant == keep]
    sets: dict[tuple, list[LinkedInput]] = {}
    for li in pool:
        sets.setdefault(set_key(li), []).append(li)
    if len(sets) > 1:
        return [], aside, list(sets.values())
    return pool, aside, None


def choose(inputs: list[LinkedInput], config: dict | None = None) -> None:
    """Fill ``counted_for`` (run key -> needs) on every input, and say why an input is not counted.

    ``not_counted_reason`` is set only for an input that is counted for nothing because of two candidate sources;
    other reasons (a variant, a lower role, a need nothing reads) go to ``link["not_counted"]``, and a conflict that
    concerns only some runs or needs to ``link["conflicts"]`` (run key -> need -> sentence).
    """
    primary_variant = ((config or {}).get("variants") or {}).get("primary")
    by_run: dict[str, list[LinkedInput]] = {}
    for li in inputs:
        li.counted_for = {k: [] for k in li.run_keys}
        li.not_counted_reason = None
        li.link.pop("conflicts", None)
        li.link.pop("not_counted", None)
        for k in li.run_keys:
            by_run.setdefault(k, []).append(li)
    reasons: dict[int, str] = {}
    for run_key, mine in by_run.items():
        first = ("reads", "genome", "transcriptome")
        needs = sorted({n for li in mine for n in li.needs}, key=lambda n: (first.index(n) if n in first else len(first), n))
        for need in needs:
            have = [li for li in mine if need in li.needs]
            if need == "basecaller.sample_sheet":               # sheets do not compete: their rows are put together
                for li in have:
                    li.counted_for[run_key].append(need)
                continue
            named = sorted({_sample(li) for li in have} - {ALL})
            for sample in (named or [ALL]):
                group = [li for li in have if _sample(li) in (sample, ALL)]
                counted, aside, conflict = _pick(group, need, primary_variant)
                reasons.update(aside)
                if conflict:
                    paths = sorted(s[0].candidate.path for s in conflict)
                    what = NEED_WORDS.get(need, need)
                    whom = f" of sample {sample}" if sample != ALL else ""
                    msg = (f"two candidate sources for {what}{whom}: {paths[0]} and {paths[1]}"
                           + (f" (and {len(paths) - 2} more)" if len(paths) > 2 else "")
                           + "; nothing was added up. Name the one to count: put primary = true on its [[inputs]] rule, or give the other one a variant")
                    for s in conflict:
                        for li in s:
                            li.link.setdefault("conflicts", {}).setdefault(run_key, {})[need] = msg
                    continue
                for li in counted:
                    if need in INACTIVE_NEEDS:
                        reasons.setdefault(id(li), f"no section of this version of FlowCellar reads {NEED_WORDS.get(need, need)}")
                    elif need not in li.counted_for[run_key]:
                        li.counted_for[run_key].append(need)
    for li in inputs:
        if not li.run_keys or any(li.counted_for.values()):
            continue
        conflicts = li.link.get("conflicts") or {}
        first = next((m for per in conflicts.values() for m in per.values()), None)
        if first:
            li.not_counted_reason = first
        elif id(li) in reasons:
            li.link["not_counted"] = reasons[id(li)]
        elif li.candidate.role in ROLE_NEEDS:
            li.link["not_counted"] = "another kind of file of this sample is the source of record for the reads"
