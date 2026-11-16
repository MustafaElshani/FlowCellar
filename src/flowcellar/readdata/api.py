"""What the other parts call in READ DATA. Signatures are the contract.

Reading BAM files needs pysam (the extra ``bam``). Without it everything that does not open a BAM still works:
collected files (pychopper statistics, quantification metadata, sample sheets), read summaries, linking by location
and by sample name; a BAM is then reported as not readable, with one line that names the extra.
"""
from __future__ import annotations

import time
from typing import Mapping, Sequence

from flowcellar.contracts import (
    Candidate,
    InputRule,
    LinkedInput,
    PartialRef,
    PassOptions,
    ProgressEvent,
    ReadData,
    RunIdentity,
)
from flowcellar.model import partial as partial_model

SHEET_SNIFF = "csv_columns:barcode,alias,flow_cell_id|position_id"


def default_input_rules() -> list[InputRule]:
    """The generic rules shipped with the package: the layouts MinKNOW and Dorado write themselves, and recognition
    by content. Never a lab's own folder names (those live in the library's flowcellar.toml)."""
    return [
        # MinKNOW with live basecalling, and dorado demux: <run or output folder>/bam_pass/<barcode or alias>/*.bam
        InputRule(role="dorado.bam", glob="**/bam_pass/{sample}/*.bam"),
        InputRule(role="dorado.bam", glob="**/bam_fail/{sample}/*.bam"),
        InputRule(role="dorado.bam", glob="**/bam_pass/*.bam"),
        InputRule(role="dorado.bam", glob="**/bam_fail/*.bam"),
        # dorado demux before 1.0: one file per barcode, named after the kit
        InputRule(role="dorado.bam", glob="**/SQK-*_barcode*.bam"),
        # told by content, whatever the file is called
        InputRule(role="basecaller.sample_sheet", sniff=SHEET_SNIFF),
    ]


def sniff(path: str) -> str | None:
    """The role of a small file told by its content (a basecaller sample sheet by its columns; a partial), or None.
    Reads at most the first 64 kB. Also tells pychopper statistics, quantification metadata, read summaries with
    basecalled columns, and whether a BAM is unaligned (``dorado.bam``) or aligned."""
    from flowcellar.readdata.sample_sheet import sniff as _sniff
    return _sniff(path)


def candidates_for_run(run: RunIdentity, paths: Sequence[str] = (), explicit: Mapping[str, Sequence[str]] | None = None,
                       sample_sheet: str | None = None) -> list[Candidate]:
    """For ``flowcellar run RUN_DIR [TIER2_PATH...] --input ROLE=PATH --sample-sheet CSV`` without a library: the
    read-data files under the run folder itself and under ``paths``, by the default rules and by content."""
    from flowcellar.readdata.finder import candidates_for_run as _find
    return _find(run, paths, explicit, sample_sheet)


def link(candidates: Sequence[Candidate], runs: Sequence[RunIdentity], *, config: Mapping | None = None,
         probes_dir: str | None = None) -> list[LinkedInput]:
    """Attach files to runs by the identifiers inside them first and by configured patterns second, make
    equal fingerprints one input, and choose ONE source of record per run, sample and need. Two candidates that no
    rule separates are both returned with ``not_counted_reason`` set; nothing is ever summed. Identity evidence is
    cached under ``probes_dir`` by fingerprint. ``config`` is the parsed flowcellar.toml (variants, link.rules,
    samples).

    Besides the members of the contract, ``LinkedInput.link`` may hold ``conflicts`` (run key -> need -> sentence,
    when a file is counted for one run and in conflict for another), ``not_counted`` (why a file that is no conflict
    is not counted: a variant, a need nothing reads) and ``labels`` (sample labels seen in a BAM)."""
    from flowcellar.readdata.link import link as _link
    return _link(list(candidates), list(runs), config=dict(config or {}), probes_dir=probes_dir)


def read_pass(inp: LinkedInput, partials_dir: str, known_runs: Sequence[str], opts: PassOptions = PassOptions()) -> PartialRef | None:
    """One exact pass over one input file (BAM, table), every read routed to its run; cached as a partial under
    ``partials_dir``. None for a small collected file that needs no pass.

    The partial is a pure function of the file: reads are keyed by the run their own identifiers name
    (``known_runs`` is not needed for that), reads that name no run are kept apart, and the link decides at assembly
    where those belong. A time-limited pass (``opts.deadline``) is cached with its coverage and done again next time.
    """
    from flowcellar.readdata import link as L
    from flowcellar.readdata import passes
    c = inp.candidate
    if c.role not in L.PASS_ROLES or not inp.partial_id:
        return None
    found = partial_model.find_partial(partials_dir, inp.partial_id)
    if found:
        try:
            p = partial_model.read_partial(found)
            if p.get("state") in ("ok", "truncated"):       # a file that is cut short stays so while its fingerprint holds
                return PartialRef(inp.partial_id, found, c.role, c.path, p["state"], partial_model.is_complete(p), True, 0.0)
        except Exception:       # noqa: BLE001  an unreadable cache entry is made again
            pass
    t0 = time.monotonic()
    jobs = 2 if (opts.gentle and (opts.jobs == 0 or opts.jobs > 2)) else opts.jobs

    def progress(done: int, total: int) -> None:
        dt = max(time.monotonic() - t0, 1e-6)
        opts.progress(ProgressEvent("reading", path=c.path, fraction=done / max(total, 1), mb_per_s=c.size * done / max(total, 1) / 1e6 / dt))

    p = passes.run_pass(inp, jobs=jobs, threads=opts.threads, deadline=opts.deadline, stall_timeout_s=opts.stall_timeout_s,
                        progress=progress if opts.progress else None, inprocess=L.inprocess())
    seconds = time.monotonic() - t0
    if p["state"] == "error":
        if opts.progress:
            opts.progress(ProgressEvent("error", path=c.path, message="; ".join(e.get("message", "") for e in p.get("errors") or [])[:300]))
        return PartialRef(inp.partial_id, "", c.role, c.path, "error", False, False, seconds)
    path = partial_model.write_partial(partials_dir, p)
    if opts.progress:
        opts.progress(ProgressEvent("done", path=c.path, fraction=(p.get("coverage") or {}).get("fraction"), mb_per_s=c.size / 1e6 / max(seconds, 1e-6)))
    return PartialRef(p["partial_id"], path, c.role, c.path, p["state"], partial_model.is_complete(p), False, seconds)


def pass_signature() -> str:
    """An addition to the contract: one short text that changes whenever a pass over a large file would store
    something else than before (the versions of the BAM pass and of the table pass, and the settings that shape
    them). The scan keeps it with its links: when it changed, the partial id of every file is worked out again, and
    a file read by the earlier pass is read again."""
    from flowcellar.readdata import link as L
    from flowcellar.readdata import readtable
    from flowcellar.util.fingerprint import hash_obj
    return hash_obj({"bam": dict(L.PASS_MODULES), "tables": dict(readtable.MODULES), "tables_config": readtable.config_hash()}, 16)


def link_version() -> int:
    """An addition to the contract: the version of the rules that attach files to runs. The scan keeps it in the
    signature of its saved links: when it changed, the files are linked again (nothing is opened for that)."""
    from flowcellar.readdata import link as L
    return L.LINK_VERSION


def assembly_signature() -> str:
    """An addition to the contract: one short text that changes when the sections of a run would be put together
    differently from the same links and the same partials (the versions of the section builders). The scan keeps it
    with every run: when it changed, the bundle is assembled and judged again. No large file is read for that."""
    from flowcellar.readdata import fate, pychopper, rnacomp, sections, sections_align
    from flowcellar.util.fingerprint import hash_obj
    return hash_obj({"sections": sections.VERSION, "align": sections_align.VERSION, "fate": fate.VERSION, "pychopper": pychopper.VERSION,
                     "rnacomp": getattr(rnacomp, "VERSION", "")}, 16)


def pass_kind(role: str) -> str | None:
    """An addition to the contract: which kind of large file a role is, for ``scan --read``: ``"bam"`` (BAM files and
    read summaries), ``"tables"`` (annotated read tables), or None for a small collected file that needs no pass."""
    from flowcellar.readdata import link as L
    if role in L.TABLE_PASS_ROLES:
        return "tables"
    return "bam" if role in L.PASS_ROLES else None


def read_data(run: RunIdentity, linked: Sequence[LinkedInput], partials_dir: str | None, *, tier1: dict | None = None,
              config: Mapping | None = None) -> ReadData:
    """The read-level sections of one run from its partials and its small collected files. Reads no BAM. ``tier1`` is
    the run's measured bundle so far (read fate starts at "sequenced"). Inputs whose pass has not run are listed in
    ``ReadData.waiting``."""
    from flowcellar.readdata.sections import build_read_data
    return build_read_data(run, list(linked), partials_dir, tier1=tier1, config=dict(config or {}))


def read_checks(bundle: dict, profiles: Mapping | None = None) -> list[dict]:
    """The read-data rules of a bundle in the stored rule form (``tier: "read_level"``), each evaluated only on a read
    population it applies to. ``profiles`` is the ``[profiles]`` table of the library's configuration: a
    ``[profiles.<assay group>.rules.<check id>]`` entry with ``warn``, ``fail`` or ``enabled = false`` replaces the
    built-in limit of that check. Called by ``flowcellar.engine.api.judge``."""
    from flowcellar.readdata.checks import read_checks as _checks
    return _checks(bundle, dict(profiles) if profiles else None)


def card_columns(bundle: dict) -> dict:
    """The read-data columns of the index row."""
    from flowcellar.readdata.cards import card_columns as _cols
    return _cols(bundle)


def sample_rows(bundle: dict) -> list[dict]:
    """The rows this run adds to the ``samples`` table of the index."""
    from flowcellar.readdata.cards import sample_rows as _rows
    return _rows(bundle)


def candidates_in_tree(top: str, rules: Sequence[InputRule] | None = None, *, root_id: str | None = None,
                       min_age_seconds: float = 0.0, prune: Sequence[str] = ()) -> list[Candidate]:
    """An addition to the contract, for callers without a library walk: the files below ``top`` that ``rules``
    (default: ``default_input_rules()``) name by pattern, resolved by directed descent from their anchors."""
    from flowcellar.readdata import patterns
    found, _ = patterns.find(top, list(rules) if rules is not None else default_input_rules(), root_id=root_id,
                             min_age_seconds=min_age_seconds, prune=tuple(prune))
    return found
