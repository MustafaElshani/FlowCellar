"""The types that cross the boundaries between the four parts (ENGINE, LIBRARY, READ DATA, VIEWER).

The functions that take and return these types are in
``flowcellar.engine.api``, ``flowcellar.readdata.api``, ``flowcellar.library.api`` and ``flowcellar.report.api``.
Everything here is plain data: no behaviour, no I/O.

A *card* is a plain dict ``{"columns": {...}, "facts": {...}}``: ``columns`` is the index row of one run,
``facts`` is what the status engine needs from this
run when it is history for a later one (opaque to everybody but ENGINE). LIBRARY adds ``run_key`` and ``samples``.
A *bundle*, a *status file* and a *partial* are plain dicts in the shapes of ``flowcellar.model``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Literal

RunKind = Literal["run", "run_raw_only"]
RunState = Literal["finished", "running", "interrupted", "raw_only"]
Tier = Literal["1", "2", "all"]
Card = dict


@dataclass(frozen=True)
class RunIdentity:
    """One MinKNOW run folder, identified from its small files (ENGINE ``identify``). Cheap: no large file is read."""
    run_key: str
    flow_cell_id: str
    run_short: str                          # first 8 hex of the protocol run id
    run_dir: str                            # absolute path
    kind: RunKind = "run"
    run_state: RunState = "finished"
    protocol_run_id: str | None = None      # None for a folder with raw data only
    started: str | None = None              # RFC 3339 UTC
    files: dict[str, str] = field(default_factory=dict)     # tier 1 role -> absolute path, for the roles present
    sizes: dict[str, int] = field(default_factory=dict)     # tier 1 role -> bytes
    newest_mtime: float = 0.0               # newest mtime among the files of the folder itself
    facts: dict = field(default_factory=dict)   # device_type, product_code, kit, experiment_id, sample_id, channel_count


@dataclass(frozen=True)
class RunContext:
    """What the caller knows about where a run sits. All defaults = no library (``flowcellar run -o``)."""
    library_id: str | None = None
    library_name: str | None = None
    root_id: str | None = None
    relpath: str | None = None              # of the run folder under its root
    profile: str = "standard"               # size profile of the bundle
    assay_override: str | None = None       # --assay
    stale_after_hours: float = 6.0          # running vs interrupted
    deterministic: bool = False             # pin times, leave timings out


@dataclass(frozen=True)
class Declared:
    """What a person stated on the command line (``run``) about one run."""
    flow_cell_use: int | None = None        # --flow-cell-use N; --reused = 2
    flow_cell_check_pores: int | None = None
    baseline_group: str | None = None


@dataclass
class ProgressEvent:
    what: str                               # "reading", "done", "skipped", "error"
    path: str | None = None
    fraction: float | None = None
    mb_per_s: float | None = None
    message: str | None = None


@dataclass
class PassOptions:
    """How to run a pass over a large input."""
    jobs: int = 0                           # worker processes; 0 = min(8, cpus), 2 for a gentle root
    threads: int = 4                        # htslib threads for an unindexed BAM
    deadline: float | None = None           # time.monotonic() value at which the pass must stop; None = no limit
    gentle: bool = False                    # root marked live = true: read with 2 workers
    stall_timeout_s: float = 120.0
    progress: Callable[[ProgressEvent], None] | None = None


@dataclass(frozen=True)
class PartialRef:
    """A partial on disk (``flowcellar.model.partial``), as the result of a pass or of a cache look-up."""
    partial_id: str
    path: str
    role: str
    input_path: str
    state: str                              # ok, partially_read, truncated, error
    complete: bool                          # state ok and coverage full
    cached: bool                            # True: found, nothing was read
    seconds: float = 0.0


@dataclass(frozen=True)
class InputRule:
    """One way of finding read-data files: an ``[[inputs]]`` entry of flowcellar.toml, or a built-in default."""
    role: str
    glob: str | None = None                 # with {sample} and {variant:a|b} captures; relative to a root
    sniff: str | None = None                # e.g. "csv_columns:barcode,alias,flow_cell_id|position_id"
    variant: str | None = None
    primary: bool | None = None
    exclude: tuple[str, ...] = ()
    source: str = "default"                 # "default" (shipped) or "config"


@dataclass(frozen=True)
class Candidate:
    """A file the discovery walk found and gave a role to (LIBRARY -> READ DATA). Nothing inside it was read."""
    role: str
    path: str
    size: int
    mtime_ns: int
    root_id: str | None = None
    relpath: str | None = None
    sample: str | None = None               # the {sample} capture
    variant: str | None = None
    primary: bool | None = None
    rule: str | None = None                 # the rule that matched, in words
    experiment_dir: str | None = None       # nearest ancestor holding a run folder, a sample sheet or analysis/


@dataclass
class LinkedInput:
    """A candidate attached to runs, with its evidence (READ DATA ``link``). Equal fingerprints are one input."""
    candidate: Candidate
    fingerprint: dict                       # {"method", "value"}
    run_keys: list[str]                     # runs this input may feed; empty = not linked
    link: dict                              # {"method", "confidence", "evidence", "rule", "confirmed"}
    alt_paths: list[str] = field(default_factory=list)
    shares: dict[str, float] = field(default_factory=dict)   # run_key -> share of reads, when a file holds several runs
    scope: Literal["run", "shared", "sample_all_runs"] = "run"
    sample: str | None = None               # canonical alias after the alias map
    variant: str | None = None
    needs: list[str] = field(default_factory=list)           # reads, genome, transcriptome, or the role itself
    counted_for: dict[str, list[str]] = field(default_factory=dict)   # run_key -> needs it is the source of record for
    not_counted_reason: str | None = None
    partial_id: str | None = None           # of the pass this input needs; None for small collected files


@dataclass
class ReadData:
    """Everything read data adds to one run's bundle (READ DATA ``read_data`` -> ENGINE ``assemble``)."""
    modules: dict[str, dict] = field(default_factory=dict)      # sections by id: len, qual, aln, cdna, funnel, ...
    samples: list[dict] = field(default_factory=list)           # bundle.samples
    read_level: dict = field(default_factory=lambda: {"state": "absent"})   # bundle.tiers.read_level
    inputs: list[dict] = field(default_factory=list)            # provenance.inputs entries, ids "rd1", "rd2", ...
    basecalling: list[dict] = field(default_factory=list)       # identity.basecalling
    partial_ids: list[str] = field(default_factory=list)        # partials merged
    warnings: list[dict] = field(default_factory=list)          # provenance.warnings entries
    waiting: list[str] = field(default_factory=list)            # paths of inputs whose pass has not run yet


@dataclass
class Served:
    """One HTTP answer of ``flowcellar serve`` (VIEWER ``serve_file`` -> LIBRARY serve)."""
    body: bytes
    content_type: str
    headers: dict[str, str] = field(default_factory=dict)
    status: int = 200
