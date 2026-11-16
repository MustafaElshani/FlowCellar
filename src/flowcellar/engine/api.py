"""What the other parts call in ENGINE. Signatures are the contract.

The bodies live in the modules of this package; this file only names them.
"""
from __future__ import annotations

from typing import Sequence

from flowcellar.contracts import Card, Declared, PartialRef, PassOptions, ReadData, RunContext, RunIdentity


def identify(run_dir: str) -> RunIdentity:
    """Identify a run folder from its small files (final summary, sample sheet, the head of the report, file names):
    under a second, no large file opened. A folder with raw data only: ``kind="run_raw_only"``, identity from folder
    and file names. Raises InputError (exit code 4) when ``run_dir`` is not a run folder; a folder that is named like
    one but holds only basecaller output is not a run folder."""
    from flowcellar.engine import identity
    return identity.identify(run_dir)


def quick_bundle(run: RunIdentity, ctx: RunContext = RunContext()) -> dict:
    """The measured bundle in quick detail from the small files (seconds). A stub bundle for a raw-only folder;
    ``tiers.run_health.state = "live"`` for a running run. No status member."""
    from flowcellar.engine import assemble as _a
    return _a.assemble(run, (), None, ctx)


def full_pass(run: RunIdentity, partials_dir: str, opts: PassOptions = PassOptions()) -> list[PartialRef]:
    """One streaming pass over the sequencing summary and one over the pore-scan table, each cached as a partial
    under ``partials_dir`` keyed by the input's fingerprint. A second call reads nothing (``cached=True``).
    Stops at ``opts.deadline`` and says so (``complete=False``). Returns nothing for a running run."""
    from flowcellar.engine import passes
    return passes.full_pass(run, partials_dir, opts)


def assemble(run: RunIdentity, partials: Sequence[PartialRef] = (), read_data: ReadData | None = None,
             ctx: RunContext = RunContext()) -> dict:
    """The measured bundle: quick content, plus the fragments of ``partials`` (full detail when the summary partial is
    complete), plus the sections of ``read_data``. Reads no large input; under a second. No status member."""
    from flowcellar.engine import assemble as _a
    return _a.assemble(run, partials, read_data, ctx)


def card(bundle: dict, status_file: dict | None = None) -> Card:
    """The index row of one run: identity, tier 1 numbers, sparklines and, when ``status_file`` is given, the verdict
    columns and stored sentences, as ``{"columns": {...}, "facts": {...}}`` (the engine's
    columns; ``facts`` is whatever ``judge`` needs from this run when it is history for a later one)."""
    from flowcellar.engine import card as _c
    return _c.card(bundle, status_file)


def judge(bundle: dict, history: Sequence[Card], *, declared: Declared | None = None,
          profile: dict | None = None, baseline_groups: Sequence[dict] = (), flow_cell: dict | None = None,
          in_library: bool = True) -> dict:
    """Judge one bundle and return its status file (``flowcellar.model.status.new_status_file``).

    ``history`` holds the cards of the OTHER runs the library knows (any order); the engine keeps those that
    started before this run and are comparable, so a run is judged as of its start. ``profile`` is the
    ``[profiles]`` table of the configuration; ``flow_cell`` the ``[flow_cell]`` table. With ``in_library=False``
    (``flowcellar run -o``) history is empty and flow cell use is unknown unless ``declared`` states it.
    Read-data checks are obtained from ``flowcellar.readdata.api.read_checks(bundle)``.
    """
    from flowcellar.engine import status
    return status.judge(bundle, history, declared=declared, profile=profile, baseline_groups=baseline_groups,
                        flow_cell=flow_cell, in_library=in_library)


def judge_all(bundles: Sequence[dict], *, known: Sequence[Card] = (), profile: dict | None = None,
              baseline_groups: Sequence[dict] = (), flow_cell: dict | None = None) -> dict[str, tuple[dict, Card]]:
    """Judge a set of runs in order of run start, each against the runs before it: ``{run_key: (status file, card)}``.

    ``known`` holds the cards of runs that are already judged and are not judged again. A library that is filled in
    one go then shows what a library kept from the first day would show. An addition to the contract:
    LIBRARY may call it instead of its own loop over ``judge`` and ``card``.
    """
    from flowcellar.engine import card as _c
    from flowcellar.engine import status
    cards: list[Card] = list(known)
    out: dict[str, tuple[dict, Card]] = {}

    def start(b):
        pr = (b.get("identity") or {}).get("protocol_run") or {}
        return (pr.get("started") or "", b.get("bundle_id") or "")
    pre = {b.get("bundle_id"): _c.card(b, None) for b in bundles}
    for b in sorted(bundles, key=start):
        key = b.get("bundle_id")
        # the runs of this set that start later are history too (they never judge an earlier run; they are listed
        # as later loadings of the flow cell)
        others = cards + [c for k, c in pre.items() if k != key and k not in out]
        sf = status.judge(b, others, profile=profile, baseline_groups=baseline_groups, flow_cell=flow_cell, in_library=True)
        c = _c.card(b, sf)
        out[key] = (sf, c)
        cards.append(c)
    return out


def baselines_today(cards: Sequence[Card]) -> dict:
    """The ``baselines`` table of the index: per comparison group and metric as of today (never used for a verdict)."""
    from flowcellar.engine import status
    return status.baselines_today(cards)


def history_hash(card_: Card, history: Sequence[Card]) -> str | None:
    """The hash ``judge`` would store in ``status.judged.history_hash`` for this run against ``history``: lets the
    library flag runs "judged before older runs were added" without judging again."""
    from flowcellar.engine import status
    return status.history_hash(card_, history)


def flow_cell_columns(cards: Sequence[Card]) -> dict[str, dict]:
    """``{run_key: {"use", "use_of", "restart", "restart_of"}}`` from every card the library holds, as of today.

    A stored verdict is frozen at the run's start, so a status file cannot know of loadings that came later. The
    index may refresh these four columns from this function at every build (they are library knowledge, not a
    judgement). An addition to the contract."""
    from flowcellar.engine import flowcell_use
    facts = [c["facts"] for c in cards or () if isinstance(c, dict) and isinstance(c.get("facts"), dict) and c["facts"].get("run_key")]
    uses = flowcell_use.annotate([dict(f) for f in facts], None)
    return {k: {"use": u["index"], "use_of": u["of"], "restart": u["restart"], "restart_of": u["restart_of"]} for k, u in uses.items()}
