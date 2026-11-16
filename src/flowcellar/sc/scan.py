"""What a scan of a library does with the analyses its walk found (called once by ``flowcellar.library.scan``).

Order of work, so that the hub shows every analysis after the first minute:

1. each anchor file the walk reported (the ``metrics_summary.csv`` of an ``outs/`` folder, the report of a
   workflow output folder) becomes its analysis folders; a folder whose files have the same names, sizes and
   times as last time is not opened;
2. new analyses are read without the one slow step (the pass over the workflow's per-read table), linked to the
   runs of the library, judged and stored; the index is written;
3. the per-read tables are read, one analysis at a time, inside what is left of the scan's time box; an
   analysis that had to wait says so and is finished by the next scan.

Nothing is written outside the library, and a folder that cannot be read costs one line, not the scan.
"""
from __future__ import annotations

import os
import time

from flowcellar.model.status import VERDICT_WORD
from flowcellar.sc import PIPELINE_WORD, ROLE_PREFIX, detect
from flowcellar.sc import collect as collect_mod
from flowcellar.sc import store as store_mod
from flowcellar.sc.catalogue import unit_words
from flowcellar.sc.detect import Claim
from flowcellar.util.fingerprint import hash_obj


def claims_of(anchors) -> list[Claim]:
    """The analysis folders behind the anchor files of a walk, each once."""
    seen: dict[str, Claim] = {}
    for c in anchors:
        if not str(c.role).startswith(ROLE_PREFIX):
            continue
        for claim in detect.claims_from_anchor(c.role, c.path, getattr(c, "relpath", None)):
            seen.setdefault(claim.path, claim)
    return [seen[p] for p in sorted(seen)]


def line(row: dict) -> str:
    """One line for the person: label, what it is, the verdict and its first reason."""
    word = VERDICT_WORD.get(row.get("verdict"), "Not assessed")
    what = " · ".join(x for x in ({"visium": "Visium", "visium_hd": "Visium HD", "singlecell": "single cell"}.get(row.get("assay"), ""),
                                  (PIPELINE_WORD.get(row.get("pipeline"), "") + (" " + row["version"] if row.get("version") else "")).strip()) if x)
    text = (row.get("reasons") or [row.get("headline") or ""])[0]
    return f"{str(row.get('label') or row.get('id'))[:34]:<34} {what:<34} {word:<10} {text}".rstrip()


def run(scan, anchors) -> list[dict]:
    """Bring the analyses of ``scan.lib`` up to date. ``scan`` is the running ``library.scan.Scan``; ``anchors``
    are the candidates of the walk whose role starts with ``analysis.``. Returns the rows of the index table."""
    lib, o = scan.lib, scan.o
    r = store_mod.Refresher(lib, say=scan.say, note=scan.note)
    claims = claims_of(anchors)
    state = r.state
    known_claims = state.setdefault("claims", {})
    by_id: dict[str, list[Claim]] = {}
    fresh: list[Claim] = []
    for c in claims:
        sig = hash_obj(collect_mod.signature(c), 24)
        k = known_claims.get(c.path)
        if k and k.get("sig") == sig and k.get("id"):
            by_id.setdefault(k["id"], []).append(c)           # unchanged since the last scan: not opened
        else:
            fresh.append(c)
    for aid, cs in r.identify(fresh).items():
        by_id.setdefault(aid, []).extend(cs)
        for c in cs:
            known_claims[c.path] = {"id": aid, "sig": hash_obj(collect_mod.signature(c), 24)}
    for path in [p for p in known_claims if p not in {c.path for c in claims}]:
        if not scan._unseen_is_unknown(path, None) and not os.path.isdir(path):
            known_claims.pop(path, None)
    scan.note("analyses.found", anchors=len(anchors), folders=len(claims), analyses=len(by_id), opened=len(fresh))

    no_pass = bool(o.quick or o.tier == "1" or not scan._reads_kind("analyses"))
    workers = int(o.jobs or scan.cfg["scan"].get("jobs") or 0)

    def options(read_pass: bool) -> collect_mod.Options:
        deadline = None
        if scan.deadline is not None:
            deadline = time.time() + max(0.0, scan.deadline - _now(scan))
        return collect_mod.Options(read_pass=read_pass, workers=workers, deadline=deadline)

    runs = store_mod.run_facts(scan.cards, scan.states, scan.cfg.get("roots") or [])
    profiles = scan.cfg.get("profiles")
    # ---- step 2: everything that is cheap, for every analysis ---------------------------------------------
    quick_opts, full_opts = options(False), options(True)
    measured: dict[str, dict] = {}
    slow: list[str] = []
    for aid in sorted(by_id):
        best = store_mod._best(by_id[aid])
        has_table = "read_summary" in best.files
        full_key = store_mod._cache_key(aid, best, full_opts)
        entry = state["analyses"].get(aid) or {}
        full_ready = entry.get("collected") == full_key and os.path.isfile(r.store.cache_path(full_key))
        if full_ready or not has_table:
            m = r.collect(aid, by_id[aid], full_opts)            # read before (nothing is opened), or nothing slow to read
        else:
            if not no_pass:
                slow.append(aid)
            m = r.collect(aid, by_id[aid], quick_opts)
        if m is not None:
            measured[aid] = m
        scan._lease.touch_if_due()
    if measured:
        r.finish(measured, runs, profiles=profiles)
    r.mark_unseen(set(by_id), unknown=lambda p: scan._unseen_is_unknown(p, _root_of(scan, p)))
    for aid in r.new:
        row = (r.cards.get(aid) or {}).get("row")
        if row:
            scan.say("analysis.new", "analysis  " + line(row), analysis=aid, verdict=row.get("verdict"))
    scan._analyses = r.rows()
    if r.changed:
        scan.hub_dirty.update(r.changed)
        scan._write_index(scanning=True, force=True)
        scan._write_hub()

    # ---- step 3: the per-read tables, inside the time box --------------------------------------------------
    waiting = 0
    for i, aid in enumerate(slow):
        if scan._out_of_time():
            waiting = len(slow) - i
            break
        best = store_mod._best(by_id[aid])
        label = ((r.cards.get(aid) or {}).get("row") or {}).get("label") or aid
        size = os.path.getsize(best.files["read_summary"]) if os.path.isfile(best.files["read_summary"]) else 0
        scan.say("analysis.reading", f"analysis  {label}: reading its per-read table ({size / 1e9:.1f} GB), {i + 1} of {len(slow)}",
                 analysis=aid, bytes=size)
        t0 = time.monotonic()
        m = r.collect(aid, by_id[aid], options(True))
        scan._lease.touch_if_due()
        if m is None:
            continue
        r.finish({aid: m}, runs, profiles=profiles)
        entry = r.state["analyses"].get(aid) or {}
        if entry.get("read_pass") == "waiting":
            waiting += 1
        row = (r.cards.get(aid) or {}).get("row") or {}
        scan.say("analysis.read", "analysis  " + line(row) + f"   ({time.monotonic() - t0:.0f} s)", analysis=aid)
        scan._analyses = r.rows()
        scan.hub_dirty.update(r.changed)
        scan._write_index(scanning=True, force=True)
        scan._write_hub()
    if no_pass:
        waiting = sum(1 for aid in by_id if "read_summary" in store_mod._best(by_id[aid]).files
                      and (r.state["analyses"].get(aid) or {}).get("read_pass") == "waiting")
    r.store.write_state(r.state)
    rows = r.rows()
    scan._analyses = rows
    if rows or claims:
        n = len(rows)
        units = sorted({unit_words(x.get("assay"))["plural"] for x in rows})
        parts = [f"{n} analys{'is' if n == 1 else 'es'} ({', '.join(units)})" if units else f"{n} analyses"]
        if r.new:
            parts.append(f"{len(r.new)} new")
        if waiting:
            parts.append(f"{waiting} with a per-read table still to be read (scan again)")
        if r.failed:
            parts.append(f"{len(r.failed)} could not be read")
        scan.say("analyses.done", "analyses: " + " · ".join(parts), analyses=n, new=len(r.new), waiting=waiting, failed=len(r.failed))
    return rows


def _now(scan) -> float:
    from flowcellar.library import scan as scan_mod
    return scan_mod._monotonic()


def _root_of(scan, path: str) -> str | None:
    best = None
    for root in scan.cfg.get("roots") or []:
        rp = root["path"].rstrip("/")
        if path == rp or path.startswith(rp + "/"):
            if best is None or len(rp) > len(best["path"].rstrip("/")):
                best = root
    return best["id"] if best else None
