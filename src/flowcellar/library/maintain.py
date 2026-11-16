"""Work on a library that reads no input: rebuilding the index, judging again, checking (``recompute``,
``validate``, ``doctor``)."""
from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import os
import platform
import sys
import time
from importlib import resources

from flowcellar import LAYOUT_VERSION, SCHEMA_VERSION, __version__
from flowcellar.engine import api as engine_api
from flowcellar.library import index as index_mod
from flowcellar.library.cards import card_for, judge_settings
from flowcellar.library.lease import Lease, read_holder
from flowcellar.library.store import LibraryStore, is_library
from flowcellar.model import bundle as bundle_model
from flowcellar.model import schema as schema_model
from flowcellar.model import status as status_model
from flowcellar.report import api as report_api
from flowcellar.util import log as log_mod
from flowcellar.util.errors import FlowCellarError, NotBuiltError, ValidationError
from flowcellar.util.timeutil import now_z


def rebuild_index_and_hub(store: LibraryStore, *, hub: bool = True) -> dict:
    """Write the index from the cards on disk, then the hub. Keeps the notices and progress of the last scan."""
    index = index_mod.index_of_library(store)
    store.write_index(index)
    if hub:
        try:
            report_api.write_hub(store, sidecars=bool(store.config["hub"]["sidecars"]), reports=bool(store.config["hub"]["reports"]))
        except FlowCellarError as e:
            log_mod.event("hub.error", f"the hub was not written: {e.message}", level=30)
    return index


# ------------------------------------------------------------------------------------------------ recompute
def _strip(sf: dict) -> dict:
    d = json.loads(json.dumps(sf))
    st = d.get("status") or {}
    st.pop("evaluated", None)
    if isinstance(st.get("judged"), dict):
        st["judged"].pop("at", None)
    return d


def recompute_status(store: LibraryStore, keys: list[str], *, rebaseline: bool = False, wait: bool = True) -> dict:
    """Judge runs again with the current limits; no input is read and no bundle is written.

    Each run is judged against the library's runs that started before it. Without ``rebaseline`` a run whose
    history differs from the one it was judged against (older runs were added since) is LEFT AS IT IS and
    reported: its verdict stays frozen until the person asks for it with ``--rebaseline``.
    Returns ``{"judged": [...], "changed": [{"run_key", "was", "now"}], "held": [...], "failed": {...}}``.
    """
    out = {"judged": [], "changed": [], "held": [], "failed": {}}
    with Lease(store.layout.lease, wait=wait):
        cards = store.cards()

        def start(k: str) -> str:
            return str((cards.get(k) or {}).get("columns", {}).get("started") or "")
        for key in sorted(keys, key=lambda k: (start(k), k)):                  # in order of run start
            try:
                bundle = store.read_measured(key)
                old = store.read_status_file(key)
                history = [c for k, c in cards.items() if k != key]
                if not rebaseline and old is not None and key in cards:
                    stored = ((old.get("status") or {}).get("judged") or {}).get("history_hash")
                    try:
                        now_hash = engine_api.history_hash(cards[key], history)
                    except NotBuiltError:
                        now_hash = stored
                    if stored is not None and now_hash is not None and stored != now_hash:
                        out["held"].append(key)
                        continue
                sf = engine_api.judge(bundle, history, in_library=True, **judge_settings(store.config))
                sf["content_sha256"] = bundle.get("content_sha256")
                rel = sf.setdefault("relations", {})
                for member in ("paths", "projects"):                           # library knowledge is kept as it is
                    rel[member] = ((old or {}).get("relations") or {}).get(member) or []
                was = ((old or {}).get("status") or {}).get("verdict")
                now = (sf.get("status") or {}).get("verdict")
                if old is None or _strip(old) != _strip(sf):
                    store.write_status_file(key, sf)
                    if was != now:
                        out["changed"].append({"run_key": key, "was": was, "now": now})
                st = store.read_state(key) or {}
                card = card_for(bundle, sf, extra={"bundle_bytes": st.get("bundle_bytes"), "bundle_gzip_bytes": st.get("bundle_gzip_bytes"),
                                                   "source_state": st.get("source_state", "ok")})
                if cards.get(key) != card:
                    store.write_card(key, card)
                cards[key] = card
                out["judged"].append(key)
            except NotBuiltError:
                raise
            except Exception as e:  # noqa: BLE001  one run must not stop the others
                out["failed"][key] = getattr(e, "message", None) or str(e)
        _refresh_history_notice(store, cards)
        rebuild_index_and_hub(store)
    return out


def _refresh_history_notice(store: LibraryStore, cards: dict) -> None:
    """The notice "N runs were judged before older runs were added" of the last scan, counted again after judging:
    it must not go on asking for the command that has just been run."""
    kept = store.read_aux(store.layout.notices, None)
    if not isinstance(kept, dict) or not any(n.get("code") == "judged_before_older_runs" for n in kept.get("notices") or []):
        return
    stale = []
    try:
        for key, card in cards.items():
            stored = (card.get("columns") or {}).get("history_hash")
            if stored is None:
                continue
            now = engine_api.history_hash(card, [c for k, c in cards.items() if k != key])
            if now is not None and now != stored:
                stale.append(key)
    except NotBuiltError:
        return
    notices = [n for n in kept["notices"] if n.get("code") != "judged_before_older_runs"]
    if stale:
        notices.append({"code": "judged_before_older_runs", "level": "info", "root": None, "since": None,
                        "message": f"{len(stale)} run{'s were' if len(stale) != 1 else ' was'} judged before older runs "
                                   "were added to the library. Their verdicts are kept; "
                                   f"'flowcellar recompute --status-only --rebaseline -L {store.path}' judges them again.",
                        "runs": sorted(stale)[:50]})
    kept["notices"] = notices
    store.write_aux(store.layout.notices, kept)


# ------------------------------------------------------------------------------------------------ validate
def check_bundle(bundle: dict, *, where: str, status_file: dict | None = None) -> list[str]:
    """Problems of one bundle: the schema, and the invariants that hold for every bundle whatever its content."""
    problems: list[str] = []
    try:
        bundle_model.check_version(bundle, where)
    except ValidationError as e:
        return [e.message]
    problems += schema_model.schema_errors(bundle)
    if bundle.get("kind") == "run" and bundle.get("bundle_id") != (bundle.get("identity") or {}).get("run_key"):
        problems.append("/bundle_id: differs from identity.run_key")
    stored = bundle.get("content_sha256")
    if stored and stored != bundle_model.content_sha256(bundle):
        problems.append("/content_sha256: does not match the content (the file was changed after it was written)")
    parts = ((bundle.get("provenance") or {}).get("partials")) or []
    if parts != sorted(parts) or len(set(parts)) != len(parts):
        problems.append("/provenance/partials: must be sorted and without repeats")
    st = bundle.get("status")
    if st:
        tiers = st.get("by_tier") or {}
        if "run_health" in tiers and "read_level" in tiers:
            want = status_model.combine_verdicts(tiers["run_health"], tiers["read_level"])
            if st.get("verdict") != want:
                problems.append(f"/status/verdict: is {st.get('verdict')!r} but its two parts give {want!r}")
        for i, r in enumerate(st.get("rules") or []):
            if r.get("outcome") != "triggered" and r.get("severity") not in (None, "none"):
                problems.append(f"/status/rules/{i}: severity must be none unless the rule triggered")
        for part in ("run_health", "read_level"):
            if tiers.get(part) in ("pass", "warn", "fail"):
                worst = status_model.worst_of_rules(st.get("rules") or [], part)
                if status_model._RANK[worst] > status_model._RANK[tiers[part]]:
                    problems.append(f"/status/by_tier/{part}: is {tiers[part]!r} but a triggered rule says {worst!r}")
    if status_file is not None and status_file.get("content_sha256") not in (None, stored):
        problems.append("status file: was written for other content than this bundle holds (scan the library again)")
    if bundle.get("kind") == "analysis":            # the add-on knows what else must hold for its bundles
        from flowcellar.sc import check as sc_check
        problems += sc_check.problems(bundle)
    return problems


def validate_paths(paths: list[str]) -> list[dict]:
    """``[{"path", "run_key", "ok", "problems"}]`` for bundle files and for every bundle of a library folder."""
    out = []
    for path in paths:
        if os.path.isdir(path):
            if not is_library(path):
                out.append({"path": path, "run_key": None, "ok": False, "problems": ["a folder that is not a FlowCellar library"]})
                continue
            store = LibraryStore.open(path)
            for key in store.run_keys():
                where = store.layout.bundle(key)
                try:
                    sf = store.read_status_file(key)
                    problems = check_bundle(store.composed(key), where=where, status_file=sf)
                    if sf is None:
                        problems.append("no status file: the run was never judged (scan the library again)")
                except (FlowCellarError, OSError, ValueError) as e:
                    if isinstance(e, FlowCellarError) and e.exit_code != ValidationError.exit_code:
                        raise
                    problems = [getattr(e, "message", None) or str(e)]
                out.append({"path": where, "run_key": key, "ok": not problems, "problems": problems})
            out.extend(_validate_analyses(path))
        else:
            try:
                bundle = bundle_model.read_bundle(path)
                problems = check_bundle(bundle, where=path)
                key = bundle.get("bundle_id")
            except (ValidationError, OSError, ValueError) as e:
                problems, key = [getattr(e, "message", None) or f"cannot be read ({e})"], None
            out.append({"path": path, "run_key": key, "ok": not problems, "problems": problems})
    return out


def _validate_analyses(library_path: str) -> list[dict]:
    """The same checks for the single-cell and spatial analyses of a library, when it holds any."""
    try:
        if not any(n.startswith("an_") for n in os.listdir(os.path.join(library_path, "bundles"))):
            return []
    except OSError:
        return []
    from flowcellar.sc import store as sc_store
    astore = sc_store.AnalysisStore(library_path)
    out = []
    for aid in astore.ids():
        where = astore.bundle_path(aid)
        try:
            sf = astore.read_status_file(aid)
            problems = check_bundle(astore.composed(aid), where=where, status_file=sf)
            if sf is None:
                problems.append("no status file: the analysis was never judged (scan the library again)")
        except (FlowCellarError, OSError, ValueError) as e:
            if isinstance(e, FlowCellarError) and e.exit_code != ValidationError.exit_code:
                raise
            problems = [getattr(e, "message", None) or str(e)]
        out.append({"path": where, "run_key": aid, "ok": not problems, "problems": problems})
    return out


# ------------------------------------------------------------------------------------------------ doctor
def _fs_type(path: str) -> str | None:
    """File-system type of the mount that holds ``path`` (Linux: /proc/mounts)."""
    try:
        best, kind = "", None
        real = os.path.realpath(path)
        with open("/proc/mounts") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 3 and (real == parts[1] or real.startswith(parts[1].rstrip("/") + "/")) and len(parts[1]) >= len(best):
                    best, kind = parts[1], parts[2]
        return kind
    except OSError:
        return None


def doctor(library_dir: str | None = None) -> dict:
    """Facts about the installation and, when given, a library and its roots. Reads no run data."""
    from flowcellar.discovery import check_roots
    rep: dict = {"flowcellar": __version__, "schema": SCHEMA_VERSION, "layout": LAYOUT_VERSION,
                 "python": platform.python_version(), "platform": f"{platform.system()} {platform.machine()}",
                 "executable": sys.executable, "packages": {}, "checks": []}

    def check(name: str, ok: bool | None, message: str) -> None:
        rep["checks"].append({"name": name, "ok": ok, "message": message})

    for pkg, need in (("numpy", "required"), ("pyarrow", "required"), ("click", "required"),
                      ("pysam", "extra 'bam': reading BAM files"), ("jsonschema", "extra 'validate': flowcellar validate")):
        try:
            importlib.import_module(pkg)
            try:
                rep["packages"][pkg] = importlib.metadata.version(pkg)
            except importlib.metadata.PackageNotFoundError:
                rep["packages"][pkg] = "installed"
        except Exception:  # noqa: BLE001
            rep["packages"][pkg] = None
            check(f"package {pkg}", False if need == "required" else None,
                  f"not installed ({need})" + ("" if need == "required" else f'; pip install "flowcellar[{need.split(chr(39))[1]}]"'))
    # the extra of the single-cell and spatial add-on is looked up, not imported: a doctor's visit loads no HDF5 library
    try:
        rep["packages"]["h5py"] = importlib.metadata.version("h5py")
    except importlib.metadata.PackageNotFoundError:
        rep["packages"]["h5py"] = None
        check("package h5py", None, "not installed (extra 'sc': the HDF5 files of a Space Ranger analysis); pip install \"flowcellar[sc]\"")
    try:
        man = report_api.viewer_manifest()
        files = man.get("files") if isinstance(man.get("files"), dict) else {}
        wrong = []
        for name, meta in files.items():                  # the hashes are computed again, not taken on trust
            try:
                data = resources.files("flowcellar").joinpath("viewer_dist", name).read_bytes()
            except OSError:
                wrong.append(f"{name} is missing")
                continue
            if isinstance(meta, dict) and meta.get("sha256") and hashlib.sha256(data).hexdigest() != meta["sha256"]:
                wrong.append(f"{name} does not match its hash in the manifest")
        rep["viewer"] = {"built": True, "files": len(files), "problems": wrong}
        check("viewer", not wrong, f"the built viewer is in the package ({len(files)} files, hashes match)" if not wrong
              else "the built viewer is damaged: " + "; ".join(wrong) + " (build it again)")
    except FlowCellarError as e:
        rep["viewer"] = {"built": False}
        check("viewer", False, e.message)
    if library_dir is None:
        return rep
    path = os.path.abspath(os.path.expanduser(library_dir))
    lib: dict = {"path": path}
    rep["library"] = lib
    try:
        store = LibraryStore.open(path)
    except FlowCellarError as e:
        check("library", False, e.message)
        return rep
    lib.update(id=store.id, name=store.name, created=store.created, file_system=_fs_type(path),
               writable=os.access(path, os.W_OK), runs=len(store.run_keys()))
    check("library", True, f"{store.name}: {lib['runs']} runs, layout {store.marker.get('layout_version')}, "
                           f"file system {lib['file_system'] or 'unknown'}")
    if not lib["writable"]:
        check("library writable", False, "this user cannot write to the library folder")
    for w in store.loaded.warnings:
        check("configuration", None, w)
    holder = read_holder(store.layout.lease)
    lib["lease"] = holder.words() if holder else None
    if holder:
        check("lease", None, f"a scan holds the library: {holder.words()}")
    sizes = []
    for k in store.run_keys():
        try:
            sizes.append(os.path.getsize(store.layout.bundle(k)))
        except OSError:
            pass
    if sizes:
        lib["bundle_gzip_bytes"] = {"largest": max(sizes), "median": sorted(sizes)[len(sizes) // 2], "total": sum(sizes)}
        check("bundle sizes", max(sizes) <= (1 << 20), f"largest bundle {max(sizes) / 1e3:.0f} kB compressed, median "
                                                      f"{sorted(sizes)[len(sizes) // 2] / 1e3:.0f} kB (a bundle above 1 MB is unusual)")
    missing_status = [k for k in store.run_keys() if not os.path.isfile(store.layout.status_file(k))]
    if missing_status:
        check("status files", False, f"{len(missing_status)} runs have no status file; scan the library again")
    recorded = (store.read_aux(store.layout.roots_state, {}) or {}).get("roots", {})
    d = store.config["discovery"]
    t0 = time.monotonic()
    outcomes = check_roots(store.config["roots"], recorded={k: v.get("first_level_names") or [] for k, v in recorded.items()},
                           check_seconds=float(d["root_check_seconds"]), stall_seconds=float(d["stall_timeout_seconds"]), walk=False)
    lib["roots"] = []
    for r in store.config["roots"]:
        o = outcomes[r["id"]]
        lib["roots"].append({"id": r["id"], "path": r["path"], "state": o.state, "reason": o.reason,
                             "listing_seconds": round(o.check_seconds, 3), "first_level": len(o.first_level),
                             "file_system": _fs_type(r["path"]) if o.state == "reachable" else None,
                             "last_scanned": recorded.get(r["id"], {}).get("last_scanned")})
        check(f"root {r['id']}", o.state == "reachable",
              f"{r['path']}: " + (f"answers in {o.check_seconds:.2f} s, {len(o.first_level)} folders at the first level, file system "
                                  f"{_fs_type(r['path']) or 'unknown'}" if o.state == "reachable" else str(o.reason)))
    lib["roots_checked_seconds"] = round(time.monotonic() - t0, 2)
    lib["checked"] = now_z()
    return rep
