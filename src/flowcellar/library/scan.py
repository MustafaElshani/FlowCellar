"""``flowcellar scan``: build or update a library from any number of roots.

Order of work, chosen so that a person can open the hub after the first minutes and trust what it says::

    lease -> check and walk every root (child processes, deadline) -> identify what is new or changed
      -> link read-data files to runs
      -> QUICK bundle for every new run; judged in order of run start; one line per run     "hub ready"
      -> full tier 1 passes, newest runs first, inside the time box; the hub is rewritten as results land
      -> read-data passes, newest runs first, inside what is left of the time box
      -> judge once more in order of run start; final index and hub; lease released

What makes the second scan cheap: a run whose MinKNOW files have the same names, sizes and times as last time
is not opened at all; pass results are cached as partials keyed by the fingerprint of their input; links are
kept while the set of files is the same. A scan with nothing new lists folders, opens no input, and writes
only its log, the index and the hub's copy of the index.

What never happens: a root that does not answer is skipped and nothing known under it changes; a run whose
folder disappeared keeps its bundle; a verdict is not touched because a later run arrived.
"""
from __future__ import annotations

import dataclasses
import json
import os
import signal
import socket
import threading
import time
from dataclasses import dataclass, field

from flowcellar import __version__
from flowcellar.config import loader
from flowcellar.contracts import Candidate, LinkedInput, PartialRef, PassOptions, RunContext, RunIdentity
from flowcellar.discovery import REACHABLE, Discovery, discover
from flowcellar.engine import api as engine_api
from flowcellar.library import index as index_mod
from flowcellar.library.cards import card_for, judge_settings
from flowcellar.library.lease import Lease
from flowcellar.library.projects import merge_projects, project_name
from flowcellar.library.store import LibraryStore, is_library
from flowcellar.model import partial as partial_model
from flowcellar.model.status import VERDICT_WORD
from flowcellar.readdata import api as readdata_api
from flowcellar.report import api as report_api
from flowcellar.util import log as log_mod
from flowcellar.util.errors import ExitCode, FlowCellarError, LibraryError, NotBuiltError, UsageError
from flowcellar.util.fingerprint import hash_obj
from flowcellar.util.timeutil import now_z

PASS_ROLES = ("minknow.sequencing_summary", "minknow.pore_scan")      # the tier 1 inputs that need a pass
# pass states of a read-data file with nothing left to read: read whole, a small collected file, or a file that is
# itself cut off (counted as far as it goes; it is read again only when the file changes)
READ_STATES = ("done", "none", "truncated")
INDEX_EVERY_SECONDS = 60.0
#: the kinds of large read-data files ``scan --read`` can choose between, in the order a scan reads them
READ_KINDS = ("bam", "tables", "analyses")

_monotonic = time.monotonic          # tests replace this to make the time box deterministic


@dataclass
class ScanOptions:
    roots: tuple[str, ...] = ()               # folders given on the command line (new ones are remembered)
    name: str | None = None
    tier: str = "all"                         # "1", "2" or "all"
    only: tuple[str, ...] = ()                # run keys or flow cell ids
    since: str | None = None                  # YYYY-MM-DD
    time_box_minutes: float | None = None     # None: the configuration; 0: no limit
    quick: bool = False                       # verdicts from the small files only: no pass over a large file
    read: tuple[str, ...] = ()                # of the read data, read only these kinds in this scan (READ_KINDS); () = all
    wait: bool | None = None                  # None: wait when stderr is not a terminal
    dry_run: bool = False
    retry_failed: bool = False
    force: bool = False
    no_index: bool = False
    no_hub: bool = False
    jobs: int = 0
    threads: int = 0
    config_file: str | None = None
    lease_poll_seconds: float = 5.0
    lease_stale_seconds: float = 180.0


@dataclass
class ScanResult:
    exit_code: int = 0
    library: str = ""
    created: bool = False
    runs: int = 0
    new: list[str] = field(default_factory=list)
    changed: list[str] = field(default_factory=list)       # bundle rewritten
    moved: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    failed: list[str] = field(default_factory=list)
    offline_roots: list[str] = field(default_factory=list)
    waiting_tier1: int = 0
    waiting_read_data: int = 0
    log: str | None = None
    planned: dict = field(default_factory=dict)             # --dry-run: what would be done


@dataclass
class _Run:
    """One run as this scan sees it."""
    key: str
    folders: list = field(default_factory=list)             # [(root id or None, FoundFolder)], the primary first
    ident: RunIdentity | None = None                        # set when identify ran or when it is needed
    identified: bool = False
    state: dict | None = None                               # the stored scan state; None for a new run
    new: bool = False
    need_quick: bool = False
    need_pass: bool = False
    started: str | None = None


def _sig(files) -> str | None:
    return hash_obj([list(f) for f in files], 24) if files else None


def _short(path: str) -> str:
    home = os.path.expanduser("~")
    return "~" + path[len(home):] if home and home != "/" and (path == home or path.startswith(home + os.sep)) else path


def _gb(n: float) -> str:
    return f"{n / 1e9:,.0f} GB" if n >= 10e9 else f"{n / 1e9:,.1f} GB" if n >= 1e8 else f"{n / 1e6:,.0f} MB"


def _span(minutes: float) -> str:
    if minutes < 1:
        return "less than a minute"
    minutes = int(round(minutes))
    if minutes < 60:
        return f"{minutes} min"
    h, m = divmod(minutes, 60)
    return f"{h} h {m:02d} min" if m else f"{h} h"


def _about(minutes: float) -> str:
    return _span(minutes) if minutes < 1 else "about " + _span(minutes)


def _ident_to_state(ident: RunIdentity) -> dict:
    rel = {role: os.path.relpath(p, ident.run_dir) for role, p in ident.files.items()}
    return {"flow_cell_id": ident.flow_cell_id, "run_short": ident.run_short, "protocol_run_id": ident.protocol_run_id,
            "started": ident.started, "files_rel": rel, "sizes": dict(ident.sizes), "newest_mtime": ident.newest_mtime,
            "facts": ident.facts}


def _ident_from_state(key: str, st: dict, run_dir: str) -> RunIdentity:
    i = st.get("identity") or {}
    fc, _, short = key.partition("_")
    return RunIdentity(run_key=key, flow_cell_id=i.get("flow_cell_id") or fc, run_short=i.get("run_short") or short,
                       run_dir=run_dir, kind=st.get("kind", "run"), run_state=st.get("run_state", "finished"),
                       protocol_run_id=i.get("protocol_run_id"), started=i.get("started"),
                       files={r: os.path.join(run_dir, p) for r, p in (i.get("files_rel") or {}).items()},
                       sizes=dict(i.get("sizes") or {}), newest_mtime=float(i.get("newest_mtime") or 0.0),
                       facts=dict(i.get("facts") or {}))


def _linked_to_json(li: LinkedInput) -> dict:
    return dataclasses.asdict(li)


def _linked_from_json(d: dict) -> LinkedInput:
    d = dict(d)
    d["candidate"] = Candidate(**d["candidate"])
    return LinkedInput(**d)


def _strip_times(status_file: dict) -> dict:
    """A status file without the moments it was written: two judgements of the same facts compare equal."""
    d = json.loads(json.dumps(status_file))
    st = d.get("status") or {}
    st.pop("evaluated", None)
    if isinstance(st.get("judged"), dict):
        st["judged"].pop("at", None)
    for p in (d.get("relations") or {}).get("paths") or []:
        p.pop("seen", None)
    return d


class Terminated(KeyboardInterrupt):
    """SIGTERM during a scan. A KeyboardInterrupt, so every way out that Ctrl-C takes is taken."""


class Scan:
    """One scan of one library. ``Scan(dir, options).run()`` returns a ScanResult; nothing here prints except
    through ``util.log`` (stderr) and the scan's own log file."""

    def __init__(self, library_dir: str, options: ScanOptions | None = None):
        self.dir = os.path.abspath(os.path.expanduser(library_dir))
        self.o = options or ScanOptions()
        self.result = ScanResult(library=self.dir)
        self.lib: LibraryStore | None = None
        self.cfg: dict = {}
        self.runs: dict[str, _Run] = {}
        self.cards: dict[str, dict] = {}
        self.states: dict[str, dict] = {}
        self.roots_state: dict = {"roots": {}}
        self.notices: list[dict] = []
        self.linked: list[LinkedInput] = []
        self.inputs_for: dict[str, list[LinkedInput]] = {}
        self.passes: dict[str, dict] = {}                   # fingerprint value -> the state of its read pass
        self.changed: set[str] = set()                      # runs whose bundle, paths or links changed in this scan
        self.hub_dirty: set[str] = set()
        self.failed: dict[str, str] = {}
        self.offline: dict[str, str] = {}                   # root id -> reason
        self.unreadable: list[str] = []                     # folders under a reachable root that could not be listed
        self.not_walked: dict[str, str] = {}                # known run -> why its folder was not looked at
        self.found: Discovery | None = None
        self.single = False                                 # add_single: no root was walked
        self.gone_inputs: set[str] = set()                  # read-data files that are no longer there
        self.t0 = _monotonic()
        self.deadline: float | None = None
        self.started_z = now_z()
        self._log_fh = None
        self._last_index = 0.0
        self._last_line = self.t0
        self._hub_warned = False
        self.t1_bytes_done = 0
        self.t1_seconds = 0.0
        self.t1_waiting_bytes = 0
        self.t1_total = 0
        self.t1_done = 0
        self.rd_bytes_done = 0
        self.rd_seconds = 0.0
        self.counts = {"identify": 0, "quick": 0, "full_pass": 0, "read_pass": 0, "assemble": 0, "judge": 0, "link": 0,
                       "bundles_written": 0, "status_written": 0, "cards_written": 0}
        self._analyses: list[dict] | None = None            # rows of the single-cell and spatial analyses; None = not looked at

    # ================================================================================================ talking
    def say(self, event: str, message: str, *, level: int = 20, **fields) -> None:
        """One line for the person (stderr) and the same event in the scan's log file. While one run is added by
        hand only warnings are shown: the command that asked for it prints its own lines."""
        if not self.single or level >= 30:
            log_mod.event(event, message, level=level, **fields)
        self.note(event, message=message, **fields)

    def note(self, event: str, **fields) -> None:
        """One line in ``logs/scan-*.jsonl`` only: a decision or an error."""
        if self._log_fh is None:
            return
        self._log_fh.write(json.dumps({"t": now_z(), "event": event, **fields}, ensure_ascii=False, default=str) + "\n")
        self._log_fh.flush()

    def _open_log(self) -> None:
        os.makedirs(self.lib.layout.logs, exist_ok=True)
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        host = "".join(ch for ch in socket.gethostname() if ch.isalnum() or ch in "-_.")[:40] or "host"
        path = os.path.join(self.lib.layout.logs, f"scan-{stamp}-{host}-{os.getpid()}.jsonl")
        n = 1
        while os.path.exists(path):                        # a second scan by this process within the same second
            n += 1
            path = os.path.join(self.lib.layout.logs, f"scan-{stamp}-{host}-{os.getpid()}-{n}.jsonl")
        self._log_fh = open(path, "a", encoding="utf-8")
        self.result.log = path

    def _fail(self, key: str, what: str, err: BaseException) -> None:
        msg = f"{what}: {getattr(err, 'message', None) or err}"
        self.failed[key] = msg
        self.say("run.error", f"{key}   could not be read ({msg})", level=30, run_key=key, error=msg,
                 kind=type(err).__name__)

    # ================================================================================================ opening
    def _new_roots(self, existing: list[dict]) -> list[dict]:
        """The command-line folders the library does not know yet, as ``[[roots]]`` entries."""
        have = {r["path"] for r in existing}
        paths = []
        for p in self.o.roots:
            ap = loader.normalise_root_path(p)
            if ap in have or ap in paths:
                continue
            inside = next((r for r in existing if ap.startswith(r["path"].rstrip("/") + "/")), None)
            if inside is not None:
                self.say("root.covered", f"{_short(ap)} lies inside the root '{inside['id']}': it is scanned already",
                         path=ap, root=inside["id"])
                continue
            paths.append(ap)
        ids = loader.root_ids_for(paths, {r["id"] for r in existing})
        return [{"id": i, "path": p, "live": False, "project_depth": 1, "project_pattern": None} for i, p in zip(ids, paths)]

    def _open(self) -> None:
        exists = os.path.isdir(self.dir) and is_library(self.dir)
        if not exists:
            if os.path.isdir(self.dir) and [n for n in os.listdir(self.dir) if n != loader.CONFIG_NAME]:
                raise LibraryError(f"{self.dir}: exists and is not a FlowCellar library",
                                   hint="give -L a new or empty folder, or the folder of an existing library")
            had_config = os.path.isfile(os.path.join(self.dir, loader.CONFIG_NAME))
            draft = loader.load(os.path.join(self.dir, loader.CONFIG_NAME)).config["roots"] if os.path.isdir(self.dir) else []
            new = self._new_roots(draft)
            if not new and not draft and not self.o.config_file:
                raise UsageError(f"{self.dir} is not a library yet and no folder to scan was given",
                                 hint="flowcellar scan /folder/with/runs ... -L " + _short(self.dir))
            if self.o.dry_run:
                loaded = loader.load(self.o.config_file) if self.o.config_file else loader.merge({})
                loaded.config["roots"] = loaded.config["roots"] + [r for r in new if r["path"] not in
                                                                   {x["path"] for x in loaded.config["roots"]}]
                self.cfg = loaded.config
                self.result.created = True
                return
            self.lib = LibraryStore.create(self.dir, name=self.o.name, roots=new, config_file=self.o.config_file)
            if had_config and new:                         # a flowcellar.toml a person wrote first: the roots are added to it
                loader.append_roots(self.lib.layout.config, new, now_z())
                self.lib.reload_config(self.o.config_file)
            self.result.created = True
            n = len(self.lib.config["roots"])
            self.say("library.created", f"created library {_short(self.dir)} ({n} root{'s' if n != 1 else ''}; "
                     "edit flowcellar.toml to add file patterns)", library=self.dir, roots=n)
        else:
            self.lib = LibraryStore.open(self.dir, config_file=self.o.config_file)
            new = self._new_roots(self.lib.config["roots"])
            if new and self.o.dry_run:
                self.lib.config["roots"] = self.lib.config["roots"] + new
            elif new and self.o.config_file:
                self.lib.config["roots"] = self.lib.config["roots"] + new
                self.say("root.not_remembered", "the new folders are scanned now but not remembered: --config names "
                         "another configuration file; add them there", level=30)
            elif new:
                loader.append_roots(self.lib.layout.config, new, now_z())
                self.lib.reload_config(self.o.config_file)
                for r in new:
                    self.say("root.added", f"added root '{r['id']}': {_short(r['path'])}", root=r["id"], path=r["path"])
        self.cfg = self.lib.config
        for w in self.lib.loaded.warnings:
            self.say("config.warning", w, level=30)

    # ================================================================================================ the scan
    def run(self) -> ScanResult:
        self._open()
        if not self.cfg["roots"]:
            raise UsageError("this library has no roots yet", hint="flowcellar scan /folder/with/runs -L " + _short(self.dir))
        if self.o.dry_run:
            return self._dry_run()
        wait = self.o.wait if self.o.wait is not None else not log_mod.is_tty()
        lease = Lease(self.lib.layout.lease, wait=wait, poll_seconds=self.o.lease_poll_seconds,
                      stale_seconds=self.o.lease_stale_seconds, on_wait=self._on_wait)
        try:
            lease.acquire(heartbeat=False)
        except OSError as e:
            raise LibraryError(f"{self.dir}: not writable ({e.strerror or e})") from None
        self._lease = lease
        self.t0 = _monotonic()                              # the time box starts now, not while waiting for the lease
        restore = self._catch_sigterm()
        try:
            self._open_log()
            self.note("scan.start", library=self.dir, version=__version__, options=dataclasses.asdict(self.o),
                      host=socket.gethostname(), pid=os.getpid())
            if lease.taken_over is not None:
                msg = f"The lease of an earlier scan ({lease.taken_over.words()}) was stale and has been taken over."
                self.notices.append({"code": "stale_lease_taken", "level": "info", "message": msg, "root": None,
                                     "since": self.started_z})
                self.say("lease.stale", msg, level=30)
            self._load()
            self._discover()
            lease.start_heartbeat()
            self._reconcile()
            self._link(saved_only=True)                     # what earlier scans linked: no file is opened
            self._refresh_projects()
            self._quick_phase()
            if self.changed or self.hub_dirty or any(r.need_pass for r in self.runs.values()) or not os.path.isfile(self.lib.layout.index):
                self._write_index(scanning=True, force=True)
                self._write_hub(first=True)
            if not self._no_read_data():
                # linking opens every read-data file it has not seen before (minutes on a first scan of a large
                # archive), so it comes after the verdicts from the small files and the first hub
                self._link()
                self._refresh_projects()
            self._tier1_phase()
            self._tier2_phase()
            self._final_judge()
            self._analyses_phase()
            self._finish()
        except KeyboardInterrupt:
            self._settle_after_interrupt()
            raise
        finally:
            restore()
            lease.release()
            if self._log_fh:
                self._log_fh.close()
                self._log_fh = None
        return self.result

    def _on_wait(self, holder) -> None:
        """Said once, on stderr, whoever waits for the lease (a scan, or ``run -L``)."""
        msg = "another scan of this library is running" + (f" ({holder.words()})" if holder else "") + "; waiting for it"
        log_mod.event("lease.wait", msg)
        self.note("lease.wait", message=msg)

    def _catch_sigterm(self):
        """A plain ``kill`` takes the same way out as Ctrl-C: the index, the hub and the lease are settled. Returns
        the function that puts the earlier handler back."""
        if threading.current_thread() is not threading.main_thread():
            return lambda: None

        def on_term(signum, _frame):
            raise Terminated()
        try:
            before = signal.signal(signal.SIGTERM, on_term)
        except (ValueError, OSError):
            return lambda: None

        def restore() -> None:
            try:
                signal.signal(signal.SIGTERM, before)
            except (ValueError, OSError, TypeError):
                pass
        return restore

    def _settle_after_interrupt(self) -> None:
        """Ctrl-C or SIGTERM: what was finished is kept, and the index and the hub's copy of it say that no scan is
        running any more (the page would otherwise keep reloading an index that says "scanning")."""
        self.note("scan.interrupted")
        if self.lib is None or self.found is None:
            return                                          # nothing was looked at yet: the index of the last scan stands
        try:
            self._save_links()
        except Exception:  # noqa: BLE001
            pass
        try:
            self._write_index(scanning=False, force=True)
            self._write_hub()
        except Exception:  # noqa: BLE001  the way out must not fail
            pass

    # ------------------------------------------------------------------------------------------ what is known
    def _load(self) -> None:
        self.cards = self.lib.cards()
        self.states = self.lib.states()
        self.roots_state = self.lib.read_aux(self.lib.layout.roots_state, {"roots": {}}) or {"roots": {}}
        self.roots_state.setdefault("roots", {})
        links = self.lib.read_aux(self.lib.layout.links, {}) or {}
        self.passes = dict(links.get("passes") or {})
        self._links_saved = links
        self.notices.extend(n for n in self.lib.loaded.notices)
        box = self.o.time_box_minutes if self.o.time_box_minutes is not None else self.cfg["scan"]["time_box_minutes"]
        self.time_box = float(box or 0)
        self.deadline = self.t0 + self.time_box * 60 if self.time_box > 0 else None

    def _root(self, root_id: str | None) -> dict | None:
        return next((r for r in self.cfg["roots"] if r["id"] == root_id), None)

    def _abs(self, p: dict) -> str:
        root = self._root(p.get("root"))
        if root is None:
            return p.get("path") or p.get("relpath") or ""
        return root["path"] if p.get("relpath") in ("", ".", None) else os.path.join(root["path"], p["relpath"])

    # ------------------------------------------------------------------------------------------ roots and walk
    def _discover(self) -> None:
        recorded = {rid: s.get("first_level_names") or [] for rid, s in self.roots_state["roots"].items()}
        sniffed = self.lib.read_aux(self.lib.layout.sniffed, {}) or {}
        before = json.dumps(sniffed, sort_keys=True)
        self.found = discover(self.cfg, recorded=recorded, sniff_cache=sniffed, tick=self._lease.touch_if_due)
        if json.dumps(sniffed, sort_keys=True) != before:
            self.lib.write_aux(self.lib.layout.sniffed, sniffed)
        now = now_z()
        for r in self.cfg["roots"]:
            o = self.found.roots[r["id"]]
            st = self.roots_state["roots"].setdefault(r["id"], {"first_scanned": None})
            st["path"] = r["path"]
            if o.state == REACHABLE:
                st.update(state=REACHABLE, reason=None, offline_since=None, last_scanned=now,
                          first_level_names=o.first_level, listings=o.walk.counts.get("listings") if o.walk else None,
                          seconds=round(o.seconds, 2))
                st["first_scanned"] = st.get("first_scanned") or now
                self.note("root.walked", root=r["id"], path=r["path"], seconds=round(o.seconds, 3),
                          counts=o.walk.counts if o.walk else None, runs=len(o.walk.runs) if o.walk else 0)
            else:
                st.update(state=o.state, reason=o.reason, offline_since=st.get("offline_since") or now)
                self.offline[r["id"]] = o.reason or o.state
                word = "offline" if o.state == "offline" else "skipped"
                msg = (f"Root '{r['id']}' ({r['path']}) is {word} since {st['offline_since']}: {o.reason}. "
                       "Nothing known under it was changed.")
                self.notices.append({"code": "root_offline", "level": "warn", "message": msg, "root": r["id"],
                                     "since": st["offline_since"]})
                self.say("root.offline", f"root '{r['id']}' skipped: {o.reason} ({_short(r['path'])}); nothing under it changes",
                         level=30, root=r["id"], state=o.state, reason=o.reason)
        for rid, sub in self.found.over_budget:
            budget = self.cfg["discovery"]["subtree_listing_budget"]
            msg = (f"{sub}: more than {budget:,} folders without a run folder or an anchor; add it to prune or raise "
                   "subtree_listing_budget in flowcellar.toml. Runs below it were not looked for in this scan.")
            self.notices.append({"code": "subtree_over_budget", "level": "warn", "message": msg, "root": rid, "since": now})
            self.say("walk.over_budget", msg, level=30, root=rid, path=sub)
        unread = list(getattr(self.found, "unreadable", ()) or ())
        if unread:
            self.unreadable = [p for _, p in unread]
            n = len(unread)
            shown = ", ".join(p for _, p in unread[:3]) + (f" and {n - 3} more" if n > 3 else "")
            rids = {rid for rid, _ in unread}
            msg = (f"{n} folder{'s' if n != 1 else ''} could not be read (no permission, or a read error): {shown}. "
                   "Runs and files below were not looked for; what the library already holds from there is kept as it was.")
            self.notices.append({"code": "folders_unreadable", "level": "warn", "message": msg,
                                 "root": next(iter(rids)) if len(rids) == 1 else None, "since": now})
            self.say("walk.unreadable", msg, level=30, paths=self.unreadable[:50], n=n)
        for p in self.found.problems:
            self.say("config.pattern", f"input pattern not used: {p}", level=30)
        self.note("walk.done", seconds=round(self.found.seconds, 3), **self.found.counts())

    def _unseen_is_unknown(self, path: str, root_id: str | None) -> bool:
        """True when a path that was not found in this scan may still be there (its root or subtree was not walked,
        or a folder above it could not be read)."""
        if self.single or (root_id is not None and root_id in self.offline):
            return True
        if any(path == sub or path.startswith(sub.rstrip("/") + "/") for sub in self.unreadable):
            return True
        return any(path == sub or path.startswith(sub.rstrip("/") + "/") for _, sub in self.found.over_budget)

    def _why_not_walked(self, path: str, root_id: str | None) -> str | None:
        """A known run folder that the walk did not report. None: it is really gone. Otherwise the reason in words:
        the folder is still there (the walk does not go there), or a folder above it cannot be read."""
        root = self._root(root_id)
        try:
            os.stat(path)
        except FileNotFoundError:
            # gone, unless the nearest folder above that still exists cannot be listed
            d = os.path.dirname(path)
            top = root["path"] if root else os.path.dirname(d)
            while d and d != os.path.dirname(d):
                if os.path.isdir(d):
                    return None if os.access(d, os.R_OK | os.X_OK) else f"{d} cannot be read"
                if d == top:
                    break
                d = os.path.dirname(d)
            return None
        except PermissionError:
            return "a folder above it cannot be read"
        except OSError as e:
            return f"it could not be looked at ({e.strerror or e})"
        if not os.path.isdir(path):
            return None
        if root is not None:
            from flowcellar.model import runfolder
            rel = os.path.relpath(path, root["path"])
            for name in rel.split(os.sep)[:-1]:
                if runfolder.is_pruned(name, tuple(self.cfg["discovery"]["prune"])):
                    return f"it lies under a folder named '{name}', which a scan does not enter"
            return "a scan did not reach it (it lies deeper than a scan looks, or beside another run's files)"
        return "it lies under no folder this library scans"

    # ------------------------------------------------------------------------------------------ new, moved, gone
    def _identify(self, path: str) -> RunIdentity | None:
        self.counts["identify"] += 1
        self.note("identify", path=path)
        try:
            ident = engine_api.identify(path)
            if ident.run_state in ("running", "interrupted") and ident.newest_mtime:
                # running or interrupted is a matter of the configured stale_after_hours
                fresh = time.time() - ident.newest_mtime < float(self.cfg["discovery"]["stale_after_hours"]) * 3600
                ident = dataclasses.replace(ident, run_state="running" if fresh else "interrupted")
            return ident
        except FlowCellarError as e:
            if isinstance(e, NotBuiltError):
                raise
            self.failed[path] = e.message
            self.say("run.error", f"{_short(path)}: not read as a run ({e.message})", level=30, path=path, error=e.message)
            return None

    def _in_scope(self, key: str, started: str | None) -> bool:
        if self.o.only:
            fc = key.split("_", 1)[0]
            if key not in self.o.only and fc not in self.o.only:
                return False
        if self.o.since and (not started or started[:10] < self.o.since):
            return False
        return True

    def _reconcile(self) -> None:
        tombs = self.lib.tombstones()
        min_age = float(self.cfg["discovery"]["min_age_seconds"])
        by_path: dict[str, str] = {}
        by_sig: dict[str, str] = {}
        for k, st in self.states.items():
            for p in st.get("paths") or []:
                by_path[self._abs(p)] = k
                if p.get("sig") and st.get("kind") == "run":
                    by_sig[p["sig"]] = k
        seen: dict[str, _Run] = {}
        now = time.time()
        for root_id, folder in self.found.runs:
            sig = _sig(folder.files)
            key = None
            known = by_path.get(folder.path)
            if not self.o.force and known:
                st = self.states[known]
                same = any(self._abs(p) == folder.path and p.get("sig") == sig for p in st.get("paths") or [])
                if same and st.get("generator") == __version__ and st.get("run_state") != "running" \
                        and (sig or folder.kind == st.get("kind")):
                    key = known
            if key is None and not self.o.force and sig and by_sig.get(sig) and folder.kind == "run":
                st = self.states[by_sig[sig]]
                if st.get("generator") == __version__ and st.get("run_state") != "running":
                    key = by_sig[sig]                      # the same files at another place: moved or mirrored
            ident = None
            if key is None:
                ident = self._identify(folder.path)
                if ident is None:
                    continue
                key = ident.run_key
                root = self._root(root_id)
                if ident.newest_mtime and now - ident.newest_mtime < min_age and not (root or {}).get("live") \
                        and key not in self.states:
                    self.say("run.settling", f"{key}   files changed in the last {int(min_age)} s: left for the next scan",
                             run_key=key, path=folder.path)
                    continue
            if key in tombs:
                self.note("run.tombstone", run_key=key, path=folder.path)
                continue
            r = seen.setdefault(key, _Run(key=key, state=self.states.get(key)))
            r.folders.append((root_id, folder))
            if ident is not None and (r.ident is None or (r.ident.kind != "run" and ident.kind == "run")):
                r.ident, r.identified = ident, True

        root_order = [r["id"] for r in self.cfg["roots"]]
        for key, r in seen.items():
            st = r.state
            prev_primary = self._abs(st["paths"][0]) if st and st.get("paths") else None
            # the folder that is read: one with reports, the one read before, else the most complete copy, then
            # the order of the roots in the configuration
            r.folders.sort(key=lambda rf: (rf[1].kind != "run", rf[1].path != prev_primary, -len(rf[1].files),
                                           -sum(size for _, size, _ in rf[1].files),
                                           root_order.index(rf[0]) if rf[0] in root_order else 99, rf[1].relpath))
            if r.ident is not None and r.ident.run_dir != r.folders[0][1].path:
                again = self._identify(r.folders[0][1].path)           # identify ran on a mirror: use the primary
                r.ident = again or r.ident
            r.started = (r.ident.started if r.ident else None) or ((st or {}).get("identity") or {}).get("started")
            if st is None:
                r.new = True
                if not self._in_scope(key, r.started):
                    self.note("run.out_of_scope", run_key=key)
                    continue
                r.need_quick = True
                r.need_pass = r.ident.kind == "run"
            else:
                touched = r.identified or self.o.force
                t1 = st.get("t1", "none")
                if touched:
                    # its files changed (or FlowCellar did): a run that has full detail keeps its bundle until the
                    # reads are counted again; anything else gets a fresh quick bundle now
                    kind = r.ident.kind if r.ident else st.get("kind")
                    state_now = r.ident.run_state if r.ident else st.get("run_state")
                    r.need_quick = not (t1 == "done" and kind == "run" and state_now == st.get("run_state")
                                        and state_now != "running")
                    r.need_pass = kind == "run"
                    if kind == "run" and t1 == "done" and self._in_scope(key, r.started):
                        st["t1"] = "waiting"
                        self.lib.write_state(key, st)
                else:
                    r.need_pass = st.get("kind") == "run" and (t1 in ("waiting", "partial")
                                                              or (t1 == "error" and self.o.retry_failed))
                if not self._in_scope(key, r.started):
                    r.need_quick = r.need_pass = False
                self._update_paths(r)
            self.runs[key] = r

        for want in self.o.only:
            if not any(k == want or k.split("_", 1)[0] == want for k in set(self.runs) | set(self.states)):
                self.say("scan.only", f"--only {want}: no run of this library or under its roots has this run key or flow cell id",
                         level=30)
        # a run that has a bundle but was never judged (an earlier scan failed there) is judged in this one
        for key in self.states:
            if key not in self.cards and key not in tombs and self.lib.has_run(key):
                self.changed.add(key)

        # runs the library holds that were not seen in this scan
        for key, st in self.states.items():
            if key in seen or key in tombs:
                continue
            paths = st.get("paths") or []
            if not paths:
                continue
            # a folder under no root of the configuration (added by hand, or its root was taken out of
            # flowcellar.toml) is simply looked for where it was
            keep = []
            for p in paths:
                ap = self._abs(p)
                if self._unseen_is_unknown(ap, p.get("root")) or (self._root(p.get("root")) is None and os.path.isdir(ap)):
                    keep.append(p)
                    continue
                why = self._why_not_walked(ap, p.get("root"))     # one stat: the folder may be there and only not walked
                if why is not None:
                    keep.append(p)
                    self.not_walked.setdefault(key, why)
            if len(keep) == len(paths):
                continue                                   # everything it had is where a scan could not or does not look
            if keep:
                st["paths"] = keep
                self._paths_changed(key, st)
            elif st.get("source_state") != "missing":
                st["source_state"], st["missing_since"] = "missing", now_z()
                self.result.missing.append(key)
                self.say("run.missing", f"{key}   its folder is no longer there; the run stays in the library",
                         run_key=key, was=[self._abs(p) for p in paths])
                self._paths_changed(key, st)
        if self.not_walked:
            keys = sorted(self.not_walked)
            n = len(keys)
            why = self.not_walked[keys[0]]
            same = all(w == why for w in self.not_walked.values())
            msg = (f"{n} run{'s' if n != 1 else ''} of this library {'were' if n != 1 else 'was'} not looked at by the last scan "
                   f"({', '.join(keys[:4])}{', ...' if n > 4 else ''}): "
                   + (why if same else f"{why}; the others for similar reasons")
                   + f". {'They stay' if n != 1 else 'It stays'} in the library as {'they were' if n != 1 else 'it was'}.")
            self.notices.append({"code": "runs_not_walked", "level": "info", "message": msg, "root": None, "since": now_z(),
                                 "runs": keys[:50]})
            self.say("run.not_walked", msg, runs=keys[:50], n=n)

    def _update_paths(self, r: _Run) -> None:
        """A known run seen again: note a new place, a second place, or that it is back."""
        st = r.state
        now = now_z()
        old = {(p.get("root"), p.get("relpath")): p for p in st.get("paths") or []}
        new = []
        for root_id, f in r.folders:
            prev = old.get((root_id, f.relpath), {})
            new.append({"root": root_id, "relpath": f.relpath, "path": f.path, "sig": _sig(f.files),
                        "first_seen": prev.get("first_seen") or now, "last_seen": now})
        kept = [p for p in st.get("paths") or [] if (p.get("root"), p.get("relpath")) not in {(x["root"], x["relpath"]) for x in new}
                and self._unseen_is_unknown(self._abs(p), p.get("root"))]
        where_old = [(p.get("root"), p.get("relpath")) for p in st.get("paths") or []]
        where_new = [(p["root"], p["relpath"]) for p in new + kept]
        sig_changed = any(old.get((p["root"], p["relpath"]), {}).get("sig") != p["sig"] for p in new)
        was_missing = st.get("source_state") == "missing"
        if where_old != where_new or was_missing or sig_changed:
            moved = set(where_old) != set(where_new)
            st["paths"] = new + kept
            st["source_state"], st["missing_since"] = "ok", None
            if moved:
                self.result.moved.append(r.key)
                more = set(where_new) > set(where_old)
                self.say("run.moved", f"{r.key}   " + ("also found at" if more else "found at a new place:")
                         + f" {_short(new[-1]['path'] if more else new[0]['path'])}; nothing is computed again",
                         run_key=r.key, paths=[p["path"] for p in new])
            if moved or was_missing:
                self._paths_changed(r.key, st)             # the status file's relations and the card follow
            else:
                self.lib.write_state(r.key, st)

    def _paths_changed(self, key: str, st: dict) -> None:
        """Paths are library knowledge: update the state, the relations of the status file and the card. No bundle
        is touched and nothing is judged."""
        self.lib.write_state(key, st)
        sf = self.lib.read_status_file(key)
        card = self.cards.get(key)
        if sf is None or card is None:
            return
        rel = sf.setdefault("relations", {})
        rel["paths"] = self._relation_paths(st)
        rel["projects"] = self._projects(key, st)
        self.lib.write_status_file(key, sf)
        self.counts["status_written"] += 1
        cols = card["columns"]
        roots = []
        for p in rel["paths"]:
            if p["root"] not in roots:
                roots.append(p["root"])
        first = (st.get("paths") or [{}])[0]
        cols.update(roots=roots, project=rel["projects"][0]["name"] if rel["projects"] else None,
                    project_ids=[p["id"] for p in rel["projects"]], source_state=st.get("source_state", "ok"),
                    folder=os.path.basename((first.get("relpath") or first.get("path") or "").rstrip("/")) or cols.get("folder"))
        card["projects"] = rel["projects"]
        self.lib.write_card(key, card)
        self.counts["cards_written"] += 1
        self.hub_dirty.add(key)

    def _refresh_projects(self) -> None:
        """Projects follow the configuration: when ``project_depth`` or a link changed, the runs it concerns get
        their new project. Library knowledge only: no bundle is touched and nothing is judged."""
        for key, st in self.states.items():
            card = self.cards.get(key)
            if card is not None and (card.get("projects") or []) != self._projects(key, st):
                self.note("run.project", run_key=key)
                self._paths_changed(key, st)

    def _relation_paths(self, st: dict) -> list[dict]:
        return [{"root": p["root"], "relpath": p.get("relpath") or ".", "seen": p.get("last_seen") or now_z()}
                for p in st.get("paths") or [] if p.get("root")]

    def _projects(self, key: str, st: dict) -> list[dict]:
        found = [(project_name(self._root(p.get("root")), p.get("relpath")), p.get("root")) for p in st.get("paths") or []]
        skip = loader.not_project_names(self.cfg)           # a pipeline's own folder at project depth names no project
        for li in self.inputs_for.get(key, []):
            if li.counted_for.get(key):
                c = li.candidate
                found.append((project_name(self._root(c.root_id), c.relpath, skip), c.root_id))
        return merge_projects(found, [r["id"] for r in self.cfg["roots"]])

    # ------------------------------------------------------------------------------------------ read-data links
    def _all_identities(self) -> list[RunIdentity]:
        out = []
        for key in sorted(set(self.states) | set(self.runs)):
            ident = self._ident(key)
            if ident is not None:
                out.append(ident)
        return out

    def _ident(self, key: str) -> RunIdentity | None:
        r = self.runs.get(key)
        if r is not None and r.ident is not None:
            return r.ident
        st = self.states.get(key)
        if st is None:
            return None
        run_dir = r.folders[0][1].path if r is not None and r.folders else (self._abs(st["paths"][0]) if st.get("paths") else "")
        ident = _ident_from_state(key, st, run_dir)
        if r is not None:
            r.ident = ident
        return ident

    def _link(self, saved_only: bool = False) -> None:
        """Attach read-data files to runs. ``saved_only``: take the links of the last scan as they are."""
        saved_only = saved_only or self._no_read_data()
        if saved_only:
            self.linked = [_linked_from_json(d) for d in self._links_saved.get("linked") or []]
        else:
            cands = list(self.found.candidates)
            here = {c.path for c in cands}
            # kept as they were: files under a root that could not be looked at, and files that are gone but were
            # read before (an input that disappeared keeps its place in the bundle)
            kept = []
            for li in (_linked_from_json(d) for d in self._links_saved.get("linked") or []):
                if li.candidate.path in here:
                    continue
                if self._unseen_is_unknown(li.candidate.path, li.candidate.root_id):
                    kept.append(li)
                elif self._pass_state(li) in READ_STATES:
                    kept.append(li)
                    self.gone_inputs.add(li.candidate.path)
            keys = sorted((set(self.states) | set(self.runs)) - self.lib.tombstones())
            sig = hash_obj({"files": sorted([c.path, c.size, c.mtime_ns, c.role, c.sample, c.variant] for c in cands),
                            "runs": keys, "config": loader.shaping_hash(self.cfg), "version": __version__,
                            "passes": readdata_api.pass_signature(), "linking": readdata_api.link_version()}, 32)
            if sig == self._links_saved.get("signature") and not self.o.force:
                self.linked = [_linked_from_json(d) for d in self._links_saved.get("linked") or []]
                self.note("link.kept", files=len(cands))
            else:
                self.counts["link"] += 1
                self.note("link", files=len(cands), runs=len(keys))
                if cands:
                    self.say("scan.linking", f"linking {len(cands):,} read-data file{'s' if len(cands) != 1 else ''} to runs "
                             "(a file that was never seen before is opened once to see which run it belongs to)", files=len(cands))
                fresh = readdata_api.link(cands, self._all_identities(), config=self.cfg, probes_dir=self.lib.layout.probes) if cands else []
                have = {li.candidate.path for li in fresh}
                self.linked = list(fresh) + [li for li in kept if li.candidate.path not in have]
                self._links_signature = sig
                self._save_links(sig)
        if self.gone_inputs:
            self.note("inputs.gone", paths=sorted(self.gone_inputs)[:50], n=len(self.gone_inputs))
        self.inputs_for = {}
        for li in self.linked:
            for k in li.run_keys:
                self.inputs_for.setdefault(k, []).append(li)
        for rid in [r["id"] for r in self.cfg["roots"]]:
            if rid in self.offline:
                continue
            mine = [li for li in self.linked if li.candidate.root_id == rid]
            before = self.roots_state["roots"].setdefault(rid, {}).get("read_data_files") or {}
            # without a new look at the files the count of the last full look stands
            found = before.get("found", len(mine)) if saved_only else sum(1 for c in self.found.candidates if c.root_id == rid)
            self.roots_state["roots"][rid]["read_data_files"] = {
                "found": found, "read": sum(1 for li in mine if self._pass_state(li) in READ_STATES)}
        if saved_only and not self._no_read_data():
            return                                          # the real link follows in this scan and brings the notice
        both = [li for li in self.linked if li.not_counted_reason and li.run_keys]
        if both:
            names = ", ".join(sorted({os.path.basename(li.candidate.path) for li in both})[:4])
            runs = sorted({k for li in both for k in li.run_keys if self._shows_read_data(k)})
            self.notices.append({"code": "two_candidate_sources", "level": "warn", "root": None, "since": None, "runs": runs[:50],
                                 "message": f"{len(both)} files of {len(runs)} run{'s' if len(runs) != 1 else ''} are not counted because "
                                            "two candidate sources describe the same reads and no rule chooses between them "
                                            f"({names}{', ...' if len(both) > 4 else ''}). To choose, put primary = true on the "
                                            "[[inputs]] entry of the files to count, in the library's flowcellar.toml, and scan again."})

    def _save_links(self, sig: str | None = None) -> None:
        sig = sig or getattr(self, "_links_signature", None) or self._links_saved.get("signature")
        new = {"signature": sig, "generator": __version__, "linked": [_linked_to_json(li) for li in self.linked],
               "passes": self.passes}
        if json.dumps(new, sort_keys=True, default=str) != json.dumps(self._links_saved, sort_keys=True, default=str):
            self.lib.write_aux(self.lib.layout.links, new)
        self._links_saved = json.loads(json.dumps(new, default=str))

    def _pass_state(self, li: LinkedInput) -> str:
        """Where the pass over one input stands. A small collected file (a statistics table, a sample sheet) has
        no partial id: it needs no pass and is read when the bundle is put together, so it never waits."""
        if not li.partial_id:
            return "none"
        if self._read_by_earlier_pass(li):
            return "waiting"
        return (self.passes.get(li.fingerprint.get("value", "")) or {}).get("pass", "waiting")

    def _read_by_earlier_pass(self, li: LinkedInput) -> bool:
        """True when what the library holds of this file was read by an earlier version of the pass (its partial has
        another id than the pass of today would write): the file waits to be read again, and comes first."""
        p = self.passes.get(li.fingerprint.get("value", "")) or {}
        return bool(li.partial_id and p.get("partial_id") and p["partial_id"] != li.partial_id)

    def _reads_kind(self, kind: str | None) -> bool:
        """``scan --read``: whether this scan reads large read-data files of this kind (all kinds without it)."""
        return not self.o.read or kind in self.o.read

    def _files_to_read(self) -> list[LinkedInput]:
        """The linked inputs that need a pass (BAMs, read summaries). Small tables are read when a bundle is put
        together and are not counted as files to read."""
        return [li for li in self.linked if li.partial_id and self._wants_pass(li)]

    def _shows_read_data(self, key: str) -> bool:
        """True for a run of this library whose bundle can hold read data (a run with MinKNOW reports)."""
        st = self.states.get(key)
        if st is not None:
            return st.get("kind", "run") == "run"
        r = self.runs.get(key)
        return bool(r is not None and r.ident is not None and r.ident.kind == "run")

    def _wants_pass(self, li: LinkedInput) -> bool:
        """An input is read only when it is the source of record of something for a run this library holds. This
        one rule gives every count of "files to read" and "files waiting" (the scan's lines, the progress on Home)
        and decides what the scan reads: a file whose only runs are not in the library is never opened, and never
        counted as waiting."""
        if not li.run_keys or li.not_counted_reason:
            return False
        if not li.counted_for:
            return any(self._shows_read_data(k) for k in li.run_keys)
        return any(needs and self._shows_read_data(k) for k, needs in li.counted_for.items())

    def _rd_hash(self, key: str) -> str | None:
        if not self.inputs_for.get(key):
            return None
        view = [[li.fingerprint.get("value"), li.candidate.role, li.sample, li.variant, sorted(li.counted_for.get(key, [])),
                 li.not_counted_reason, self._pass_state(li) if self._wants_pass(li) else "none"]
                for li in self.inputs_for.get(key, [])]
        # with the versions of the section builders: after an upgrade that puts the same files together differently
        # the bundle is assembled again (from the partials: nothing is read)
        return hash_obj([readdata_api.assembly_signature(), sorted(view, key=json.dumps)], 24)

    # ------------------------------------------------------------------------------------------ one run
    def _context(self, key: str) -> RunContext:
        r = self.runs.get(key)
        st = self.states.get(key) or {}
        if r is not None and r.folders:
            root_id, rel = r.folders[0][0], r.folders[0][1].relpath
        else:
            p = (st.get("paths") or [{}])[0]
            root_id, rel = p.get("root"), p.get("relpath")
        return RunContext(library_id=self.lib.id, library_name=self.lib.name, root_id=root_id, relpath=rel,
                          stale_after_hours=float(self.cfg["discovery"]["stale_after_hours"]))

    def _pass_options(self, key: str) -> PassOptions:
        r = self.runs.get(key)
        root = self._root(r.folders[0][0]) if r is not None and r.folders else None
        live = bool((root or {}).get("live"))
        jobs = self.o.jobs or int(self.cfg["scan"]["jobs"])
        if live:
            jobs = min(jobs or int(self.cfg["scan"]["live_workers"]), int(self.cfg["scan"]["live_workers"]))
        return PassOptions(jobs=jobs, threads=self.o.threads or int(self.cfg["scan"]["threads"]), deadline=None,
                           gentle=live, stall_timeout_s=float(self.cfg["discovery"]["stall_timeout_seconds"]))

    def _new_state(self, r: _Run) -> dict:
        now = now_z()
        return {"run_key": r.key, "kind": r.ident.kind, "run_state": r.ident.run_state, "generator": __version__,
                "identity": _ident_to_state(r.ident),
                "paths": [{"root": rid, "relpath": f.relpath, "path": f.path, "sig": _sig(f.files), "first_seen": now,
                           "last_seen": now} for rid, f in r.folders],
                "source_state": "ok", "missing_since": None, "inputs": [], "partials_planned": [], "partials_merged": [],
                "t1": "waiting" if r.ident.kind == "run" else "none", "rd_hash": None, "content_sha256": None,
                "bundle_bytes": None, "bundle_gzip_bytes": None, "assembled": None, "judged": None, "errors": []}

    def _refresh_state(self, r: _Run) -> dict:
        """The state of a run after identify ran on it in this scan."""
        st = self.states.get(r.key)
        if st is None:
            st = self._new_state(r)
        else:
            st.update(kind=r.ident.kind, run_state=r.ident.run_state, generator=__version__, identity=_ident_to_state(r.ident))
            if st.get("t1") == "none" and r.ident.kind == "run":
                st["t1"] = "waiting"
        self.states[r.key] = st
        return st

    def _tier1_inputs(self, ident: RunIdentity, refs: list[PartialRef], st: dict, reset: bool = False) -> list[dict]:
        """The MinKNOW files of a run as state entries. ``reset``: the passes are to be done again."""
        by_role = {ref.role: ref for ref in refs}
        old = {} if reset else {i.get("role"): i for i in st.get("inputs") or [] if i.get("tier") == 1}
        r = self.runs.get(ident.run_key)
        times = {name: m for name, _, m in (r.folders[0][1].files if r is not None and r.folders else ())}
        out = []
        for role, path in sorted(ident.files.items()):
            ref = by_role.get(role)
            prev = old.get(role, {})
            needs = role in PASS_ROLES or ref is not None
            entry = {"tier": 1, "role": role, "path": path, "size": ident.sizes.get(role),
                     "mtime_ns": times.get(os.path.relpath(path, ident.run_dir)),
                     "fingerprint": None,
                     "partial_id": (ref.partial_id or prev.get("partial_id")) if ref else prev.get("partial_id") if needs else None,
                     "pass": ("done" if ref.complete else "error" if ref.state == "error" else "partial") if ref
                     else (prev.get("pass", "waiting") if needs else "none"),
                     "attempts": int(prev.get("attempts") or 0) + (1 if ref and not ref.cached else 0),
                     "error": None if not ref or ref.complete else "could not be opened" if not ref.partial_id else ref.state}
            out.append(entry)
        return out

    def _read_inputs(self, key: str) -> list[dict]:
        out = []
        for li in self.inputs_for.get(key, []):
            c = li.candidate
            p = self.passes.get(li.fingerprint.get("value", "")) or {}
            if self._read_by_earlier_pass(li):
                p = {}                                      # what an earlier version of the pass read does not count
            out.append({"tier": 2, "role": c.role, "path": c.path, "size": c.size, "mtime_ns": c.mtime_ns,
                        "state": "missing" if c.path in self.gone_inputs else "ok",
                        "fingerprint": li.fingerprint, "partial_id": p.get("partial_id") or li.partial_id,
                        "pass": p.get("pass", "waiting") if self._wants_pass(li) else "none",
                        "attempts": int(p.get("attempts") or 0), "error": p.get("error")})
        return out

    def _t1_refs(self, st: dict) -> list[PartialRef]:
        """The tier 1 partials a run has, from its state: nothing is read but the small partial files."""
        out = []
        for i in st.get("inputs") or []:
            pid = i.get("partial_id")
            if i.get("tier") != 1 or not pid:
                continue
            path = partial_model.find_partial(self.lib.layout.partials, pid)
            if path is None:
                continue
            try:
                p = partial_model.read_partial(path)
            except Exception:  # noqa: BLE001  a damaged cache file is as good as none
                continue
            out.append(PartialRef(partial_id=pid, path=path, role=i["role"], input_path=i["path"], state=p.get("state", "ok"),
                                  complete=partial_model.is_complete(p), cached=True))
        return out

    def _store(self, key: str, bundle: dict, st: dict) -> bool:
        """Write a measured bundle when its content changed. True when the file was written."""
        from flowcellar.model import bundle as bundle_model
        sha = bundle_model.content_sha256(bundle)
        st["partials_merged"] = sorted((bundle.get("provenance") or {}).get("partials") or [])
        st["assembled"] = now_z()
        if sha == st.get("content_sha256") and self.lib.has_run(key):
            return False
        info = self.lib.write_measured(key, bundle)
        self.counts["bundles_written"] += 1
        st.update(content_sha256=info["content_sha256"], bundle_bytes=info["bytes"], bundle_gzip_bytes=info["gzip_bytes"])
        if info["warn"]:
            self.say("bundle.large", f"{key}   its bundle is {info['gzip_bytes'] / 1e6:.1f} MB compressed (more than 1 MB)",
                     level=30, run_key=key)
        self.changed.add(key)
        self.hub_dirty.add(key)
        if key not in self.result.changed:
            self.result.changed.append(key)
        self.note("bundle.written", run_key=key, bytes=info["bytes"], gzip_bytes=info["gzip_bytes"], sha256=info["content_sha256"])
        return True

    def _quick(self, r: _Run) -> bool:
        key = r.key
        try:
            st = self._refresh_state(r)
            self.counts["quick"] += 1
            t = _monotonic()
            bundle = engine_api.quick_bundle(r.ident, self._context(key))
            st["inputs"] = self._tier1_inputs(r.ident, [], st, reset=True) + self._read_inputs(key)
            if st.get("t1") in ("done", "partial", "truncated", "error") and r.ident.kind == "run":
                st["t1"] = "waiting"                     # quick content replaced what was there: count again
            self._store(key, bundle, st)
            self.lib.write_state(key, st)
            self.note("quick", run_key=key, seconds=round(_monotonic() - t, 3))
            self.changed.add(key)
            return True
        except NotBuiltError:
            raise
        except Exception as e:  # noqa: BLE001  one run must not stop the scan
            self._fail(key, "quick report", e)
            return False

    def _assemble(self, key: str) -> bool:
        """Put a run's bundle together again from its partials and its read data. Reads no large input."""
        st = self.states[key]
        ident = self._ident(key)
        try:
            self.counts["assemble"] += 1
            inputs = self.inputs_for.get(key, [])
            rd = None
            if inputs:
                tier1 = self.lib.read_measured(key) if self.lib.has_run(key) else None
                rd = readdata_api.read_data(ident, inputs, self.lib.layout.partials, tier1=tier1, config=self.cfg)
            bundle = engine_api.assemble(ident, self._t1_refs(st), rd, self._context(key))
            st["rd_hash"] = self._rd_hash(key)
            st["inputs"] = [i for i in st.get("inputs") or [] if i.get("tier") == 1] + self._read_inputs(key)
            st["partials_planned"] = sorted({i["partial_id"] for i in st["inputs"] if i.get("partial_id") and i.get("pass") != "none"})
            self._store(key, bundle, st)
            self.lib.write_state(key, st)
            return True
        except NotBuiltError:
            raise
        except Exception as e:  # noqa: BLE001
            self._fail(key, "assembling the bundle", e)
            return False

    def _judge(self, key: str) -> None:
        """Judge one run against the cards of the others and write its status file and card when they changed."""
        st = self.states[key]
        try:
            bundle = self.lib.read_measured(key)
            history = [c for k, c in self.cards.items() if k != key]
            self.counts["judge"] += 1
            sf = engine_api.judge(bundle, history, in_library=True, **judge_settings(self.cfg))
            sf["content_sha256"] = bundle.get("content_sha256")
            rel = sf.setdefault("relations", {})
            rel["paths"] = self._relation_paths(st)
            rel["projects"] = self._projects(key, st)
            old = self.lib.read_status_file(key)
            if old is None or _strip_times(old) != _strip_times(sf):
                self.lib.write_status_file(key, sf)
                self.counts["status_written"] += 1
                st["judged"] = now_z()
                self.lib.write_state(key, st)
                self.hub_dirty.add(key)
                self.note("judged", run_key=key, verdict=(sf.get("status") or {}).get("verdict"))
            else:
                sf = old
            card = card_for(bundle, sf, extra={"bundle_bytes": st.get("bundle_bytes"), "bundle_gzip_bytes": st.get("bundle_gzip_bytes"),
                                               "source_state": st.get("source_state", "ok")})
            if self.cards.get(key) != card:
                self.lib.write_card(key, card)
                self.counts["cards_written"] += 1
                self.hub_dirty.add(key)
            self.cards[key] = card
        except NotBuiltError:
            raise
        except Exception as e:  # noqa: BLE001
            self._fail(key, "judging", e)

    def _by_start(self, keys) -> list[str]:
        def start(k: str) -> str:
            r = self.runs.get(k)
            return (r.started if r else None) or ((self.states.get(k) or {}).get("identity") or {}).get("started") or ""
        return sorted(keys, key=lambda k: (start(k), k))

    # ------------------------------------------------------------------------------------------ the phases
    def _line(self, key: str) -> str:
        c = (self.cards.get(key) or {}).get("columns") or {}
        word = VERDICT_WORD.get(c.get("verdict"), "Not assessed") + (" so far" if c.get("provisional") else "")
        reasons = c.get("reasons") or []
        if reasons:
            text = reasons[0]
        elif c.get("provisional") and c.get("pending"):
            text = f"Run health: {c.get('checks_done') or 0} checks done, {c['pending']} wait for the read count"
        else:
            text = c.get("headline") or ""
        return f"{key:<19} {word:<18} {text}".rstrip()

    def _quick_phase(self) -> None:
        todo = [r for r in self.runs.values() if r.need_quick]
        new = [r for r in self.runs.values() if r.new and r.need_quick]
        waiting = [r for r in self.runs.values() if r.need_pass]
        self.t1_total = len(waiting)
        self.t1_waiting_bytes = sum(self._pass_bytes(r.key) for r in waiting)
        if new or self.result.created:
            with_reports = sum(1 for r in new if r.ident.kind == "run")
            raw = len(new) - with_reports
            to_count = sum(self._pass_bytes(r.key) for r in new if r.need_pass)
            tail = f", {_gb(to_count)} of read tables to count" if to_count and self.o.tier != "2" else ""
            self.say("scan.found", f"found {len(new)} {'new ' if not self.result.created else ''}run{'s' if len(new) != 1 else ''} "
                     f"({with_reports} with MinKNOW reports, {raw} with raw data only){tail}",
                     new=len(new), with_reports=with_reports, raw_only=raw, bytes_to_count=to_count)
        done = []
        for r in sorted(todo, key=lambda r: (r.started or "", r.key), reverse=True):       # newest first
            if r.ident is None:
                r.ident = self._identify(r.folders[0][1].path)
                if r.ident is None:
                    continue
            if self._quick(r):
                done.append(r.key)
                if r.new:
                    self.result.new.append(r.key)
        for key in self._by_start(done):                                                    # judged oldest first
            self._judge(key)
        for key in sorted(done, key=lambda k: (self.runs[k].started or "", k), reverse=True):                # newest first
            if key in self.cards:
                self.say("run.quick", self._line(key), run_key=key, verdict=self.cards[key]["columns"].get("verdict"))
        self.result.runs = len(self.cards)

    def _pass_bytes(self, key: str) -> int:
        ident = self._ident(key)
        return sum(int(ident.sizes.get(role) or 0) for role in PASS_ROLES) if ident else 0

    def _no_read_data(self) -> bool:
        """True when this scan neither links nor reads read-data files (``--tier 1`` or ``--quick``)."""
        return self.o.tier == "1" or self.o.quick

    def _out_of_time(self) -> bool:
        return self.deadline is not None and _monotonic() >= self.deadline

    def _tier1_phase(self) -> None:
        todo = [r for r in self.runs.values() if r.need_pass]
        if self.o.tier == "2" or not todo:
            return
        if self.o.quick:
            self.say("scan.quick", f"--quick: the reads of {len(todo)} run{'s' if len(todo) != 1 else ''} are "
                     "counted by the next scan without --quick", waiting=len(todo))
            return
        todo.sort(key=lambda r: (r.started or "", r.key), reverse=True)                     # newest runs first
        announced = False
        for r in todo:
            key = r.key
            if self._out_of_time():
                if not announced:
                    left = sum(1 for x in todo if self.states.get(x.key, {}).get("t1") != "done")
                    self.say("scan.time_box", f"time box of {_span(self.time_box)} reached: the reads of {left} "
                             f"run{'s' if left != 1 else ''} are counted by the next scan", waiting=left)
                    announced = True
                continue
            if key not in self.states or not self.lib.has_run(key):
                continue                                    # its quick bundle failed
            if r.identified:
                self._refresh_state(r)
            st = self.states[key]
            ident = self._ident(key)
            size = self._pass_bytes(key)
            try:
                if self.o.force:
                    for i in st.get("inputs") or []:
                        path = partial_model.find_partial(self.lib.layout.partials, i.get("partial_id") or "") if i.get("tier") == 1 else None
                        if path:
                            os.unlink(path)
                self.counts["full_pass"] += 1
                self.note("pass.start", run_key=key, bytes=size)
                t = _monotonic()
                refs = list(engine_api.full_pass(ident, self.lib.layout.partials, self._pass_options(key)))
                dt = _monotonic() - t
                read = any(not ref.cached for ref in refs)
                self.note("pass.done", run_key=key, seconds=round(dt, 3), read=read,
                          partials=[{"role": ref.role, "id": ref.partial_id, "cached": ref.cached, "complete": ref.complete,
                                     "state": ref.state} for ref in refs])
                if read:
                    self.t1_bytes_done += size
                    self.t1_seconds += dt
                # the pass must answer for every large file the run has: one it could not open (no permission, the
                # prepare step failed or stalled) comes back as an error or not at all, and is never "counted"
                wanted = [role for role in PASS_ROLES if ident.files.get(role)]
                got = {ref.role for ref in refs}
                unopened = [role for role in wanted if role not in got] \
                    + [ref.role for ref in refs if ref.state == "error" and not ref.partial_id]
                if unopened and ident.run_state != "running":
                    have = {ref.role for ref in refs}
                    refs += [PartialRef(partial_id="", path="", role=role, input_path=ident.files[role], state="error",
                                        complete=False, cached=False) for role in unopened if role not in have]
                st["inputs"] = self._tier1_inputs(ident, refs, st) + [i for i in st.get("inputs") or [] if i.get("tier") == 2]
                bad = [ref for ref in refs if ref.state == "error" and ref.partial_id]
                cut = [ref for ref in refs if not ref.complete and ref.state == "truncated"]
                # "truncated": the file itself is cut off. It was counted as far as it goes and there is nothing
                # more to read until the file changes, so the run does not wait for a later scan.
                if ident.run_state == "running":
                    st["t1"] = "done"                        # its files still change: it is looked at again at every scan
                elif unopened:
                    st["t1"] = "waiting"                     # not counted: the next scan tries again
                else:
                    st["t1"] = "error" if bad else "done" if all(ref.complete for ref in refs) \
                        else "truncated" if len(cut) == sum(1 for ref in refs if not ref.complete) else "partial"
                if st["t1"] == "truncated" and read:
                    self.say("run.truncated", f"{key}   " + ", ".join(os.path.basename(ref.input_path) for ref in cut)
                             + ": the file is cut off; its reads were counted as far as it goes", level=30, run_key=key)
                if unopened and ident.run_state != "running":
                    names = ", ".join(os.path.basename(ident.files[role]) for role in unopened)
                    st.setdefault("errors", []).append({"when": now_z(), "what": "pass", "roles": unopened,
                                                        "error": "could not be opened"})
                    st["errors"] = st["errors"][-20:]
                    self.failed[key] = f"{names} could not be opened; its reads are not counted, the next scan tries again"
                    self.say("run.error", f"{key}   {self.failed[key]}", level=30, run_key=key, roles=unopened)
                elif bad:
                    st.setdefault("errors", []).append({"when": now_z(), "what": "pass", "roles": [b.role for b in bad]})
                    self.failed[key] = "a pass over " + ", ".join(b.role for b in bad) + " ended in an error"
                    self.say("run.error", f"{key}   {self.failed[key]}", level=30, run_key=key)
                self.lib.write_state(key, st)
                if self._assemble(key):
                    self._judge(key)
                    if st["t1"] != "waiting":
                        self.t1_done += 1
                        self.t1_waiting_bytes = max(0, self.t1_waiting_bytes - size)
                    self._progress_line()
                    self._write_index(scanning=True)
            except NotBuiltError:
                raise
            except Exception as e:  # noqa: BLE001
                st["t1"] = "error"
                st.setdefault("errors", []).append({"when": now_z(), "what": "pass", "error": str(e)})
                self.lib.write_state(key, st)
                self._fail(key, "counting the reads", e)

    def _progress_line(self) -> None:
        left = self.t1_total - self.t1_done
        now = _monotonic()
        if left <= 0 or now - self._last_line < 15.0:
            return
        self._last_line = now
        rate = self.t1_bytes_done / self.t1_seconds if self.t1_seconds > 0 and self.t1_bytes_done else None
        eta = f", {_about(self.t1_waiting_bytes / rate / 60)} left ({rate / 1e6:,.0f} MB/s)" if rate else ""
        self.say("scan.counting", f"counting reads: {self.t1_done} of {self.t1_total} runs done{eta}",
                 done=self.t1_done, total=self.t1_total)

    def _reading_line(self) -> None:
        now = _monotonic()
        if now - self._last_line < 15.0:
            return
        self._last_line = now
        wanted = self._files_to_read()
        read = sum(1 for li in wanted if self._pass_state(li) in READ_STATES)
        rate = self.rd_bytes_done / self.rd_seconds if self.rd_seconds > 0 and self.rd_bytes_done else None
        self.say("scan.reading", f"read data: {read} of {len(wanted)} files read" + (f" ({rate / 1e6:,.0f} MB/s)" if rate else ""),
                 read=read, total=len(wanted))

    def _tier2_phase(self) -> None:
        if self._no_read_data():
            stale = [k for k in self.states if self.states[k].get("rd_hash") != self._rd_hash(k)]
        else:
            order = self._by_start([k for k in self.inputs_for if k in self.states])[::-1]   # newest runs first
            known = sorted(self.states)
            announced = False
            tried: set[str] = set()                          # no file is read twice in one scan
            # Two rounds, each newest runs first. What an earlier version of the pass had read is read again first:
            # a number the library showed before an upgrade is back before anything new is counted.
            again = [li for li in self.linked if self._read_by_earlier_pass(li) and self._wants_pass(li)
                     and self._reads_kind(readdata_api.pass_kind(li.candidate.role))]
            if again:
                size = sum(li.candidate.size for li in {li.fingerprint.get("value", ""): li for li in again}.values())
                n = len({li.fingerprint.get("value", "") for li in again})
                self.say("scan.read_again", f"read data: {n} file{'s' if n != 1 else ''} ({_gb(size)}) that an earlier version "
                         "of FlowCellar read are read again first, because this version counts more from each file", files=n, bytes=size)
            for earlier in ((True, False) if again else (False,)):
                for key in order:
                    if not self._in_scope(key, None if not self.o.since else ((self.states[key].get("identity") or {}).get("started"))):
                        continue
                    read_now = False
                    for li in self.inputs_for[key]:
                        fp = li.fingerprint.get("value", "")
                        p = self.passes.get(fp) or {}
                        if not self._wants_pass(li) or fp in tried:
                            continue
                        if not self._reads_kind(readdata_api.pass_kind(li.candidate.role)):
                            continue
                        was_earlier = self._read_by_earlier_pass(li)
                        if earlier and not was_earlier:
                            continue
                        if self._pass_state(li) in READ_STATES and not self.o.force:
                            continue
                        if p.get("pass") in ("error", "stalled") and not was_earlier and not (self.o.retry_failed or self.o.force):
                            continue
                        if self._out_of_time():
                            if not announced:
                                self.say("scan.time_box", f"time box of {_span(self.time_box)} reached: the remaining read-data "
                                         "files are read by the next scan")
                                announced = True
                            break
                        tried.add(fp)
                        self._read_one(li, known)
                        read_now = True
                        self._reading_line()
                    if read_now and self.lib.has_run(key) and self.states[key].get("rd_hash") != self._rd_hash(key):
                        # the run's files are read (or the time is up): its bundle, verdict and page follow now, not at
                        # the end of the scan, so that the hub fills in while the reading goes on
                        self._save_links()
                        if self._assemble(key):
                            self._judge(key)
                            self._write_index(scanning=True)
            self._save_links()
            stale = [k for k in self.states if self.states[k].get("rd_hash") != self._rd_hash(k)]
        for key in self._by_start(stale)[::-1]:
            if self.lib.has_run(key) and self._assemble(key):
                self._judge(key)
                self._write_index(scanning=True)

    def _read_one(self, li: LinkedInput, known: list[str]) -> None:
        fp = li.fingerprint.get("value", "")
        c = li.candidate
        p = self.passes.setdefault(fp, {"pass": "waiting", "attempts": 0, "error": None})
        before = p.get("partial_id")
        if before and li.partial_id and before != li.partial_id:
            # read by an earlier version of the pass: its record starts again, and its partial is of no use any more
            p.update({"pass": "waiting", "attempts": 0, "error": None})
            stale = partial_model.find_partial(self.lib.layout.partials, before)
            if stale:
                try:
                    os.unlink(stale)
                except OSError:
                    pass
        try:
            if self.o.force and (p.get("partial_id") or li.partial_id):
                path = partial_model.find_partial(self.lib.layout.partials, p.get("partial_id") or li.partial_id)
                if path:
                    os.unlink(path)
            self.counts["read_pass"] += 1
            self.note("read_pass.start", path=c.path, role=c.role, bytes=c.size)
            t = _monotonic()
            opts = PassOptions(jobs=self.o.jobs or int(self.cfg["scan"]["jobs"]),
                               threads=self.o.threads or int(self.cfg["scan"]["threads"]),
                               gentle=bool((self._root(c.root_id) or {}).get("live")),
                               stall_timeout_s=float(self.cfg["discovery"]["stall_timeout_seconds"]))
            ref = readdata_api.read_pass(li, self.lib.layout.partials, known, opts)
            dt = _monotonic() - t
            if ref is None:
                p.update({"pass": "none", "error": None, "path": c.path})
            else:
                p.update({"partial_id": ref.partial_id, "path": c.path, "seconds": round(dt, 2),
                          "pass": "done" if ref.complete else ref.state if ref.state in ("error", "truncated") else "partial",
                          "error": None if ref.complete else ref.state,
                          "attempts": int(p.get("attempts") or 0) + (0 if ref.cached else 1)})
                if not ref.cached:
                    self.rd_bytes_done += c.size
                    self.rd_seconds += dt
                if ref.state == "error":
                    self.failed[c.path] = "the pass ended in an error"
            self.note("read_pass.done", path=c.path, seconds=round(dt, 3), state=p["pass"],
                      cached=bool(ref and ref.cached), partial=ref.partial_id if ref else None)
        except NotBuiltError:
            raise
        except Exception as e:  # noqa: BLE001
            p.update({"pass": "error", "error": str(e), "attempts": int(p.get("attempts") or 0) + 1, "path": c.path})
            self.failed[c.path] = str(e)
            self.say("input.error", f"{_short(c.path)}: not read ({e})", level=30, path=c.path, error=str(e))

    def _final_judge(self) -> None:
        """Judge every run that changed in this scan once more, oldest first, so that the stored verdicts do not
        depend on the order in which the passes happened to finish."""
        for key in self._by_start(k for k in self.changed if k in self.states and self.lib.has_run(k)):
            self._judge(key)

    # ------------------------------------------------------------------------------------------ analyses
    def _analysis_rows(self) -> list[dict]:
        """The rows of the index's analyses table: what this scan made of them, else what the library holds."""
        if self._analyses is None:
            self._analyses = index_mod.analysis_rows(self.lib)
        return self._analyses

    def _analyses_phase(self) -> None:
        """Single-cell and spatial analyses the walk found (their anchor files): read, linked to the runs, judged.
        The add-on is loaded only when there is such a folder or the library already holds an analysis; one that
        cannot be read costs a line, never the scan."""
        anchors = list(getattr(self.found, "analyses", ()) or ())
        if self.single or self.o.only or (not anchors and not self._analysis_rows()):
            return
        try:
            from flowcellar.sc import scan as sc_scan
            sc_scan.run(self, anchors)
        except NotBuiltError:
            raise
        except Exception as e:  # noqa: BLE001
            self.say("analyses.error", f"the analyses were not brought up to date ({getattr(e, 'message', None) or e})",
                     level=30, kind=type(e).__name__)

    # ------------------------------------------------------------------------------------------ index and hub
    def _waiting(self) -> tuple[int, int, int, int]:
        t1 = [k for k, st in self.states.items() if st.get("kind") == "run" and st.get("t1") in ("waiting", "partial")
              and st.get("source_state") != "missing"]
        t1_bytes = sum(sum(int(((st.get("identity") or {}).get("sizes") or {}).get(role) or 0) for role in PASS_ROLES)
                       for k, st in self.states.items() if k in t1)
        rd = [li for li in self.linked if self._wants_pass(li) and self._pass_state(li) in ("waiting", "partial")]
        return len(t1), t1_bytes, len(rd), sum(li.candidate.size for li in rd)

    def _progress(self, scanning: bool) -> dict:
        t1_runs, t1_bytes, rd_files, rd_bytes = self._waiting()
        wanted = self._files_to_read()
        read = sum(1 for li in wanted if self._pass_state(li) in READ_STATES)
        rate1 = self.t1_bytes_done / self.t1_seconds if self.t1_seconds > 0 and self.t1_bytes_done else None
        rate2 = self.rd_bytes_done / self.rd_seconds if self.rd_seconds > 0 and self.rd_bytes_done else None
        eta1 = round(t1_bytes / rate1 / 60, 1) if rate1 and t1_bytes else None
        eta2 = round(rd_bytes / rate2 / 60, 1) if rate2 and rd_bytes else None
        counted = sum(1 for st in self.states.values() if st.get("kind") == "run" and st.get("t1") in ("done", "truncated"))
        with_reports = sum(1 for st in self.states.values() if st.get("kind") == "run")
        parts = []
        if t1_runs:
            parts.append(f"Counting reads: {counted} of {with_reports} runs done"
                         + (f", {_about(eta1)} left" if eta1 and scanning else "") + ".")
        if rd_files:
            parts.append(f"Read data: {read} of {len(wanted)} files read"
                         + (f"; {_about(eta2)} of reading left" if eta2 and scanning else "") + "; newest runs first.")
        if not scanning and (t1_runs or rd_files):
            parts.append("The next scan goes on from here.")
        return {"state": "scanning" if scanning else "idle", "started": self.started_z if scanning else None,
                "host": socket.gethostname() if scanning else None,
                "time_box_minutes": self.time_box if scanning else None,
                "tier1": {"runs_waiting": t1_runs, "bytes_waiting": t1_bytes, "bytes_done": self.t1_bytes_done,
                          "mb_per_s": round(rate1 / 1e6, 1) if rate1 else None, "eta_minutes": eta1},
                "read_data": {"files_total": len(wanted), "files_read": read, "files_waiting": rd_files,
                              "bytes_waiting": rd_bytes, "eta_minutes": eta2},
                "message": " ".join(parts) or None}

    def _late_notices(self) -> list[dict]:
        out = []
        stale = []
        try:
            for key, card in self.cards.items():
                stored = card["columns"].get("history_hash")
                if stored is None:
                    continue
                now = engine_api.history_hash(card, [c for k, c in self.cards.items() if k != key])
                if now is not None and now != stored:
                    stale.append(key)
        except NotBuiltError:
            stale = []
        if stale:
            out.append({"code": "judged_before_older_runs", "level": "info", "root": None, "since": None,
                        "message": f"{len(stale)} run{'s were' if len(stale) != 1 else ' was'} judged before older runs "
                                   "were added to the library. Their verdicts are kept; "
                                   f"'flowcellar recompute --status-only --rebaseline -L {self.dir}' judges them again.",
                        "runs": sorted(stale)[:50]})
        bad_inputs = [p for p in self.failed if os.sep in p]
        if bad_inputs:
            out.append({"code": "inputs_failed", "level": "warn", "root": None, "since": None,
                        "message": f"{len(bad_inputs)} file{'s' if len(bad_inputs) != 1 else ''} or folder"
                                   f"{'s' if len(bad_inputs) != 1 else ''} could not be read in the last scan; "
                                   "'flowcellar scan --retry-failed' tries again."})
        return out

    def _index_notices(self, baselines: dict | None) -> list[dict]:
        out = [n for n in self.notices]
        for g in (baselines or {}).get("groups") or []:
            if any(isinstance(m, dict) and m.get("spread_guard") for m in (g.get("metrics") or {}).values()):
                out.append({"code": "baseline_spread_guard", "level": "info", "root": None, "since": None,
                            "message": f"Comparison group '{g.get('group')}': its runs differ so much that lab-typical "
                                       "values are not used for verdicts. Splitting the group in flowcellar.toml "
                                       "([[baselines.groups]]) may help."})
        return out

    def _write_index(self, *, scanning: bool, force: bool = False) -> None:
        if self.o.no_index:
            return
        now = _monotonic()
        if not force and now - self._last_index < INDEX_EVERY_SECONDS:
            return
        self._last_index = now
        try:
            baselines = engine_api.baselines_today(list(self.cards.values()))
        except NotBuiltError:
            baselines = {"groups": []}
        notices = self._index_notices(baselines) + ([] if scanning else self._late_notices())
        progress = self._progress(scanning)
        scanned = now_z()
        self.lib.write_aux(self.lib.layout.roots_state, self.roots_state)
        self.lib.write_aux(self.lib.layout.notices, {"scanned": scanned, "notices": notices, "progress": progress})
        roots = []
        for r in self.cfg["roots"]:
            s = self.roots_state["roots"].get(r["id"], {})
            roots.append({**r, "state": s.get("state", "reachable"), "reason": s.get("reason"),
                          "offline_since": s.get("offline_since"), "last_scanned": s.get("last_scanned"),
                          "read_data_files": s.get("read_data_files") or {"found": 0, "read": 0}})
        index = index_mod.build_index(
            list(self.cards.values()),
            library={"id": self.lib.id, "name": self.lib.name, "path": self.lib.path, "created": self.lib.created, "scanned": scanned},
            roots=roots, notices=notices, progress=progress, baselines=baselines, analyses=self._analysis_rows())
        self.lib.write_index(index)
        self.note("index.written", runs=index["runs"]["n"], scanning=scanning)
        if not force:
            self._write_hub()

    def _write_hub(self, first: bool = False, final: bool = False) -> None:
        if self.o.no_hub or self.o.no_index:
            return
        only = None if (first or final) else sorted(self.hub_dirty)
        try:
            out = report_api.write_hub(self.lib, sidecars=bool(self.cfg["hub"]["sidecars"]), reports=bool(self.cfg["hub"]["reports"]),
                                       only=only)
            self.hub_dirty.clear()
            self.note("hub.written", written=len(out.get("written") or []), removed=len(out.get("removed") or []))
            if first:
                self.say("hub.ready", f"hub ready: {_short(os.path.join(self.lib.layout.hub, 'index.html'))}     "
                         + ("open it now; it fills in as the counting goes on" if any(r.need_pass for r in self.runs.values())
                            and self.o.tier != "2" and not self.o.quick else "open it in a browser, or run: flowcellar serve "
                            + _short(self.dir) + " --open"), hub=self.lib.layout.hub)
        except FlowCellarError as e:
            if not self._hub_warned:
                self.say("hub.error", f"the hub was not written: {e.message}", level=30)
                self._hub_warned = True

    def _finish(self) -> None:
        self._save_links()
        for rid in [r["id"] for r in self.cfg["roots"]]:
            if rid in self.offline:
                continue
            mine = [li for li in self.linked if li.candidate.root_id == rid]
            rs = self.roots_state["roots"].setdefault(rid, {})
            rs["read_data_files"] = {"found": (rs.get("read_data_files") or {}).get("found", len(mine)),
                                     "read": sum(1 for li in mine if self._pass_state(li) in READ_STATES)}
        self._write_index(scanning=False, force=True)
        self._write_hub(final=True)
        t1_runs, _, rd_files, _ = self._waiting()
        res = self.result
        res.runs, res.waiting_tier1, res.waiting_read_data = len(self.cards), t1_runs, rd_files
        res.failed, res.offline_roots = sorted(self.failed), sorted(self.offline)
        tally = index_mod._tally(c["columns"].get("verdict") for c in self.cards.values())
        words = [f"{tally['pass']} Pass", f"{tally['warn']} Warning", f"{tally['fail']} Fail", f"{tally['incomplete']} Incomplete"]
        if tally["running"]:
            words.append(f"{tally['running']} Running")
        if tally["unknown"]:
            words.append(f"{tally['unknown']} Not assessed")
        parts = [f"{len(self.cards)} run{'s' if len(self.cards) != 1 else ''}", ", ".join(words)]
        if not (res.new or res.changed or res.moved or res.missing):
            parts.insert(0, "nothing new")
        if t1_runs:
            parts.append(f"the reads of {t1_runs} run{'s' if t1_runs != 1 else ''} are still to be counted (scan again)")
        cut = sum(1 for st in self.states.values() if st.get("t1") == "truncated" and st.get("source_state") != "missing")
        if cut:
            parts.append(f"{cut} run{'s have' if cut != 1 else ' has'} a read table that is cut off (counted as far as it goes)")
        wanted = [li for li in self.linked if self._wants_pass(li)]
        found_n = len(self.found.candidates)
        if not wanted and self._no_read_data():
            parts.append("read data not looked at (" + ("--quick" if self.o.quick else "--tier 1") + ")")
        elif not wanted and not found_n:
            parts.append("read data not linked yet (no file patterns)" if not self.cfg["inputs"] else "no read-data files found")
        elif not wanted:
            parts.append(f"read data: {found_n:,} file{'s' if found_n != 1 else ''} found, none linked to a run")
        else:
            big = self._files_to_read()
            read = sum(1 for li in big if self._pass_state(li) in READ_STATES)
            tables = len(wanted) - len(big)
            parts.append((f"read data: {read} of {len(big)} files read" + (f", {rd_files} waiting" if rd_files else "") if big
                          else "read data: no file needs a pass")
                         + (f" ({tables:,} small table{'s' if tables != 1 else ''} linked)" if tables else ""))
        if self.offline:
            parts.append("skipped: " + ", ".join(f"root '{r}'" for r in sorted(self.offline)))
        if self.unreadable:
            n = len(self.unreadable)
            parts.append(f"{n} folder{'s' if n != 1 else ''} could not be read and {'were' if n != 1 else 'was'} skipped")
        if self.failed:
            parts.append(f"{len(self.failed)} could not be read (see {_short(res.log)})" if res.log else f"{len(self.failed)} could not be read")
        res.exit_code = int(ExitCode.PARTIAL) if (self.failed or self.offline or self.unreadable) else int(ExitCode.OK)
        self.say("scan.done", "done: " + " · ".join(parts), seconds=round(_monotonic() - self.t0, 2), counts=self.counts,
                 exit_code=res.exit_code, runs=len(self.cards), verdicts=tally, waiting_tier1=t1_runs,
                 waiting_read_data=rd_files)

    # ------------------------------------------------------------------------------------------ one run by path
    def add_single(self, run_dir: str, *, quick: bool = False) -> str:
        """``flowcellar run RUN_DIR -L DIR``: add or refresh one run folder under the scan lease. No root is walked,
        so nothing else the library knows changes state. Returns the run key."""
        from flowcellar.discovery.walk import FoundFolder
        self.lib = LibraryStore.open(self.dir, config_file=self.o.config_file)
        self.cfg = self.lib.config
        self.single = True
        run_dir = os.path.abspath(run_dir)
        wait = self.o.wait if self.o.wait is not None else True
        lease = Lease(self.lib.layout.lease, wait=wait, poll_seconds=self.o.lease_poll_seconds,
                      stale_seconds=self.o.lease_stale_seconds, on_wait=self._on_wait)
        lease.acquire()
        self.t0 = _monotonic()
        restore = self._catch_sigterm()
        try:
            self._open_log()
            self.note("add.start", library=self.dir, run_dir=run_dir, version=__version__)
            self._load()
            # no root is walked and nothing is linked afresh here: what the last scan said about them still holds
            kept = (self.lib.read_aux(self.lib.layout.notices, {}) or {}).get("notices") or []
            self.notices.extend(n for n in kept if n.get("code") in ("root_offline", "subtree_over_budget", "two_candidate_sources"))
            self.found = Discovery(roots={})
            ident = engine_api.identify(run_dir)
            key = ident.run_key
            tomb = self.lib.layout.tombstone(key)
            if os.path.isfile(tomb):
                os.unlink(tomb)                             # adding a run by hand takes back an earlier rm
            root = max((r for r in self.cfg["roots"] if run_dir == r["path"] or run_dir.startswith(r["path"].rstrip("/") + "/")),
                       key=lambda r: len(r["path"]), default=None)
            rel = (os.path.relpath(run_dir, root["path"]) if run_dir != root["path"] else "") if root else None
            folder = FoundFolder(ident.kind, run_dir, rel if rel is not None else run_dir)
            st = self.states.get(key)
            r = _Run(key=key, folders=[(root["id"] if root else None, folder)], ident=ident, identified=True, state=st,
                     new=st is None, need_quick=True, need_pass=ident.kind == "run" and not quick, started=ident.started)
            self.runs[key] = r
            if st is not None:
                have = [self._abs(p) for p in st.get("paths") or []]
                if run_dir not in have:
                    now = now_z()
                    st.setdefault("paths", []).insert(0, {"root": root["id"] if root else None, "relpath": rel, "path": run_dir,
                                                          "sig": None, "first_seen": now, "last_seen": now})
                    st["source_state"], st["missing_since"] = "ok", None
            self.linked = [_linked_from_json(d) for d in self._links_saved.get("linked") or []]
            try:
                cands = list(readdata_api.candidates_for_run(ident))
                fresh = readdata_api.link(cands, [ident], config=self.cfg, probes_dir=self.lib.layout.probes) if cands else []
            except FlowCellarError as e:
                fresh = []
                self.say("input.error", f"read data of {key} not linked ({e.message})", level=30)
            have_paths = {li.candidate.path for li in fresh}
            self.linked = [li for li in self.linked if li.candidate.path not in have_paths] + list(fresh)
            for li in self.linked:
                for k in li.run_keys:
                    self.inputs_for.setdefault(k, []).append(li)
            self.o.only = (key,)
            self._quick_phase()
            self._write_index(scanning=True, force=True)
            self._tier1_phase()
            if not quick:
                self._tier2_phase()
            self._final_judge()
            self.found.candidates = [li.candidate for li in self.linked]
            self._finish()
            if key in self.failed:
                raise FlowCellarError(f"{key}: {self.failed[key]}")
            return key
        except KeyboardInterrupt:
            self._settle_after_interrupt()
            raise
        finally:
            restore()
            lease.release()
            if self._log_fh:
                self._log_fh.close()
                self._log_fh = None

    # ------------------------------------------------------------------------------------------ --dry-run
    def _dry_run(self) -> ScanResult:
        """Say what a scan would do. Lists folders and reads the small identity files; writes nothing."""
        recorded = {}
        states: dict[str, dict] = {}
        tombs: set[str] = set()
        if self.lib is not None:
            rs = self.lib.read_aux(self.lib.layout.roots_state, {"roots": {}}) or {"roots": {}}
            recorded = {rid: s.get("first_level_names") or [] for rid, s in rs.get("roots", {}).items()}
            states = self.lib.states()
            tombs = self.lib.tombstones()
        found = discover(self.cfg, recorded=recorded)
        by_path = {}
        for k, st in states.items():
            for p in st.get("paths") or []:
                root = self._root(p.get("root"))
                by_path[os.path.join(root["path"], p["relpath"]) if root and p.get("relpath") not in ("", ".", None)
                        else (root["path"] if root else p.get("path"))] = k
        new, known, to_count, raw = [], [], 0, 0
        for _root_id, folder in found.runs:
            key = by_path.get(folder.path)
            if key is None:
                try:
                    ident = engine_api.identify(folder.path)
                except NotBuiltError:
                    raise
                except FlowCellarError as e:
                    self.say("run.error", f"{_short(folder.path)}: not read as a run ({e.message})", level=30)
                    continue
                key = ident.run_key
                if key in tombs or not self._in_scope(key, ident.started):
                    continue
                if key in states or key in new:
                    known.append(key)
                    continue
                new.append(key)
                raw += ident.kind != "run"
                to_count += sum(int(ident.sizes.get(role) or 0) for role in PASS_ROLES)
            else:
                known.append(key)
        waiting = [k for k, st in states.items() if st.get("t1") in ("waiting", "partial")]
        gone = [k for k in states if k not in known and k not in tombs
                and not any(p.get("root") in {r for r, o in found.roots.items() if o.state != REACHABLE} for p in states[k].get("paths") or [])]
        for rid, o in found.roots.items():
            if o.state != REACHABLE:
                self.say("root.offline", f"root '{rid}' would be skipped: {o.reason}", level=30)
        where = "would create library " + _short(self.dir) if self.result.created else "library " + _short(self.dir)
        self.say("dry_run", f"{where}: {len(new)} new run{'s' if len(new) != 1 else ''} ({len(new) - raw} with MinKNOW "
                 f"reports, {raw} with raw data only), {_gb(to_count)} of read tables to count; "
                 f"{len(waiting)} known run{'s' if len(waiting) != 1 else ''} still waiting for the read count; "
                 f"{len(gone)} known run{'s' if len(gone) != 1 else ''} not found any more; "
                 f"{len(found.candidates):,} read-data file{'s' if len(found.candidates) != 1 else ''} found. Nothing was written.")
        for key in sorted(new):
            self.say("dry_run.new", f"  would add {key}")
        self.result.planned = {"new": sorted(new), "known": sorted(set(known)), "waiting": sorted(waiting), "gone": sorted(gone),
                               "bytes_to_count": to_count, "read_data_files": len(found.candidates),
                               "offline_roots": sorted(r for r, o in found.roots.items() if o.state != REACHABLE)}
        self.result.exit_code = int(ExitCode.OK)
        return self.result


def scan(library_dir: str, options: ScanOptions | None = None) -> ScanResult:
    return Scan(library_dir, options).run()
