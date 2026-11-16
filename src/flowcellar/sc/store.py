"""Analyses inside a library: where their files live, the table of the index, and what a scan does with them.

    bundles/<analysis id>.fcb.json.gz     the measured content of one analysis (kind "analysis")
    status/<analysis id>.status.json      its verdict
    cards/<analysis id>.card.json         its row of the index
    state/analyses.json                   bookkeeping: the folders of each analysis, what was read, what failed
    state/analysis_links.json             flow cell runs a person stated for an analysis
    partials/an/<key>.fcc.json.gz         cache of what was read from an analysis folder; safe to delete

An analysis id (``an_spaceranger_0a1b2c3d4e5f``) never looks like a run key, so the files sit beside those of
the runs and the run code does not see them. Copies of one output folder have one id and are one analysis with
several paths. Everything here is cheap except ``collect``, which is cached by a signature of the files it
read: a second scan with nothing new opens no file of an analysis.
"""
from __future__ import annotations

import os
import time

from flowcellar import __version__
from flowcellar.model import bundle as bundle_model
from flowcellar.model import status as status_model
from flowcellar.sc import (
    ANALYSIS_ID_RE,
    MODULE_ID,
    MODULE_VERSION,
    PIPELINE_SPACERANGER,
    PIPELINE_WF,
    is_analysis_id,
    judge,
    link,
)
from flowcellar.sc import collect as collect_mod
from flowcellar.sc.detect import Claim
from flowcellar.util import atomic
from flowcellar.util.errors import FlowCellarError, InputError
from flowcellar.util.fingerprint import hash_obj
from flowcellar.util.timeutil import now_z

CARD_SUFFIX = ".card.json"
STATE_FORMAT = "flowcellar.analyses"
#: the summary values an index row carries (the scalars the hub table and the run page show; no map, no curve)
ROW_KEYS = ("barcodes_cells", "barcodes_detected", "barcodes_genes_median", "barcodes_umis_median", "barcodes_reads_per_cell_mean",
            "spatial_units_under_tissue", "spatial_genes_per_unit", "spatial_umis_per_unit", "spatial_reads_per_unit",
            "spatial_fraction_reads_under_tissue", "spatial_umis_per_mm2", "saturation_value", "sr_saturation",
            "reads_valid_barcode_fraction", "sr_valid_barcodes", "reads_full_length_fraction", "sr_usable_reads_fraction",
            "reads_gene_assigned_fraction", "mito_fraction_median", "sr_antisense", "reads_input", "sr_reads",
            "barcodes_unique_molecules", "unit_stat", "chemistry", "kit")


def check_id(analysis_id: str) -> str:
    """An analysis id is the only thing joined into a library path here: refuse anything else."""
    if not is_analysis_id(analysis_id):
        raise InputError(f"not an analysis id: {analysis_id!r}", hint="an analysis id looks like an_spaceranger_0a1b2c3d4e5f")
    return analysis_id


class AnalysisStore:
    """The analysis files of one library folder. Reading for everybody; writing for the scan and the command."""

    def __init__(self, library_path: str):
        self.path = os.path.abspath(library_path)
        j = lambda *p: os.path.join(self.path, *p)  # noqa: E731
        self.bundles, self.status, self.cards_dir = j("bundles"), j("status"), j("cards")
        self.state_file, self.links_file = j("state", "analyses.json"), j("state", "analysis_links.json")
        self.cache_dir = j("partials", "an")

    # ---- paths --------------------------------------------------------------------------------------------
    def bundle_path(self, aid: str) -> str:
        return os.path.join(self.bundles, check_id(aid) + bundle_model.BUNDLE_SUFFIX)

    def status_path(self, aid: str) -> str:
        return os.path.join(self.status, check_id(aid) + status_model.STATUS_SUFFIX)

    def card_path(self, aid: str) -> str:
        return os.path.join(self.cards_dir, check_id(aid) + CARD_SUFFIX)

    def cache_path(self, key: str) -> str:
        return os.path.join(self.cache_dir, key + ".fcc.json.gz")

    # ---- reading ------------------------------------------------------------------------------------------
    def ids(self) -> list[str]:
        try:
            names = os.listdir(self.bundles)
        except OSError:
            return []
        suf = bundle_model.BUNDLE_SUFFIX
        return sorted(n[: -len(suf)] for n in names if n.endswith(suf) and ANALYSIS_ID_RE.match(n[: -len(suf)]))

    def has(self, aid: str) -> bool:
        return is_analysis_id(aid) and os.path.isfile(self.bundle_path(aid))

    def read_measured(self, aid: str) -> dict:
        path = self.bundle_path(aid)
        if not os.path.isfile(path):
            raise InputError(f"{aid}: not in this library")
        return bundle_model.read_bundle(path)

    def read_status_file(self, aid: str) -> dict | None:
        path = self.status_path(aid)
        return status_model.read_status(path) if os.path.isfile(path) else None

    def composed(self, aid: str) -> dict:
        sf = self.read_status_file(aid) or {}
        return compose(self.read_measured(aid), sf.get("status"))

    def composed_gzip(self, aid: str) -> bytes:
        return bundle_model.to_gzip_bytes(self.composed(aid))

    def card(self, aid: str) -> dict | None:
        path = self.card_path(aid)
        return atomic.read_json(path) if os.path.isfile(path) else None

    def cards(self) -> dict[str, dict]:
        out = {}
        try:
            names = os.listdir(self.cards_dir)
        except OSError:
            return out
        for n in sorted(names):
            aid = n[: -len(CARD_SUFFIX)] if n.endswith(CARD_SUFFIX) else ""
            if is_analysis_id(aid):
                try:
                    c = atomic.read_json(os.path.join(self.cards_dir, n))
                except (OSError, ValueError):
                    continue
                if isinstance(c, dict) and "row" in c:
                    out[aid] = c
        return out

    def state(self) -> dict:
        try:
            d = atomic.read_json(self.state_file) if os.path.isfile(self.state_file) else None
        except (OSError, ValueError):
            d = None
        if not isinstance(d, dict) or d.get("format") != STATE_FORMAT:
            d = {"format": STATE_FORMAT, "analyses": {}, "anchors": {}}
        d.setdefault("analyses", {})
        d.setdefault("anchors", {})
        return d

    def stated_links(self) -> dict[str, list]:
        try:
            d = atomic.read_json(self.links_file) if os.path.isfile(self.links_file) else {}
        except (OSError, ValueError):
            d = {}
        return d if isinstance(d, dict) else {}

    # ---- writing ------------------------------------------------------------------------------------------
    def write_measured(self, aid: str, bundle: dict) -> dict:
        measured, _ = bundle_model.split(bundle)
        path = self.bundle_path(aid)
        prov = measured.setdefault("provenance", {})
        if os.path.isfile(path):
            try:
                old = bundle_model.read_bundle(path).get("provenance") or {}
                prov["created"] = old.get("created") or prov.get("created")
                prov["revision"] = int(old.get("revision") or 1) + 1
            except Exception:  # noqa: BLE001  a damaged old file is simply replaced
                prov.setdefault("revision", 1)
        else:
            prov.setdefault("revision", 1)
        measured["content_sha256"] = bundle_model.content_sha256(measured)
        raw = atomic.dumps(bundle_model.order_members(measured)).encode("utf-8")
        gz = atomic.gzip_bytes(raw)
        atomic.atomic_write_bytes(path, gz)
        return {"content_sha256": measured["content_sha256"], "bytes": len(raw), "gzip_bytes": len(gz)}

    def write_status_file(self, aid: str, sf: dict) -> None:
        status_model.write_status(self.status_path(aid), sf)

    def write_card(self, aid: str, card: dict) -> None:
        atomic.atomic_write_json(self.card_path(aid), card)

    def write_state(self, state: dict) -> None:
        atomic.atomic_write_json(self.state_file, state, indent=1)

    def write_stated_links(self, links: dict) -> None:
        atomic.atomic_write_json(self.links_file, links, indent=1)

    def remove(self, aid: str) -> None:
        for p in (self.bundle_path(aid), self.status_path(aid), self.card_path(aid)):
            try:
                os.unlink(p)
            except FileNotFoundError:
                pass


def compose(measured: dict, status: dict | None) -> dict:
    """An analysis bundle as it leaves a library or the command: measured content and its verdict."""
    return bundle_model.compose(measured, status=status, relations={"analyses": [], "superseded_by": None})


def has_analyses(library_path: str) -> bool:
    """True when a library holds any analysis. Listing only: nothing of the add-on is needed to ask."""
    try:
        return any(n.startswith("an_") for n in os.listdir(os.path.join(library_path, "cards")))
    except OSError:
        return False


# ---------------------------------------------------------------------------------------------------------
# the index row of one analysis
# ---------------------------------------------------------------------------------------------------------

def _section(bundle: dict) -> dict:
    return ((bundle.get("modules") or {}).get(MODULE_ID)) or {}


def link_run_key(ln: dict) -> str | None:
    """The run key a link points at, whether or not that run is in the library."""
    if ln.get("run_key"):
        return ln["run_key"]
    if ln.get("flow_cell_id") and ln.get("protocol_run_id_prefix"):
        return f"{ln['flow_cell_id']}_{ln['protocol_run_id_prefix'][:8]}"
    return None


def _cut(text: str, limit: int = 200) -> str:
    text = str(text or "")
    return text if len(text) <= limit else text[: limit - 1].rstrip() + "\u2026"


def _folder(path) -> str | None:
    """The last two folder names of an analysis (``outs`` left out): what tells two analyses of one sample apart."""
    parts = [x for x in str(path or "").replace("\\", "/").split("/") if x]
    if parts and parts[-1] == "outs":
        parts.pop()
    return "/".join(parts[-2:]) or None


def card_for(bundle: dict, status_file: dict | None, *, extra: dict | None = None) -> dict:
    """``{"analysis_id", "row", "facts"}``: the index row of one analysis and what the library needs of it when
    it pairs and links other analyses (its sample names, its molecules, where its sibling should be)."""
    ana = (bundle.get("identity") or {}).get("analysis") or {}
    sec = _section(bundle)
    s, params = sec.get("summary") or {}, sec.get("params") or {}
    st = (status_file or {}).get("status") or {}
    rs = judge.reasons(st)
    inputs = (bundle.get("provenance") or {}).get("inputs") or []
    notes = sec.get("notes") or []
    alarms = [n for n in notes if n.get("code") == "pipeline_alarm" and n.get("level") in ("warning", "error")]
    from_report = sorted(k for k, v in (sec.get("value_sources") or {}).items() if (v or {}).get("kind") in ("html_table", "embedded_json"))
    unit = ana.get("unit") or {}
    row = {
        "id": bundle["bundle_id"], "group": ana.get("group_key"), "label": ana.get("label"), "variant": ana.get("variant"),
        "pipeline": ana.get("pipeline"), "version": ana.get("pipeline_version"), "assay": ana.get("assay"),
        "unit": unit.get("name"), "unit_um": unit.get("size_um"), "slide": ana.get("slide_serial_capture_area"),
        "reference": ana.get("reference"), "siblings": list(ana.get("sibling_analyses") or []),
        "date": next((i.get("mtime") for i in inputs if i.get("mtime")), None),
        "verdict": st.get("verdict") or "unknown", "headline": st.get("headline") or "",
        # the sentences of the first two reasons, cut: the page of the analysis has them whole
        "reasons": [_cut(r["message"]) for r in rs[:2]], "reason_levels": [r["severity"] for r in rs[:2]],
        "reason_checks": [r["check"] for r in rs[:2]], "reasons_n": len(rs),
        "alarms": len(alarms), "from_report": len(from_report), "partial": sec.get("state") == "partial",
        "links": [{"run": link_run_key(ln), "fc": ln.get("flow_cell_id"), "share": ln.get("share"), "method": ln.get("method"),
                   "confidence": ln.get("confidence"), "resolved": bool(ln.get("resolved"))} for ln in ana.get("links") or []],
        "source_state": "ok", "folder": _folder(ana.get("path")),
    }
    for k in ROW_KEYS:
        if s.get(k) is not None:
            row[k] = s[k]
    row.update(extra or {})
    facts = {"names": [n for n in (params.get("sample_names") or []) if n] + ([ana.get("label")] if ana.get("label") else []),
             "molecules": s.get("barcodes_unique_molecules"), "sibling_hint": params.get("sibling_hint"),
             "sketch": params.get("read_id_sketch"), "path": ana.get("path")}
    return {"analysis_id": bundle["bundle_id"], "row": row, "facts": facts}


def index_rows(library_path: str) -> list[dict]:
    """The ``analyses`` table of the index from the cards on disk, newest analysis first."""
    rows = [c["row"] for c in AnalysisStore(library_path).cards().values()]
    rows.sort(key=lambda r: (str(r.get("date") or ""), str(r.get("id"))), reverse=True)
    return rows


# ---------------------------------------------------------------------------------------------------------
# runs of the library, as the ladder needs them
# ---------------------------------------------------------------------------------------------------------

def run_facts(cards: dict, states: dict, roots: list[dict]) -> list[link.RunFact]:
    """``RunFact`` for every run a library knows, from its cards (columns ``id``, ``fc``, ``run_name``) and its
    scan states (the folders the run was found in)."""
    by_id = {r.get("id"): r for r in roots}
    out = []
    for key, card in cards.items():
        cols = card.get("columns") or {}
        fc = cols.get("fc") or key.split("_")[0]
        paths = []
        for p in ((states.get(key) or {}).get("paths")) or []:
            root = by_id.get(p.get("root"))
            if root is not None:
                paths.append(root["path"] if p.get("relpath") in ("", ".", None) else os.path.join(root["path"], p["relpath"]))
            elif p.get("path"):
                paths.append(p["path"])
        out.append(link.RunFact(run_key=key, flow_cell_id=fc, run_short=key.split("_")[-1], sample_id=cols.get("run_name"),
                                paths=tuple(paths)))
    return out


def run_facts_of_library(lib) -> list[link.RunFact]:
    return run_facts(lib.cards(), lib.states(), lib.config.get("roots") or [])


# ---------------------------------------------------------------------------------------------------------
# refresh: claims -> bundles, links, verdicts, cards
# ---------------------------------------------------------------------------------------------------------

def _best(claims: list[Claim]) -> Claim:
    """The copy to read when one analysis lies in several folders: the one with the most files, then the
    shortest path."""
    return sorted(claims, key=lambda c: (-len(c.files), len(c.path), c.path))[0]


def _cache_key(aid: str, claim: Claim, opts: collect_mod.Options) -> str:
    return hash_obj({"id": aid, "signature": collect_mod.signature(claim), "module": MODULE_VERSION, "options": opts.shaping(),
                     "generator": __version__}, 40)


def _pairs(cards: dict[str, dict]) -> dict[str, str]:
    """Which analysis is the sibling of which (the workflow analysis and the Space Ranger output of the same
    reads). From what the workflow's own files say: the path of the Space Ranger BAM in its parameters, a Space
    Ranger folder inside its output folder, or a sibling a person named. A parameter path that does not exist
    here still pairs when exactly one Space Ranger analysis of the library ends in the same three folders (the
    same tree copied to another disk)."""
    sr_by_path: dict[str, str] = {}
    for aid, c in cards.items():
        if c["row"].get("pipeline") == PIPELINE_SPACERANGER:
            for p in c.get("paths") or [c["facts"].get("path")]:
                if p:
                    sr_by_path[os.path.normpath(p)] = aid
    pairs: dict[str, str] = {}
    for aid, c in cards.items():
        if c["row"].get("pipeline") != PIPELINE_WF:
            continue
        hint = c["facts"].get("sibling_hint") or {}
        cands = [hint.get("stated"), hint.get("outs")] + list(hint.get("outs_below") or [])
        hit = next((sr_by_path[os.path.normpath(p)] for p in cands if p and os.path.normpath(p) in sr_by_path), None)
        if hit is None and hint.get("tail"):
            tail = "/" + hint["tail"].strip("/")
            same = sorted({a for p, a in sr_by_path.items() if p.endswith(tail)})
            hit = same[0] if len(same) == 1 else None
        if hit is not None and hit not in pairs.values():
            pairs[aid] = hit
    out = dict(pairs)
    out.update({v: k for k, v in pairs.items()})
    return out


class Refresher:
    """Brings the analyses of a library up to date from a list of claims. Used by a scan (all claims the walk
    found) and by ``flowcellar analysis DIR -L`` (the claims of one folder)."""

    def __init__(self, lib, *, say=None, note=None):
        self.lib = lib
        self.store = AnalysisStore(lib.path)
        self.say = say or (lambda event, message, **kw: None)
        self.note = note or (lambda event, **kw: None)
        self.state = self.store.state()
        self.cards = self.store.cards()
        self.failed: dict[str, str] = {}
        self.changed: set[str] = set()
        self.new: list[str] = []
        self.waiting = 0                               # analyses whose per-read table is still to be read

    # ---- step 1: which analysis is each claim ------------------------------------------------------------
    def identify(self, claims: list[Claim]) -> dict[str, list[Claim]]:
        by_id: dict[str, list[Claim]] = {}
        for c in claims:
            try:
                aid = collect_mod.quick_id(c)
            except (FlowCellarError, OSError, ValueError, KeyError) as e:
                self.failed[c.path] = getattr(e, "message", None) or str(e)
                self.say("analysis.error", f"{c.path}   not readable as an analysis ({self.failed[c.path]})", level=30, path=c.path)
                continue
            by_id.setdefault(aid, []).append(c)
        return by_id

    # ---- step 2: read what is new or changed ---------------------------------------------------------------
    def collect(self, aid: str, claims: list[Claim], opts: collect_mod.Options, *, force: bool = False) -> dict | None:
        claim = _best(claims)
        key = _cache_key(aid, claim, opts)
        entry = self.state["analyses"].setdefault(aid, {})
        is_new = not self.store.has(aid)
        entry.update(paths=sorted({c.path for c in claims}), claim=claim.to_json(), seen=now_z(), source_state="ok")
        cache = self.store.cache_path(key)
        if not force and entry.get("collected") == key and os.path.isfile(cache):
            try:
                return atomic.read_json(cache)
            except (OSError, ValueError):
                pass
        t0 = time.monotonic()
        try:
            measured = collect_mod.collect(claim, opts)
        except (FlowCellarError, OSError, ValueError, KeyError, MemoryError) as e:
            msg = getattr(e, "message", None) or f"{type(e).__name__}: {e}"
            self.failed[aid] = msg
            entry["error"] = msg
            self.say("analysis.error", f"{aid}   could not be read ({msg})", level=30, analysis=aid, path=claim.path)
            return None
        entry.pop("error", None)
        notes = (_section(measured).get("notes")) or []
        partial_pass = any(n.get("code") == "read_pass_partial" for n in notes)
        skipped_pass = any(n.get("code") == "read_pass_skipped" for n in notes) and not opts.read_pass and "read_summary" in claim.files
        entry["read_pass"] = "waiting" if (partial_pass or skipped_pass) else "done" if "read_summary" in claim.files else "none"
        if not partial_pass:
            os.makedirs(os.path.dirname(cache), exist_ok=True)
            atomic.atomic_write_json(cache, measured)
            entry["collected"] = key
        if is_new:
            self.new.append(aid)
        self.note("analysis.collected", analysis=aid, path=claim.path, seconds=round(time.monotonic() - t0, 2),
                  read_pass=entry["read_pass"])
        return measured

    # ---- step 3: pair, link, judge, store ------------------------------------------------------------------
    def finish(self, measured: dict[str, dict], runs: list[link.RunFact], *, profiles: dict | None = None,
               stated: dict[str, list] | None = None) -> None:
        """``measured`` holds the collected bundle of every analysis that is to be (re)written; analyses of the
        library that are not in it take part through their cards (pairing, inherited links)."""
        stated = dict(self.store.stated_links() if stated is None else stated)
        prelim = dict(self.cards)
        for aid, b in measured.items():
            c = card_for(b, None)
            c["paths"] = (self.state["analyses"].get(aid) or {}).get("paths") or [c["facts"].get("path")]
            prelim[aid] = c
        pairs = _pairs(prelim)
        # the workflow's analyses first: a Space Ranger output takes the links of its sibling
        order = sorted(measured, key=lambda a: (prelim[a]["row"].get("pipeline") != PIPELINE_WF, a))
        final_links: dict[str, list] = {}
        for aid in order:
            b = measured[aid]
            ana = b["identity"]["analysis"]
            sib = pairs.get(aid)
            sib_links = None
            if sib and ana.get("pipeline") == PIPELINE_SPACERANGER:
                sib_links = final_links.get(sib)
                if sib_links is None and sib in self.cards:
                    try:
                        sib_links = (self.store.read_measured(sib)["identity"]["analysis"].get("links")) or None
                    except (FlowCellarError, OSError, ValueError):
                        sib_links = None
            params = _section(b).get("params") or {}
            intrinsic = [ln for ln in ana.get("links") or [] if ln.get("method") in ("read_tag", "bam_header_rg")]
            try:
                links = link.resolve(intrinsic=intrinsic, runs=runs, sibling_links=sib_links, sketch=params.get("read_id_sketch"),
                                     names=prelim[aid]["facts"]["names"], analysis_path=ana.get("path") or "",
                                     stated=[tuple(x) for x in stated.get(aid) or []] or None)
            except FlowCellarError as e:
                links = link.join(intrinsic, runs)
                self.say("analysis.link", f"{aid}   {e.message}", level=30, analysis=aid)
            final_links[aid] = links
            ana["links"] = links
            ana["sibling_analyses"] = [sib] if sib else []
            if sib and ana.get("pipeline") == PIPELINE_WF:
                # the workflow does not know its slide: it takes the group of the Space Ranger output of the same reads
                ana["group_key"] = prelim[sib]["row"].get("group") or ana.get("group_key")
            b["provenance"]["library"] = {"id": self.lib.id, "name": self.lib.name}
            self._store(aid, b, profiles, prelim.get(sib) if sib else None)
        self.store.write_state(self.state)

    def _store(self, aid: str, bundle: dict, profiles: dict | None, sibling_card: dict | None) -> None:
        entry = self.state["analyses"].setdefault(aid, {})
        sha = bundle_model.content_sha256(bundle)
        if sha != entry.get("content_sha256") or not self.store.has(aid):
            info = self.store.write_measured(aid, bundle)
            entry.update(content_sha256=info["content_sha256"], bundle_bytes=info["bytes"], bundle_gzip_bytes=info["gzip_bytes"],
                         written=now_z())
            self.changed.add(aid)
            self.note("analysis.written", analysis=aid, bytes=info["bytes"], gzip_bytes=info["gzip_bytes"])
        measured = dict(bundle, content_sha256=entry["content_sha256"])
        sib_summary = {"barcodes_unique_molecules": (sibling_card or {}).get("facts", {}).get("molecules")} if sibling_card else None
        st = judge.judge(measured, profiles=profiles, in_library=True, sibling=sib_summary)
        sf = status_model.new_status_file(aid, entry["content_sha256"], st, {"analyses": [], "superseded_by": None})
        old = self.store.read_status_file(aid)
        if old is None or _strip_times(old) != _strip_times(sf):
            self.store.write_status_file(aid, sf)
            self.changed.add(aid)
        else:
            sf = old
        card = card_for(measured, sf, extra={"bundle_bytes": entry.get("bundle_bytes"), "bundle_gzip_bytes": entry.get("bundle_gzip_bytes"),
                                             "source_state": entry.get("source_state", "ok"), "copies": len(entry.get("paths") or [])})
        card["paths"] = entry.get("paths") or []
        if self.cards.get(aid) != card:
            self.store.write_card(aid, card)
            self.changed.add(aid)
        self.cards[aid] = card

    # ---- analyses that were not found again ----------------------------------------------------------------
    def mark_unseen(self, seen: set[str], *, unknown=lambda path: True) -> None:
        """An analysis whose folder is gone keeps its bundle; its row says so. ``unknown(path)`` is true when the
        folder may still be there (its root was not looked at)."""
        for aid, entry in self.state["analyses"].items():
            if aid in seen or not self.store.has(aid):
                continue
            paths = entry.get("paths") or []
            state = "ok" if (not paths or any(unknown(p) for p in paths) or any(os.path.isdir(p) for p in paths)) else "missing"
            if entry.get("source_state", "ok") != state:
                entry["source_state"] = state
                card = self.cards.get(aid)
                if card:
                    card["row"]["source_state"] = state
                    self.store.write_card(aid, card)
                    self.changed.add(aid)
        self.store.write_state(self.state)

    def rows(self) -> list[dict]:
        rows = [c["row"] for c in self.cards.values()]
        rows.sort(key=lambda r: (str(r.get("date") or ""), str(r.get("id"))), reverse=True)
        return rows


def _strip_times(sf: dict) -> dict:
    """A status file without the moments it was made at: what "did the judgement change" compares."""
    out = dict(sf)
    st = dict(out.get("status") or {})
    st.pop("evaluated", None)
    judged = dict(st.get("judged") or {})
    judged.pop("at", None)
    st["judged"] = judged
    out["status"] = st
    return out


def relink(lib, *, say=None) -> list[dict]:
    """Resolve the links of every analysis again against the runs the library holds now, without reading any
    analysis folder (a run that arrived later resolves the links that named it). Returns the index rows."""
    r = Refresher(lib, say=say)
    measured = {}
    for aid in r.store.ids():
        try:
            measured[aid] = r.store.read_measured(aid)
        except (FlowCellarError, OSError, ValueError):
            continue
    if measured:
        r.finish(measured, run_facts_of_library(lib), profiles=lib.config.get("profiles"))
    return r.rows()
