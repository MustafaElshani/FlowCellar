"""The hub's search grammar in Python, for ``flowcellar ls --where`` and ``recompute --where``.

    query := term (space term)*          all terms must match
    term  := ["-"] field ":" value | word | "quoted words"
    value := item ("," item)*            any of the items
    item  := text | number[unit] | cmp number[unit] | low ".." high | date

The viewer has the same grammar in ``viewer/src/hub/query.js``; the shared vectors in ``tests/vectors/query``
hold the two equal. Rows are the run columns of the index (``library.index.rows``), each optionally with
``samples``: the names and aliases of its samples, which ``rows_with_samples`` fills in. Besides the named
fields, any run column of the index can be searched by its own id (``--where 'pores_last:>500'``).
A term that cannot be read is shown struck through and ignored by the hub; for a command it is an error,
because a command must not act on more runs than was meant.
"""
from __future__ import annotations

import datetime as dt
import re
import unicodedata
from typing import Callable

from flowcellar.util.errors import UsageError
from flowcellar.util.timeutil import parse_ts

UNIT = {"bases": {"gb": 1e9, "mb": 1e6, "kb": 1e3, "b": 1}, "count": {"m": 1e6, "k": 1e3}, "pct": {"%": 0.01},
        "hours": {"h": 1, "d": 24}, "num": {"k": 1e3, "m": 1e6}}
_TOKEN = re.compile(r'(-?)(?:([a-z_][a-z0-9_]*):)?(?:"([^"]*)"?|(\S+))', re.IGNORECASE)
_NUM = re.compile(r"^(<=|>=|<|>)?\s*(-?\d+(?:\.\d+)?)\s*([a-z%]*)$", re.IGNORECASE)
_DAY = re.compile(r"^(\d{4})(?:-(\d{1,2}))?(?:-(\d{1,2}))?$")
_SPAN = re.compile(r"^(\d+)([dwm])$", re.IGNORECASE)


def fold(s) -> str:
    """Lower case without accents: how text is compared."""
    s = unicodedata.normalize("NFD", "" if s is None else str(s))
    return "".join(ch for ch in s if not "̀" <= ch <= "ͯ").lower()


def tokenize(q: str) -> list[dict]:
    out = []
    for m in _TOKEN.finditer(q or ""):
        v = m.group(3) if m.group(3) is not None else m.group(4)
        if v is None or (v == "" and not m.group(2)):
            continue
        out.append({"neg": m.group(1) == "-", "field": m.group(2).lower() if m.group(2) else None, "value": v, "raw": m.group(0)})
    return out


def _parse_num(s: str, unit: str):
    m = _NUM.match(s)
    if not m:
        return None
    table, u = UNIT.get(unit) or UNIT["num"], m.group(3).lower()
    mul = table.get(u) if u else (0.01 if unit == "pct" else 1)
    if mul is None:
        return None
    dec = len(m.group(2).split(".")[1]) if "." in m.group(2) else 0
    return m.group(1), float(m.group(2)) * mul, 0.5 * 10 ** (-dec) * mul


def _num_test(item: str, unit: str) -> Callable[[float], bool] | None:
    r = item.split("..")
    if len(r) == 2:
        a, b = _parse_num(r[0], unit), _parse_num(r[1], unit)
        if not a or not b:
            return None
        return lambda x: a[1] <= x <= b[1]
    p = _parse_num(item, unit)
    if not p:
        return None
    cmp, v, half = p
    if cmp == "<":
        return lambda x: x < v
    if cmp == "<=":
        return lambda x: x <= v
    if cmp == ">":
        return lambda x: x > v
    if cmp == ">=":
        return lambda x: x >= v
    return lambda x: v - half <= x < v + half          # a bare number matches what rounds to it


def _day_range(s: str):
    m = _DAY.match(s)
    if not m:
        return None
    y, mo, d = int(m.group(1)), int(m.group(2) or 1), int(m.group(3) or 1)
    try:
        a = dt.datetime(y, mo, d, tzinfo=dt.UTC)
        if m.group(3):
            b = a + dt.timedelta(days=1)
        elif m.group(2):
            b = dt.datetime(y + (mo == 12), 1 if mo == 12 else mo + 1, 1, tzinfo=dt.UTC)
        else:
            b = dt.datetime(y + 1, 1, 1, tzinfo=dt.UTC)
    except ValueError:
        return None
    return a.timestamp() * 1000, b.timestamp() * 1000


def _lev(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def _gb(k: str):
    return lambda r: None if r.get(k) is None else r[k] * 1e9


# The same tables as viewer/src/hub/query.js: the two are held equal by the vectors in tests/vectors/query.
_STATUS_WORDS = {"warn": ["warn", "warning"], "pass": ["pass"], "fail": ["fail"], "incomplete": ["incomplete"],
                 "running": ["running"], "unknown": ["unknown", "not_assessed"]}
_PART_WORDS = {**_STATUS_WORDS, "not_available": ["not_available", "none", "not", "no", "na"]}


def _mods(r) -> list:
    return r.get("modifiers") or []


def _has_report(r) -> bool:
    return r.get("detail_t1") not in (None, "none")


def projects_of(r: dict) -> list[str]:
    """The ids and names of the projects of a run, folded. A ``project:`` value that is exactly one of these (for
    any run) means that project and no other; any other value is looked for inside the names."""
    vals = [r.get("project")] + list(r.get("project_ids") or []) + list(r.get("project_names") or [])
    return [fold(v) for v in vals if v not in (None, "")]


FIELDS: dict[str, dict] = {
    "status": {"type": "enum", "get": lambda r: r.get("verdict"), "words": _STATUS_WORDS},
    "health": {"type": "enum", "get": lambda r: r.get("health"), "words": _PART_WORDS},
    "data": {"type": "enum", "get": lambda r: r.get("data"), "words": _PART_WORDS},
    "fc": {"type": "text", "get": lambda r: r.get("fc"), "match": "prefix"},
    "id": {"type": "text", "get": lambda r: r.get("id"), "match": "prefix"},
    "project": {"type": "project", "get": lambda r: " ".join([str(r.get("project") or "")] + list(r.get("project_ids") or []))},
    "root": {"type": "text", "get": lambda r: " ".join(r.get("roots") or [])},
    "sample": {"type": "text", "get": lambda r: " ".join(r.get("samples") or [])},
    "name": {"type": "text", "get": lambda r: r.get("run_name")},
    "exp": {"type": "text", "get": lambda r: r.get("experiment")},
    "assay": {"type": "text", "get": lambda r: f"{r.get('assay_group') or ''} {r.get('assay') or ''}"},
    "kit": {"type": "text", "get": lambda r: r.get("kit")},
    "fctype": {"type": "text", "get": lambda r: r.get("product")},
    "device": {"type": "text", "get": lambda r: r.get("device_type")},
    "pos": {"type": "text", "get": lambda r: r.get("position")},
    "minknow": {"type": "text", "get": lambda r: r.get("minknow"), "match": "prefix"},
    "detail": {"type": "text", "get": lambda r: r.get("detail_t1")},
    "check": {"type": "text", "get": lambda r: " ".join(r.get("reason_checks") or [])},
    "yield": {"type": "num", "unit": "bases", "get": _gb("gb")},
    "yield24": {"type": "num", "unit": "bases", "get": _gb("gb_24h")},
    "yield48": {"type": "num", "unit": "bases", "get": _gb("gb_48h")},
    "yield72": {"type": "num", "unit": "bases", "get": _gb("gb_72h")},
    "passed": {"type": "num", "unit": "bases", "get": _gb("t2_passed_gb")},
    "reads": {"type": "num", "unit": "count", "get": lambda r: r.get("reads")},
    "pores": {"type": "num", "unit": "num", "get": lambda r: r.get("pores")},
    "hours": {"type": "num", "unit": "hours", "get": lambda r: r.get("hours")},
    "n50": {"type": "num", "unit": "bases", "get": lambda r: r.get("t2_n50") if r.get("t2_n50") is not None else r.get("n50_est")},
    "q": {"type": "num", "unit": "num", "get": lambda r: r.get("t2_q_med")},
    "occ": {"type": "num", "unit": "pct", "get": lambda r: r.get("occ_first")},
    "usable": {"type": "num", "unit": "pct", "get": lambda r: r.get("t2_usable")},
    "samples": {"type": "num", "unit": "num", "get": lambda r: r.get("t2_samples")},
    "use": {"type": "num", "unit": "num", "get": lambda r: r.get("use")},
    "is": {"type": "flag", "flags": {
        "reused": lambda r: (r.get("use") or 0) > 1,
        "first-use": lambda r: r.get("use") == 1 and _has_report(r),
        "restart": lambda r: bool(r.get("restart")),
        # checks of this run still wait for the read count; a run whose read table is cut off is not "counting"
        "counting": lambda r: r.get("detail_t1") == "quick" and (r.get("pending") or 0) > 0
        and "reads_partly_counted" not in _mods(r) and "reads_not_counted" not in _mods(r),
        "quick": lambda r: r.get("detail_t1") == "quick",
        "cut": lambda r: "reads_partly_counted" in _mods(r) or "reads_not_counted" in _mods(r),
        "conflict": lambda r: (r.get("t2_conflicts") or 0) > 0,
        "waiting": lambda r: (r.get("t2_files_waiting") or 0) > 0,
        "full": lambda r: r.get("detail_t1") == "full",
        "early": lambda r: "baseline_provisional" in _mods(r),
        "raw": lambda r: r.get("cause") == "raw_only",
        "running": lambda r: r.get("run_state") == "running" or r.get("verdict") == "running",
        "offline": lambda r: r.get("source_state") in ("offline", "missing"),
        "stale": lambda r: "judged_before_older_runs" in _mods(r)}},
    "has": {"type": "flag", "flags": {
        "reads": lambda r: r.get("read_data") not in (None, "absent"),
        "history": lambda r: "no_history" not in _mods(r) and _has_report(r),
        "report": _has_report,
        "project": lambda r: bool(r.get("project_ids") or r.get("project"))}},
}
ALIASES = {"verdict": "status", "flowcell": "fc", "flow_cell": "fc", "run": "name", "experiment": "exp", "product": "fctype",
           "position": "pos", "read_data": "data", "occupancy": "occ"}
DATE_FIELDS = ("after", "before", "on", "year", "last")
# what a word without a field looks through
_HAY = ("id", "fc", "run_name", "experiment", "project", "kit", "product", "assay", "assay_group", "device_type", "position",
        "folder", "headline")


def _hay(r: dict) -> str:
    parts = [r.get(k) for k in _HAY] + list(r.get("samples") or []) + list(r.get("project_ids") or []) + list(r.get("roots") or [])
    return fold(" ".join(str(p) for p in parts if p))


def _started_ms(r: dict):
    t = parse_ts(r.get("started")) if r.get("started") else None
    return None if t is None else t * 1000


def compile_term(t: dict, now_ms: float, columns=(), projects=()) -> dict:
    """``{"fn": row -> bool}`` or ``{"error": text}`` for one token. ``projects``: every project id and name of the
    rows (folded), for the exact match of ``project:``."""
    if not t["field"]:
        w = fold(t["value"])
        return {"fn": lambda r: w in _hay(r)}
    f, flip = ALIASES.get(t["field"], t["field"]), False
    if f == "no":
        f, flip = "has", True
    if f in DATE_FIELDS:
        lo, hi = float("-inf"), float("inf")
        if f == "last":
            m = _SPAN.match(t["value"])
            if not m:
                return {"error": "expected a span such as 14d, 8w or 6m"}
            lo = now_ms - int(m.group(1)) * {"d": 1, "w": 7, "m": 30.44}[m.group(2).lower()] * 864e5
        else:
            r_ = _day_range(t["value"])
            if not r_:
                return {"error": "expected a date such as 2026, 2026-11 or 2026-11-04"}
            if f == "after":
                lo = r_[0]
            elif f == "before":
                hi = r_[0]
            else:
                lo, hi = r_

        def in_range(r):
            ms = _started_ms(r)
            return ms is not None and lo <= ms < hi      # a run without a start time matches no date
        return {"field": f, "fn": in_range}
    spec = FIELDS.get(f)
    if spec is None and f in columns:                      # any run column of the index, by its id
        def get(r, _f=f):
            return r.get(_f)
        spec = {"type": "auto", "get": get}
    if spec is None:
        best, bd = None, 3
        for k in list(FIELDS) + list(ALIASES) + list(DATE_FIELDS) + ["no"]:
            d = _lev(f, k)
            if d < bd:
                best, bd = k, d
        return {"error": f"'{t['field']}' is not a field." + (f" Did you mean {best}?" if best else "")}
    items = [x for x in t["value"].split(",") if x]
    if not items:
        return {"error": "no value"}
    kind = spec["type"]
    if kind == "auto":
        tests = [_num_test(it, "num") for it in items]
        want = [fold(x) for x in items]

        def auto(r):
            v = spec["get"](r)
            if isinstance(v, bool) or v is None:
                return v is not None and fold(str(v).lower()) in want
            if isinstance(v, (int, float)):
                return any(fn(v) for fn in tests if fn)
            return any(w in fold(" ".join(map(str, v)) if isinstance(v, list) else v) for w in want)
        return {"field": f, "fn": auto}
    if kind == "num":
        tests = [_num_test(it, spec["unit"]) for it in items]
        if any(x is None for x in tests):
            return {"error": "could not read the number"}

        def num(r):
            x = spec["get"](r)
            return isinstance(x, (int, float)) and x == x and any(fn(x) for fn in tests)
        return {"field": f, "fn": num}
    if kind == "enum":
        want = [fold(x) for x in items]
        # a value that is no verdict is an error, as for is: and has: ("wran" must not read as "no warnings")
        known = list(spec["words"]) + [w for ws in spec["words"].values() for w in ws]
        bad = next((x for x in want if x not in known), None)
        if bad is not None:
            best, bd = None, 3
            for k in known:
                d = _lev(bad, k)
                if d < bd:
                    best, bd = k, d
            return {"error": f"'{bad}' is not a value of {f}; one of {', '.join(spec['words'])}." + (f" Did you mean {best}?" if best else "")}

        def enum(r):
            v = spec["get"](r)
            words = spec["words"].get(v) or [fold(v)]
            return any(x in words or fold(v) == x for x in want)
        return {"field": f, "fn": enum}
    if kind == "project":
        want = [fold(x) for x in items]
        names = set(projects)

        def project(r):
            mine = projects_of(r)
            text = " ".join(mine)
            return any((w in mine) if w in names else (w in text) for w in want)
        return {"field": f, "fn": project}
    if kind == "text":
        want = [fold(x) for x in items]
        prefix = spec.get("match") == "prefix"

        def text(r):
            v = fold(spec["get"](r))
            return any(v.startswith(w) if prefix else w in v for w in want)
        return {"field": f, "fn": text}
    fns = [spec["flags"].get(fold(it)) for it in items]
    if any(x is None for x in fns):
        return {"error": "unknown value; one of " + ", ".join(spec["flags"])}
    return {"field": f, "fn": lambda r: any(fn(r) for fn in fns) != flip}


def parse(query: str, *, now_ms: float | None = None, columns=(), projects=()) -> tuple[list[dict], list[str]]:
    """(terms, errors). Each term: ``{"raw", "neg", "fn"}``; errors name the token and what is wrong with it."""
    if now_ms is None:
        now_ms = dt.datetime.now(dt.UTC).timestamp() * 1000
    terms, errors = [], []
    for t in tokenize(query):
        c = compile_term(t, now_ms, columns, projects)
        if "fn" in c:
            terms.append({"raw": t["raw"], "neg": t["neg"], "fn": c["fn"]})
        else:
            errors.append(f"{t['raw']}: {c['error']}")
    return terms, errors


def matches(row: dict, terms: list[dict]) -> bool:
    return all(bool(t["fn"](row)) != t["neg"] for t in terms)


def filter_rows(rows: list[dict], query: str | None, *, now: str | None = None, strict: bool = True) -> list[dict]:
    """The rows that match ``query`` (all of them for an empty query). ``now`` is an RFC 3339 time for ``last:``.
    With ``strict`` a term that cannot be read is a usage error; otherwise it is ignored, as the hub does."""
    if not query or not query.strip():
        return list(rows)
    now_ms = parse_ts(now) * 1000 if now else None
    columns = set().union(*(r.keys() for r in rows)) if rows else set()
    terms, errors = parse(query, now_ms=now_ms, columns=columns, projects={p for r in rows for p in projects_of(r)})
    if errors and strict:
        raise UsageError("cannot read the search: " + "; ".join(errors),
                         hint="examples: status:warn   fc:TST00001   after:2026-03 kit:PCB114   reads:>1M")
    return [r for r in rows if matches(r, terms)]


def rows_with_samples(index: dict) -> list[dict]:
    """The runs of an index as rows, each with the names of its samples (what free-text search looks through)."""
    from flowcellar.library.index import rows as index_rows
    out = index_rows(index)
    by_run: dict[str, list[str]] = {}
    for s in index.get("samples") or []:
        for n in (s.get("name"), s.get("alias")):
            if n and n not in by_run.setdefault(s.get("run"), []):
                by_run[s.get("run")].append(n)
    for r in out:
        r["samples"] = by_run.get(r.get("id"), [])
    return out


def run_vector(vector: dict) -> list[dict]:
    """One shared test vector (``tests/vectors/query``): ``[{"ids", "inactive"}]``, one entry per query, the way
    the hub behaves: a term that cannot be read is left out and named, and the others still filter."""
    args = vector["args"]
    now_ms = parse_ts(args["now"]) * 1000
    columns = set().union(*(r.keys() for r in args["rows"])) if args["rows"] else set()
    projects = {p for r in args["rows"] for p in projects_of(r)}
    out = []
    for q in args["queries"]:
        terms, inactive = [], []
        for t in tokenize(q):
            c = compile_term(t, now_ms, columns, projects)
            if "fn" in c:
                terms.append({"raw": t["raw"], "neg": t["neg"], "fn": c["fn"]})
            else:
                inactive.append(t["raw"])
        out.append({"ids": [r["id"] for r in args["rows"] if matches(r, terms)], "inactive": inactive})
    return out
