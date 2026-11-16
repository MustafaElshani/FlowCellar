"""Merging the stored form of the aggregate kinds, and rolling a cube up.

All functions return a new object and leave their arguments untouched. Merging is associative and
commutative for every kind that is mergeable (``curves`` and ``grids`` are not: they raise NotMergeable).
"""
from __future__ import annotations

import copy

from flowcellar.model import axes as _axes

MERGE_RULES = ("sum", "min", "max", "mean", "first", "last", "same")


class NotMergeable(ValueError):
    """Two aggregates that cannot be added: different axes outside one family, a curve, a grid."""


# ---------------------------------------------------------------------------------------- scalars
def _merge_value(a, b, rule: str):
    """One value from two under a merge rule. None is "not known": a sum with an unknown part is unknown; for
    every other rule the known side is kept."""
    if rule == "sum":
        return None if a is None or b is None else a + b
    if a is None:
        return b
    if b is None:
        return a
    if rule == "min":
        return min(a, b)
    if rule == "max":
        return max(a, b)
    if rule == "first":
        return a
    if rule == "last":
        return b
    if rule == "same":
        return a if a == b else None
    raise NotMergeable(f"rule {rule!r} needs a weight and cannot be applied to two values")


def merge_scalars(a: dict, b: dict, rules: dict | None = None) -> dict:
    """Add; other rules per id in ``rules`` (the section's ``scalar_merge``). An id present on one side only is
    kept as it is. A sum with an unknown (None) part is unknown."""
    rules = rules or {}
    out = dict(a)
    for k, v in b.items():
        out[k] = _merge_value(a[k], v, rules.get(k, "sum")) if k in a else v
    return out


# ------------------------------------------------------------------------------------- histograms
def _to_axis(h: dict, axis: str) -> dict:
    """Re-bin a stored histogram exactly onto a coarser axis of its family."""
    if h["axis"] == axis:
        return h
    f = _axes.coarsen_factor(h["axis"], axis)
    out = {k: v for k, v in h.items() if k not in ("n", "w", "start", "axis")}
    out["axis"] = axis
    start = h["start"] // f
    for key in ("n", "w"):
        if key not in h:
            continue
        vals = [0] * (h["start"] - start * f) + list(h[key])
        vals += [0] * (-len(vals) % f)
        out[key] = [sum(vals[i:i + f]) for i in range(0, len(vals), f)]
    out["start"] = start
    out.setdefault("n", [])
    return out


def _coarser_common(a: str, b: str) -> str:
    if a == b:
        return a
    if b in _axes.AXES[a].get("coarsens_to", []):
        return b
    if a in _axes.AXES[b].get("coarsens_to", []):
        return a
    raise NotMergeable(f"axes {a} and {b} are not levels of one family")


def _add_offset(sa: int, va: list, sb: int, vb: list) -> tuple[int, list]:
    if not va:
        return (sb, list(vb)) if vb else (0, [])
    if not vb:
        return sa, list(va)
    s = min(sa, sb)
    e = max(sa + len(va), sb + len(vb))
    out = [0] * (e - s)
    for i, v in enumerate(va):
        out[sa - s + i] += v
    for i, v in enumerate(vb):
        out[sb - s + i] += v
    return s, out


def merge_hist(a: dict, b: dict) -> dict:
    """Add bin by bin. Two levels of one axis family merge exactly on the coarser one."""
    axis = _coarser_common(a["axis"], b["axis"])
    a, b = _to_axis(a, axis), _to_axis(b, axis)
    start, n = _add_offset(a["start"], a.get("n", []), b["start"], b.get("n", []))
    out: dict = {"axis": axis, "start": start, "n": n}
    if "w" in a or "w" in b:
        if ("w" in a) != ("w" in b) and (a.get("n") and b.get("n")):
            raise NotMergeable("one histogram has weights and the other has none")
        _, out["w"] = _add_offset(a["start"], a.get("w", []), b["start"], b.get("w", []))
        if len(out["w"]) != len(n):
            out["w"] = out["w"] + [0] * (len(n) - len(out["w"]))
        if a.get("w_name") or b.get("w_name"):
            out["w_name"] = a.get("w_name") or b.get("w_name")
    for k in ("under", "over", "missing", "w_under", "w_over"):
        if k in a or k in b:
            out[k] = a.get(k, 0) + b.get(k, 0)
    for k, fn in (("min", min), ("max", max)):
        vals = [x for x in (a.get(k), b.get(k)) if x is not None]
        if vals:
            out[k] = fn(vals)
    # trim empty ends so that the result does not depend on the order of merging
    while out["n"] and out["n"][0] == 0 and (not out.get("w") or out["w"][0] == 0):
        out["n"].pop(0); out["start"] += 1
        if "w" in out:
            out["w"].pop(0)
    while out["n"] and out["n"][-1] == 0 and (not out.get("w") or out["w"][-1] == 0):
        out["n"].pop()
        if "w" in out:
            out["w"].pop()
    if not out["n"]:
        out["start"] = 0
    return out


def merge_hist2d(a: dict, b: dict) -> dict:
    """Add cell by cell (sparse ``coo`` encoding; cells stay sorted by (i, j))."""
    if (a["x"], a["y"]) != (b["x"], b["y"]):
        raise NotMergeable(f"2-D histograms on different axes: {a['x']} x {a['y']} and {b['x']} x {b['y']}")
    if a.get("enc", "coo") != "coo" or b.get("enc", "coo") != "coo":
        raise NotMergeable("only the coo encoding is merged")
    has_w = "w" in a or "w" in b
    cells: dict[tuple[int, int], list[int]] = {}
    for h in (a, b):
        w = h.get("w") or [0] * len(h["n"])
        for i, j, n, ww in zip(h.get("i", []), h.get("j", []), h["n"], w):
            c = cells.setdefault((i, j), [0, 0])
            c[0] += n
            c[1] += ww
    keys = sorted(k for k, v in cells.items() if v[0] or v[1])
    out: dict = {"x": a["x"], "y": a["y"], "enc": "coo", "i": [k[0] for k in keys], "j": [k[1] for k in keys],
                 "n": [cells[k][0] for k in keys]}
    if has_w:
        out["w"] = [cells[k][1] for k in keys]
        if a.get("w_name") or b.get("w_name"):
            out["w_name"] = a.get("w_name") or b.get("w_name")
    for k in ("outside", "w_outside", "missing"):
        if k in a or k in b:
            out[k] = a.get(k, 0) + b.get(k, 0)
    return out


# ----------------------------------------------------------------------------------------- series
def merge_series(a: dict, b: dict) -> dict:
    """Increments add bin by bin. A level (listed in ``levels``) follows its rule in ``merge``; ``mean`` needs the
    array named by ``weights``."""
    if a["axis"] != b["axis"]:
        raise NotMergeable(f"series on different axes: {a['axis']} and {b['axis']}")
    levels = set(a.get("levels", [])) | set(b.get("levels", []))
    rules = {**(b.get("merge") or {}), **(a.get("merge") or {})}
    s = min(a["start"], b["start"])
    la = max((len(v) for v in a["values"].values()), default=0)
    lb = max((len(v) for v in b["values"].values()), default=0)
    n = max(a["start"] + la, b["start"] + lb) - s

    def spread(src: dict, name: str, fill):
        vals = src["values"].get(name)
        out = [fill] * n
        if vals is not None:
            off = src["start"] - s
            out[off: off + len(vals)] = vals
        return out
    wname = a.get("weights") or b.get("weights")
    values: dict = {}
    for name in list(a["values"]) + [k for k in b["values"] if k not in a["values"]]:
        if name not in levels:
            va, vb = spread(a, name, 0), spread(b, name, 0)
            values[name] = [x + y for x, y in zip(va, vb)]
            continue
        rule = rules.get(name)
        va, vb = spread(a, name, None), spread(b, name, None)
        if rule == "mean":
            if not wname:
                raise NotMergeable(f"level {name!r} has rule mean and the series names no weights")
            wa, wb = spread(a, wname, 0), spread(b, wname, 0)
            values[name] = [y if x is None else x if y is None else ((x * p + y * q) / (p + q) if p + q else None)
                            for x, y, p, q in zip(va, vb, wa, wb)]
        elif rule in MERGE_RULES and rule != "sum":
            values[name] = [_merge_value(x, y, rule) for x, y in zip(va, vb)]
        else:
            raise NotMergeable(f"level {name!r} has no merge rule")
    out: dict = {"axis": a["axis"], "start": s}
    if a.get("unit") or b.get("unit"):
        out["unit"] = a.get("unit") or b.get("unit")
    out["values"] = values
    if levels:
        out["levels"] = sorted(levels)
    if rules:
        out["merge"] = rules
    if wname:
        out["weights"] = wname
    return out


# ----------------------------------------------------------------------------------------- counts
def merge_counts(a: dict, b: dict) -> dict:
    """Union of the categories in order of first appearance, counts added. A category that is not known on one
    side (listed in ``missing``) stays not known: unknown plus a number is unknown."""
    cats = list(a["cats"]) + [c for c in b["cats"] if c not in a["cats"]]
    missing = [c for c in cats if c in (a.get("missing") or []) or c in (b.get("missing") or [])]
    ia = {c: k for k, c in enumerate(a["cats"])}
    ib = {c: k for k, c in enumerate(b["cats"])}

    def add(key):
        return [(a[key][ia[c]] if c in ia else 0) + (b[key][ib[c]] if c in ib else 0) for c in cats]
    out: dict = {"cats": cats, "n": add("n")}
    if "w" in a and "w" in b:
        out["w"] = add("w")
        if a.get("w_name") or b.get("w_name"):
            out["w_name"] = a.get("w_name") or b.get("w_name")
    if a.get("other") is not None or b.get("other") is not None:
        out["other"] = int(a.get("other") or 0) + int(b.get("other") or 0)     # what was folded away adds up too
    if a.get("role") or b.get("role"):
        out["role"] = a.get("role") or b.get("role")
    if missing:
        out["missing"] = missing
    if a.get("sources") or b.get("sources"):
        out["sources"] = {**(b.get("sources") or {}), **(a.get("sources") or {})}
    return out


# --------------------------------------------------------------------------------------- channels
def merge_channels(a: dict, b: dict) -> dict:
    """Per-channel arrays of one layout: add, or the per-value rule in ``merge`` (max for the last read time)."""
    if a.get("layout") != b.get("layout") or a.get("wells") != b.get("wells"):
        raise NotMergeable("channel maps of different layouts")
    rules = {**(b.get("merge") or {}), **(a.get("merge") or {})}
    values: dict = {}
    for name in list(a["values"]) + [k for k in b["values"] if k not in a["values"]]:
        va, vb = a["values"].get(name), b["values"].get(name)
        if va is None or vb is None:
            values[name] = list(va if vb is None else vb)
            continue
        if len(va) != len(vb):
            raise NotMergeable(f"channel value {name!r} has {len(va)} and {len(vb)} entries")
        rule = rules.get(name, "sum")
        values[name] = [_merge_value(x, y, rule) for x, y in zip(va, vb)]
    out = {k: v for k, v in a.items() if k not in ("values", "merge")}
    out["values"] = values
    if rules:
        out["merge"] = rules
    return out


# ----------------------------------------------------------------------------------------- tables
def merge_table(a: dict, b: dict) -> dict:
    """``merge``: ``concat`` (rows of both), ``by_key`` (rows with an equal key are added column by column where
    numeric), ``replace`` (the later one wins)."""
    how = a.get("merge", "concat")
    if how == "replace":
        return copy.deepcopy(b)
    if [c["id"] for c in a["columns"]] != [c["id"] for c in b["columns"]]:
        raise NotMergeable("tables with different columns")
    out = {k: copy.deepcopy(v) for k, v in a.items() if k != "rows"}
    if how == "concat":
        out["rows"] = [list(r) for r in a["rows"]] + [list(r) for r in b["rows"]]
        if "total_rows" in a or "total_rows" in b:
            out["total_rows"] = a.get("total_rows", len(a["rows"])) + b.get("total_rows", len(b["rows"]))
        return out
    if how == "by_key":
        ids = [c["id"] for c in a["columns"]]
        kidx = [ids.index(k) for k in a.get("key", [])]
        rows: dict[tuple, list] = {}
        for r in list(a["rows"]) + list(b["rows"]):
            key = tuple(r[i] for i in kidx)
            if key not in rows:
                rows[key] = list(r)
            else:
                cur = rows[key]
                for i, v in enumerate(r):
                    if i not in kidx and isinstance(v, (int, float)) and not isinstance(v, bool) \
                            and isinstance(cur[i], (int, float)) and not isinstance(cur[i], bool):
                        cur[i] = cur[i] + v
        out["rows"] = list(rows.values())
        return out
    raise NotMergeable(f"table merge rule {how!r}")


_MERGERS = {"hists": merge_hist, "hists2d": merge_hist2d, "series": merge_series, "counts": merge_counts,
            "channels": merge_channels, "tables": merge_table}


def merge_kind(kind: str, a: dict, b: dict):
    """Merge two aggregates of one kind (not ``scalars``, which need the section's rules)."""
    if kind not in _MERGERS:
        raise NotMergeable(f"aggregates of kind {kind!r} are not mergeable")
    return _MERGERS[kind](a, b)


# ------------------------------------------------------------------------------------------ cells
def merge_cell(a: dict, b: dict, scalar_merge: dict | None = None, *, at: list | None = None) -> dict:
    """Two cells (or two ``run_level`` objects) added member by member. ``at`` defaults to that of ``a``."""
    out: dict = {}
    if "at" in a or at is not None:
        out["at"] = list(at if at is not None else a["at"])
    if "scalars" in a or "scalars" in b:
        out["scalars"] = merge_scalars(a.get("scalars", {}), b.get("scalars", {}), scalar_merge)
    for kind in ("hists", "hists2d", "series", "counts", "channels", "tables"):
        if kind not in a and kind not in b:
            continue
        ka, kb = a.get(kind, {}), b.get(kind, {})
        out[kind] = {name: (merge_kind(kind, ka[name], kb[name]) if name in ka and name in kb
                            else copy.deepcopy(ka.get(name, kb.get(name))))
                     for name in list(ka) + [n for n in kb if n not in ka]}
    for kind in ("curves", "grids"):
        if kind in a and kind in b:
            raise NotMergeable(f"cells that both hold {kind} cannot be merged")
        if kind in a or kind in b:
            out[kind] = copy.deepcopy(a.get(kind, b.get(kind)))
    return out


def merge_cells(cells_a: list[dict], cells_b: list[dict], scalar_merge: dict | None = None) -> list[dict]:
    """Two cell lists of one section (same dims): cells with an equal ``at`` are merged, the others kept."""
    out: dict[tuple, dict] = {}
    for c in list(cells_a) + list(cells_b):
        key = tuple(c["at"])
        out[key] = merge_cell(out[key], c, scalar_merge) if key in out else copy.deepcopy(c)
    return list(out.values())


def rollup(cells: list[dict], dims: list[str], keep: list[str], scalar_merge: dict | None = None,
           where: dict | None = None) -> list[dict]:
    """Sum the cube over every dimension that is not in ``keep``.

    ``where`` restricts first: ``{"filter": ["pass"]}``. ``rollup(cells, dims, [])`` is the run total;
    ``rollup(cells, ["sample", "filter"], ["sample"])`` the per-sample view. Cells never overlap, so the sum is exact.
    """
    idx = [dims.index(d) for d in keep]
    rolled = len(idx) != len(dims)
    out: dict[tuple, dict] = {}
    for c in cells:
        if where and any(c["at"][dims.index(d)] not in allowed for d, allowed in where.items()):
            continue
        key = tuple(c["at"][i] for i in idx)
        if rolled and ("curves" in c or "grids" in c):
            # curves and grids are derived per cell and cannot be added: a view over fewer dimensions has none,
            # however many of its cells held one
            c = {k: v for k, v in c.items() if k not in ("curves", "grids")}
        if key in out:
            out[key] = merge_cell(out[key], c, scalar_merge, at=list(key))
        else:
            cc = copy.deepcopy(c)
            cc["at"] = list(key)
            out[key] = cc
    return list(out.values())
