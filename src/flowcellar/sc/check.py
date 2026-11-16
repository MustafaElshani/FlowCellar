"""Invariants of an analysis bundle that the schema cannot say: what ``flowcellar validate`` checks besides it."""
from __future__ import annotations

import base64

from flowcellar.sc import MODULE_ID, SECTION_CAP_BYTES, is_analysis_id
from flowcellar.util import atomic


def problems(bundle: dict) -> list[str]:
    out: list[str] = []
    if not is_analysis_id(bundle.get("bundle_id")):
        out.append(f"/bundle_id: {bundle.get('bundle_id')!r} is not an analysis id")
    ana = (bundle.get("identity") or {}).get("analysis") or {}
    if bundle.get("bundle_id") and ana.get("id") and not str(bundle["bundle_id"]).endswith(str(ana["id"])[:12]):
        out.append("/bundle_id: does not end in the first 12 characters of identity.analysis.id")
    sec = (bundle.get("modules") or {}).get(MODULE_ID)
    if not isinstance(sec, dict):
        return out + ["/modules/sc: an analysis bundle has no section 'sc'"]
    size = len(atomic.dumps(sec))
    if size > SECTION_CAP_BYTES + 4096:
        out.append(f"/modules/sc: the section is {size / 1000:.0f} kB; the limit is {SECTION_CAP_BYTES // 1000} kB")
    cells = sec.get("cells") or []
    if len(cells) != 1 or cells[0].get("at") != []:
        out.append("/modules/sc/cells: one cell without dimensions is expected")
        return out
    cell = cells[0]
    axes = bundle.get("axes") or {}
    for name, h in (cell.get("hists") or {}).items():
        if h.get("axis") not in axes:
            out.append(f"/modules/sc/cells/0/hists/{name}: its axis {h.get('axis')!r} is not defined in the bundle")
        if "w" in h and len(h["w"]) != len(h.get("n") or []):
            out.append(f"/modules/sc/cells/0/hists/{name}: n and w differ in length")
    for name, g in (cell.get("grids") or {}).items():
        want = int(g.get("rows", 0)) * int(g.get("cols", 0))
        for ch, c in (g.get("channels") or {}).items():
            try:
                n = len(base64.b64decode(c.get("data") or ""))
            except ValueError:
                n = -1
            if n != want:
                out.append(f"/modules/sc/cells/0/grids/{name}/channels/{ch}: holds {n} cells, the grid has {want}")
    for name, c in (cell.get("counts") or {}).items():
        if len(c.get("cats") or []) != len(c.get("n") or []):
            out.append(f"/modules/sc/cells/0/counts/{name}: cats and n differ in length")
    for name, c in (cell.get("curves") or {}).items():
        for series, y in (c.get("y") or {}).items():
            if len(y) != len(c.get("x") or []):
                out.append(f"/modules/sc/cells/0/curves/{name}: series {series} and x differ in length")
    shares = [ln.get("share") for ln in ana.get("links") or [] if isinstance(ln.get("share"), (int, float))]
    if shares and sum(shares) > 1.01:
        out.append("/identity/analysis/links: the shares of reads add up to more than 1")
    known = {i.get("id") for i in (bundle.get("provenance") or {}).get("inputs") or []}
    for key, src in (sec.get("value_sources") or {}).items():
        if src.get("input") is not None and src["input"] not in known:
            out.append(f"/modules/sc/value_sources/{key}: names an input the bundle does not list")
    return out
