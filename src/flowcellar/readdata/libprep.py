"""What the library preparation and annotation steps of a pipeline wrote as tables, as sections of a run:

* ``cdna.pychopper``: every category of the pychopper statistics tables (``pychopper``);
* ``rna.biotype`` and ``rna.genes``: the annotated read tables, through their partials (``readtable``, ``rnacomp``).

``sections.build_read_data`` calls ``add_sections`` once, after the classes of ``cdna`` are collected. Everything
here follows the rule of the rest of read data: one source of record per run, sample and need (``record.choose``
has already decided which files are counted), nothing is summed over two candidates, and a file that could not be
read is named.
"""
from __future__ import annotations

from flowcellar.contracts import ReadData
from flowcellar.readdata import pychopper, readtable, rnacomp
from flowcellar.util.errors import FlowCellarError
from flowcellar.util.text import clean_text

PYCHOPPER_ROLE = "pychopper.stats"
TABLE_ROLE = readtable.ROLE


def _weak(li) -> bool:
    return li.link.get("confidence") in ("medium", "low") and not li.link.get("confirmed")


def _several(li) -> int | None:
    return len(li.run_keys) if len(li.run_keys) > 1 else None


def add_sections(rd: ReadData, ins: list, samples, *, config: dict, barcoded: bool, conflicts: dict[str, list[str]]) -> None:
    """Add the sections to ``rd.modules``. ``ins`` are the inputs of the run as ``sections`` holds them (id, role,
    li, counted, partial, waiting, state, samples, used_by, notes, records); ``samples`` its sample table."""
    from flowcellar.readdata.sections import project_of

    def key_of(x) -> str:
        skey = samples.key_for(None, x.li.sample, barcoded=barcoded, source="file_name", project=project_of(x.li, config))
        x.samples.add(skey)
        return skey

    # ---- pychopper: the small tables are read whole, here
    entries = []
    for x in ins:
        if x.role != PYCHOPPER_ROLE or PYCHOPPER_ROLE not in x.counted:
            continue
        try:
            stats = pychopper.read_stats(x.li.candidate.path)
        except (OSError, FlowCellarError) as e:
            if not x.read_error:
                x.state, x.read_error = "error", str(e)
            continue
        skey = key_of(x)
        x.used_by.add(pychopper.SECTION)
        if any(e["key"] == skey for e in entries):          # two tables that name one sample: the first one stands, never a sum
            continue
        entries.append({"key": skey, "name": samples.name(skey), "input": x.id, "stats": stats, "several": _several(x.li), "weak": _weak(x.li)})
    sec = pychopper.build_section(entries, order=samples.order(), expected=samples.expected_keys(), conflicts=conflicts.get(PYCHOPPER_ROLE, []))
    if sec is not None:
        rd.modules[pychopper.SECTION] = sec

    # ---- annotated read tables: through the partial of the pass over each
    read, waiting = [], []
    for x in ins:
        if x.role != TABLE_ROLE or TABLE_ROLE not in x.counted:
            continue
        size = x.li.candidate.size
        if x.partial is None:
            if x.waiting:
                waiting.append({"input": x.id, "size": size})
            continue
        st = x.partial.get("state")
        if st not in ("ok", "truncated", "partially_read"):
            continue
        skey = key_of(x)
        facts = x.partial.get("facts") or {}
        x.records = facts.get("primary")
        several = _several(x.li)
        if several:
            x.notes.append({"level": "info", "code": "covers_several_runs", "value": several,
                            "message": clean_text(f"This table has no per-read run identifiers and covers the sample on {several} flow cells "
                                                  f"({', '.join(sorted(x.li.run_keys))}); it is attached unsplit and must be counted once.", 600)})
        read.append({"key": skey, "name": samples.name(skey), "input": x.id, "partial": x.partial, "state": st, "size": size,
                     "several": several, "weak": _weak(x.li)})
        x.used_by.update((rnacomp.SECTION_BIOTYPE, rnacomp.SECTION_GENES))
    # two tables of one sample (the same sample under two input rules) would be two entries: never add them up
    seen: dict[str, dict] = {}
    doubled: list[str] = []
    for e in read:
        if e["key"] in seen:
            doubled.append(e["name"])
        else:
            seen[e["key"]] = e
    table_conflicts = list(conflicts.get(TABLE_ROLE, []))
    if doubled:
        table_conflicts.append(f"more than one annotated read table is counted for {', '.join(sorted(set(doubled)))}; only the first one is shown")
    rd.modules.update(rnacomp.build_sections(list(seen.values()), order=samples.order(), expected=samples.expected_keys(),
                                             waiting=waiting, conflicts=table_conflicts))

