"""Helpers of the read-data tests (imported by name: this folder is on the path while its tests run)."""
from __future__ import annotations

import os

from flowcellar.contracts import Candidate, PassOptions, RunIdentity
from flowcellar.readdata import api
from flowcellar.testing import readdata as T


def ident(fr, **facts) -> RunIdentity:
    """The identity of a fixture run folder, as the engine would report it."""
    return RunIdentity(run_key=fr.run_key, flow_cell_id=fr.flow_cell_id, run_short=fr.protocol_run_id[:8], run_dir=fr.run_dir,
                       protocol_run_id=fr.protocol_run_id, started=fr.spec.start, files=dict(fr.files),
                       facts={"experiment_id": fr.spec.experiment, **facts})


def ident_of(ref: T.RunRef, run_dir: str = "") -> RunIdentity:
    return RunIdentity(run_key=ref.run_key, flow_cell_id=ref.flow_cell_id, run_short=ref.run8, run_dir=run_dir,
                       protocol_run_id=ref.protocol_run_id, facts={"experiment_id": ref.experiment})


def cand(role: str, path: str, *, sample=None, variant=None, primary=None, rule="test rule", exp=None, root_id=None, relpath=None) -> Candidate:
    st = os.stat(path)
    return Candidate(role=role, path=path, size=st.st_size, mtime_ns=st.st_mtime_ns, root_id=root_id, relpath=relpath, sample=sample,
                     variant=variant, primary=primary, rule=rule, experiment_dir=exp)


def tier1_of(fr, *, group="cdna", barcoded=True, family="minion") -> dict:
    """The part of a run's measured bundle that read data looks at."""
    return {"identity": {"run_key": fr.run_key, "assay": {"group": group, "barcoded": barcoded, "chemistry": "r10.4.1"},
                         "flow_cell": {"family": family}},
            "samples": [],
            "modules": {"yield": {"summary": {"yield_reads": fr.expected["reads"]}},
                        "run": {"summary": {"run_completion_frac": 1.0, "run_duration_h": 8.0}}}}


def bundle_of(tier1: dict, rd) -> dict:
    """A bundle-shaped dict holding what read_checks, card_columns and sample_rows look at."""
    idn = dict(tier1["identity"])
    idn["basecalling"] = rd.basecalling
    return {"bundle_id": idn["run_key"], "identity": idn, "modules": {**tier1["modules"], **rd.modules}, "samples": rd.samples,
            "tiers": {"read_level": rd.read_level}, "provenance": {"inputs": rd.inputs}}


def do_passes(linked, partials_dir, known=(), **kw) -> dict:
    """read_pass for every input that is the source of record of something; path -> PartialRef."""
    out = {}
    for li in linked:
        if li.run_keys and not li.not_counted_reason and any(li.counted_for.values()):
            ref = api.read_pass(li, str(partials_dir), list(known), PassOptions(**kw))
            if ref is not None:
                out[li.candidate.path] = ref
    return out


def project_candidates(proj: T.ProjectFixture, *, exp=None, roles=None) -> list[Candidate]:
    """The candidates a library scan with a lab-like configuration would hand over for a fixture project."""
    exp = exp or proj.project_dir
    out = [cand("basecaller.sample_sheet", proj.sample_sheet, rule="told by its content", exp=exp)]
    for role, by_alias in proj.files.items():
        if roles and role not in roles:
            continue
        for alias, path in by_alias.items():
            out.append(cand(role, path, sample=alias, rule=f"configured rule {role}", exp=exp))
    return out
