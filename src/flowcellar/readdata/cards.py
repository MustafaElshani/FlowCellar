"""What read data adds to the index: the columns of a run's row and the rows of the samples table
. Computed from the bundle alone."""
from __future__ import annotations

from flowcellar.readdata import sections as S

UNCLASSIFIED = "unclassified"


def _summary(bundle: dict, module: str) -> dict:
    return ((bundle.get("modules") or {}).get(module) or {}).get("summary") or {}


def card_columns(bundle: dict) -> dict:
    """``read_data``, ``t2_population``, ``t2_samples``, ``t2_samples_with_data``, ``t2_reads``, ``t2_passed_gb``,
    ``t2_n50``, ``t2_q_med``, ``t2_ident_med``, ``t2_mapped_frac``, ``t2_usable``, ``coverage``, ``t2_files_waiting``,
    ``t2_conflicts`` (files of this run that are not counted because two candidate sources describe the same reads),
    ``t2_full_length`` (full-length share of the reads pychopper classified, over the samples of the run).
    A number that FlowCellar has not counted itself is None ("not available yet"), never an estimate."""
    rl = ((bundle.get("tiers") or {}).get("read_level") or {})
    mods = bundle.get("modules") or {}
    state = rl.get("state", "absent")
    waiting = sum(1 for i in ((bundle.get("provenance") or {}).get("inputs") or [])
                  if any(n.get("code") == "waiting" for n in i.get("notes") or []))
    conflicts = sum(1 for i in ((bundle.get("provenance") or {}).get("inputs") or [])
                    if str(i.get("not_counted_reason") or "").startswith("two candidate sources"))
    out = {"read_data": state, "t2_population": None, "t2_samples": None, "t2_samples_with_data": None, "t2_reads": None,
           "t2_passed_gb": None, "t2_n50": None, "t2_q_med": None, "t2_ident_med": None, "t2_mapped_frac": None, "t2_usable": None,
           "coverage": None, "t2_files_waiting": waiting, "t2_conflicts": conflicts}
    out["t2_full_length"] = None                  # full-length share of the pychopper tables (library preparation)
    if state == "absent":
        return out
    out["t2_population"] = (rl.get("population") or {}).get("stage")
    out["t2_samples"] = rl.get("samples_expected")
    out["t2_samples_with_data"] = rl.get("samples_with_data")
    out["t2_reads"] = rl.get("reads_covered")
    ln = mods.get("len") or {}
    counted = ln.get("state") == "ok" or (ln.get("state") == "partial" and ln.get("partial_reason") == "samples_missing")
    if ln.get("cells") and counted:               # numbers of a partly counted file are not totals: not shown in the hub
        passed = [c for c in ln["cells"] if c["at"][1] == "pass"] or ln["cells"]
        out["t2_passed_gb"] = round(sum((c.get("scalars") or {}).get("bases", 0) for c in passed) / 1e9, 4)
        out["t2_n50"] = (ln.get("summary") or {}).get("len_n50")
        out["t2_q_med"] = _summary(bundle, "qual").get("qual_median")
    al = mods.get("aln") or {}
    if al.get("cells") and (al.get("state") == "ok" or al.get("partial_reason") == "samples_missing"):
        out["t2_ident_med"] = (al.get("summary") or {}).get("aln_identity_median")
        out["t2_mapped_frac"] = (al.get("summary") or {}).get("aln_mapped_frac")
    out["t2_usable"] = _summary(bundle, "funnel").get("funnel_usable_frac")
    cd = mods.get("cdna.pychopper") or mods.get("cdna") or {}
    if cd.get("cells") and (cd.get("state") == "ok" or cd.get("partial_reason") == "samples_missing"):
        out["t2_full_length"] = (cd.get("summary") or {}).get("cdna_full_length_frac")
    exp, cov = rl.get("reads_expected"), rl.get("reads_covered")
    if exp and cov is not None:
        out["coverage"] = round(min(cov / exp, 1.0), 6)
    return out


def sample_rows(bundle: dict) -> list[dict]:
    """One row per sample of the run (not the pseudo-sample of reads without a barcode): ``sample_uid``, ``name``,
    ``alias``, ``barcode``, ``project``, ``run``, ``key``, and ``reads``, ``bases``, ``n50``, ``q_med`` (None until read)."""
    mods = bundle.get("modules") or {}
    run_key = (bundle.get("identity") or {}).get("run_key") or bundle.get("bundle_id")
    ln, ql = mods.get("len") or {}, mods.get("qual") or {}
    # a sample whose file is read has its numbers even while other samples wait; a partly read file gives no total
    ok = ln.get("state") == "ok" or ln.get("partial_reason") in ("samples_missing", "waiting")
    rows = []
    for s in bundle.get("samples") or []:
        if s.get("key") == UNCLASSIFIED:
            continue
        row = {"sample_uid": s.get("sample_uid"), "name": s.get("name"), "alias": s.get("alias"), "barcode": s.get("barcode"),
               "project": s.get("project"), "run": run_key, "key": s.get("key"), "reads": None, "bases": None, "n50": None, "q_med": None}
        if ok:
            mine = [c for c in ln.get("cells") or [] if c["at"][0] == s.get("key")]
            if mine:
                tot = S._total(mine, pass_only=True)
                sc = (tot or {}).get("scalars") or {}
                row["reads"], row["bases"] = sc.get("reads"), sc.get("bases")
                row["n50"] = S.len_summary(tot)["len_n50"]
                q = [c for c in ql.get("cells") or [] if c["at"][0] == s.get("key")]
                row["q_med"] = S.qual_median(S._total(q, pass_only=True)) if q else None
        rows.append(row)
    return rows
