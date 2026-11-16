"""Fake library-preparation tables with known true values: a complete pychopper statistics table, and an annotated
read table (one row per read given to a transcript). Invented names only.

``write_pychopper_stats_full`` writes every category the real tables hold, in their order, and returns what a
correct reader must find. ``write_annotated_table`` writes rows drawn from a seed and returns the true counts, worked
out with plain loops over the same rows (no histogram code of the package is used for the expected values).
"""
from __future__ import annotations

import os
import random
from dataclasses import dataclass, field

DEFAULT_CONFIGS = ("SSP,-VNP", "VNP,-SSP", "SSP,-VNP,-VNP", "VNP,VNP,-SSP", "SSP,SSP,-VNP", "VNP,-SSP,-SSP", "SSP,-SSP", "VNP,-VNP",
                   "SSP", "VNP", "-SSP", "-VNP", "SSP,SSP", "VNP,VNP", "-VNP,-VNP", "SSP,-VNP,SSP,-VNP")
AUTOTUNE_Q = tuple([1e-05] + [round(5.0 * k / 29, 13) for k in range(1, 30)])


@dataclass
class PychopperSpec:
    """The true numbers of one table. Everything else follows from them."""
    primers_found: int = 9000
    plus: int = 4400                               # full-length reads on the + strand
    rescued: dict[int, int] = field(default_factory=lambda: {2: 40, 3: 6, 5: 1})      # pieces per fused read -> reads
    unusable_hits: dict[int, int] = field(default_factory=lambda: {0: 120, 2: 30, 3: 200, 4: 9, 25: 1})   # hits -> unusable reads
    qc_fail: int = 310
    len_fail: int = 55
    umi: tuple[int, int] | None = (8600, 8450)     # Umi_detected, Umi_detected_final; None = UMI detection off
    cutoff_index: int | None = 1                   # which of the tuning cut-offs was chosen; None = no tuning rows
    configs: tuple[str, ...] = DEFAULT_CONFIGS
    seed: int = 5


def pychopper_truth(spec: PychopperSpec) -> dict:
    """Everything a table written from ``spec`` holds, as plain numbers."""
    rnd = random.Random(spec.seed)
    rescued_reads = sum(spec.rescued.values())
    pieces = sum(k * v for k, v in spec.rescued.items())
    unusable = sum(spec.unusable_hits.values())
    pass_reads = spec.primers_found + rescued_reads + unusable
    rescue_plus = pieces // 3
    # rescued PIECES by the primer hits on their read: a read cut into k pieces has 2k hits, a few have one more
    rescue_hits: dict[int, int] = {}
    for k, v in spec.rescued.items():
        extra = v // 5
        rescue_hits[2 * k] = rescue_hits.get(2 * k, 0) + (v - extra) * k
        if extra:
            rescue_hits[2 * k + 1] = rescue_hits.get(2 * k + 1, 0) + extra * k
    # usable reads by the bases that were not kept: most around 75, a few at 0 and far out
    usable = spec.primers_found + rescued_reads
    trimmed: dict[int, int] = {}
    left = usable
    for bases, share in ((0, 0.01), (41, 0.04), (72, 0.30), (73, 0.25), (81, 0.20), (90, 0.10), (299, 0.004), (300, 0.003), (1234, 0.002)):
        n = min(left, int(usable * share))
        if n:
            trimmed[bases] = n
            left -= n
    trimmed[75] = trimmed.get(75, 0) + left
    # reads with at least one hit, by configuration
    with_hits = pass_reads - spec.unusable_hits.get(0, 0)
    configs: dict[str, int] = {}
    left = with_hits
    for i, name in enumerate(spec.configs):
        n = min(left, max(1, int(with_hits * (0.44 if i < 2 else 0.07 * 0.6 ** (i - 1)))))
        if n:
            configs[name] = n
            left -= n
    if left:
        configs[spec.configs[0]] = configs.get(spec.configs[0], 0) + left
    autotune = []
    if spec.cutoff_index is not None:
        for i, q in enumerate(AUTOTUNE_Q):
            autotune.append((q, 9000 - 37 * abs(i - spec.cutoff_index) - rnd.randrange(5) * (i != spec.cutoff_index)))
    sc = {"pass_reads": pass_reads, "qc_fail": spec.qc_fail, "len_fail": spec.len_fail, "primers_found": spec.primers_found,
          "rescued_reads": rescued_reads, "rescue_segments": pieces, "unusable": unusable, "oriented": spec.primers_found + pieces,
          "strand_plus": spec.plus, "strand_minus": spec.primers_found - spec.plus,
          "rescue_strand_plus": rescue_plus, "rescue_strand_minus": pieces - rescue_plus,
          "unusable_no_hit": spec.unusable_hits.get(0, 0), "unusable_multi_hit": sum(v for k, v in spec.unusable_hits.items() if k >= 3),
          "hits_reads": with_hits, "trimmed_reads": usable, "trimmed_bases": sum(k * v for k, v in trimmed.items()),
          "umi_detected": spec.umi[0] if spec.umi else None, "umi_detected_final": spec.umi[1] if spec.umi else None,
          "primer_configs": len(configs),
          "input": pass_reads + spec.qc_fail, "written": spec.primers_found + pieces - spec.len_fail}
    return {"scalars": sc, "cutoff_q": AUTOTUNE_Q[spec.cutoff_index] if spec.cutoff_index is not None else None,
            "segments": dict(sorted(spec.rescued.items())), "unusable_hits": dict(sorted(spec.unusable_hits.items())),
            "rescue_hits": dict(sorted(rescue_hits.items())), "trimmed": dict(sorted(trimmed.items())),
            "configs": sorted(configs.items(), key=lambda kv: (-kv[1], kv[0])), "autotune": autotune}


def write_pychopper_stats_full(path: str, spec: PychopperSpec | None = None, *, whole_numbers: bool = False, crlf: bool = False,
                               extra_category: bool = False) -> dict:
    """Write the table in the layout of the real files (the UMI rows of ``ReadStats`` come after ``Hits``; values are
    written as ``123.0``) and return ``pychopper_truth(spec)``. ``whole_numbers`` writes ``123``; ``extra_category``
    adds a category no reader knows."""
    spec = spec or PychopperSpec()
    t = pychopper_truth(spec)
    sc = t["scalars"]
    rows: list[tuple[str, str, float]] = [
        ("ReadStats", "PassReads", sc["pass_reads"]), ("ReadStats", "LenFail", sc["len_fail"]), ("ReadStats", "QcFail", sc["qc_fail"]),
        ("Classification", "Primers_found", sc["primers_found"]), ("Classification", "Rescue", sc["rescue_segments"]),
        ("Classification", "Unusable", sc["unusable"]),
        ("Strand", "+", sc["strand_plus"]), ("Strand", "-", sc["strand_minus"]),
        ("RescueStrand", "+", sc["rescue_strand_plus"]), ("RescueStrand", "-", sc["rescue_strand_minus"])]
    rows += [("UnclassHitNr", str(k), v) for k, v in t["unusable_hits"].items()]
    rows += [("RescueHitNr", str(k), v) for k, v in t["rescue_hits"].items()]
    rows += [("RescueSegmentNr", str(k), v) for k, v in t["segments"].items()]
    rows += [("Unusable", str(k), v) for k, v in t["trimmed"].items()]
    rows += [("Hits", k, v) for k, v in t["configs"]]
    if spec.umi:
        rows += [("ReadStats", "Umi_detected", spec.umi[0]), ("ReadStats", "Umi_detected_final", spec.umi[1])]
    if extra_category:
        rows += [("FutureCategory", "something", 7), ("FutureCategory", "else", 8)]
    lines = ["Category\tName\tValue"]
    for c, n, v in rows:
        lines.append(f"{c}\t{n}\t{int(v) if whole_numbers else float(v)}")
    for q, n in t["autotune"]:
        lines.append(f"AutotuneSample\t{q!r}\t{float(n)}")
    if t["cutoff_q"] is not None:
        lines.append(f"Parameter\tCutoff(q)\t{t['cutoff_q']!r}")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    nl = "\r\n" if crlf else "\n"
    with open(path, "w", newline="") as fh:
        fh.write(nl.join(lines) + nl)
    return t


# ------------------------------------------------------------------------------------ annotated read table
BIOTYPES = (("protein_coding", 0.74), ("retained_intron", 0.08), ("nonsense_mediated_decay", 0.06), ("lncRNA", 0.05), ("Mt_rRNA", 0.04),
            ("processed_pseudogene", 0.015), ("rRNA", 0.008), ("Mt_tRNA", 0.004), ("snRNA", 0.002), ("TEC", 0.001))
TABLE_COLUMNS = ("read_id", "transcript_id", "gene_id", "gene_name", "length", "transcript_type", "mapping_quality", "FLAG", "ePolyA_tail_length")


def _genes(n_genes: int) -> list[tuple[str, str, str]]:
    """(gene id, gene name, biotype): invented ids; a few genes carry mitochondrial names."""
    out = []
    names = [b for b, _ in BIOTYPES]
    for i in range(n_genes):
        bt = names[0] if i % 3 else names[(i // 3) % len(names)]
        if bt in ("Mt_rRNA", "Mt_tRNA") or i in (1, 7):
            name = f"MT-TEST{i}"
        else:
            name = f"GENE{i}"
        out.append((f"FAKEG{i:011d}.1", name, bt))
    return out


def make_table_rows(n: int, *, seed: int = 3, n_genes: int = 240, tail: str = "dorado", not_primary: int = 0, twice: int = 0) -> list[dict]:
    """``n`` primary rows (and ``not_primary`` secondary or supplementary ones) sorted by transcript, as the real
    tables are. ``tail``: ``"dorado"`` (-1 and 0 for no estimate), ``"na"`` (``NA`` and empty cells for no value),
    ``"mixed"`` (all of them), or ``"none"`` (no tail value at all). ``twice`` reads appear under one id in two rows."""
    rnd = random.Random(seed)
    genes = _genes(n_genes)
    weights = [1.0 / (i + 1) ** 1.1 for i in range(n_genes)]
    rows = []
    for i in range(n + not_primary):
        g = rnd.choices(range(n_genes), weights)[0]
        gid, gname, bt = genes[g]
        tx = f"FAKET{g:011d}.{1 + rnd.randrange(3)}"
        r = rnd.random()
        if tail == "none":
            value = None
        elif r < 0.90:
            base = 22 if bt.startswith("Mt_") else 95
            value = str(max(1, int(rnd.gauss(base, base / 3))) if rnd.random() > 0.004 else 1000 + rnd.randrange(300))
        elif tail == "dorado":
            value = "-1" if r < 0.96 else "0"
        elif tail == "na":
            value = "NA" if r < 0.96 else ""
        else:
            value = ("-1", "0", "NA", "")[int((r - 0.90) * 40) % 4]
        flag = 0 if rnd.random() < 0.97 else 16
        if i >= n:
            flag |= 256 if (i - n) % 2 == 0 else 2048
        rows.append({"read_id": f"{rnd.getrandbits(32):08x}-0000-4000-8000-{i:012x}", "transcript_id": tx, "gene_id": gid, "gene_name": gname,
                     "length": 400 + (g * 37) % 5000, "transcript_type": bt, "mapping_quality": rnd.choice((0, 1, 60)), "FLAG": flag,
                     "ePolyA_tail_length": value})
    primary = [r for r in rows if not r["FLAG"] & 0x900]
    for k in range(min(twice, len(primary) // 2)):
        primary[2 * k + 1]["read_id"] = primary[2 * k]["read_id"]
    rows.sort(key=lambda r: (r["transcript_id"], r["read_id"]))
    return rows


def table_truth(rows: list[dict], *, has_tail_column: bool = True) -> dict:
    """The true counts of a row list, by plain loops: per biotype ``{"reads", "tail_estimated", "tail_zero",
    "tail_not_found", "tail_no_value", "tail_sum", "tails": sorted estimated lengths}``, and the gene numbers."""
    bio: dict[str, dict] = {}
    genes: dict[str, list] = {}
    transcripts = set()
    n_primary = secondary = supplementary = 0
    for r in rows:
        if r["FLAG"] & 0x100:
            secondary += 1
            continue
        if r["FLAG"] & 0x800:
            supplementary += 1
            continue
        n_primary += 1
        b = bio.setdefault(r["transcript_type"], {"reads": 0, "tail_estimated": 0, "tail_zero": 0, "tail_not_found": 0, "tail_no_value": 0,
                                                  "tail_sum": 0, "tails": []})
        b["reads"] += 1
        v = r["ePolyA_tail_length"] if has_tail_column else None
        if v in (None, "", "NA"):
            b["tail_no_value"] += 1
        elif int(v) > 0:
            b["tail_estimated"] += 1
            b["tail_sum"] += int(v)
            b["tails"].append(int(v))
        elif int(v) == 0:
            b["tail_zero"] += 1
        else:
            b["tail_not_found"] += 1
        g = genes.setdefault(r["gene_id"], [0, r["gene_name"]])
        g[0] += 1
        transcripts.add(r["transcript_id"])
    for b in bio.values():
        b["tails"].sort()
    ranked = sorted(genes.values(), key=lambda a: (-a[0], a[1]))
    counts = [a[0] for a in ranked]
    return {"rows": len(rows), "primary": n_primary, "secondary": secondary, "supplementary": supplementary, "biotypes": bio,
            "genes_ge1": len(counts), "genes_ge10": sum(1 for c in counts if c >= 10), "genes_ge100": sum(1 for c in counts if c >= 100),
            "transcripts_ge1": len(transcripts), "reads_top1": sum(counts[:1]), "reads_top10": sum(counts[:10]),
            "reads_top100": sum(counts[:100]), "reads_mito_genes": sum(a[0] for a in ranked if a[1].startswith("MT-")),
            "top_genes": [(a[1], a[0]) for a in ranked[:10]], "gene_counts": counts}


def write_annotated_table(path: str, rows: list[dict], *, tail_column: bool = True, final_newline: bool = True,
                          columns: tuple[str, ...] | None = None) -> dict:
    """Write the rows as an annotated read table and return ``table_truth``. Without ``tail_column`` the table has the
    eight columns of the older tables."""
    cols = list(columns or TABLE_COLUMNS)
    if not tail_column:
        cols = [c for c in cols if c != "ePolyA_tail_length"]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    lines = [",".join(cols)]
    for r in rows:
        lines.append(",".join("" if r[c] is None else str(r[c]) for c in cols))
    with open(path, "w", newline="") as fh:
        fh.write("\n".join(lines) + ("\n" if final_newline else ""))
    return table_truth(rows, has_tail_column=tail_column)


def lower_median(values: list[int]) -> float | None:
    """The lower of the two middle values (the rule the tail quantiles follow)."""
    return float(sorted(values)[(len(values) + 1) // 2 - 1]) if values else None


def lower_quantile(values: list[int], q: float) -> float | None:
    """The value at rank ceil(q x n) of the sorted values."""
    if not values:
        return None
    import math
    return float(sorted(values)[max(1, math.ceil(q * len(values))) - 1])
