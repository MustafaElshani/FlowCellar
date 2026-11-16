"""An invented laboratory on disk: the example data of the README and of ``docs/quick-start.md``.

    python -m flowcellar.testing.demo OUT_DIR [--scale 0.01] [--jobs 4] [--seed 7]

writes two folders to scan and a library folder that holds only its configuration::

    OUT_DIR/minknow/      what a sequencer writes: <experiment>/<sample>/<run folder>; 15 runs of 13 flow cells in all
    OUT_DIR/projects/     what a pipeline leaves beside them: basecalled and aligned BAM files, pychopper statistics,
                          annotated read tables, quantification metadata, single-cell and spatial analysis folders
    OUT_DIR/library/      flowcellar.toml: the file patterns of this example, and its output limits

    flowcellar scan OUT_DIR/minknow OUT_DIR/projects -L OUT_DIR/library

Everything is invented and derived from the seed: equal arguments give the same runs, reads and names. The files are
real in form (a scan reads them as it reads a laboratory's own) and small in size: with ``--scale 0.01`` every run
holds a hundredth of the reads of a full-size run, and the library's ``flowcellar.toml`` scales the output limits by
the same factor so that the runs are judged as full-size runs would be. Counts, lengths, qualities and the time axis
of a run agree between its MinKNOW folder and its read-data files, because both are written from one table of reads.

The laboratory (four projects, June to October 2026):

* ``demo_cdna_timecourse``: PCR-cDNA with six barcoded samples on MinION flow cells. Eight runs, of which one was
  stopped after half an hour and started again on the same flow cell, and one was loaded with too little library.
  Three runs have read data: one pool that went on two flow cells, whose reads were merged per sample before
  pychopper and alignment, so that those files cover both flow cells; and the last pool on one flow cell, with every
  kind of file (basecalled reads, alignments, pychopper statistics, annotated read tables, salmon).
* ``demo_genome``: genomic DNA, one sample per flow cell; the first flow cell was washed and used a second time.
  The newest run has two basecalling passes (one with a modified-base model) and an alignment.
* ``demo_direct_rna``: direct RNA on an RNA flow cell, with an alignment to the transcripts.
* ``demo_spatial``: cDNA of a Visium and a Visium HD slide on PromethION flow cells, each analysed with Space Ranger
  and with wf-single-cell.

BAM files need ``pysam`` and the analysis folders ``h5py``; without one of them the files that need it are left out
and the example says so.
"""
from __future__ import annotations

import json
import os
import sys
import time
from dataclasses import dataclass, field
from functools import cache, lru_cache

import numpy as np

from flowcellar.testing import fixtures as fx

MODEL_DNA = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
MODEL_MODS = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2"
MODEL_RNA = "rna004_130bps_sup@v5.2.0"
KIT_CDNA = "SQK-PCB114-24"
ORIGIN_S = 252.0                     # the fixture run folders start reading data this long after the protocol starts
ID_CHUNK = 1 << 18                   # read ids are made in chunks, so that a file can ask for a slice of them
#: an invented genome of seven chromosomes and a mitochondrial genome, 114 Mb
GENOME = (("chr1", 21_400_000), ("chr2", 18_900_000), ("chr3", 17_200_000), ("chr4", 15_600_000), ("chr5", 13_100_000),
          ("chr6", 11_800_000), ("chrX", 16_300_000), ("chrM", 16_600))
BIOTYPES = (("protein_coding", 0.80), ("lncRNA", 0.09), ("retained_intron", 0.04), ("nonsense_mediated_decay", 0.03),
            ("processed_pseudogene", 0.025), ("snoRNA", 0.01), ("TEC", 0.005))
SAMPLES = ("ctrl_1", "ctrl_2", "ctrl_3", "treated_1", "treated_2", "treated_3")           # the pool that went on two flow cells
SAMPLES_F = ("ctrl_4", "ctrl_5", "ctrl_6", "treated_4", "treated_5", "treated_6")          # the pool after it: one flow cell
#: output limits of full-size runs that the example scales (check id -> assay group -> (warn, fail)); the catalogue's own
FULL_SIZE_LIMITS = {
    "yield_reads_at_72h": {"cdna": (5_000_000, 2_500_000), "drna": (800_000, 400_000)},
    "yield_estimated_bases_at_72h": {"cdna": (6_000_000_000, 3_000_000_000), "gdna": (10_000_000_000, 5_000_000_000),
                                     "drna": (1_500_000_000, 700_000_000)},
    "yield_basecalled_bases": {"cdna": (4_000_000_000, 2_000_000_000), "gdna": (8_000_000_000, 4_000_000_000),
                               "drna": (1_000_000_000, 500_000_000)},
}


# --------------------------------------------------------------------------------------------- the runs
@dataclass(frozen=True)
class DemoRun:
    """One invented run. ``reads`` is the number of a full-size run; ``--scale`` multiplies it."""
    label: str
    flow_cell_id: str
    assay: str                              # cdna, gdna, drna
    experiment: str
    sample: str
    start: str
    reads: int
    seed: int
    hours: float = 72.2                     # MinKNOW goes on for some minutes after the planned time
    planned_hours: float | None = 72.0
    family: str = "minion"
    kit: str = KIT_CDNA
    root: str = "minknow"
    state: str = "PROTOCOL_COMPLETED"
    pore_frac: float = 0.74                 # wells with a pore at the first scan
    occupancy: float | None = None          # share of working pores busy with a strand at the start; None = drawn for the run
    tau_h: float = 95.0                     # how long pores last
    aliases: tuple[str, ...] = ()           # the samples of a barcoded run, barcode01 onwards
    read_data: bool = False
    dead_block: float = 0.0                 # share of channels in one block that never gave a read

    @property
    def model(self) -> str:
        return MODEL_RNA if self.assay == "drna" else MODEL_DNA


def demo_runs() -> list[DemoRun]:
    cd = dict(assay="cdna", experiment="demo_cdna_timecourse", aliases=SAMPLES)
    return [
        DemoRun("c1", "TST00101", sample="pool_a", start="2026-06-09T09:40:00Z", reads=7_300_000, seed=101, root="projects", **cd),
        DemoRun("c2", "TST00102", sample="pool_b", start="2026-06-23T10:05:00Z", reads=6_600_000, seed=102, pore_frac=0.70, root="projects", **cd),
        DemoRun("c3a", "TST00103", sample="pool_c", start="2026-07-07T09:15:00Z", reads=7_000_000, seed=103, hours=0.6,
                state="PROTOCOL_FINISHED_WITH_ERROR", **cd),
        DemoRun("c3b", "TST00103", sample="pool_c_restart", start="2026-07-07T10:20:00Z", reads=7_000_000, seed=104, **cd),
        DemoRun("c4", "TST00104", sample="pool_d", start="2026-07-28T11:30:00Z", reads=4_700_000, seed=105, occupancy=0.61, **cd),
        DemoRun("c5", "TST00105", sample="pool_e", start="2026-08-18T09:55:00Z", reads=7_800_000, seed=106, pore_frac=0.78, read_data=True, **cd),
        DemoRun("c6", "TST00106", sample="pool_e_fc2", start="2026-08-25T10:10:00Z", reads=6_400_000, seed=107, pore_frac=0.71, read_data=True, **cd),
        DemoRun("c7", "TST00107", sample="pool_f", start="2026-09-15T09:30:00Z", reads=7_100_000, seed=108, read_data=True,
                **{**cd, "aliases": SAMPLES_F}),
        DemoRun("g1", "TST00201", assay="gdna", experiment="demo_genome", sample="isolate_1", start="2026-08-04T13:20:00Z", reads=1_150_000,
                seed=201, hours=48.2, planned_hours=48.0, kit="SQK-LSK114", tau_h=62.0),
        DemoRun("g2", "TST00201", assay="gdna", experiment="demo_genome", sample="isolate_2", start="2026-08-11T09:45:00Z", reads=520_000,
                seed=202, hours=48.2, planned_hours=48.0, kit="SQK-LSK114", tau_h=50.0, pore_frac=0.31, occupancy=0.90),
        DemoRun("g3", "TST00202", assay="gdna", experiment="demo_genome", sample="isolate_3", start="2026-10-13T10:00:00Z", reads=1_420_000,
                seed=203, kit="SQK-LSK114", tau_h=66.0, read_data=True),
        DemoRun("r1", "TST00301", assay="drna", experiment="demo_direct_rna", sample="total_rna_1", start="2026-09-22T14:10:00Z",
                reads=1_650_000, seed=301, kit="SQK-RNA004", tau_h=60.0, occupancy=0.90, read_data=True),
        DemoRun("s1", "TST00401", assay="cdna", experiment="demo_spatial", sample="visium_a", start="2026-10-20T09:00:00Z", reads=74_000_000,
                seed=401, family="promethion", kit="SQK-PCS114", pore_frac=0.71, tau_h=105.0),
        DemoRun("s2", "TST00402", assay="cdna", experiment="demo_spatial", sample="visium_hd_b", start="2026-10-27T09:20:00Z", reads=81_000_000,
                seed=402, family="promethion", kit="SQK-PCS114", pore_frac=0.74, tau_h=110.0, dead_block=0.03),
        DemoRun("s3", "TST00403", assay="cdna", experiment="demo_spatial", sample="visium_hd_b_fc2", start="2026-10-27T09:35:00Z",
                reads=47_000_000, seed=403, family="promethion", kit="SQK-PCS114", pore_frac=0.385, tau_h=90.0),
    ]


def character(run: DemoRun) -> dict:
    """What differs from one run to the next without anybody choosing it: how well the library loaded, how long the
    pores of this flow cell lasted, how long its molecules were, how good its reads."""
    rng = np.random.default_rng(run.seed + 3000)
    return {"occupancy": run.occupancy if run.occupancy is not None else float(rng.uniform(0.885, 0.965)),
            "tau_h": run.tau_h * float(rng.uniform(0.82, 1.2)), "length": float(rng.uniform(0.9, 1.1)), "quality": float(rng.normal(0.0, 0.45)),
            "decline": float(rng.uniform(0.02, 0.07))}


def _fixture_spec(run: DemoRun) -> fx.RunSpec:
    extra = {}
    if run.assay == "drna":
        extra = dict(product="FLO-MIN004RA", protocol_script="sequencing_MIN004RA_RNA_e8_2_130K", speed_bps=130)
    if run.family == "promethion":
        extra = dict(position="P2S-00000-B" if run.label == "s3" else "P2S-00000-A")
    return fx.RunSpec(flow_cell_id=run.flow_cell_id, family=run.family, kit=run.kit, experiment=run.experiment, sample=run.sample,
                      start=run.start, hours=run.hours, planned_hours=run.planned_hours, reads=40, seed=run.seed,
                      first_scan_pore_frac=run.pore_frac, occupancy=character(run)["occupancy"], state=run.state, host="demohost-01", **extra)


@cache
def _ids_of(run: DemoRun) -> tuple[str, str]:
    """(protocol run id, acquisition id) exactly as ``fixtures.make_run`` draws them from the seed."""
    rng = np.random.default_rng(run.seed)
    return fx._uuid(rng), fx._hex(rng, 40)


def run_key(run: DemoRun) -> str:
    return f"{run.flow_cell_id}_{_ids_of(run)[0][:8]}"


def _epoch(run: DemoRun) -> float:
    return fx._epoch(run.start)


# ------------------------------------------------------------------------------------- the flow cell over time
def flow_cell_model(run: DemoRun) -> dict:
    """Shares of the channels in each state per minute, and the wells at every pore scan: pores that die off faster
    late in the run, a step back up at each scan when another well takes over, occupancy that sinks slowly."""
    rng = np.random.default_rng(run.seed + 5000)
    channels = fx._FAMILIES[run.family][0]
    minutes = int(round(run.hours * 60))
    t_h = np.arange(minutes) / 60.0
    period = 1.5
    n_scans = int(run.hours / period) + 1
    scan_h = np.arange(n_scans) * period

    me = character(run)
    tau = me["tau_h"]

    def survival(h):
        return np.exp(-np.power(np.maximum(h, 0.0) / tau, 1.35))
    p0 = run.pore_frac
    # ---- wells at each scan: a pore dies for good; a share of the living ones is held in reserve
    total_wells = channels * 4
    state = rng.choice(len(fx.SCAN_CATEGORIES), size=total_wells, p=fx._scan_p(p0))
    life_h = tau * np.power(-np.log(np.maximum(rng.random(total_wells), 1e-12)), 1 / 1.35)
    scans = []
    for k in range(n_scans):
        alive = (state <= 1) & (life_h > scan_h[k])
        st_k = np.where(alive, 0, np.where(state <= 1, rng.choice([2, 3, 4, 6], size=total_wells, p=[0.25, 0.2, 0.45, 0.1]), state))
        # per channel the first living well sequences; the others wait (MinKNOW lists a growing share as reserve)
        grid = (st_k == 0).reshape(channels, 4)
        extra = grid & (np.cumsum(grid, axis=1) > 1)
        reserve = extra & (rng.random((channels, 4)) < min(0.9, 0.08 + 0.5 * scan_h[k] / max(run.hours, 1.0)))
        st_k = st_k.copy()
        st_k[reserve.ravel()] = 1
        scans.append(st_k.astype(np.int8))
    with_pore = np.array([((s <= 1).reshape(channels, 4)).any(axis=1).mean() for s in scans])
    # ---- per minute: channels that hold a working pore step up at each scan and sink until the next
    k_of = np.minimum((t_h / period).astype(int), n_scans - 1)
    u = (t_h - scan_h[k_of]) / period
    lost = 1.0 - survival(t_h) / survival(0)
    active = with_pore[k_of] * (1.0 - (0.005 + 0.04 * lost) * u)
    wobble = np.convolve(rng.normal(0, 1, minutes + 60), np.ones(60) / 60, mode="valid")[:minutes]
    occ = np.clip(me["occupancy"] - me["decline"] * (t_h / max(run.hours, 1.0)) + 0.012 * wobble, 0.05, 0.995)
    if me["occupancy"] < 0.75:                                   # too little library: occupancy sinks further as it is used up
        occ = np.clip(occ - 0.10 * (t_h / max(run.hours, 1.0)), 0.05, 0.995)
    recovering = 0.022 + 0.035 * lost + 0.004 * np.abs(wobble)
    adapter = {"cdna": 0.09, "gdna": 0.025, "drna": 0.07}[run.assay]
    shares = np.zeros((len(fx.CHANNEL_STATES), minutes))
    working = np.clip(active - recovering - 0.004, 0.0, None)
    seq = working * occ
    idx = {s: i for i, s in enumerate(fx.CHANNEL_STATES)}
    shares[idx["strand"]] = seq * (1 - adapter)
    shares[idx["adapter"]] = seq * adapter
    shares[idx["pore"]] = working - seq
    shares[idx["unavailable"]] = recovering
    shares[idx["unblocking"]] = 0.004
    shares[idx["zero"]] = 0.006 + 0.05 * lost
    shares[idx["saturated"]] = 0.003 + 0.02 * lost
    shares[idx["multiple"]] = 0.006
    shares[idx["unclassified"]] = 0.002
    shares[idx["no_pore"]] = np.clip(1.0 - shares.sum(0), 0.0, None)
    return {"channels": channels, "minutes": minutes, "n_scans": n_scans, "scan_s": (scan_h * 3600).astype(int), "scans": scans,
            "shares": shares, "seq_share": seq, "period": period}


# ------------------------------------------------------------------------------------------ the reference
@lru_cache(maxsize=2)
def reference(seed: int = 5) -> dict:
    """Invented genes and transcripts on the invented genome: where each gene lies, how strongly it is expressed, the
    length and type of each transcript. Transcript names follow the layout of a reference of transcripts
    (``transcript|gene|-|-|transcript name|gene name|length|type|``)."""
    rng = np.random.default_rng(seed)
    n_genes, n_mito = 9000, 13
    nuclear = np.array([ln for name, ln in GENOME if name != "chrM"], dtype=np.float64)
    chrom = rng.choice(len(nuclear), size=n_genes, p=nuclear / nuclear.sum())
    chrom[:n_mito] = len(GENOME) - 1
    span = np.clip(rng.lognormal(np.log(9000), 0.9, n_genes), 1500, 120_000).astype(np.int64)
    span[:n_mito] = 1100
    start = np.array([rng.integers(10_000, GENOME[c][1] - span[i] - 10_000) if c != len(GENOME) - 1 else 300 + 1200 * i
                      for i, c in enumerate(chrom)], dtype=np.int64)
    names = [b for b, _ in BIOTYPES]
    bt = rng.choice(len(names), size=n_genes, p=np.array([p for _, p in BIOTYPES]) / sum(p for _, p in BIOTYPES))
    biotype = [names[i] for i in bt]
    gname = [f"g{i:04d}" for i in range(n_genes)]
    expr = rng.lognormal(0.0, 2.0, n_genes)
    for i in range(n_mito):                                      # mitochondrial genes: named MT-, two of them rRNA
        gname[i] = f"MT-g{i + 1:02d}"
        biotype[i] = "Mt_rRNA" if i < 2 else "protein_coding"
    expr[:n_mito] = np.sort(expr)[-200:].mean() * np.array([2.2, 1.6] + [0.7] * (n_mito - 2))
    for i in range(n_mito, n_mito + 3):                          # three nuclear rRNA genes
        biotype[i] = "rRNA"
    expr[n_mito:n_mito + 3] = np.sort(expr)[-200:].mean() * 0.5
    n_tx = 1 + np.minimum(rng.poisson(1.3, n_genes), 4)
    n_tx[:n_mito] = 1
    gene_of = np.repeat(np.arange(n_genes), n_tx)
    k_in_gene = np.concatenate([np.arange(k) for k in n_tx])
    tlen = np.clip(rng.lognormal(np.log(1250), 0.55, len(gene_of)), 280, 9000).astype(np.int64)
    tlen[gene_of < n_mito] = rng.integers(500, 1600, int((gene_of < n_mito).sum()))
    tlen = np.minimum(tlen, span[gene_of] - 100)
    w_tx = expr[gene_of] * np.power(0.32, k_in_gene) * rng.uniform(0.6, 1.4, len(gene_of))
    tx_type = [biotype[g] if k == 0 or biotype[g] != "protein_coding" else ("protein_coding", "retained_intron", "nonsense_mediated_decay",
                                                                           "protein_coding", "protein_coding")[k]
               for g, k in zip(gene_of, k_in_gene)]
    tx_id = [f"DEMOT{i:08d}.1" for i in range(len(gene_of))]
    gene_id = [f"DEMOG{g:08d}.1" for g in range(n_genes)]
    tx_name = [f"{tx_id[i]}|{gene_id[g]}|-|-|{gname[g]}-{201 + k}|{gname[g]}|{tlen[i]}|{tx_type[i]}|"
               for i, (g, k) in enumerate(zip(gene_of, k_in_gene))]
    return {"gene_chrom": chrom, "gene_start": start, "gene_span": span, "gene_strand": rng.integers(0, 2, n_genes), "gene_name": gname,
            "gene_id": gene_id, "gene_of": gene_of, "tlen": tlen, "tx_weight": w_tx / w_tx.sum(), "tx_type": tx_type, "tx_id": tx_id,
            "tx_name": tx_name, "n_tx": n_tx, "mito_tx": gene_of < n_mito}


# ---------------------------------------------------------------------------------------- the reads of a run
def n_reads(run: DemoRun, scale: float) -> int:
    n = run.reads * scale
    if run.hours < 2:                                             # the run that was stopped: what half an hour gives
        n *= run.hours / 72.0 * 1.6
    return max(200, int(round(n)))


def reads_of(run: DemoRun, scale: float) -> dict:
    """The reads of one run as arrays: when each started, how long it took, its channel, how it ended, its length and
    quality, the sample it was given to, the transcript it came from and what became of it. Both the MinKNOW files
    and the read-data files of the run are written from this one table."""
    rng = np.random.default_rng(run.seed + 9000)
    fc = flow_cell_model(run)
    me = character(run)
    n = n_reads(run, scale)
    channels, minutes = fc["channels"], fc["minutes"]
    run_s = minutes * 60
    t_frac = np.arange(minutes) / max(minutes - 1, 1)
    v0 = 131.0 if run.assay == "drna" else 412.0
    speed_min = v0 * (1.0 - 0.055 * t_frac)
    rate = fc["seq_share"] * speed_min
    cdf = np.cumsum(rate) / rate.sum()
    start = np.sort((np.searchsorted(cdf, rng.random(n)) + rng.random(n)) * 60.0)
    start = np.minimum(start, run_s - 30.0)
    tt = start / run_s
    out: dict = {"n": n, "start_s": start}
    # ---- which sample, and a quality that sinks a little during the run
    q_top = {"cdna": 25.8, "gdna": 26.8, "drna": 23.0}[run.assay]
    q = q_top + me["quality"] - rng.gamma(2.4, 2.35, n) - 0.8 * tt + rng.normal(0.0, 1.1, n)
    low = rng.random(n) < {"cdna": 0.055, "gdna": 0.035, "drna": 0.09}[run.assay]
    q[low] = rng.normal(7.4, 1.4, int(low.sum()))
    alias = np.full(n, -1, dtype=np.int16)
    if run.aliases:
        share = np.array([0.19, 0.165, 0.15, 0.18, 0.16, 0.155][:len(run.aliases)])
        share = share * rng.uniform(0.9, 1.1, len(share))
        unclassified = np.where(low, 0.45, 0.062)
        alias = rng.choice(len(share), size=n, p=share / share.sum()).astype(np.int16)
        lost = rng.random(n) < unclassified
        alias[lost] = -1
        q[lost & ~low] -= 3.0
    out["alias"] = alias
    # ---- what the molecule was
    if run.assay == "gdna":
        length = np.clip(rng.lognormal(np.log(6400 * me["length"]), 0.95, n), 150, 240_000)
        length[low] *= 0.45
        tx = part = tail = None
        extra = np.zeros(n)
    else:
        ref = reference()
        tx = np.searchsorted(np.cumsum(ref["tx_weight"]), rng.random(n)).clip(0, len(ref["tlen"]) - 1)
        tlen = ref["tlen"][tx]
        p_full = (0.74 if run.assay == "cdna" else 0.62) * np.exp(-tlen / (6500.0 if run.assay == "cdna" else 5200.0))
        whole = rng.random(n) < p_full
        part = np.where(whole, tlen - rng.integers(0, 12, n), (tlen * rng.uniform(0.12, 0.98, n) * min(1.0, me["length"])).astype(np.int64))
        part = np.maximum(part, 60)
        mito = ref["mito_tx"][tx]
        tail = np.clip(rng.gamma(4.2, 23.0 if run.assay == "cdna" else 27.0, n), 8, 420)
        tail[mito] = np.clip(rng.gamma(3.0, 13.0, int(mito.sum())), 5, 160)
        tail = tail.astype(np.int64)
        extra = np.full(n, 138.0) if run.assay == "cdna" else np.full(n, 62.0)       # primers and adapters the basecaller leaves on
        length = part + tail + extra
    length = np.maximum(length, 40).astype(np.int64)
    out.update(q=np.clip(q, 2.5, 45.0).astype(np.float32), length=length, tx=tx, part=part, tail=tail)
    # ---- how it went through the pore
    reason = rng.choice(4, size=n, p=(0.948, 0.022, 0.018, 0.012))
    speed = np.clip(v0 * (1.0 - 0.055 * tt) * rng.normal(1.0, 0.045, n), v0 * 0.55, v0 * 1.3)
    est = np.maximum(np.rint(length * rng.normal(1.0, 0.035, n)), 20).astype(np.int64)
    est[reason == 2] = np.maximum(est[reason == 2] // 6, 20)
    ratio = 0.6
    events = np.maximum(np.rint(est / ratio), 1).astype(np.int64)
    duration = np.round(events * ratio / speed + rng.uniform(0.05, 0.4, n), 4)
    duration = np.minimum(duration, run_s - start - 0.5).clip(0.2)
    # how much each channel gives belongs to the flow cell: a second run on it finds the same channels dead, and more
    cell = np.random.default_rng(int(run.flow_cell_id[3:]))
    w = cell.gamma(4.0, 1.0, channels)
    w[cell.random(channels) < 0.035] = 0.0
    if run.pore_frac < 0.5:
        w[rng.random(channels) < 0.3 * (0.5 - run.pore_frac) / 0.2] = 0.0
    if run.dead_block:
        a = int(channels * 0.4)
        w[a:a + int(channels * run.dead_block)] = 0.0
    out.update(reason=reason.astype(np.int8), speed=speed, events=events, duration=duration, ratio=ratio,
               channel=rng.choice(channels, size=n, p=w / w.sum()).astype(np.int32) + 1,
               mux=rng.choice([1, 2, 3, 4], size=n, p=[0.42, 0.28, 0.19, 0.11]).astype(np.int8))
    # ---- what became of it after sequencing
    called = (reason != 2) & (rng.random(n) < 0.992)                # the basecaller wrote it
    passed = called & (out["q"] >= 10.0)
    r = rng.random(n)
    if run.assay == "cdna":
        # pychopper on the reads of a sample that passed: 0 not given to it, 1 primers found, 2 rescued, 3 unusable,
        # 4 below its quality limit, 5 too short after trimming
        fate = np.zeros(n, dtype=np.int8)
        given = passed & (alias >= 0)
        fate[given] = np.select([r[given] < 0.004, r[given] < 0.004 + 0.865, r[given] < 0.004 + 0.865 + 0.009], [4, 1, 2], default=3)
        fate[(fate == 1) & (part < 90)] = 5
        out["fate"] = fate
        usable = (fate == 1) | (fate == 2)
    else:
        usable = passed
    g = rng.random(n)
    out.update(called=called, passed=passed, usable=usable,
               on_genome=usable & (g < {"cdna": 0.987, "gdna": 0.991, "drna": 0.975}[run.assay]),
               on_tx=usable & (rng.random(n) < {"cdna": 0.945, "gdna": 0.0, "drna": 0.955}[run.assay]))
    return out


def read_ids(run: DemoRun, lo: int, hi: int) -> np.ndarray:
    """The ids of reads ``lo`` to ``hi`` of a run as a bytes array, without making the ids of the whole run."""
    parts = []
    for c in range(lo // ID_CHUNK, (max(hi, lo + 1) - 1) // ID_CHUNK + 1):
        rng = np.random.default_rng([run.seed, 77, c])
        hexc = np.frombuffer(b"0123456789abcdef", np.uint8)[rng.integers(0, 16, size=(ID_CHUNK, 32), dtype=np.uint8)]
        b = np.empty((ID_CHUNK, 36), np.uint8)
        b[:, [8, 13, 18, 23]] = 45
        b[:, 14], b[:, 19] = 52, 97
        b[:, 0:8], b[:, 9:13], b[:, 15:18], b[:, 20:23], b[:, 24:36] = hexc[:, 0:8], hexc[:, 8:12], hexc[:, 13:16], hexc[:, 17:20], hexc[:, 20:32]
        a, z = max(lo, c * ID_CHUNK) - c * ID_CHUNK, min(hi, (c + 1) * ID_CHUNK) - c * ID_CHUNK
        parts.append(b[a:z].copy().view("S36").ravel())
    return np.concatenate(parts) if parts else np.empty(0, "S36")


# --------------------------------------------------------------------------------------- the MinKNOW folder
def _csv(path: str, columns: dict, *, delimiter: str = ",") -> None:
    import pyarrow as pa
    import pyarrow.csv as pacsv
    # the header is written by hand: MinKNOW does not quote its column names, the CSV writer always does
    with open(path, "wb") as fh:
        fh.write((delimiter.join(columns) + "\n").encode("ascii"))
        pacsv.write_csv(pa.table(columns), fh, write_options=pacsv.WriteOptions(delimiter=delimiter, quoting_style="none", include_header=False))


def write_run_folder(root: str, run: DemoRun, scale: float) -> fx.FixtureRun:
    """One MinKNOW run folder. The small files and the layout come from ``fixtures.make_run``; the tables that carry
    the shape of the run (reads, output per minute, channel states, pore scans) are then written from the model, and
    the report is made to agree with them."""
    import pyarrow as pa
    import pyarrow.csv as pacsv
    fr = fx.make_run(root, _fixture_spec(run))
    fc, rd = flow_cell_model(run), reads_of(run, scale)
    spec = fr.spec
    channels, minutes, n = fc["channels"], fc["minutes"], rd["n"]
    fcid, prot_id, acq_id = run.flow_cell_id, fr.protocol_run_id, fr.acquisition_id
    run8, acq8 = prot_id[:8], acq_id[:8]
    # ---- reads
    path = fr.files["minknow.sequencing_summary"]
    names = ["filename_fastq", "filename_fast5", "filename_pod5", "filename_bam", "parent_read_id", "read_id", "run_id", "channel", "mux",
             "minknow_events", "start_time", "duration", "pore_type", "experiment_id", "sample_id", "end_reason"]
    reasons = pa.array(fx.END_REASONS)
    schema = pa.schema([(c, pa.string()) for c in names[:7]] + [("channel", pa.int32()), ("mux", pa.int8()), ("minknow_events", pa.int64()),
                                                               ("start_time", pa.float64()), ("duration", pa.float64())]
                       + [(c, pa.string()) for c in names[12:]])
    with open(path, "wb") as sink, pacsv.CSVWriter(sink, schema, write_options=pacsv.WriteOptions(
            delimiter="\t", quoting_style="none", include_header=False)) as w:
        sink.write(("\t".join(names) + "\n").encode("ascii"))
        for lo in range(0, n, 1 << 20):
            hi = min(n, lo + (1 << 20))
            m = hi - lo
            ids = pa.array(read_ids(run, lo, hi)).cast(pa.string())
            dash = pa.repeat(pa.scalar("-"), m)
            pod5 = pa.array([f"{fcid}_{run8}_{acq8}_{k}.pod5" for k in range(lo // 4000, (hi - 1) // 4000 + 1)]).take(
                pa.array(np.arange(lo, hi) // 4000 - lo // 4000))
            w.write_table(pa.table([dash, dash, pod5, dash, ids, ids, pa.repeat(pa.scalar(prot_id), m), pa.array(rd["channel"][lo:hi]),
                                    pa.array(rd["mux"][lo:hi]), pa.array(rd["events"][lo:hi]), pa.array(np.round(rd["start_s"][lo:hi], 6)),
                                    pa.array(rd["duration"][lo:hi]), pa.repeat(pa.scalar("not_set"), m), pa.repeat(pa.scalar(run.experiment), m),
                                    pa.repeat(pa.scalar(run.sample), m), reasons.take(pa.array(rd["reason"][lo:hi]))], names=names))
    # ---- output per minute, cumulative
    end_min = np.minimum(((rd["start_s"] + rd["duration"]) // 60).astype(np.int64), minutes - 1)
    cum_r = np.cumsum(np.bincount(end_min, minlength=minutes))
    cum_e = np.cumsum(np.bincount(end_min, weights=rd["events"], minlength=minutes).astype(np.int64))
    cum_s = np.cumsum(np.bincount(end_min, weights=np.rint(rd["duration"] * spec.sample_rate), minlength=minutes).astype(np.int64))
    cum_b = np.rint(cum_e * rd["ratio"]).astype(np.int64)
    zero = np.zeros(minutes, np.int64)
    _csv(fr.files["minknow.throughput"], {
        "Experiment Time (minutes)": np.arange(minutes), "Reads": cum_r, "Basecalled Reads Passed": zero, "Basecalled Reads Failed": zero,
        "Basecalled Reads Skipped": zero, "Selected Raw Samples": cum_s, "Selected Events": cum_e, "Estimated Bases": cum_b,
        "Basecalled Bases": zero, "Basecalled Samples": zero})
    # ---- channel states per minute: every channel gives 60 x sample rate samples to every minute
    per_min = channels * 60 * spec.sample_rate
    samples = np.zeros((len(fx.CHANNEL_STATES), minutes), np.int64)
    for m in range(minutes):
        samples[:, m] = fx._largest_remainder(fc["shares"][:, m], per_min)
    order = np.argsort(fx.CHANNEL_STATES, kind="stable")
    _csv(fr.files["minknow.pore_activity"], {
        "Channel State": np.repeat(np.array(fx.CHANNEL_STATES)[order], minutes), "Experiment Time (minutes)": np.tile(np.arange(minutes), len(order)),
        "State Time (samples)": samples[order].ravel()})
    # ---- pore scans: one row per channel, well and scan
    rng = np.random.default_rng(run.seed + 6000)
    n_scans, wells = fc["n_scans"], 4
    total = channels * wells
    st = np.concatenate(fc["scans"])
    k = np.repeat(np.arange(n_scans), total)
    good = st <= 1
    med = np.where(good, np.round(rng.normal(195.0, 9.0, len(st)) - k * 0.4, 3), np.nan)
    sd = np.where(good, np.round(np.abs(rng.normal(2.2, 0.4, len(st))), 3), np.nan)
    volt = -180 - 5 * (k // 4)
    _csv(fr.files["minknow.pore_scan"], {
        "channel": np.tile(np.repeat(np.arange(1, channels + 1), wells), n_scans), "well": np.tile(np.tile(np.arange(1, wells + 1), channels), n_scans),
        "flow_cell_id": pa.repeat(pa.scalar(fcid), len(st)), "mux_scan_assessment": pa.array(fx.SCAN_CATEGORIES).take(pa.array(st)),
        "pore_median": pa.array(med, mask=np.isnan(med)), "pore_sd": pa.array(sd, mask=np.isnan(sd)), "scan_number": k + 1,
        "scan_offset_voltage_mv": 10 + k, "scan_running_voltage_mv": volt, "seconds_since_start_of_run": fc["scan_s"][k] + 3,
        "total_good_samples": np.where(good, 50000, 0)})
    # ---- the report: counts of the scans, the output, the messages, temperature and voltage
    with open(fr.files["minknow.report_json"], encoding="utf-8") as fh:
        rep = json.load(fh)
    counts = [{c: int(v) for c, v in zip(fx.SCAN_CATEGORIES, np.bincount(s, minlength=len(fx.SCAN_CATEGORIES)))} for s in fc["scans"]]
    pores = [c["single_pore"] + c["reserved_pore"] for c in counts]
    ch_pore = [int(((s <= 1).reshape(channels, wells)).any(axis=1).sum()) for s in fc["scans"]]
    t_start = _epoch(run)
    t0 = t_start + ORIGIN_S
    msgs = [m for m in rep["user_messages"] if m.get("identifier") not in ("sequencing.mux_scan", "mux_scan_result")]
    for i in range(n_scans):
        at = float(fc["scan_s"][i])
        msgs.append({"time": fx._iso(t0 + at - 200), "severity": "MESSAGE_SEVERITY_INFO", "identifier": "sequencing.mux_scan",
                     "user_message": "Performing Pore Scan"})
        msgs.append({"time": fx._iso(t0 + at), "severity": "MESSAGE_SEVERITY_INFO", "identifier": "mux_scan_result",
                     "user_message": f"Pore scan for flow cell {fcid} has found a total of {pores[i]} pores. "
                                     f"{ch_pore[i]} pores available for sequencing immediately",
                     "extra_data": {"flow_cell_id": fcid, "num_pores": str(ch_pore[i]), "total_pores": str(pores[i])}})
    if run.state != "PROTOCOL_COMPLETED":
        msgs.append({"time": fx._iso(t0 + minutes * 60 - 5), "severity": "MESSAGE_SEVERITY_WARNING", "identifier": "protocol_stopped_by_user",
                     "user_message": "The protocol was stopped by the user"})
    rep["user_messages"] = sorted(msgs, key=lambda m: m["time"])
    seq = rep["acquisitions"][1]
    ari = seq["acquisition_run_info"]
    ari["yield_summary"] = {"read_count": str(n), "selected_raw_samples": str(int(cum_s[-1])), "selected_events": str(int(cum_e[-1])),
                            "estimated_selected_bases": str(int(cum_b[-1]))}
    ari["writer_summary"] = {"bytes_to_write_produced": str(int(cum_s[-1]) * 2), "bytes_to_write_completed": str(int(cum_s[-1]) * 2)}
    ari["bream_info"]["mux_scan_results"] = [{"counts": counts[i], "mux_scan_timestamp": str(int(fc["scan_s"][i]))} for i in range(n_scans)]
    volts = -180 - 5 * (np.arange(n_scans) // 4)
    seq["bias_voltage"] = [{"bias_voltages": [{}] + [{"acquisition_index": str(int(fc["scan_s"][i]) * spec.sample_rate), "bias_voltage": int(volts[i]),
                                                       "time_seconds": str(int(fc["scan_s"][i]))} for i in range(n_scans)]}]
    # the fixture's last temperature sample reads low (the sample taken as a run stops); here the run goes on past it
    for t in seq["temperature"][0]["temperatures"]:
        for device in ("minion", "promethion"):
            for key in ("heatsink_temperature", "flowcell_temperature"):
                if key in t.get(device, {}):
                    t[device][key] = round(t["target_temperature"]["minimum"] + float(rng.normal(0.0, 0.02)), 4)
    with open(fr.files["minknow.report_json"], "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1)
    # ---- the sample sheet MinKNOW keeps in the run folder names the samples of a barcoded run
    if run.aliases:
        stamp = spec.start[:10].replace("-", "") + "_" + spec.start[11:16].replace(":", "")
        position, product = spec.position or fx._FAMILIES[run.family][4], spec.product or fx._FAMILIES[run.family][2]
        rows = ["protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,flow_cell_product_code,kit,alias,barcode"]
        rows += [f"{prot_id},{position},{fcid},{run.sample},{run.experiment},{product},{run.kit},{a},barcode{i + 1:02d}" for i, a in enumerate(run.aliases)]
        fx._write(os.path.join(fr.run_dir, f"sample_sheet_{fcid}_{stamp}_{run8}.csv"), "\n".join(rows) + "\n")
    with open(fr.files["minknow.final_summary"], encoding="utf-8") as fh:
        text = fh.read().replace("pod5_files_in_final_dest=1\n", f"pod5_files_in_final_dest={n // 4000 + 1}\n")
    if run.state != "PROTOCOL_COMPLETED":
        os.remove(fr.files["minknow.final_summary"])            # a run that was stopped by hand before MinKNOW wrote its summary
    else:
        fx._write(fr.files["minknow.final_summary"], text)
    _set_times(fr.run_dir, t0 + minutes * 60 + 90)
    return fr


def _set_times(path: str, epoch: float) -> None:
    """Give a folder and everything in it the time it would have had: a scan leaves what changed in the last two
    minutes for the next scan, and a run folder that is still changing counts as running."""
    epoch = min(epoch, time.time() - 600)
    for d, _dirs, files in os.walk(path):
        for name in files:
            os.utime(os.path.join(d, name), (epoch, epoch))
        os.utime(d, (epoch, epoch))


# ------------------------------------------------------------------------------------------- BAM files
_POOL: dict = {}


def _pools() -> tuple[str, str]:
    """A long invented sequence and a matching string of base qualities; reads are slices of them."""
    if not _POOL:
        rng = np.random.default_rng(424242)
        size = 1 << 22
        _POOL["seq"] = np.frombuffer(b"ACGT", np.uint8)[rng.integers(0, 4, size, dtype=np.uint8)].tobytes().decode("ascii")
        # base qualities in stretches of some tens of bases, as a basecaller writes them (and so that the files stay small)
        runs = rng.geometric(1 / 40.0, size // 20)
        _POOL["qual"] = (np.repeat(np.clip(rng.normal(27, 7, len(runs)), 3, 50).astype(np.uint8), runs)[:size] + 33).tobytes().decode("ascii")
    return _POOL["seq"], _POOL["qual"]


def _pysam():
    import pysam
    pysam.set_verbosity(0)
    return pysam


def _to_bam(sam: str, bam: str, *, sort: bool, threads: int = 2) -> None:
    pysam = _pysam()
    if os.path.getsize(sam) > (2 << 30):                          # a large file of a large example: more hands
        threads = 6
    if sort:
        pysam.sort("-@", str(threads), "-m", "512M", "-o", bam, sam, catch_stdout=False)
        pysam.index(bam)
    else:
        pysam.view("-b", "-@", str(threads), "-o", bam, sam, catch_stdout=False)
    os.remove(sam)


def _rg(run: DemoRun, model: str) -> str:
    return f"{_ids_of(run)[0]}_{model}"


def _dorado_header(run: DemoRun, model: str, barcodes: list[str], *, demux: bool) -> str:
    prot = _ids_of(run)[0]
    cl = f'dorado "basecaller" "/models/{model.split("@v5.2.0_")[0] if "_5mCG" in model else model}" "./pod5" "-r"'
    if "_5mCG" in model:
        cl += ' "--modified-bases" "5mCG_5hmCG"'
    if run.aliases:
        cl += f' "--kit-name" "{run.kit}" "--no-trim"'
    cl += ' "--min-qscore" "10"'
    if run.assay != "gdna":
        cl += ' "--estimate-poly-a"'
    lines = ["@HD\tVN:1.6\tSO:unknown", f"@PG\tID:basecaller\tPN:dorado\tVN:1.3.1+0000000\tCL:{cl}"]
    if demux:
        lines.append(f'@PG\tID:demux\tPP:basecaller\tPN:dorado\tVN:1.3.1+0000000\tCL:dorado "demux" "--kit-name" "{run.kit}" "--no-trim" "--output-dir" "./out" "calls.bam"')
    ds = f"runid={prot} basecall_model={model} experiment_id={run.experiment} model_stride=6"
    base = f"PU:{run.flow_cell_id}\tPM:demohost-01\tDT:{run.start}\tPL:ONT\tDS:{ds}\tLB:{run.sample}"
    lines.append(f"@RG\tID:{_rg(run, model)}\t{base}")
    for b in barcodes:
        lines.append(f"@RG\tID:{_rg(run, model)}_{run.kit}_{b}\t{base}\tSM:{b}\tal:{b}")
    return "\n".join(lines) + "\n"


def _read_columns(run: DemoRun, rd: dict, sel: np.ndarray, model: str, *, trimmed: bool, seed: int) -> dict:
    """What every record of the chosen reads carries, as lists of text ready to be joined into SAM lines."""
    rng = np.random.default_rng([run.seed, seed])
    prot, acq = _ids_of(run)
    n = len(sel)
    ids = [x.decode("ascii") for x in read_ids(run, 0, rd["n"])[sel]] if n else []
    length = rd["length"][sel].copy()
    lead = np.zeros(n, np.int64)
    if trimmed:                                                    # pychopper cut the primers off and wrote where the read lay
        lead = rng.integers(48, 96, n)
        length = np.maximum(rd["part"][sel] + rd["tail"][sel], 40)
    seq_pool, qual_pool = _pools()
    off = rng.integers(0, len(seq_pool) - int(length.max(initial=1)) - 1, n) if n else np.zeros(0, np.int64)
    us = np.rint((_epoch(run) + ORIGIN_S + rd["start_s"][sel]) * 1e6).astype("int64").astype("datetime64[us]")
    st = [s + "+00:00" for s in np.datetime_as_string(us, unit="ms")] if n else []
    labels = [run.aliases[a] if a >= 0 else "unclassified" for a in rd["alias"][sel]] if run.aliases else [None] * n
    rg = _rg(run, model)
    tags = []
    q, du, ch, mux = rd["q"][sel], rd["duration"][sel], rd["channel"][sel], rd["mux"][sel]
    tail = rd["tail"][sel] if rd["tail"] is not None else None
    no_tail = rng.random(n) < (0.065 if run.assay == "cdna" else 0.05)
    pod5 = sel // 4000
    for i in range(n):
        t = f"qs:f:{q[i]:.3f}\tdu:f:{du[i]:.4f}\tns:i:{int(du[i] * 5000)}\tts:i:{10 + i % 40}\tmx:i:{mux[i]}\tch:i:{ch[i]}\tst:Z:{st[i]}\trn:i:{100 + i % 9000}" \
            f"\tfn:Z:{run.flow_cell_id}_{prot[:8]}_{acq[:8]}_{pod5[i]}.pod5\tsv:Z:pa\tdx:i:0\tRG:Z:{rg}"
        if labels[i]:
            t += f"_{labels[i]}\tBC:Z:{labels[i]}"
        if tail is not None:
            t += f"\tpt:i:{-1 if no_tail[i] else tail[i] + int(rng.integers(-6, 7))}"
        tags.append(t)
    return {"n": n, "ids": ids, "length": length, "off": off, "tags": tags, "lead": lead, "q": q, "labels": labels}


def write_dorado_bam(path: str, run: DemoRun, scale: float, *, alias: int | None, passed: bool, model: str | None = None) -> int:
    """An unaligned BAM as the basecaller writes it: the reads of one sample (or of the whole run) that passed or
    failed the quality filter. Returns the number of records."""
    rd = reads_of(run, scale)
    model = model or run.model
    keep = rd["called"] & (rd["passed"] == passed)
    if alias is not None:
        keep &= rd["alias"] == alias
    sel = np.nonzero(keep)[0]
    c = _read_columns(run, rd, sel, model, trimmed=False, seed=11 + (alias or 0) * 2 + int(passed))
    seq_pool, qual_pool = _pools()
    barcodes = [f"barcode{alias + 1:02d}"] if (alias is not None and alias >= 0) else []
    mods = "\tMM:Z:C+m?,3,7,12;\tML:B:C,224,18,241" if "_5mCG" in model else ""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    sam = path + ".sam"
    with open(sam, "w", encoding="ascii") as fh:
        fh.write(_dorado_header(run, model, barcodes, demux=bool(run.aliases)))
        L, off, ids, tags = c["length"], c["off"], c["ids"], c["tags"]
        for i in range(c["n"]):
            a, z = off[i], off[i] + L[i]
            fh.write(f"{ids[i]}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq_pool[a:z]}\t{qual_pool[a:z]}\t{tags[i]}{mods}\n")
    _to_bam(sam, path, sort=False)
    return c["n"]


def _errors(rng, q: np.ndarray, cols: np.ndarray, assay: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mismatches, inserted and deleted bases of each read: the measured quality follows the basecaller's own with
    some scatter, a little below it."""
    qe = q * (0.86 if assay == "drna" else 0.93) + rng.normal(0.0, 1.3, len(q))
    e = np.power(10.0, -np.clip(qe, 3.0, 40.0) / 10.0)
    return rng.binomial(cols, np.minimum(e * 0.40, 0.3)), rng.binomial(cols, np.minimum(e * 0.22, 0.2)), rng.binomial(cols, np.minimum(e * 0.38, 0.3))


def _usable(rd: dict, alias: int | None) -> np.ndarray:
    """The reads an alignment of one sample holds. A read pychopper rescued was cut in two pieces: it is there twice."""
    keep = rd["usable"] if alias is None else (rd["usable"] & (rd["alias"] == alias))
    sel = np.nonzero(keep)[0]
    if "fate" in rd:
        sel = np.sort(np.concatenate([sel, sel[rd["fate"][sel] == 2]]))
    return sel


def _aligned_header(contigs, minimap_cl: str) -> str:
    return ("@HD\tVN:1.6\tSO:unsorted\n" + "".join(f"@SQ\tSN:{n}\tLN:{ln}\n" for n, ln in contigs)
            + f"@PG\tID:minimap2\tPN:minimap2\tVN:2.31-r0000\tCL:{minimap_cl}\n")


def write_genome_bam(path: str, parts: list[tuple[DemoRun, int | None]], scale: float, *, sample: str) -> int:
    """A sorted alignment to the genome of the usable reads of one sample, which may come from several runs. cDNA and
    RNA reads are spliced and pile up on their genes; genomic reads cover the genome with slow waves and a few gaps."""
    ref = reference()
    seq_pool, qual_pool = _pools()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    sam = path + ".sam"
    assay = parts[0][0].assay
    cl = (f"minimap2 -ax splice -uf -t 8 --MD -y genome.fa analysis/fastq/{sample}_full_length.fastq.gz" if assay == "cdna"
          else f"minimap2 -ax {'splice -uf -k14' if assay == 'drna' else 'lr:hq'} -t 8 -y genome.fa analysis/fastq/{sample}.fastq.gz")
    nuclear = np.array([ln for name, ln in GENOME[:-1]], dtype=np.float64)
    written = 0
    with open(sam, "w", encoding="ascii") as fh:
        fh.write(_aligned_header(GENOME, cl))
        for run, alias in parts:
            rd = reads_of(run, scale)
            sel = _usable(rd, alias)
            rng = np.random.default_rng([run.seed, 31, alias or 0])
            c = _read_columns(run, rd, sel, run.model, trimmed=assay == "cdna", seed=40 + (alias or 0))
            n, L = c["n"], c["length"]
            mapped = rd["on_genome"][sel]
            soft_l = rng.integers(0, 9, n)
            soft_r = rng.integers(0, 7, n)
            if assay == "gdna":
                chrom = rng.choice(len(nuclear), size=n, p=nuclear / nuclear.sum())
                mito = (rng.random(n) < 0.004) & (L < 14_000)       # a read longer than the mitochondrial genome does not lie on it
                chrom[mito] = len(GENOME) - 1
                clen = np.array([GENOME[k][1] for k in chrom])
                # coverage is uneven the way a genome is: it drifts a little along each chromosome, a few stretches take
                # hardly any read and a few take more than their share. By rejection: a read drawn there moves on.
                pos = (rng.random(n) * np.maximum(clen - L - 10, 1)).astype(np.int64)
                shape = np.random.default_rng(9173)                   # the same genome for every sample
                waves = [(shape.uniform(2.0, 31.0, 7), shape.uniform(0.0, 6.283, 7), shape.uniform(0.012, 0.04, 7)) for _ in GENOME]
                gaps = [(shape.uniform(0.05, 0.95, 3), shape.uniform(0.003, 0.012, 3), shape.uniform(0.02, 0.35, 3)) for _ in GENOME]
                peaks = [(shape.uniform(0.05, 0.95), shape.uniform(0.002, 0.006)) for _ in GENOME]
                for _ in range(4):
                    x = pos / clen
                    weight = np.full(n, 0.70)
                    for k in range(len(GENOME) - 1):
                        here = chrom == k
                        xs = x[here]
                        f, ph, amp = waves[k]
                        w = 0.70 + (amp[:, None] * np.sin(6.283 * f[:, None] * xs[None, :] + ph[:, None])).sum(axis=0)
                        for centre, half, depth in zip(*gaps[k]):
                            w = np.where(np.abs(xs - centre) < half, depth * 0.70, w)
                        w = np.where(np.abs(xs - peaks[k][0]) < peaks[k][1], 1.0, w)
                        weight[here] = w
                    again = (rng.random(n) > weight) & ~mito
                    pos[again] = (rng.random(int(again.sum())) * np.maximum(clen[again] - L[again] - 10, 1)).astype(np.int64)
                spans = None
                supp = mapped & (rng.random(n) < 0.014) & (L > 2000)
            else:
                g = ref["gene_of"][rd["tx"][sel]]
                chrom = ref["gene_chrom"][g]
                body = rd["part"][sel]
                soft_r = soft_r + rd["tail"][sel]                 # the tail is not in the genome: clipped
                if assay == "drna":
                    soft_l = soft_l + 30
                room = np.maximum(ref["gene_span"][g] - body, 0)
                spans = (room * rng.uniform(0.55, 1.0, n)).astype(np.int64)     # the introns the read jumps
                pos = ref["gene_start"][g] + ((room - spans) * rng.random(n)).astype(np.int64)
                L = np.maximum(L, soft_l + soft_r + 30)
                supp = mapped & (rng.random(n) < 0.008) & (L > 1200)
            aligned = np.maximum(L - soft_l - soft_r, 20)
            mism, ins, dele = _errors(rng, c["q"], aligned, assay)
            ins = np.minimum(ins, aligned // 4)
            m = aligned - ins
            # no alignment runs past the end of its chromosome, and a second piece of a read has room there too
            ends = np.array([ln for _, ln in GENOME], dtype=np.int64)[chrom]
            pos = np.clip(np.minimum(pos, ends - (m + dele + (spans if spans is not None else 0)) - 2), 0, None)
            supp &= pos + 41_000 < ends
            reverse = rng.random(n) < 0.5
            mapq = rng.choice([60, 47, 22, 3, 0], size=n, p=[0.89, 0.03, 0.03, 0.02, 0.03])
            n_intr = np.minimum(rng.poisson(3.2, n) + 1, 9) if spans is not None else None
            ids, off, tags, lead = c["ids"], c["off"], c["tags"], c["lead"]
            for i in range(n):
                a, z = off[i], off[i] + L[i]
                name = f"{lead[i]}:{lead[i] + L[i]}|{ids[i]}" if assay == "cdna" else ids[i]
                if not mapped[i]:
                    fh.write(f"{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq_pool[a:z]}\t{qual_pool[a:z]}\t{tags[i]}\n")
                    written += 1
                    continue
                mi = int(m[i])
                if spans is None or spans[i] < 40 or mi < 60:
                    body_cigar = f"{mi // 2}M" + (f"{ins[i]}I" if ins[i] else "") + (f"{dele[i]}D" if dele[i] else "") + f"{mi - mi // 2}M"
                else:
                    k = int(n_intr[i])
                    cut = np.sort(rng.integers(20, mi - 20, k)) if mi > 60 else np.array([mi // 2])
                    blocks = np.diff(np.concatenate([[0], cut, [mi]]))
                    gaps = np.maximum((spans[i] * rng.dirichlet(np.ones(len(cut)))).astype(np.int64), 30)
                    body_cigar = ""
                    for j, b in enumerate(blocks):
                        if b <= 0:
                            continue
                        body_cigar += f"{b}M"
                        if j == 0 and ins[i]:
                            body_cigar += f"{ins[i]}I1M"
                        if j == 0 and dele[i]:
                            body_cigar += f"{dele[i]}D1M"
                        if j < len(gaps):
                            body_cigar += f"{gaps[j]}N"
                    if ins[i] or dele[i]:
                        # the two single matches written after the insertion and the deletion come off the last block
                        tail_m = (1 if ins[i] else 0) + (1 if dele[i] else 0)
                        last = body_cigar.rfind("M")
                        prev = max(body_cigar.rfind(x, 0, last) for x in "NID")
                        body_cigar = body_cigar[:prev + 1] + f"{int(body_cigar[prev + 1:last]) - tail_m}M" + body_cigar[last + 1:]
                cigar = (f"{soft_l[i]}S" if soft_l[i] else "") + body_cigar + (f"{soft_r[i]}S" if soft_r[i] else "")
                nm = int(mism[i] + ins[i] + dele[i])
                flag = 16 if reverse[i] else 0
                extra = f"NM:i:{nm}\tms:i:{mi}\tAS:i:{max(mi - 3 * nm, 1)}\tnn:i:0\ttp:A:P" + ("\tts:A:+" if spans is not None else "")
                if supp[i]:
                    extra += f"\tSA:Z:{GENOME[int(chrom[i])][0]},{int(pos[i]) + 40_000},+,{L[i] - 400}S400M,60,4;"
                fh.write(f"{name}\t{flag}\t{GENOME[int(chrom[i])][0]}\t{int(pos[i]) + 1}\t{mapq[i]}\t{cigar}\t*\t0\t0\t{seq_pool[a:z]}\t{qual_pool[a:z]}"
                         f"\t{extra}\t{tags[i]}\n")
                written += 1
                if spans is not None and mapq[i] < 10:              # a read that fits two places: the other one, without a sequence
                    o = int((chrom[i] + 1) % (len(GENOME) - 1))
                    fh.write(f"{name}\t{256 | flag}\t{GENOME[o][0]}\t{int(pos[i]) % (GENOME[o][1] - 200_000) + 1}\t0\t{cigar}\t*\t0\t0\t*\t*"
                             f"\tNM:i:{nm + 9}\ttp:A:S\t{tags[i]}\n")
                    written += 1
                if supp[i]:
                    fh.write(f"{name}\t2048\t{GENOME[int(chrom[i])][0]}\t{int(pos[i]) + 40_001}\t60\t{L[i] - 400}H400M\t*\t0\t0\t{seq_pool[a:a + 400]}"
                             f"\t{qual_pool[a:a + 400]}\tNM:i:4\ttp:A:P\tSA:Z:{GENOME[int(chrom[i])][0]},{int(pos[i]) + 1},+,{cigar},60,{nm};\t{tags[i]}\n")
                    written += 1
    _to_bam(sam, path, sort=True)
    return written


def write_transcriptome_bam(path: str, parts: list[tuple[DemoRun, int | None]], scale: float, *, sample: str) -> dict:
    """A sorted alignment to the transcripts: reads end near the 3' end of their transcript and reach as far towards
    its 5' end as they are long; the other transcripts of the gene appear as secondary records without a sequence.
    Returns the counts the tables of the sample are written from."""
    ref = reference()
    seq_pool, qual_pool = _pools()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    sam = path + ".sam"
    assay = parts[0][0].assay
    contigs = list(zip(ref["tx_name"], (int(x) for x in ref["tlen"])))
    first_tx = np.concatenate([[0], np.cumsum(ref["n_tx"])[:-1]])
    records = primary = mapped_n = 0
    with open(sam, "w", encoding="ascii") as fh:
        fh.write(_aligned_header(contigs, f"minimap2 -ax map-ont -uf -t 8 -y transcripts.fa analysis/fastq/{sample}"
                                          f"{'_full_length' if assay == 'cdna' else ''}.fastq.gz"))
        for run, alias in parts:
            rd = reads_of(run, scale)
            sel = _usable(rd, alias)
            rng = np.random.default_rng([run.seed, 32, alias or 0])
            c = _read_columns(run, rd, sel, run.model, trimmed=assay == "cdna", seed=50 + (alias or 0))
            n = c["n"]
            mapped = rd["on_tx"][sel]
            tx = rd["tx"][sel]
            tlen = ref["tlen"][tx]
            body = np.minimum(rd["part"][sel], tlen)
            back = np.minimum(rng.integers(0, 14, n), tlen - body)
            soft_l = rng.integers(0, 6, n) + (30 if assay == "drna" else 0)
            soft_r = rd["tail"][sel] + rng.integers(0, 5, n)
            L = body + soft_l + soft_r
            mism, ins, dele = _errors(rng, c["q"], body, assay)
            ins = np.minimum(ins, body // 4)
            # the alignment covers the matched and the deleted bases of the transcript and ends where the read does
            dele = np.minimum(dele, tlen - back - (body - ins))
            start = tlen - back - (body - ins) - dele
            reverse = rng.random(n) < (0.03 if assay == "cdna" else 0.004)
            gene = ref["gene_of"][tx]
            others = ref["n_tx"][gene] - 1
            mapq = np.where(others == 0, 60, rng.choice([60, 30, 12, 1, 0], size=n, p=[0.35, 0.2, 0.15, 0.1, 0.2]))
            n_sec = np.minimum(rng.poisson(0.7 * others), 3)
            ids, off, tags, lead = c["ids"], c["off"], c["tags"], c["lead"]
            for i in range(n):
                a, z = off[i], off[i] + L[i]
                name = f"{lead[i]}:{lead[i] + L[i]}|{ids[i]}" if assay == "cdna" else ids[i]
                primary += 1
                if not mapped[i]:
                    fh.write(f"{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq_pool[a:z]}\t{qual_pool[a:z]}\t{tags[i]}\n")
                    records += 1
                    continue
                mapped_n += 1
                mi = int(body[i] - ins[i])
                d = int(dele[i])
                cigar = (f"{soft_l[i]}S" if soft_l[i] else "") + f"{mi // 2}M" + (f"{ins[i]}I" if ins[i] else "") + (f"{d}D" if d else "") \
                    + f"{mi - mi // 2}M" + (f"{soft_r[i]}S" if soft_r[i] else "")
                nm = int(mism[i] + ins[i] + d)
                fh.write(f"{name}\t{16 if reverse[i] else 0}\t{ref['tx_name'][tx[i]]}\t{int(start[i]) + 1}\t{mapq[i]}\t{cigar}\t*\t0\t0\t{seq_pool[a:z]}"
                         f"\t{qual_pool[a:z]}\tNM:i:{nm}\tms:i:{mi}\tAS:i:{max(mi - 3 * nm, 1)}\ttp:A:P\t{tags[i]}\n")
                records += 1
                for s in range(int(n_sec[i])):
                    o = int(first_tx[gene[i]] + (tx[i] - first_tx[gene[i]] + 1 + s) % (others[i] + 1))
                    span = int(min(body[i], ref["tlen"][o]))
                    fh.write(f"{name}\t256\t{ref['tx_name'][o]}\t{int(ref['tlen'][o]) - span + 1}\t0\t{L[i] - span}S{span}M\t*\t0\t0\t*\t*"
                             f"\tNM:i:{nm + 6}\ttp:A:S\t{tags[i]}\n")
                    records += 1
    _to_bam(sam, path, sort=True)
    return {"records": records, "primary": primary, "mapped": mapped_n}


# ----------------------------------------------------------------------------------------------- tables
def write_sample_tables(project_dir: str, sample: str, parts: list[tuple[DemoRun, int]], scale: float, *, quant: bool = False) -> dict:
    """The small files a cDNA pipeline leaves per sample, made from the merged reads of the sample: pychopper
    statistics, the annotated read table (one row per read with transcript, gene, type and poly(A) tail) and, with
    ``quant``, the metadata of a salmon run. Their numbers agree with the BAM files of the sample."""
    from flowcellar.testing import libprep
    ref = reference()
    ana = os.path.join(project_dir, "analysis")
    fate = np.concatenate([reads_of(r, scale)["fate"][(reads_of(r, scale)["alias"] == a)] for r, a in parts])
    rng = np.random.default_rng([parts[0][0].seed, 61, parts[0][1]])
    counts = np.bincount(fate, minlength=6)
    u = int(counts[3])
    spec = libprep.PychopperSpec(
        primers_found=int(counts[1] + counts[5]), plus=int((counts[1] + counts[5]) * rng.uniform(0.49, 0.51)),
        rescued={2: int(counts[2])} if counts[2] else {},
        unusable_hits={k: v for k, v in ((0, int(u * 0.42)), (1, int(u * 0.31)), (3, int(u * 0.17)), (4, u - int(u * 0.42) - int(u * 0.31) - int(u * 0.17))) if v},
        qc_fail=int(counts[4]), len_fail=int(counts[5]), umi=None, cutoff_index=int(rng.integers(1, 6)), seed=int(parts[0][0].seed + parts[0][1]))
    py = libprep.write_pychopper_stats_full(os.path.join(ana, "pychopper_reports", f"{sample}_pychopper_stats.tsv"), spec)
    # ---- the annotated read table: the reads on a transcript, sorted by transcript as the real tables are
    rows = []
    mapped_total = 0
    for run, alias in parts:
        rd = reads_of(run, scale)
        sel = np.nonzero(rd["on_tx"] & (rd["alias"] == alias))[0]
        mapped_total += len(sel)
        r2 = np.random.default_rng([run.seed, 62, alias])
        ids = read_ids(run, 0, rd["n"])[sel]
        tx = rd["tx"][sel]
        tail = np.where(r2.random(len(sel)) < 0.07, -1, rd["tail"][sel] + r2.integers(-5, 6, len(sel)))
        flag = np.where(r2.random(len(sel)) < 0.03, 16, 0)
        mq = r2.choice([60, 60, 60, 30, 1, 0], size=len(sel))
        g = ref["gene_of"][tx]
        for i in range(len(sel)):
            t = int(tx[i])
            rows.append((ref["tx_id"][t], ids[i].decode("ascii"), ref["gene_id"][g[i]], ref["gene_name"][g[i]], int(ref["tlen"][t]), ref["tx_type"][t],
                         int(mq[i]), int(flag[i]), int(tail[i])))
    rows.sort()
    path = os.path.join(ana, "polya", f"{sample}_readid_annotated.csv")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="ascii", newline="") as fh:
        fh.write(",".join(libprep.TABLE_COLUMNS) + "\n")
        for t, rid, gid, gname, ln, typ, mq, flag, tail in rows:
            fh.write(f"{rid},{t},{gid},{gname},{ln},{typ},{mq},{flag},{tail}\n")
    if quant:                                                      # the metadata of the salmon run that counted the reads of the sample
        from flowcellar.testing import readdata
        readdata.write_salmon_meta(os.path.join(ana, "counts", "salmon", f"{sample}_salmon"), py["scalars"]["written"], int(mapped_total * 0.972))
    return {"pychopper": py["scalars"], "table_rows": len(rows), "on_transcript": mapped_total}


# ---------------------------------------------------------------------------------------- the whole example
CONFIG_HEAD = """\
# The configuration of the FlowCellar example library (written by python -m flowcellar.testing.demo).
# `flowcellar scan` adds the folders it is given as [[roots]] entries; everything else in this file is the part a
# laboratory writes once: where its pipeline leaves the files that are not inside a run folder.

[library]
name = "FlowCellar example library"

# ---- read data: patterns relative to each scanned folder ------------------------------------------------------
# Basecalled BAM files in bam_pass/ and bam_fail/ folders are found without a pattern. The genomic project keeps two
# basecalling passes per sample, in folders named after the pass: the one named in [variants] is counted.
[[inputs]]
role = "dorado.bam"
glob = "*/basecalled/{variant:sup|sup_5mC}/{sample}.bam"

[[inputs]]
role = "align.genome_bam"
glob = "*/analysis/alignments/{sample}_genome_sorted.bam"

[[inputs]]
role = "align.transcriptome_bam"
glob = "*/analysis/alignments/{sample}_transcriptome_sorted.bam"

[[inputs]]
role = "pychopper.stats"
glob = "*/analysis/pychopper_reports/{sample}_pychopper_stats.tsv"

[[inputs]]
role = "polya.read_table"
glob = "*/analysis/polya/{sample}_readid_annotated.csv"

[[inputs]]
role = "quant.salmon"
glob = "*/analysis/counts/salmon/{sample}_salmon/aux_info/meta_info.json"

[[inputs]]
role = "quant.oarfish"
glob = "*/analysis/counts/oarfish/{sample}.meta_info.json"

[[inputs]]
role = "analysis.spaceranger"
glob = "*/analysis/*/spaceranger/outs/metrics_summary.csv"

[[inputs]]
role = "analysis.wf_single_cell"
glob = "*/analysis/*/wf-single-cell/wf-single-cell-report.html"

[[link.rules]]                    # tables without identifiers inside follow their sample
name = "experiment-folder"
roles = ["pychopper.stats", "polya.read_table", "quant.salmon", "quant.oarfish"]
via = "sample"
within = "experiment"

[variants]
primary = "sup"
"""

CONFIG_LIMITS = """
# ---- output limits of this example ----------------------------------------------------------------------------
# The example runs hold {percent} of the reads of full-size runs, so that they are written and scanned in minutes.
# The limits below are the built-in output limits multiplied by the same factor: the runs are judged as full-size
# runs would be. A library of real runs does not need this block.
"""


def config_text(scale: float) -> str:
    text = CONFIG_HEAD
    if scale < 1.0:
        pct = f"{scale * 100:g} %"
        text += CONFIG_LIMITS.format(percent=pct)
        for group in ("cdna", "gdna", "drna"):
            for check, by_group in FULL_SIZE_LIMITS.items():
                if group in by_group:
                    warn, fail = by_group[group]
                    text += f"[profiles.{group}.rules.{check}]\nwarn = {int(round(warn * scale))}\nfail = {int(round(fail * scale))}\n"
    return text


@dataclass
class DemoTree:
    out_dir: str
    scale: float
    roots: dict[str, str] = field(default_factory=dict)
    library: str = ""
    runs: dict[str, str] = field(default_factory=dict)            # label -> run folder
    run_keys: dict[str, str] = field(default_factory=dict)
    files: int = 0
    skipped: list[str] = field(default_factory=list)


def _task(args: tuple):
    """One file (or one run folder) of the example; runs in a worker process."""
    kind = args[0]
    if kind == "run":
        _, root, run, scale, move_to = args
        fr = write_run_folder(root, run, scale)
        if move_to:                                               # an older run, moved from the sequencer's folder to the project share
            import shutil
            dest = os.path.join(move_to, run.sample, os.path.basename(fr.run_dir))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.move(fr.run_dir, dest)
            _set_times(os.path.dirname(dest), _end(run) + 90)
            return ("run", run.label, dest)
        return ("run", run.label, fr.run_dir)
    if kind == "dorado":
        _, path, run, scale, alias, passed, model, when = args
        write_dorado_bam(path, run, scale, alias=alias, passed=passed, model=model)
    elif kind == "genome":
        _, path, parts, scale, sample, when = args
        write_genome_bam(path, parts, scale, sample=sample)
    elif kind == "transcriptome":
        _, path, parts, scale, sample, when = args
        counts = write_transcriptome_bam(path, parts, scale, sample=sample)
        if parts[0][0].assay == "drna":                           # the quantification that followed: how many reads it counted
            from flowcellar.testing import readdata
            meta = os.path.join(os.path.dirname(os.path.dirname(path)), "counts", "oarfish", f"{sample}.meta_info.json")
            readdata.write_oarfish_meta(meta, int(counts["mapped"] * 0.968), discard_score=int(counts["mapped"] * 0.021),
                                        discard_supp=int(counts["mapped"] * 0.011))
    elif kind == "tables":
        _, project_dir, sample, parts, scale, when = args
        write_sample_tables(project_dir, sample, parts, scale, quant=len(parts) == 1)
        return ("tables", sample, project_dir)
    elif kind == "analysis":
        _, base, sample, parts, scale = args
        from flowcellar.testing import demo_spatial
        demo_spatial.write_analysis(base, sample, parts, _ids_of, scale=scale)
        return ("analysis", sample, base)
    return (kind, args[1], None)


def _end(run: DemoRun) -> float:
    return _epoch(run) + ORIGIN_S + run.hours * 3600


def make_demo(out_dir: str | os.PathLike, *, scale: float = 0.01, jobs: int = 4, seed: int = 7, progress=None) -> DemoTree:
    """Write the example under ``out_dir`` (which must not hold an earlier example) and return where everything is."""
    from concurrent.futures import ProcessPoolExecutor, as_completed
    out = os.path.abspath(os.fspath(out_dir))
    say = progress or (lambda text: None)
    tree = DemoTree(out, scale, roots={"minknow": os.path.join(out, "minknow"), "projects": os.path.join(out, "projects")},
                    library=os.path.join(out, "library"))
    for d in tree.roots.values():
        if os.path.isdir(d) and os.listdir(d):
            raise FileExistsError(f"{d} is not empty: give a new folder, or remove the earlier example first")
    runs = demo_runs()
    by = {r.label: r for r in runs}
    try:
        _pysam()
        bam = True
    except ImportError:
        bam = False
        tree.skipped.append("BAM files (pysam is not installed: pip install 'flowcellar[bam]')")
    tasks: list[tuple] = []
    for r in runs:
        # the two oldest runs were moved to the project share: one project under two scanned folders
        move_to = None if r.root == "minknow" else os.path.join(tree.roots["projects"], r.experiment, "runs")
        tasks.append(("run", tree.roots["minknow"] if r.root == "minknow" else os.path.join(out, ".moving"), r, scale, move_to))
    proj = tree.roots["projects"]
    done_at: dict[str, float] = {}
    if bam:
        # ---- cDNA: the six samples of one pool on one flow cell, and the six of the last pool on two flow cells
        cd = os.path.join(proj, "demo_cdna_timecourse")
        pair = [by["c5"], by["c6"]]
        for r in pair + [by["c7"]]:
            base = os.path.join(cd, "basecalled", f"{r.start[:10].replace('-', '')}_{r.flow_cell_id}")
            prot, acq = _ids_of(r)
            for a in list(range(len(r.aliases))) + [-1]:
                label = r.aliases[a] if a >= 0 else "unclassified"
                for passed in (True, False):
                    name = f"{r.flow_cell_id}_{'pass' if passed else 'fail'}_{label}_{prot[:8]}_{acq[:8]}_0.bam"
                    tasks.append(("dorado", os.path.join(base, "bam_pass" if passed else "bam_fail", label, name), r, scale, a, passed, None, _end(r) + 86400))
        for group in (pair, [by["c7"]]):                           # two flow cells, merged per sample before pychopper; one flow cell
            when = _end(group[-1]) + 3 * 86400
            for a, s in enumerate(group[0].aliases):
                parts = [(r, a) for r in group]
                tasks.append(("genome", os.path.join(cd, "analysis", "alignments", f"{s}_genome_sorted.bam"), parts, scale, s, when))
                tasks.append(("transcriptome", os.path.join(cd, "analysis", "alignments", f"{s}_transcriptome_sorted.bam"), parts, scale, s, when))
                tasks.append(("tables", cd, s, parts, scale, when))
        # ---- genomic DNA: two basecalling passes and an alignment of the newest run
        g = by["g3"]
        gd = os.path.join(proj, "demo_genome")
        for variant, model in (("sup", MODEL_DNA), ("sup_5mC", MODEL_MODS)):
            tasks.append(("dorado", os.path.join(gd, "basecalled", variant, f"{g.sample}.bam"), g, scale, None, True, model, _end(g) + 86400))
        tasks.append(("genome", os.path.join(gd, "analysis", "alignments", f"{g.sample}_genome_sorted.bam"), [(g, None)], scale, g.sample, _end(g) + 2 * 86400))
        # ---- direct RNA: basecalled reads and an alignment to the transcripts
        d = by["r1"]
        dd = os.path.join(proj, "demo_direct_rna")
        prot, acq = _ids_of(d)
        for passed in (True, False):
            name = f"{d.flow_cell_id}_{'pass' if passed else 'fail'}_{prot[:8]}_{acq[:8]}_0.bam"
            tasks.append(("dorado", os.path.join(dd, "basecalled", "bam_pass" if passed else "bam_fail", name), d, scale, None, passed, None, _end(d) + 86400))
        tasks.append(("transcriptome", os.path.join(dd, "analysis", "alignments", f"{d.sample}_transcriptome_sorted.bam"), [(d, None)], scale, d.sample,
                      _end(d) + 2 * 86400))
    # ---- single-cell and spatial analyses: a Visium slide from one flow cell, a Visium HD slide from two
    try:
        import h5py  # noqa: F401
        if not bam:
            raise ImportError
        base = os.path.join(proj, "demo_spatial", "analysis")
        tasks.append(("analysis", base, "visium_a", [(by["s1"], 100)], scale))
        tasks.append(("analysis", base, "visium_hd_b", [(by["s2"], 61), (by["s3"], 39)], scale))
    except ImportError:
        tree.skipped.append("the single-cell and spatial analysis folders (h5py and pysam are needed: pip install 'flowcellar[all]')")
    t0 = time.time()
    results = []
    # new processes, not copies of this one: a copy of a process that already runs threads is not safe to work in
    import multiprocessing
    with ProcessPoolExecutor(max_workers=max(1, jobs), mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = {pool.submit(_task, t): t for t in tasks}
        for k, f in enumerate(as_completed(futures), 1):
            res = f.result()
            results.append(res)
            t = futures[f]
            if t[0] != "run" and isinstance(t[-1], float) and t[0] != "tables":
                done_at[t[1]] = t[-1]
            if k % 10 == 0 or k == len(tasks):
                say(f"  {k} of {len(tasks)} parts written ({time.time() - t0:.0f} s)")
    for res in results:
        if res[0] == "run":
            tree.runs[res[1]] = res[2]
            tree.run_keys[res[1]] = run_key(by[res[1]])
    # ---- the sample sheets the basecaller was given: which sample sat on which barcode of each flow cell
    if bam:
        for r in (by["c5"], by["c6"], by["c7"]):
            rows = ["experiment_id,kit,flow_cell_id,barcode,sample_id,alias"]
            rows += [f"{r.experiment},{r.kit},{r.flow_cell_id},barcode{i + 1:02d},{r.sample},{a}" for i, a in enumerate(r.aliases)]
            fx._write(os.path.join(proj, "demo_cdna_timecourse", "basecalled", f"{r.start[:10].replace('-', '')}_{r.flow_cell_id}", "samples.csv"),
                      "\n".join(rows) + "\n")
    import shutil
    shutil.rmtree(os.path.join(out, ".moving"), ignore_errors=True)
    # ---- file times: analyses a few days after their runs
    for path, when in done_at.items():
        for p in (path, path + ".bai"):
            if os.path.exists(p):
                os.utime(p, (min(when, time.time() - 600),) * 2)
    newest = max(_end(r) for r in runs)
    for d, _dirs, files in os.walk(proj):
        if os.sep + "runs" + os.sep in d + os.sep:
            continue
        for name in files:
            p = os.path.join(d, name)
            if p not in done_at and not (p.endswith(".bai") and p[:-4] in done_at):
                os.utime(p, (min(newest + 4 * 86400, time.time() - 600),) * 2)
    os.makedirs(tree.library, exist_ok=True)
    with open(os.path.join(tree.library, "flowcellar.toml"), "w", encoding="utf-8") as fh:
        fh.write(config_text(scale))
    tree.files = sum(len(files) for _d, _s, files in os.walk(out))
    return tree


def parser():
    import argparse
    p = argparse.ArgumentParser(prog="python -m flowcellar.testing.demo",
                                description="Write the FlowCellar example: invented run folders, read-data files and analysis folders, "
                                            "and a library folder that holds only its configuration.")
    p.add_argument("out_dir")
    p.add_argument("--scale", type=float, default=0.01,
                   help="size of the runs as a share of full-size runs (default 0.01: about a minute and 1.5 GB; 0.1: some minutes and 12 GB, "
                        "with smoother charts)")
    p.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1), help="worker processes (default: up to 8)")
    p.add_argument("--seed", type=int, default=7, help=argparse.SUPPRESS)
    return p


def main(argv: list[str] | None = None) -> int:
    p = parser()
    a = p.parse_args(argv)
    if not 0.0005 <= a.scale <= 1.0:
        p.error("--scale must lie between 0.0005 and 1")
    t0 = time.time()
    print(f"writing the FlowCellar example to {os.path.abspath(a.out_dir)} (scale {a.scale:g}, {a.jobs} worker processes)")
    try:
        tree = make_demo(a.out_dir, scale=a.scale, jobs=a.jobs, seed=a.seed, progress=print)
    except FileExistsError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    size = sum(os.path.getsize(os.path.join(d, f)) for d, _s, files in os.walk(tree.out_dir) for f in files)
    print(f"done in {time.time() - t0:.0f} s: {len(tree.runs)} run folders, {tree.files} files, {size / 1e9:.2f} GB")
    for s in tree.skipped:
        print(f"left out: {s}")
    rel = os.path.relpath(tree.out_dir)
    print("next:")
    print(f"  flowcellar scan {rel}/minknow {rel}/projects -L {rel}/library")
    print(f"  flowcellar serve {rel}/library --open")
    return 0


if __name__ == "__main__":
    sys.exit(main())
