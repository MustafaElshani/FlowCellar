"""Synthetic read-data files for tests: small BAMs written with pysam, pipeline tables, sample sheets.

Everything is invented and derived from a seed. Every writer returns what a correct reader must find, counted while
the file is written, so tests never recompute an expectation with the code under test.

    from flowcellar.testing.readdata import fake_run_ref, make_reads, write_dorado_bam, write_aligned_bam, make_project

The hazards of real files are reproduced on purpose: a first tag wrapped as ``CO:Z:<tag>:<type>:<value>``,
``ts`` present twice with two types, read names prefixed by a cDNA pre-processor, secondary records without a
sequence, supplementary records, read groups whose value carries an alias the header does not list, a file that
holds reads of two flow cells, a file without any identifier, a file cut short.

    python -m flowcellar.testing.readdata OUT_DIR        # one fake project beside a fake run folder
"""
from __future__ import annotations

import csv
import io
import json
import os
import sys
from array import array
from dataclasses import dataclass, field

import numpy as np

MODEL = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
KIT = "SQK-PCB114-24"
CONTIGS = (("chrT1", 2_000_000), ("chrT2", 1_000_000), ("chrM", 16_500))
ORIGIN_AFTER_START_S = 252.0        # the fixture run folders start reading data this long after the protocol starts
RUN_SPAN_S = 30 * 3600              # the start times of invented reads are spread over this much run time


def _pysam():
    import pysam
    pysam.set_verbosity(0)
    return pysam


def _hex(rng, n: int) -> str:
    return "".join(rng.choice(list("0123456789abcdef"), n))


def _uuid(rng) -> str:
    h = _hex(rng, 32)
    return f"{h[:8]}-{h[8:12]}-4{h[13:16]}-a{h[17:20]}-{h[20:32]}"


@dataclass(frozen=True)
class RunRef:
    """The identifiers of one run as files carry them."""
    flow_cell_id: str
    protocol_run_id: str
    acquisition_id: str
    experiment: str = "demo_experiment"
    start: str = "2026-03-02T09:00:00Z"

    @property
    def run8(self) -> str:
        return self.protocol_run_id[:8]

    @property
    def acq8(self) -> str:
        return self.acquisition_id[:8]

    @property
    def run_key(self) -> str:
        return f"{self.flow_cell_id}_{self.run8}"


def run_ref(fixture_run) -> RunRef:
    """The identifiers of a run folder made by ``flowcellar.testing.fixtures.make_run``."""
    return RunRef(fixture_run.flow_cell_id, fixture_run.protocol_run_id, fixture_run.acquisition_id,
                  fixture_run.spec.experiment, fixture_run.spec.start)


def fake_run_ref(flow_cell_id: str = "TST00009", seed: int = 9) -> RunRef:
    """Identifiers without a run folder (a run the library does not know)."""
    rng = np.random.default_rng(seed)
    return RunRef(flow_cell_id, _uuid(rng), _hex(rng, 40))


@dataclass
class Read:
    id: str
    run: RunRef
    length: int
    q: float                        # mean quality as a float32 stores it
    pod5: int
    channel: int
    alias: str | None = None        # the sample the demultiplexer gave the read to
    barcode: str | None = None
    parent: str | None = None       # id of the read this one was split from
    seq: str = field(default="", repr=False)
    start_s: float = 0.0            # start of the read, seconds after the run's time origin
    speed: float = 380.0            # bases per second through the pore
    pt: int | None = 60             # poly(A) estimate as the basecaller tags it: a length, 0, -1 (none found), None = no tag

    @property
    def du(self) -> float:
        """The duration as the ``du`` tag stores it (a 32-bit float)."""
        return float(np.float32(round(self.length / self.speed, 4)))

    @property
    def epoch(self) -> float:
        """Start time of the read on the clock, to the microsecond the ``st`` tag keeps."""
        import datetime as dt
        t0 = dt.datetime.strptime(self.run.start, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.UTC)
        return (t0 + dt.timedelta(microseconds=round((ORIGIN_AFTER_START_S + self.start_s) * 1e6))).timestamp()

    @property
    def st(self) -> str:
        import datetime as dt
        t0 = dt.datetime.strptime(self.run.start, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.UTC)
        t = t0 + dt.timedelta(microseconds=round((ORIGIN_AFTER_START_S + self.start_s) * 1e6))
        return t.strftime("%Y-%m-%dT%H:%M:%S.%f") + "+00:00"


def make_reads(n: int, run: RunRef, *, seed: int = 1, alias: str | None = None, barcode: str | None = None,
               mean_len: float = 900.0, q_mean: float = 15.0, split_every: int = 0) -> list[Read]:
    """``n`` invented reads of one run and one sample. With ``split_every`` = k, every k-th read is one of two pieces
    of a read the basecaller split (the two pieces share a parent)."""
    rng = np.random.default_rng(seed)
    lengths = np.clip(rng.lognormal(np.log(mean_len), 0.6, n), 12, 20000).astype(np.int64)
    quals = np.clip(rng.normal(q_mean, 4.0, n), 2.0, 45.0).astype(np.float32)
    out: list[Read] = []
    parent = None
    for i in range(n):
        seq = "".join(rng.choice(list("ACGT"), int(lengths[i])))
        r = Read(_uuid(rng), run, int(lengths[i]), float(quals[i]), i // 500, int(rng.integers(1, 513)), alias, barcode, None, seq)
        # set on purpose, without touching the random stream: when it started, how fast it went, its poly(A) tail
        r.start_s = float((i * 7919 + seed * 101) % RUN_SPAN_S) + (i % 1000) / 1000.0
        r.speed = 250.0 + (i * 37) % 300
        r.pt = -1 if i % 11 == 0 else 0 if i % 13 == 0 else None if i % 17 == 0 else 20 + r.length % 150
        if split_every and i % split_every == 0 and i + 1 < n:
            parent = _uuid(rng)
            r.parent = parent
        elif split_every and i % split_every == 1 and parent:
            r.parent, parent = parent, None
        out.append(r)
    return out


# ------------------------------------------------------------------------------------- expectations
def _blank() -> dict:
    return {"reads": 0, "bases": 0, "mapped": 0, "unmapped": 0, "secondary": 0, "supplementary": 0, "nm": 0, "cols": 0,
            "perfect": 0, "with_sa": 0, "q_missing": 0, "split_children": 0, "split_parents": 0, "lengths": [], "quals": [],
            "mapq": [], "identity_q": [], "_parents": set(),
            # the tags of the basecaller, one entry per read: (ten-minute clock bin or None, length, quality or None,
            # duration or None, poly(A) or None)
            "tagged": [],
            # the primary mapped records, one entry per record: a dict with length, q, m, ins, dele, skip, soft, hard,
            # mism, reverse, contig, start, blocks [(start, end) of the aligned stretches], t10
            "aligned": []}


def clock_bin(epoch: float, step: int = 600) -> int:
    return int(epoch // step)


class _Expect:
    def __init__(self):
        self.groups: dict[tuple, dict] = {}
        self.records = 0

    def g(self, key, label) -> dict:
        return self.groups.setdefault((key, label), _blank())

    def done(self) -> dict:
        for g in self.groups.values():
            g["split_parents"] = len(g.pop("_parents"))
        return {"records": self.records, "primary": sum(g["reads"] for g in self.groups.values()),
                "secondary": sum(g["secondary"] for g in self.groups.values()),
                "supplementary": sum(g["supplementary"] for g in self.groups.values()), "groups": self.groups}


# ------------------------------------------------------------------------------------------- BAMs
def _dorado_tags(a, r: Read, *, rg: str | None, bc: str | None, fn: bool = True, ints: str = "i") -> None:
    a.set_tag("qs", r.q, "f")
    a.set_tag("du", r.du, "f")
    a.set_tag("ns", r.length * 12, ints)
    a.set_tag("ts", 10, ints)
    a.set_tag("mx", 1 + r.channel % 4, ints)
    a.set_tag("ch", r.channel, ints)
    a.set_tag("st", r.st, "Z")
    a.set_tag("rn", 100 + r.pod5, ints)
    if fn:
        a.set_tag("fn", f"{r.run.flow_cell_id}_{r.run.run8}_{r.run.acq8}_{r.pod5}.pod5", "Z")
    a.set_tag("sv", "pa", "Z")
    a.set_tag("dx", 0, ints)
    if rg:
        a.set_tag("RG", rg, "Z")
    if r.pt is not None:
        a.set_tag("pt", r.pt, ints)
    if r.parent:
        a.set_tag("pi", r.parent, "Z")
    if bc:
        a.set_tag("BC", bc, "Z")


def _rg_id(run: RunRef, model: str) -> str:
    return f"{run.protocol_run_id}_{model}"


def write_dorado_bam(path: str, reads: list[Read], *, header_rg: bool = True, min_qscore: float | None = 7.0,
                     no_trim: bool = True, both_ends: bool = True, demux: bool = True, kit: str = KIT, model: str = MODEL,
                     bc_tag: bool = True, fn_tag: bool = True, rg_tag: bool = True, truncate_frac: float | None = None) -> dict:
    """An unaligned BAM as Dorado writes it. Read groups of the reads end in the alias, which the header lists only
    under the kit's barcode name (as Dorado does when a sample sheet renames barcodes)."""
    pysam = _pysam()
    runs = list({r.run.run_key: r.run for r in reads}.values())
    hd: dict = {"HD": {"VN": "1.6", "SO": "unknown"}}
    cl = ['dorado "basecaller" "/models/%s" "./pod5" "-r" "--kit-name" "%s"' % (model, kit)]
    if no_trim:
        cl.append('"--no-trim"')
    if min_qscore is not None:
        cl.append('"--min-qscore" "%g"' % min_qscore)
    if both_ends:
        cl.append('"--barcode-both-ends"')
    cl.append('"--estimate-poly-a"')
    hd["PG"] = [{"ID": "basecaller", "PN": "dorado", "VN": "1.3.1+0000000", "CL": " ".join(cl)}]
    if demux:
        hd["PG"].append({"ID": "demux", "PP": "basecaller", "PN": "dorado", "VN": "1.3.1+0000000",
                         "CL": 'dorado "demux" "--kit-name" "%s" "--output-dir" "./out"%s "/tmp/calls.bam"' % (kit, ' "--no-trim"' if no_trim else "")})
    if header_rg:
        hd["RG"] = []
        for run in runs:
            ds = f"runid={run.protocol_run_id} basecall_model={model} experiment_id={run.experiment} model_stride=6"
            base = {"ID": _rg_id(run, model), "PU": run.flow_cell_id, "PM": "fakehost-01", "DT": run.start, "PL": "ONT", "DS": ds, "LB": "demo_library"}
            hd["RG"].append(base)
            for bcode in sorted({r.barcode for r in reads if r.barcode and r.run is run}):
                hd["RG"].append({**base, "ID": f"{_rg_id(run, model)}_{kit}_{bcode}", "SM": bcode, "al": bcode})
    exp = _Expect()
    with pysam.AlignmentFile(path, "wb", header=hd) as out:
        for r in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name = r.id
            a.flag = 4
            a.query_sequence = r.seq
            a.query_qualities = array("B", [20]) * r.length
            label = r.alias or r.barcode
            rg = None
            if rg_tag:
                rg = _rg_id(r.run, model) + (f"_{label}" if label else "")
            _dorado_tags(a, r, rg=rg, bc=(label if bc_tag else None), fn=fn_tag)
            out.write(a)
            exp.records += 1
            if fn_tag or (rg_tag and header_rg) or (header_rg and len(runs) == 1):
                key = r.run.run_key
            else:
                key = "?_" + r.run.run8 if rg_tag else None      # a read group without a header gives the run id only
            g = exp.g(key, label if (bc_tag or rg_tag) else None)
            g["reads"] += 1; g["bases"] += r.length; g["unmapped"] += 1
            g["lengths"].append(r.length); g["quals"].append(r.q)
            g["tagged"].append((clock_bin(r.epoch), r.length, r.q, r.du, r.pt))
            if r.parent:
                g["split_children"] += 1; g["_parents"].add(r.parent)
    res = exp.done()
    if truncate_frac is not None:
        size = os.path.getsize(path)
        with open(path, "r+b") as fh:
            fh.truncate(int(size * truncate_frac))
        res["truncated"] = True
    return res


def write_aligned_bam(path: str, reads: list[Read], *, seed: int = 3, mapped_frac: float = 0.96, secondary_rate: float = 0.5,
                      supplementary_frac: float = 0.03, pychopper_names: bool = True, wrap: str | None = "mixed",
                      dorado_tags: bool = True, duplicate_ts: bool = True, index: str | None = "bai", model: str = MODEL,
                      full_length_input: bool = True, keep_unmapped: bool = True, fn_tag: bool = True, rg_tag: bool = True,
                      base_qualities: bool = True, contigs=CONTIGS) -> dict:
    """A coordinate-sorted BAM as minimap2 and samtools write it from a FASTQ that kept the basecaller's tags.

    ``wrap``: which tag arrives wrapped in ``CO`` (``"BC"``, ``"qs"``, ``"mixed"`` = both kinds of record, None = no
    wrapping). Secondary records carry no sequence; supplementary records are hard-clipped.
    """
    pysam = _pysam()
    rng = np.random.default_rng(seed)
    hd = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": n, "LN": ln} for n, ln in contigs],
          "PG": [{"ID": "minimap2", "PN": "minimap2", "VN": "2.31-r0000",
                  "CL": "minimap2 -ax splice -uf -t 8 --MD -y ref.fa analysis/fastq/sample%s.fastq.gz" % ("_full_length" if full_length_input else "")},
                 {"ID": "samtools", "PN": "samtools", "PP": "minimap2", "VN": "1.21", "CL": "samtools sort -o out.bam in.bam"}]}
    header = pysam.AlignmentHeader.from_dict(hd)
    exp = _Expect()
    recs: list[tuple] = []          # (sort key, segment)

    def wrapped_tag(k: int) -> str | None:
        if wrap == "BC" or (wrap == "mixed" and k % 5 != 0):
            return "BC"
        return "qs" if wrap in ("qs", "mixed") else None

    def tags(a, r: Read, k: int) -> None:
        label = r.alias or r.barcode
        wrapped = wrapped_tag(k)
        if wrapped == "BC" and label:
            a.set_tag("CO", f"BC:Z:{label}", "Z")
        elif wrapped == "qs":
            a.set_tag("CO", f"qs:f:{r.q!r}", "Z")
        if not dorado_tags:
            return
        if wrapped != "qs":
            a.set_tag("qs", r.q, "f")
        a.set_tag("du", r.du, "f")
        a.set_tag("st", r.st, "Z")
        if duplicate_ts:
            a.set_tag("ts", 10, "C", replace=False)          # the basecaller's ts beside the aligner's ts:A
        a.set_tag("ch", r.channel, "S")
        if fn_tag:
            a.set_tag("fn", f"{r.run.flow_cell_id}_{r.run.run8}_{r.run.acq8}_{r.pod5}.pod5", "Z")
        if rg_tag:
            a.set_tag("RG", _rg_id(r.run, model) + (f"_{KIT}_{r.barcode}" if r.barcode else ""), "Z")
        if wrapped != "BC" and label:
            a.set_tag("BC", label, "Z")
        if r.pt is not None:
            a.set_tag("pt", r.pt, "i")

    for k, r in enumerate(reads):
        name = f"{50 + k % 90}:{50 + k % 90 + r.length}|{r.id}" if pychopper_names else r.id
        wrapped = wrapped_tag(k)
        label = (r.alias or r.barcode) if (wrapped == "BC" or dorado_tags) else None
        if dorado_tags and fn_tag:
            key = r.run.run_key
        else:
            key = "?_" + r.run.run8 if (dorado_tags and rg_tag) else None    # no header read groups: the run id only
        g = exp.g(key, label)
        q_known = dorado_tags or wrapped == "qs"
        mapped = rng.random() < mapped_frac
        if not mapped and not keep_unmapped:
            continue
        a = pysam.AlignedSegment(header)
        a.query_name = name
        a.query_sequence = r.seq
        if base_qualities:
            a.query_qualities = array("B", [20]) * r.length
        g["reads"] += 1; g["bases"] += r.length; g["lengths"].append(r.length)
        if q_known:
            g["quals"].append(r.q)
        elif base_qualities:
            g["quals"].append(20.0)                         # no qs tag: the reader takes the mean of the base qualities
            g["q_computed"] = g.get("q_computed", 0) + 1
        else:
            g["q_missing"] += 1
        q_seen = r.q if q_known else 20.0 if base_qualities else None
        t10 = clock_bin(r.epoch) if dorado_tags else None
        g["tagged"].append((t10, r.length, q_seen, r.du if dorado_tags else None, r.pt if dorado_tags else None))
        if not mapped:
            a.flag = 4
            tags(a, r, k)
            recs.append(((1 << 40, 0, k), a))
            g["unmapped"] += 1
            exp.records += 1
            continue
        rid = int(rng.choice(len(contigs), p=[0.6, 0.3, 0.1]))
        pos = int(rng.integers(0, max(contigs[rid][1] - 60_000, 1)))
        clip = min(int(rng.integers(0, 12)), r.length - 4)
        body = r.length - clip
        ins = 2 if body > 40 and k % 3 == 0 else 0
        dele = 3 if body > 40 and k % 4 == 0 else 0
        m1 = (body - ins) // 2
        m2 = body - ins - m1
        skip = 700 if k % 5 else 0                          # four records in five are spliced
        hard = 3 if k % 9 == 0 else 0
        cigar = (([(4, clip)] if clip else []) + [(0, m1)] + ([(3, skip)] if skip else []) + ([(1, ins)] if ins else [])
                 + ([(2, dele)] if dele else []) + [(0, m2)] + ([(5, hard)] if hard else []))
        cols = m1 + m2 + ins + dele
        mism = 0 if k % 7 == 0 else int(rng.integers(0, max(cols // 25, 1) + 1))
        nm = mism + ins + dele
        reverse = bool(rng.integers(0, 2))
        a.flag = 16 if reverse else 0
        a.reference_id, a.reference_start = rid, pos
        a.mapping_quality = int(rng.choice([60, 60, 60, 30, 1, 0]))
        a.cigartuples = cigar
        a.set_tag("NM", nm, "i")
        a.set_tag("ms", 100 + k % 50, "i")
        a.set_tag("AS", 90 + k % 50, "i")
        a.set_tag("ts", "+", "A")
        a.set_tag("tp", "P", "A")
        has_supp = rng.random() < supplementary_frac and r.length > 60
        if has_supp:
            a.set_tag("SA", f"{contigs[0][0]},{pos + 5000},+,30S30M,60,1;", "Z")
        tags(a, r, k)
        recs.append(((rid, pos, k), a))
        exp.records += 1
        g["mapped"] += 1; g["nm"] += nm; g["cols"] += cols; g["perfect"] += int(nm == 0); g["with_sa"] += int(has_supp)
        g["mapq"].append(a.mapping_quality)
        g["identity_q"].append(float("inf") if nm == 0 else -10.0 * float(np.log10(nm / cols)))
        g["aligned"].append({"length": r.length, "q": q_seen, "m": m1 + m2, "ins": ins, "dele": dele, "skip": skip, "soft": clip,
                             "hard": hard, "mism": mism, "nm": nm, "cols": cols, "reverse": reverse, "contig": rid, "start": pos,
                             "blocks": [(pos, pos + m1), (pos + m1 + skip, pos + m1 + skip + dele + m2)], "t10": t10})
        n_sec = int(rng.poisson(secondary_rate))
        for s in range(n_sec):                              # secondary: no sequence, tags copied by the aligner
            b = pysam.AlignedSegment(header)
            b.query_name = name
            b.flag = 256 | (16 if s % 2 else 0)
            srid = int(rng.integers(0, len(contigs)))
            spos = int(rng.integers(0, max(contigs[srid][1] - 60_000, 1)))
            b.reference_id, b.reference_start = srid, spos
            b.mapping_quality = 0
            b.cigartuples = [(0, r.length)]
            b.set_tag("NM", 5, "i")
            b.set_tag("ts", "+", "A")
            b.set_tag("tp", "S", "A")
            tags(b, r, k)
            recs.append(((srid, spos, k + 0.5), b))
            exp.records += 1
            g["secondary"] += 1
        if has_supp:                                        # supplementary: hard-clipped piece of the same read
            c = pysam.AlignedSegment(header)
            c.query_name = name
            c.flag = 2048
            c.reference_id, c.reference_start = 0, pos + 5000
            c.mapping_quality = 60
            c.query_sequence = r.seq[:30]
            c.query_qualities = array("B", [20]) * 30
            c.cigartuples = [(5, r.length - 30), (0, 30)] if r.length > 30 else [(0, len(r.seq[:30]))]
            c.set_tag("NM", 0, "i")
            c.set_tag("SA", f"{contigs[rid][0]},{pos + 1},+,30M,60,0;", "Z")
            tags(c, r, k)
            recs.append(((0, pos + 5000, k + 0.7), c))
            exp.records += 1
            g["supplementary"] += 1
    recs.sort(key=lambda x: x[0])
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for _, a in recs:
            out.write(a)
    if index == "bai":
        pysam.index(path)
    elif index == "csi":
        pysam.index("-c", path)
    return exp.done()


def fake_transcripts(n: int = 420, seed: int = 4) -> list[tuple[str, int]]:
    """Invented transcripts named the way a reference of transcripts names them (``transcript|gene|...``): three
    transcripts per gene, lengths from 300 bases to 9 kb so that every length class is filled."""
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        ln = int(np.clip(rng.lognormal(np.log(1800), 0.9), 300, 9000))
        gene = i // 3
        out.append((f"ENST9{i:010d}.1|ENSG9{gene:010d}.1|-|-|FAKE{gene}-{201 + i % 3}|FAKE{gene}|{ln}|protein_coding|", ln))
    return out


def write_transcriptome_bam(path: str, reads: list[Read], *, transcripts=None, seed: int = 8, mapped_frac: float = 0.9,
                            secondary_rate: float = 1.2, index: str | None = "bai", model: str = MODEL) -> dict:
    """A coordinate-sorted BAM of reads aligned to a reference of transcripts, as minimap2 writes it: no skipped
    region, most reads on the transcript's own strand and ending near its 3' end, many secondary records.

    Besides the counts of ``write_aligned_bam`` the result holds, per group, ``tx``: one entry per primary mapped
    record with the transcript index, its length, the start and the end of the alignment on it and the strand; and
    ``transcripts`` (the reference) and ``genes`` (gene index per transcript)."""
    pysam = _pysam()
    rng = np.random.default_rng(seed)
    transcripts = transcripts or fake_transcripts()
    hd = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": n, "LN": ln} for n, ln in transcripts],
          "PG": [{"ID": "minimap2", "PN": "minimap2", "VN": "2.31-r0000", "CL": "minimap2 -ax map-ont -uf -t 8 -y transcripts.fa analysis/fastq/sample_full_length.fastq.gz"}]}
    header = pysam.AlignmentHeader.from_dict(hd)
    exp = _Expect()
    recs: list[tuple] = []
    weights = 1.0 / np.arange(1, len(transcripts) + 1) ** 1.1          # a few transcripts take most reads; many take none
    weights /= weights.sum()
    order = rng.permutation(len(transcripts))

    def tags(a, r: Read) -> None:
        label = r.alias or r.barcode
        a.set_tag("qs", r.q, "f")
        a.set_tag("du", r.du, "f")
        a.set_tag("st", r.st, "Z")
        a.set_tag("fn", f"{r.run.flow_cell_id}_{r.run.run8}_{r.run.acq8}_{r.pod5}.pod5", "Z")
        a.set_tag("RG", _rg_id(r.run, model) + (f"_{KIT}_{r.barcode}" if r.barcode else ""), "Z")
        if label:
            a.set_tag("BC", label, "Z")
        if r.pt is not None:
            a.set_tag("pt", r.pt, "i")

    for k, r in enumerate(reads):
        g = exp.g(r.run.run_key, r.alias or r.barcode)
        g.setdefault("tx", [])
        a = pysam.AlignedSegment(header)
        a.query_name = r.id
        a.query_sequence = r.seq
        a.query_qualities = array("B", [20]) * r.length
        g["reads"] += 1; g["bases"] += r.length; g["lengths"].append(r.length); g["quals"].append(r.q)
        g["tagged"].append((clock_bin(r.epoch), r.length, r.q, r.du, r.pt))
        if rng.random() >= mapped_frac:
            a.flag = 4
            tags(a, r)
            recs.append(((1 << 40, 0, k), a))
            g["unmapped"] += 1
            exp.records += 1
            continue
        tid = int(order[int(rng.choice(len(transcripts), p=weights))])
        tlen = transcripts[tid][1]
        span = min(r.length, tlen)                           # the read reaches the 3' end and as far 5' as it is long
        back = int(rng.integers(0, min(40, tlen - span) + 1))
        start = tlen - span - back
        clip = r.length - span
        ins = 2 if span > 60 and k % 3 == 0 else 0
        dele = 2 if ins else 0                               # as many deleted as inserted: the alignment still spans `span`
        m = span - ins
        cigar = ([(4, clip)] if clip else []) + [(0, m // 2)] + ([(1, ins), (2, dele)] if ins else []) + [(0, m - m // 2)]
        mism = 0 if k % 6 == 0 else int(rng.integers(0, max(span // 30, 1) + 1))
        nm = mism + ins + dele
        reverse = k % 25 == 0                                # one read in 25 on the other strand
        a.flag = 16 if reverse else 0
        a.reference_id, a.reference_start = tid, start
        a.mapping_quality = int(rng.choice([60, 30, 1, 0]))
        a.cigartuples = cigar
        a.set_tag("NM", nm, "i")
        tags(a, r)
        recs.append(((tid, start, k), a))
        exp.records += 1
        cols = m + ins + dele
        g["mapped"] += 1; g["nm"] += nm; g["cols"] += cols; g["perfect"] += int(nm == 0)
        g["mapq"].append(a.mapping_quality)
        g["identity_q"].append(float("inf") if nm == 0 else -10.0 * float(np.log10(nm / cols)))
        g["tx"].append({"transcript": tid, "length": tlen, "start": start, "end": start + m + dele, "reverse": reverse})
        g["aligned"].append({"length": r.length, "q": r.q, "m": m, "ins": ins, "dele": dele, "skip": 0, "soft": clip, "hard": 0,
                             "mism": mism, "nm": nm, "cols": cols, "reverse": reverse, "contig": tid, "start": start,
                             "blocks": [(start, start + m + dele)], "t10": clock_bin(r.epoch)})
        for s_ in range(int(rng.poisson(secondary_rate))):  # the other transcripts of the gene
            b = pysam.AlignedSegment(header)
            b.query_name = r.id
            b.flag = 256
            other = int(rng.integers(0, len(transcripts)))
            b.reference_id, b.reference_start = other, 0
            b.mapping_quality = 0
            b.cigartuples = [(0, min(r.length, transcripts[other][1]))]
            b.set_tag("NM", 4, "i")
            tags(b, r)
            recs.append(((other, 0, k + 0.5 + s_ / 100), b))
            exp.records += 1
            g["secondary"] += 1
    recs.sort(key=lambda x: x[0])
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for _, a in recs:
            out.write(a)
    if index == "bai":
        pysam.index(path)
    res = exp.done()
    res["transcripts"] = list(transcripts)
    res["genes"] = [i // 3 for i in range(len(transcripts))]
    return res


# ------------------------------------------------------------------------------------------ tables
def write_pychopper_stats(path: str, *, pass_reads: int, qc_fail: int, len_fail: int, primers_found: int,
                          rescued: dict[int, int], unusable: int, umi: int | None = None, cutoff: float = 0.8) -> dict:
    """A pychopper statistics table. ``rescued`` = {segments per read: reads}. The identities of the real files hold:
    Primers_found + rescued reads + Unusable = PassReads; records written = Primers_found + Rescue - LenFail."""
    rescued_reads = sum(rescued.values())
    segments = sum(k * v for k, v in rescued.items())
    if primers_found + rescued_reads + unusable != pass_reads:
        raise ValueError("classes do not add up to PassReads")
    plus = primers_found // 2
    rows = [("ReadStats", "PassReads", pass_reads), ("ReadStats", "LenFail", len_fail), ("ReadStats", "QcFail", qc_fail)]
    if umi is not None:
        rows.append(("ReadStats", "Umi_detected", umi))
    rows += [("Classification", "Primers_found", primers_found), ("Classification", "Rescue", segments), ("Classification", "Unusable", unusable),
             ("Strand", "+", plus), ("Strand", "-", primers_found - plus), ("RescueStrand", "+", segments // 2), ("RescueStrand", "-", segments - segments // 2),
             ("UnclassHitNr", "0", unusable // 3), ("UnclassHitNr", "3", unusable - unusable // 3)]
    rows += [("RescueSegmentNr", str(k), v) for k, v in sorted(rescued.items())]
    hits = {"SSP,-VNP": primers_found // 2, "VNP,-SSP": primers_found - primers_found // 2, "SSP,-VNP,-SSP": unusable // 4, "SSP,-SSP": unusable // 8}
    rows += [("Hits", k, v) for k, v in hits.items()]
    rows += [("Unusable", str(k), 3 + k) for k in range(5)]
    rows.append(("Parameter", "Cutoff(q)", cutoff))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write("Category\tName\tValue\n")
        for c, n, v in rows:
            fh.write(f"{c}\t{n}\t{float(v)}\n")
    return {"pass_reads": pass_reads, "qc_fail": qc_fail, "len_fail": len_fail, "primers_found": primers_found,
            "rescue_segments": segments, "rescued_reads": rescued_reads, "unusable": unusable, "input": pass_reads + qc_fail,
            "written": primers_found + segments - len_fail, "umi_detected": umi}


def write_salmon_meta(sample_dir: str, processed: int, mapped: int) -> str:
    """``<sample_dir>/aux_info/meta_info.json``."""
    p = os.path.join(sample_dir, "aux_info", "meta_info.json")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        json.dump({"salmon_version": "1.10.3", "mapping_type": "alignment", "num_valid_targets": 1000, "num_processed": processed,
                   "num_mapped": mapped, "percent_mapped": 100.0 * mapped / max(processed, 1), "call": "quant"}, fh, indent=4)
    return p


def write_oarfish_meta(path: str, valid_best: int, discard_score: int = 3, discard_supp: int = 11) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump({"prob_model": "logistic_coverage", "alignment_source": "from_raw_reads", "bin_width": 100,
                   "discard_table": {"discard_5p": 0, "discard_3p": 0, "discard_score": discard_score, "discard_aln_frac": 0,
                                     "discard_aln_len": 0, "discard_ori": 0, "discard_supp": discard_supp, "valid_best_aln": valid_best},
                   "num_bootstraps": 0}, fh, indent=2)
    return path


def write_kallisto_info(sample_dir: str, processed: int, pseudoaligned: int) -> str:
    p = os.path.join(sample_dir, "run_info.json")
    os.makedirs(sample_dir, exist_ok=True)
    with open(p, "w") as fh:
        json.dump({"n_targets": 1000, "n_bootstraps": 0, "n_processed": processed, "n_pseudoaligned": pseudoaligned,
                   "n_unique": pseudoaligned // 2, "p_pseudoaligned": round(100.0 * pseudoaligned / max(processed, 1), 1),
                   "kallisto_version": "0.51.1", "call": "kallisto bus --long -x bulk reads.fastq.gz"}, fh, indent="\t")
    return p


def write_sample_sheet(path: str, rows: list[dict], columns: tuple[str, ...] = ("experiment_id", "kit", "flow_cell_id", "barcode", "sample_id", "alias")) -> str:
    """A basecaller sample sheet. Dorado accepts the columns in any order; pass another order to test that."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(columns), extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow(r)
    with open(path, "w", newline="") as fh:
        fh.write(buf.getvalue())
    return path


def write_dorado_summary(path: str, reads: list[Read], *, pass_min_q: float | None = 9.0, barcode_column: bool = True) -> dict:
    """A ``dorado summary`` table (or a basecalled sequencing summary outside a run folder)."""
    exp = _Expect()
    cols = ["filename", "read_id", "run_id", "channel", "mux", "start_time", "duration", "template_start", "template_duration",
            "sequence_length_template", "mean_qscore_template"]
    if pass_min_q is not None:
        cols.append("passes_filtering")
    if barcode_column:
        cols.append("barcode")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for k, r in enumerate(reads):
            ok = pass_min_q is None or r.q >= pass_min_q
            row = [f"{r.run.flow_cell_id}_{r.run.run8}_{r.run.acq8}_{r.pod5}.pod5", r.id, r.run.protocol_run_id, r.channel, 1 + r.channel % 4,
                   round(k * 0.5, 3), round(r.length / 380.0, 4), round(k * 0.5 + 0.1, 3), round(r.length / 380.0, 4), r.length, repr(r.q)]
            if pass_min_q is not None:
                row.append("TRUE" if ok else "FALSE")
            if barcode_column:
                row.append(r.alias or r.barcode or "unclassified")
            fh.write("\t".join(str(x) for x in row) + "\n")
            exp.records += 1
            g = exp.g(r.run.run_key, (r.alias or r.barcode) if barcode_column else None)
            g["reads"] += 1; g["bases"] += r.length; g["lengths"].append(r.length); g["quals"].append(r.q)
            g.setdefault("pass", 0); g.setdefault("fail", 0)
            g["pass" if ok else "fail"] += 1
    return exp.done()


# ----------------------------------------------------------------------------------- a whole project
@dataclass
class ProjectFixture:
    """A fake project folder beside a run folder, and what a correct reader must find."""
    project_dir: str
    run: RunRef
    aliases: list[str]
    barcodes: list[str]
    sample_sheet: str
    files: dict[str, dict[str, str]]            # role -> alias -> path
    expected: dict                              # per alias: the stage counts; per file: what the writer counted


def make_project(project_dir: str, run: RunRef, *, aliases=("alpha_ctrl", "alpha_trt", "beta_ctrl"), reads_per_sample: int = 600,
                 unclassified: int = 150, seed: int = 11, dorado: bool = True, quant: bool = True, sheet_columns=None,
                 basecalled_dir: str = "demo_basecalled", index: str | None = "bai") -> ProjectFixture:
    """One project the way a lab pipeline leaves it:

    ``<project>/samples.csv``; ``<project>/<basecalled_dir>/bam_pass/<alias>/*.bam`` (and ``unclassified``);
    ``<project>/analysis/alignments/<alias>_genome_sorted.bam``; ``analysis/pychopper_reports/<alias>_pychopper_stats.tsv``;
    ``analysis/counts/{salmon,oarfish,kallisto}/...``. The numbers of one sample are consistent from file to file:
    reads in its Dorado BAM = PassReads + QcFail of its pychopper table; records written by pychopper = primary
    records of its genome BAM = reads salmon processed.
    """
    os.makedirs(project_dir, exist_ok=True)
    barcodes = [f"barcode{i + 1:02d}" for i in range(len(aliases))]
    rows = [{"experiment_id": run.experiment, "kit": KIT, "flow_cell_id": run.flow_cell_id, "barcode": b, "sample_id": a, "alias": a}
            for a, b in zip(aliases, barcodes)]
    sheet = write_sample_sheet(os.path.join(project_dir, "samples.csv"), rows, sheet_columns or ("experiment_id", "kit", "flow_cell_id", "barcode", "sample_id", "alias"))
    files: dict[str, dict[str, str]] = {"dorado.bam": {}, "align.genome_bam": {}, "pychopper.stats": {}, "quant.salmon": {}, "quant.oarfish": {}, "quant.kallisto": {}}
    exp: dict = {"samples": {}, "files": {}}
    ana = os.path.join(project_dir, "analysis")
    for i, (alias, bcode) in enumerate(zip(aliases, barcodes)):
        n = reads_per_sample + 40 * i
        reads = make_reads(n, run, seed=seed + i, alias=alias, barcode=bcode, split_every=97)
        if dorado:
            d = os.path.join(project_dir, basecalled_dir, "bam_pass", alias)
            os.makedirs(d, exist_ok=True)
            p = os.path.join(d, f"{run.flow_cell_id}_pass_{alias}_{run.run8}_{run.acq8}_0.bam")
            exp["files"][p] = write_dorado_bam(p, reads)
            files["dorado.bam"][alias] = p
        qc_fail = n // 40
        pass_reads = n - qc_fail
        unusable = pass_reads // 20
        rescued = {2: 4, 3: 1}
        primers_found = pass_reads - unusable - sum(rescued.values())
        len_fail = 6
        py = write_pychopper_stats(os.path.join(ana, "pychopper_reports", f"{alias}_pychopper_stats.tsv"), pass_reads=pass_reads, qc_fail=qc_fail,
                                   len_fail=len_fail, primers_found=primers_found, rescued=rescued, unusable=unusable, umi=primers_found // 2)
        files["pychopper.stats"][alias] = os.path.join(ana, "pychopper_reports", f"{alias}_pychopper_stats.tsv")
        written = py["written"]
        kept = reads[:written]                               # the genome BAM holds exactly the records pychopper wrote
        os.makedirs(os.path.join(ana, "alignments"), exist_ok=True)
        gp = os.path.join(ana, "alignments", f"{alias}_genome_sorted.bam")
        exp["files"][gp] = write_aligned_bam(gp, kept, seed=seed + 100 + i, index=index)
        files["align.genome_bam"][alias] = gp
        mapped = exp["files"][gp]["groups"][(run.run_key, alias)]["mapped"]
        quantified = mapped - 7
        if quant:
            files["quant.salmon"][alias] = write_salmon_meta(os.path.join(ana, "counts", "salmon", f"{alias}_salmon"), written, quantified)
            files["quant.oarfish"][alias] = write_oarfish_meta(os.path.join(ana, "counts", "oarfish", f"{alias}.meta_info.json"), quantified - 2)
            files["quant.kallisto"][alias] = write_kallisto_info(os.path.join(ana, "counts", "kallisto", "kallisto_quant", alias), written, int(written * 0.7))
        exp["samples"][alias] = {"barcode": bcode, "classified": n, "pychopper": py, "full_length": written, "genome_mapped": mapped,
                                 "quantified": quantified, "bases": sum(r.length for r in reads)}
    if dorado and unclassified:
        reads = make_reads(unclassified, run, seed=seed + 50, q_mean=9.0)
        d = os.path.join(project_dir, basecalled_dir, "bam_pass", "unclassified")
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, f"{run.flow_cell_id}_pass_unclassified_{run.run8}_{run.acq8}_0.bam")
        exp["files"][p] = write_dorado_bam(p, reads)
        files["dorado.bam"]["unclassified"] = p
        exp["unclassified"] = unclassified
    exp["basecalled"] = sum(s["classified"] for s in exp["samples"].values()) + (unclassified if dorado else 0)
    return ProjectFixture(project_dir, run, list(aliases), barcodes, sheet, files, exp)


def main(argv: list[str] | None = None) -> int:
    import argparse

    from flowcellar.testing.fixtures import RunSpec, make_run
    p = argparse.ArgumentParser(prog="python -m flowcellar.testing.readdata",
                                description="Write a fake run folder and a fake project with read-data files (invented names only).")
    p.add_argument("out_dir")
    p.add_argument("--reads", type=int, default=600, help="reads per sample")
    p.add_argument("--seed", type=int, default=11)
    a = p.parse_args(argv)
    run = make_run(os.path.join(a.out_dir, "project_one"), RunSpec(experiment="demo_experiment", sample="pool_1", seed=a.seed))
    proj = make_project(os.path.join(a.out_dir, "project_one"), run_ref(run), reads_per_sample=a.reads, seed=a.seed)
    print(json.dumps({"run_key": run.run_key, "run_dir": run.run_dir, "project_dir": proj.project_dir,
                      "samples": {k: {x: v[x] for x in ("classified", "full_length", "genome_mapped", "quantified")} for k, v in proj.expected["samples"].items()}}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
