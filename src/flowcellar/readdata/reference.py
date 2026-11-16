"""What the header of an aligned BAM says about its reference, and the layouts the read pass bins on.

Everything here is decided from the ``@SQ`` lines alone, never from a file name:

* **genome or transcriptome**: a transcriptome is told by its reference names (transcript identifiers) or, without
  them, by its shape (a thousand sequences or more, none of them a megabase long). Anything else is a genome;
* **windows along a genome**: one window length for the whole reference, chosen so that the genome is about 3,000
  windows (``TARGET_WINDOWS``), rounded up to two significant digits. Every contig that is kept starts a new window,
  so its last window is shorter. The contigs kept by name are the ``MAX_CONTIGS`` longest and the mitochondrial
  genome; the others are added up under "other contigs" and get no window;
* **transcript length classes** for the coverage along the transcript body.

Two files aligned to the same reference get the same layout and the same ``ref_id``, so their numbers add up.
"""
from __future__ import annotations

import hashlib
import math
import re

import numpy as np

TARGET_WINDOWS = 3000
MAX_CONTIGS = 64
OTHER_CONTIGS = "other contigs"
#: names a mitochondrial genome goes by, lower-cased
MITO_NAMES = frozenset(("chrm", "mt", "chrmt", "m", "mito", "chrmito", "mitochondrion", "mitochondrion_genome", "mtdna"))
#: (id, upper bound of the transcript length in bases)
TX_CLASSES = (("under_1kb", 1000), ("1_to_2kb", 2000), ("2_to_4kb", 4000), ("4kb_and_longer", None))
BODY_BINS = 100
_TX_NAME = re.compile(r"^(?:ENS[A-Z]{0,4}T\d|[NX][MR]_\d|MSTRG\.\d|STRG\.\d|TCONS_\d|transcript[:_.-])|\|ENS[A-Z]{0,4}G\d")
_GENE_FIELD = re.compile(r"^ENS[A-Z]{0,4}G\d")
NAME_SAMPLE = 2000


def classify(first_names: list[str], lengths: np.ndarray) -> tuple[str | None, str]:
    """("genome" | "transcriptome" | None, the evidence in words)."""
    n = len(lengths)
    if n == 0:
        return None, "the header lists no reference sequence"
    hits = sum(1 for s in first_names if _TX_NAME.search(s))
    if first_names and hits * 2 >= len(first_names):
        return "transcriptome", f"{hits:,} of the first {len(first_names):,} reference names are transcript identifiers"
    longest = int(lengths.max())
    if n >= 1000 and longest < 1_000_000 and float(np.median(lengths)) < 20_000:
        return "transcriptome", f"{n:,} reference sequences, the longest {longest:,} bases: too many and too short to be chromosomes"
    return "genome", f"{n:,} reference sequences, the longest {longest:,} bases"


def window_length(total: int) -> int:
    """The window length that makes ``total`` bases about ``TARGET_WINDOWS`` windows, rounded up to two digits."""
    w = max(1, math.ceil(total / TARGET_WINDOWS))
    if w < 100:
        return w
    mag = 10 ** (len(str(w)) - 2)
    return int(math.ceil(w / mag) * mag)


class Layout:
    """The reference of one BAM as the pass needs it. Built once per process from the header."""

    def __init__(self, names_of, lengths):
        """``names_of(i)`` gives the name of reference ``i`` (asked only for the few names that are needed)."""
        self.lengths = np.asarray(lengths, dtype=np.int64)
        self.n = n = int(len(self.lengths))
        self.total = int(self.lengths.sum()) if n else 0
        self.kind, self.evidence = classify([names_of(i) for i in range(min(n, NAME_SAMPLE))], self.lengths)
        h = hashlib.blake2b(digest_size=6)
        h.update(self.lengths.tobytes())
        h.update("|".join(names_of(i) for i in range(min(n, 64))).encode())
        self.ref_id = h.hexdigest()
        self.window = self.n_windows = 0
        self.kept: list[int] = []
        self.names: list[str] = []
        self.mito: int | None = None
        self.slot = self.win_first = self.tx_class = None
        self.win_first_list: list[int] = []
        if self.kind == "genome":
            self._genome(names_of)
        elif self.kind == "transcriptome":
            bounds = np.array([b for _, b in TX_CLASSES[:-1]], dtype=np.int64)
            self.tx_class = np.searchsorted(bounds, self.lengths, side="right").astype(np.int64)

    def _genome(self, names_of) -> None:
        n = self.n
        order = np.argsort(-self.lengths, kind="stable")[:MAX_CONTIGS]
        keep = set(int(i) for i in order)
        # the mitochondrial genome is short and would be lost among the unplaced contigs: look for it by name
        if n <= 100_000:
            for i in range(n):
                if int(self.lengths[i]) < 200_000 and names_of(i).lower() in MITO_NAMES:
                    keep.add(i)
                    break
        self.kept = sorted(keep)
        self.names = [names_of(i) for i in self.kept]
        self.window = window_length(self.total)
        self.slot = np.full(n, len(self.kept), dtype=np.int64)
        self.win_first = np.full(n, -1, dtype=np.int64)
        first = 0
        for k, i in enumerate(self.kept):
            self.slot[i] = k
            self.win_first[i] = first
            first += int(-(-int(self.lengths[i]) // self.window))
            if self.names[k].lower() in MITO_NAMES and int(self.lengths[i]) < 200_000:
                self.mito = k
        self.n_windows = first
        self.win_first_list = [int(v) for v in self.win_first]

    def describe(self) -> dict:
        """Plain data for the partial and the bundle (section ``params.reference``)."""
        d: dict = {"kind": self.kind, "evidence": self.evidence, "id": self.ref_id, "sequences": self.n, "total_length": self.total}
        if self.kind == "genome":
            d["window_bases"] = self.window
            d["windows"] = self.n_windows
            d["contigs"] = [{"name": self.names[k], "length": int(self.lengths[i]), "first_window": int(self.win_first[i]),
                             "windows": int(-(-int(self.lengths[i]) // self.window))} for k, i in enumerate(self.kept)]
            kept_len = int(self.lengths[self.kept].sum()) if self.kept else 0
            d["other"] = {"sequences": self.n - len(self.kept), "length": self.total - kept_len}
            d["mitochondrial"] = self.names[self.mito] if self.mito is not None else None
        elif self.kind == "transcriptome":
            d["length_classes"] = [{"id": cid, "below_bases": hi, "transcripts": int((self.tx_class == k).sum())}
                                   for k, (cid, hi) in enumerate(TX_CLASSES)]
        return d


def gene_index(names) -> tuple[np.ndarray, int] | None:
    """The gene of every transcript, as an index in order of first appearance, when the reference names carry it
    (``transcript|gene|...`` as GENCODE writes them); None otherwise."""
    if not len(names):
        return None
    first = names[0].split("|")
    if len(first) < 2 or not _GENE_FIELD.match(first[1]):
        return None
    ids: dict[str, int] = {}
    out = np.empty(len(names), dtype=np.int64)
    for i, s in enumerate(names):
        parts = s.split("|", 2)
        out[i] = ids.setdefault(parts[1] if len(parts) > 1 else s, len(ids))
    return out, len(ids)
