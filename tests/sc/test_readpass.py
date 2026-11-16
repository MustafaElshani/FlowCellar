"""The pass over the per-read table is exact, however the file is cut into ranges and blocks."""
from collections import Counter

import numpy as np
import pytest

from flowcellar.sc import axes, readpass
from flowcellar.testing import sc as fx


def _table(tmp_path, hd: bool, n: int = 6000):
    import random
    rng = random.Random(11)
    names = fx.gene_names()
    if hd:
        bcs = [f"s_002um_{r:05d}_{c:05d}-1" for r in range(0, 40, 3) for c in range(0, 40, 3)]
        rows = fx._reads(rng, bcs, names, n, 9)
    else:
        rows = fx._reads(rng, [fx._barcode(rng) for _ in range(30)], names, n, 12)
    path = tmp_path / ("hd.tsv" if hd else "spots.tsv")
    fx._write_table(str(path), rows)
    return str(path), rows


@pytest.mark.parametrize("hd", [False, True])
@pytest.mark.parametrize("cut", [dict(inprocess=True, range_bytes=1 << 30, block_bytes=1 << 30),
                                 dict(inprocess=True, range_bytes=20011, block_bytes=3001),
                                 dict(workers=3, range_bytes=90001, block_bytes=40009)])
def test_every_count_is_exact(tmp_path, hd, cut):
    path, rows = _table(tmp_path, hd)
    r = readpass.read_pass(path, **cut)
    assert r["complete"]
    assert r["rows"] == r["valid_barcode"] == len(rows)
    assert r["valid_umi"] == sum(1 for x in rows if x[2] != "-")
    assert r["gene_assigned"] == sum(1 for x in rows if x[3] != "-")
    assert r["transcript_assigned"] == sum(1 for x in rows if x[4] != "-")
    mol = Counter((x[1], x[3], x[2]) for x in rows if x[2] != "-" and x[3] != "-")
    assert r["molecules"] == len(mol) and r["molecule_reads"] == sum(mol.values())
    assert r["molecule_hist"].tolist() == np.bincount(list(mol.values())).tolist()
    unit = (lambda b: (int(b[8:13]) // 4, int(b[14:19]) // 4)) if hd else (lambda b: b)
    per = Counter(unit(x[1]) for x in rows)
    assert sorted(readpass.reads_per_unit(r).tolist()) == sorted(per.values())
    assert r["sketch"]["ids"] == sorted({x[0].replace("-", "")[:16] for x in rows if x[0].startswith(r["sketch"]["prefix"])})
    assert len(r["sketch"]["ids"]) > 20


def test_a_pass_that_ran_out_of_time_says_so(tmp_path):
    path, rows = _table(tmp_path, False)
    r = readpass.read_pass(path, inprocess=True, range_bytes=20011, block_bytes=3001, deadline=0.0)
    assert not r["complete"] and r["rows"] < len(rows)


def test_a_table_of_another_kind_is_refused(tmp_path):
    from flowcellar.util.errors import InputError
    p = tmp_path / "x.tsv"
    p.write_text("a\tb\n1\t2\n")
    with pytest.raises(InputError):
        readpass.read_pass(str(p), inprocess=True)


def test_saturation_from_the_histogram_is_the_exact_expectation():
    counts = [1] * 500 + [2] * 300 + [3] * 100 + [10] * 20 + [70] * 5 + [300] * 2
    h = axes.reads_per_molecule_hist(np.bincount(counts))
    assert sum(h["n"]) == len(counts) and sum(h["w"]) == sum(counts)
    assert axes.saturation_from_hist(h) == pytest.approx(1 - len(counts) / sum(counts), abs=1e-12)
    # half the reads: the expected number of molecules still seen, molecule by molecule
    p = 0.5
    unique = sum(1 - (1 - p) ** c for c in counts)
    assert axes.saturation_from_hist(h, p) == pytest.approx(1 - unique / (p * sum(counts)), abs=2e-4)


def test_the_two_sketch_prefixes():
    assert readpass.sketch_prefix(10 << 20) == "0000" and readpass.sketch_prefix(30 << 30) == "00000"
