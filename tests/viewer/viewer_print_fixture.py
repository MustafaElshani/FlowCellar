"""What the print check opens (viewer/verify/print.mjs), made by the engine itself from invented data.

    python tests/viewer/viewer_print_fixture.py OUT_DIR [--samples 3,13]

    OUT_DIR/run<N>/out/<run key>.flowcellar.html   a one-run report whose run page has every section: read fate, samples,
                                                  read length and quality, alignment, reads over run time, transcripts,
                                                  Pychopper and RNA composition, for N samples
    OUT_DIR/synth/library/hub/                     the hub of the synthetic library (a comparison of three runs, a flow
                                                  cell with a history)
    OUT_DIR/spatial/library/hub/                   a hub with one run and four single-cell and spatial analyses
    OUT_DIR/spatial/pair.flowcellar.html           the one-file report of a Visium HD sample analysed by both pipelines

Prints one line of JSON that names them. The spatial part needs h5py; without it the line says ``"spatial": null``.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

NAMES = ("ctrl", "treated", "knockdown", "rescue")


def make_run_report(base: Path, n: int, reads: int = 400) -> str:
    """A fixture run with ``n`` samples and every kind of read data the run page has a section for -> the report."""
    from flowcellar.testing import libprep as F
    from flowcellar.testing import readdata as T
    from flowcellar.testing.fixtures import RunSpec, make_run

    fr = make_run(base / "project", RunSpec(experiment="demo_experiment", sample=f"pool_{n}", seed=500 + n))
    ref = T.run_ref(fr)
    aliases = tuple(f"{NAMES[i % 4]}_rep{i // 4 + 1}" for i in range(n))
    proj = T.make_project(str(base / "project"), ref, aliases=aliases, reads_per_sample=reads, seed=61, quant=False)
    ana = base / "project" / "analysis"
    args = []
    for i, alias in enumerate(aliases):
        rs = T.make_reads(reads, ref, seed=80 + i, alias=alias, barcode=proj.barcodes[i])
        p = ana / "alignments" / f"{alias}_transcriptome_sorted.bam"
        T.write_transcriptome_bam(str(p), rs, seed=90 + i)
        args += ["--input", f"align.transcriptome_bam={p}"]
        spec = F.PychopperSpec(primers_found=40000 + 3100 * (i % 7), plus=19000 + 1700 * (i % 5), rescued={2: 60 + 9 * i, 3: 5 + i % 4, 6: i % 2},
                               unusable_hits={0: 300 + 40 * (i % 6), 2: 60, 3: 500 + 90 * (i % 4), 4: 20, 7: 2, 30: 1}, qc_fail=400 + 35 * i, len_fail=90 + 7 * i,
                               umi=(38000 + 2900 * (i % 7), 37500 + 2900 * (i % 7)), cutoff_index=1 if i != 2 else 0, seed=200 + i)
        F.write_pychopper_stats_full(str(ana / "pychopper_reports" / f"{alias}_pychopper_stats.tsv"), spec)
        F.write_annotated_table(str(ana / "dorado_polyA" / "annotated_readid" / f"{alias}_readid_annotated.csv"),
                                F.make_table_rows(3000 + 400 * (i % 5), seed=300 + i, not_primary=3))
    out = base / "out"
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    r = subprocess.run([sys.executable, "-m", "flowcellar", "run", fr.run_dir, str(base / "project"), *args, "--sample-sheet", proj.sample_sheet, "-o", str(out)],
                       capture_output=True, text=True, timeout=1200, env=env)
    reports = sorted(out.glob("*.flowcellar.html")) if out.is_dir() else []
    if not reports:
        raise RuntimeError("flowcellar run wrote no report:\n" + r.stdout + r.stderr)
    return str(reports[0])


def make_spatial(base: Path) -> dict | None:
    """A library with one run and four analyses (Visium and Visium HD, each by both pipelines), and a pair report."""
    try:
        import h5py  # noqa: F401
    except ImportError:
        return None
    import base64
    import gzip
    import re

    from flowcellar.library import scan as scan_mod
    from flowcellar.sc import report as sc_report
    from flowcellar.sc.store import AnalysisStore
    from flowcellar.testing import sc as fx

    tree = fx.make_spatial_tree(str(base / "data"))
    lib_dir = base / "library"
    res = scan_mod.scan(str(lib_dir), scan_mod.ScanOptions(roots=(tree.roots["seq"], tree.roots["store"]), wait=False))
    if res.exit_code != 0:
        raise RuntimeError("the scan of the spatial fixture failed")
    text = (lib_dir / "hub" / "index.js").read_text("utf-8")
    ix = json.loads(gzip.decompress(base64.b64decode(re.fullmatch(r'__FC_INDEX__\("([A-Za-z0-9+/=]*)"\);\n', text).group(1))))
    rows = {(r["pipeline"], r["assay"]): r for r in ix["analyses"]}
    store = AnalysisStore(str(lib_dir))
    wf, sr = rows[("wf-single-cell", "visium_hd")], rows[("spaceranger", "visium_hd")]
    report = sc_report.write_report([store.composed(wf["id"]), store.composed(sr["id"])], str(base / "pair.flowcellar.html"))
    return {"hub": str(lib_dir / "hub"), "report": str(report), "analysis": sr["id"], "workflow": wf["id"]}


def make(out: str, samples: tuple[int, ...] = (3, 13)) -> dict:
    import viewer_make_fixture

    base = Path(os.path.abspath(out))
    shutil.rmtree(base, ignore_errors=True)
    base.mkdir(parents=True)
    runs = {str(n): make_run_report(base / f"run{n}", n, 1200 if n <= 4 else 200) for n in samples}
    synth = viewer_make_fixture.make(str(base / "synth"))
    return {"out": str(base), "runs": runs, "hub": synth["hub"], "report": synth["reports"]["TST00001_0a1b2c3d"], "spatial": make_spatial(base / "spatial")}


if __name__ == "__main__":
    argv = sys.argv[1:]
    counts = (3, 13)
    if "--samples" in argv:
        i = argv.index("--samples")
        counts = tuple(int(x) for x in argv[i + 1].split(",") if x)
        del argv[i:i + 2]
    if len(argv) != 1:
        sys.exit(__doc__)
    print(json.dumps(make(argv[0], counts)))
