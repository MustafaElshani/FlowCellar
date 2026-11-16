"""The library-preparation tables through a whole scan: found by configured patterns, linked through the sample
sheet, the annotated tables read once and cached, the sections in a bundle that validates, the index column filled."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys

import pytest

from flowcellar.aggregates import merge
from flowcellar.model import bundle as bundle_model
from flowcellar.model import schema
from flowcellar.testing import libprep as F
from flowcellar.testing import readdata as T
from flowcellar.testing.fixtures import RunSpec, make_run

CONFIG = """\
[library]
name = "library preparation demo"

[discovery]
min_age_seconds = 0

[[inputs]]
role = "pychopper.stats"
glob = "**/analysis/pychopper_reports/{sample}_pychopper_stats.tsv"
[[inputs]]
role = "polya.read_table"
glob = "**/analysis/dorado_polyA/annotated_readid/{sample}_readid_annotated.csv"
[[inputs]]
role = "basecaller.sample_sheet"
sniff = "csv_columns:barcode,alias,flow_cell_id|position_id"

[[link.rules]]
name = "experiment-folder"
roles = ["pychopper.stats", "polya.read_table"]
via = "sample"
within = "experiment"
"""


def flowcellar(*args: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    return subprocess.run([sys.executable, "-m", "flowcellar", *args], capture_output=True, text=True, timeout=600, env=env)


@pytest.fixture(scope="module")
def scanned(tmp_path_factory):
    pytest.importorskip("pysam", reason="the fixture project writes BAM files (extra 'bam')")
    base = tmp_path_factory.mktemp("lp_scan")
    root = base / "store"
    fr = make_run(root / "project_one", RunSpec(experiment="demo_experiment", sample="pool_1", seed=81))
    proj = T.make_project(str(root / "project_one"), T.run_ref(fr), seed=83, quant=False, reads_per_sample=200)
    py, tab = {}, {}
    for i, alias in enumerate(proj.aliases):
        py[alias] = F.write_pychopper_stats_full(proj.files["pychopper.stats"][alias], F.PychopperSpec(primers_found=4000 + 500 * i, plus=2000, seed=90 + i))
        path = os.path.join(proj.project_dir, "analysis", "dorado_polyA", "annotated_readid", f"{alias}_readid_annotated.csv")
        tab[alias] = F.write_annotated_table(path, F.make_table_rows(900 + 200 * i, seed=95 + i, not_primary=2))
    lib = base / "library"
    lib.mkdir()
    (lib / "flowcellar.toml").write_text(CONFIG)
    first = flowcellar("scan", str(root), "-L", str(lib), "--no-hub")
    return {"lib": str(lib), "root": str(root), "fr": fr, "proj": proj, "py": py, "tab": tab, "first": first}


def _bundle(lib: str, run_key: str) -> dict:
    return bundle_model.read_bundle(os.path.join(lib, "bundles", bundle_model.bundle_filename(run_key)))


def test_scan_reads_the_tables_and_writes_valid_sections(scanned):
    first = scanned["first"]
    assert first.returncode == 0, first.stdout + first.stderr
    b = _bundle(scanned["lib"], scanned["fr"].run_key)
    assert schema.schema_errors(b) == []
    mods = b["modules"]
    assert {"cdna", "cdna.pychopper", "rna.biotype", "rna.genes"} <= set(mods)
    name_of = {s["key"]: s["name"] for s in b["samples"]}
    for c in mods["cdna.pychopper"]["cells"]:
        t = scanned["py"][name_of[c["at"][0]]]["scalars"]
        assert {k: c["scalars"][k] for k in t} == t
    assert mods["cdna.pychopper"]["state"] == "ok" and mods["rna.biotype"]["state"] == "ok" and mods["rna.genes"]["state"] == "ok"
    per_sample = {name_of[c["at"][0]]: c["scalars"] for c in merge.rollup(mods["rna.biotype"]["cells"], ["sample", "read_class"], ["sample"])}
    for alias, t in scanned["tab"].items():
        assert per_sample[alias]["reads"] == t["primary"]
        assert per_sample[alias]["tail_estimated"] == sum(e["tail_estimated"] for e in t["biotypes"].values())
    detected = {name_of[r[0]]: r[3] for r in mods["rna.genes"]["run_level"]["tables"]["rna_genes_detected"]["rows"]}
    assert detected == {a: t["genes_ge10"] for a, t in scanned["tab"].items()}
    # the axes the new histograms use are described in the bundle, and the inputs are named with their link
    assert {"polya_5nt", "count_1", "trimmed_bases_2"} <= set(b["axes"])
    ins = [i for i in b["provenance"]["inputs"] if i["role"] == "polya.read_table"]
    assert len(ins) == 3 and all(i["state"] == "ok" and i["link"]["method"] == "sample_name_chain" and i["relpath"].endswith("_readid_annotated.csv") for i in ins)
    assert sorted(b["provenance"]["partials"]) == b["provenance"]["partials"] and len(b["provenance"]["partials"]) >= 3
    # the index column
    with gzip.open(os.path.join(scanned["lib"], "index", "index.json.gz")) as fh:
        idx = json.load(fh)
    assert idx["runs"]["columns"]["t2_full_length"] == [mods["cdna.pychopper"]["summary"]["cdna_full_length_frac"]]


def test_second_scan_reads_no_table_again(scanned):
    lib = scanned["lib"]
    partials = sorted(os.path.join(d, f) for d, _, fs in os.walk(os.path.join(lib, "partials")) for f in fs)
    stamp = {p: os.stat(p).st_mtime_ns for p in partials}
    before = _bundle(lib, scanned["fr"].run_key)
    logs_before = set(os.listdir(os.path.join(lib, "logs")))
    second = flowcellar("scan", "-L", lib, "--no-hub")
    assert second.returncode == 0, second.stdout + second.stderr
    assert {p: os.stat(p).st_mtime_ns for p in partials} == stamp and len(partials) >= 3     # no partial was written again
    new_logs = sorted(set(os.listdir(os.path.join(lib, "logs"))) - logs_before)
    assert len(new_logs) == 1
    with open(os.path.join(lib, "logs", new_logs[0])) as fh:
        events = [json.loads(line) for line in fh if line.strip()]
    assert not [e for e in events if str(e.get("event") or "").startswith("read_pass")], "the second scan started a pass over a table"
    after = _bundle(lib, scanned["fr"].run_key)
    assert after["modules"]["rna.biotype"] == before["modules"]["rna.biotype"] and after["content_sha256"] == before["content_sha256"]


def test_validate_and_show_know_the_new_sections(scanned):
    v = flowcellar("validate", "-L", scanned["lib"], os.path.join(scanned["lib"], "bundles", bundle_model.bundle_filename(scanned["fr"].run_key)))
    assert v.returncode == 0, v.stdout + v.stderr
    s = flowcellar("show", scanned["fr"].run_key, "-L", scanned["lib"], "--section", "cdna.pychopper", "--json")
    assert s.returncode == 0, s.stdout + s.stderr
    assert json.loads(s.stdout)["summary"]["cdna_pychopper_samples"] == 3
