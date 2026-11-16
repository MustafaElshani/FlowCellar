"""The discovery walk: run folders are recognised by their files, read-data files by configured patterns, and a
root that does not answer is skipped. Listings and stat calls only."""
from __future__ import annotations

import os
import time

import pytest

from flowcellar import discovery
from flowcellar.config import loader
from flowcellar.contracts import InputRule
from flowcellar.discovery import guard
from flowcellar.discovery.patterns import PatternError, compile_rules, parse_glob
from flowcellar.discovery.walk import WalkSettings, classify, walk_root
from flowcellar.testing import trees
from flowcellar.testing.fixtures import RunSpec, make_run


def config_for(tree, **discovery_keys) -> dict:
    raw = {"roots": [{"id": "seq", "path": tree.roots["seq"]}, {"id": "store", "path": tree.roots["store"], "project_depth": 2},
                     {"id": "arch", "path": tree.roots["arch"]}]}
    if discovery_keys:
        raw["discovery"] = discovery_keys
    return loader.merge(raw).config


def touch(path, text="x", age=3600.0):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    t = time.time() - age
    os.utime(path, (t, t))
    return path


# ------------------------------------------------------------------------------------------------ classes of folders
def test_the_walk_finds_runs_raw_only_folders_and_look_alikes(tree):
    found = discovery.discover(config_for(tree))
    assert all(o.state == "reachable" for o in found.roots.values())
    runs = {f.path: (root, f) for root, f in found.runs}
    assert sorted(p for p, (_, f) in runs.items() if f.kind == "run") == sorted(r.run_dir for r in tree.runs.values())
    assert [p for p, (_, f) in runs.items() if f.kind == "run_raw_only"] == list(tree.raw_only.values())
    assert [f.path for _, f in found.lookalikes] == tree.lookalikes
    assert not set(tree.hidden) & set(runs), "a folder named 'work' is never entered"
    a1 = runs[tree.runs["a1"].run_dir]
    assert a1[0] == "seq" and a1[1].relpath.startswith("alpha_project/pool_a/")
    assert any(e.startswith("final_summary_") for e in a1[1].evidence), "the evidence names the files that decided"
    assert len(a1[1].files) == 8 and all(size > 0 and mtime > 0 for _, size, mtime in a1[1].files)
    assert any(name.startswith("other_reports/pore_scan_data_") for name, _, _ in a1[1].files)
    c = found.counts()
    assert (c["runs"], c["raw_only"], c["lookalikes"]) == (5, 1, 1) and c["pruned"] == 1 and c["listings"] > 10


def test_a_run_is_recognised_by_its_files_not_by_its_name(tmp_path):
    run = make_run(tmp_path / "root", RunSpec(reads=20, hours=1.0))
    renamed = os.path.join(os.path.dirname(run.run_dir), "my renamed folder")
    os.rename(run.run_dir, renamed)
    named_like_a_run = tmp_path / "root" / "x" / "20260101_0000_MD-000001_TST00077_aaaaaaaa"
    named_like_a_run.mkdir(parents=True)
    (named_like_a_run / "notes.txt").write_text("nothing of MinKNOW here")
    w = walk_root("r", str(tmp_path / "root"))
    assert [(f.kind, f.path) for f in w.runs] == [("run", renamed)]
    assert w.lookalikes == []


def test_the_folder_above_a_run_is_not_taken_for_raw_data(tmp_path):
    """A project folder that holds a script called fast5_to_pod5.sh and a run folder: the run inside must be found."""
    root = tmp_path / "root"
    run = make_run(root / "project_x", RunSpec(reads=20, hours=1.0))
    sample_dir = os.path.dirname(run.run_dir)
    touch(os.path.join(sample_dir, "fast5topod5.sh"))
    touch(os.path.join(sample_dir, "pod5_archive.tar.gz"))
    w = walk_root("r", str(root))
    assert [(f.kind, f.path) for f in w.runs] == [("run", run.run_dir)]
    entries = [("fast5topod5.sh", False, True), (os.path.basename(run.run_dir), True, False)]
    assert classify("some_project", entries) is None
    assert classify("20260101_0000_MD-000001_TST00077_aaaaaaaa", [("pod5", True, False)]) == "run_raw_only"
    assert classify("anything", [("pod5", True, False), ("TST00077_aaaaaaaa_bbbbbbbb_0.pod5", False, True)]) == "run_raw_only"
    assert classify("anything", [("pod5", True, False)]) is None, "raw data without a run key is not a run"


def test_symlinks_are_never_followed_and_hidden_names_are_ignored(tmp_path):
    real = make_run(tmp_path / "elsewhere", RunSpec(reads=20, hours=1.0))
    root = tmp_path / "root"
    root.mkdir()
    os.symlink(os.path.dirname(os.path.dirname(real.run_dir)), root / "link_to_runs")
    hidden = make_run(root / ".snapshots", RunSpec(reads=20, hours=1.0, seed=5, flow_cell_id="TST00005"))
    w = walk_root("r", str(root))
    assert w.runs == [] and w.counts["symlinks"] == 1 and hidden.run_dir


def test_the_configured_prune_list_is_added_to_the_built_in_one(tree):
    found = discovery.discover(config_for(tree, prune=["beta_*"]))
    keys = {os.path.basename(f.path)[-17:] for _, f in found.runs}
    assert keys == {tree.runs[k].run_key for k in ("a1", "a2", "a3")} | set(tree.raw_only)
    assert found.counts()["pruned"] == 2


def test_the_safety_depth_and_the_listing_budget(tmp_path):
    root = tmp_path / "root"
    deep = root / "top"
    for i in range(14):
        deep = deep / f"level{i}"
    run = make_run(deep, RunSpec(reads=20, hours=1.0))
    assert run.run_dir.count("/") - str(root).count("/") > 12
    w = walk_root("r", str(root))
    assert w.runs == [] and w.counts["depth_capped"] == 1, "nothing is looked for more than twelve folders down"
    trees.make_deep_tree(root, folders=60, name="wide")
    ok = make_run(root / "fine", RunSpec(reads=20, hours=1.0, seed=3, flow_cell_id="TST00003"))
    w = walk_root("r", str(root), settings=WalkSettings(budget=20))
    assert w.over_budget == [str(root / "wide")], "the subtree that was given up is named"
    assert [f.path for f in w.runs] == [ok.run_dir], "and the other first-level folders are still walked"
    assert walk_root("r", str(root), settings=WalkSettings(budget=1000)).over_budget == []


def test_a_root_inside_another_root_is_walked_once(tree):
    inner = os.path.join(tree.roots["seq"], "alpha_project")
    cfg = loader.merge({"roots": [{"id": "seq", "path": tree.roots["seq"]}, {"id": "alpha", "path": inner}]}).config
    found = discovery.discover(cfg)
    by_root = {}
    for rid, f in found.runs:
        by_root.setdefault(rid, []).append(f.path)
    assert len(by_root["alpha"]) == 2 and len(by_root["seq"]) == 2
    assert len({f.path for _, f in found.runs}) == 4


def test_a_root_that_is_itself_one_run_folder(tree):
    run = tree.runs["b1"]
    w = walk_root("one", run.run_dir)
    assert [(f.kind, f.path, f.relpath) for f in w.runs] == [("run", run.run_dir, "")]


# ------------------------------------------------------------------------------------------------ patterns
def test_globs_with_captures():
    parts = parse_glob("**/analysis/alignments/{sample}_sorted.bam")
    assert [p.kind for p in parts] == ["deep", "literal", "literal", "pattern"]
    assert parts[-1].match("liver_1_sorted.bam", {}) == {"sample": "liver_1"}
    assert parts[-1].match("liver_1_sorted.bam.bai", {}) is None
    v = parse_glob("x/{variant:main|alt.1}/{sample}.bam")
    assert v[1].match("alt.1", {}) == {"variant": "alt.1"} and v[1].match("altx1", {}) is None and v[1].match("other", {}) is None
    assert parse_glob("a/[!x]?c*.txt")[1].match("ybcd.txt", {}) == {} and parse_glob("a/[!x]?c*.txt")[1].match("xbc.txt", {}) is None
    assert parse_glob("{sample}/{sample}.bam")[1].match("b.bam", {"sample": "a"}) is None, "a capture must repeat its value"
    assert [p.kind for p in parse_glob("./a//**/**/b.txt")] == ["literal", "deep", "literal"]
    for bad in ("", "a/**", "a/b**/c", "a/{bad name}/x", "../x", "a/{v:}/x"):
        with pytest.raises(PatternError):
            parse_glob(bad)
    rules, problems = compile_rules([InputRule("a.b", glob="**/x/{sample}.bam"), InputRule("a.c", glob="a/**"),
                                     InputRule("a.d", sniff="csv_columns:a")])
    assert len(rules) == 1 and rules[0].floating and rules[0].head.text == "x" and len(problems) == 1


def _analysis_tree(tmp_path):
    root = tmp_path / "root"
    run = make_run(root / "proj_one", RunSpec(reads=20, hours=1.0))
    trees.age_tree(root)
    exp = os.path.dirname(os.path.dirname(run.run_dir))
    a = os.path.join(exp, "analysis")
    for s in ("liver", "kidney"):
        touch(os.path.join(a, "alignments", f"{s}_genome_sorted.bam"))
        touch(os.path.join(a, "alignments", f"{s}_alt_genome_sorted.bam"))
        touch(os.path.join(a, "alignments", f"{s}_genome_sorted.bam.bai"))
        touch(os.path.join(a, "counts", "salmon", f"{s}_salmon", "aux_info", "meta_info.json"))
        touch(os.path.join(exp, "called_basecalled", "main", f"{s}.bam"))
        touch(os.path.join(exp, "called_basecalled", "alt", f"{s}.bam"))
        touch(os.path.join(exp, "called_basecalled", "deep", "er", "bam_pass", s, "part_0.bam"))
    touch(os.path.join(a, "alignments", "fresh_genome_sorted.bam"), age=5)          # still being written
    touch(os.path.join(exp, "sheet.csv"), "kit,flow_cell_id,barcode,alias\nK,TST00001,barcode01,liver\n")
    touch(os.path.join(exp, "notes.csv"), "a,b\n1,2\n")
    touch(os.path.join(root, "elsewhere", "deep", "analysis", "alignments", "lung_genome_sorted.bam"))
    touch(os.path.join(root, "fixed", "place", "table.tsv"))
    return root, exp, run


RULES = [
    InputRule("align.genome_bam", glob="**/analysis/alignments/{sample}_genome_sorted.bam", exclude=("*_alt_*",), source="config"),
    InputRule("align.genome_bam", glob="**/analysis/alignments/{sample}_alt_genome_sorted.bam", variant="alt", source="config"),
    InputRule("quant.salmon", glob="**/analysis/counts/salmon/{sample}_salmon/aux_info/meta_info.json", source="config"),
    InputRule("dorado.bam", glob="**/*_basecalled/{variant:main|alt}/{sample}.bam", source="config"),
    InputRule("dorado.bam", glob="**/*_basecalled/**/bam_pass/{sample}/*.bam", source="config"),
    InputRule("pychopper.stats", glob="fixed/place/*.tsv", source="config"),
]


def test_files_are_found_from_anchors_by_directed_descent(tmp_path):
    root, exp, run = _analysis_tree(tmp_path)
    w = walk_root("r", str(root), RULES)
    got = sorted((c.role, os.path.relpath(c.path, root), c.sample, c.variant) for c in w.candidates)
    e = os.path.relpath(exp, root)
    assert got == sorted([
        ("align.genome_bam", f"{e}/analysis/alignments/liver_genome_sorted.bam", "liver", None),
        ("align.genome_bam", f"{e}/analysis/alignments/kidney_genome_sorted.bam", "kidney", None),
        ("align.genome_bam", f"{e}/analysis/alignments/liver_alt_genome_sorted.bam", "liver", "alt"),
        ("align.genome_bam", f"{e}/analysis/alignments/kidney_alt_genome_sorted.bam", "kidney", "alt"),
        ("align.genome_bam", "elsewhere/deep/analysis/alignments/lung_genome_sorted.bam", "lung", None),
        ("quant.salmon", f"{e}/analysis/counts/salmon/liver_salmon/aux_info/meta_info.json", "liver", None),
        ("quant.salmon", f"{e}/analysis/counts/salmon/kidney_salmon/aux_info/meta_info.json", "kidney", None),
        ("dorado.bam", f"{e}/called_basecalled/main/liver.bam", "liver", "main"),
        ("dorado.bam", f"{e}/called_basecalled/main/kidney.bam", "kidney", "main"),
        ("dorado.bam", f"{e}/called_basecalled/alt/liver.bam", "liver", "alt"),
        ("dorado.bam", f"{e}/called_basecalled/alt/kidney.bam", "kidney", "alt"),
        ("dorado.bam", f"{e}/called_basecalled/deep/er/bam_pass/liver/part_0.bam", "liver", None),
        ("dorado.bam", f"{e}/called_basecalled/deep/er/bam_pass/kidney/part_0.bam", "kidney", None),
        ("pychopper.stats", "fixed/place/table.tsv", None, None),
    ])
    assert w.counts["too_new"] == 1, "the file changed five seconds ago is left for the next scan"
    assert w.counts["excluded"] == 2, "without the exclude, {sample} would have captured 'liver_alt'"
    assert len(w.anchors) == 3 and [f.path for f in w.runs] == [run.run_dir]
    assert all(c.size == 1 and c.mtime_ns > 0 and c.root_id == "r" for c in w.candidates)
    liver = next(c for c in w.candidates if c.path.endswith("/liver_genome_sorted.bam"))
    assert liver.rule == "flowcellar.toml rule align.genome_bam: **/analysis/alignments/{sample}_genome_sorted.bam"
    # the anchors were not walked: only the folders a pattern can match in were listed
    assert w.counts["listings"] < 40


def test_the_rules_the_package_ships_never_hide_what_lies_below_a_folder_called_outs(tmp_path):
    # the rules for analyses (single-cell and spatial add-on) look for file names; a folder called "outs" is walked
    # like any other, so a run and a read-data file below one are still found
    root = tmp_path / "data"
    run = make_run(root / "alpha_project" / "outs", RunSpec(seed=91))
    touch(str(root / "beta_project" / "outs" / "analysis" / "alignments" / "liver_genome_sorted.bam"))
    touch(str(root / "beta_project" / "outs" / "metrics_summary.csv"))
    cfg = loader.merge({"roots": [{"id": "r", "path": str(root)}], "discovery": {"min_age_seconds": 0},
                        "inputs": [{"role": "align.genome_bam", "glob": "**/analysis/alignments/{sample}_genome_sorted.bam"}]}).config
    found = discovery.discover(cfg)
    assert [f.path for _, f in found.runs] == [run.run_dir]
    assert [os.path.relpath(c.path, root) for c in found.candidates if c.role == "align.genome_bam"] == ["beta_project/outs/analysis/alignments/liver_genome_sorted.bam"]
    assert [os.path.relpath(c.path, root) for c in found.analyses] == ["beta_project/outs/metrics_summary.csv"], \
        "the anchor file of an analysis is reported apart, and is no read-data file"
    assert all(not (r.glob or "").startswith("**/outs/") for r in discovery.input_rules(cfg))


def test_tables_are_told_by_their_content_and_every_file_knows_its_experiment_folder(tmp_path, readdata):
    root, exp, run = _analysis_tree(tmp_path)
    cfg = loader.merge({"roots": [{"id": "r", "path": str(root)}],
                        "inputs": [{"role": r.role, "glob": r.glob, "exclude": list(r.exclude)} for r in RULES]}).config
    cache: dict = {}
    found = discovery.discover(cfg, sniff_cache=cache)
    sheets = [c for c in found.candidates if c.role == "basecaller.sample_sheet"]
    assert [c.path for c in sheets] == [os.path.join(exp, "sheet.csv")] and sheets[0].rule == "told by its content"
    assert readdata.calls["sniff"] == 2 and len(cache) == 2
    by_name = {os.path.basename(c.path): c for c in found.candidates}
    assert by_name["liver_genome_sorted.bam"].experiment_dir == exp
    assert by_name["liver.bam"].experiment_dir == exp
    assert by_name["lung_genome_sorted.bam"].experiment_dir == str(root / "elsewhere" / "deep")
    assert by_name["table.tsv"].experiment_dir is None
    discovery.discover(cfg, sniff_cache=cache)
    assert readdata.calls["sniff"] == 2, "a table that did not change is not opened again"


# ------------------------------------------------------------------------------------------------ roots that do not answer
def test_a_missing_root_and_an_emptied_root_are_offline(tree):
    cfg = config_for(tree)
    names = {r["id"]: sorted(os.listdir(r["path"])) for r in cfg["roots"]}
    trees.take_offline(tree.roots["seq"], "missing")
    trees.take_offline(tree.roots["store"], "empty")
    found = discovery.discover(cfg, recorded=names)
    assert found.roots["seq"].state == "offline" and "does not exist" in found.roots["seq"].reason
    assert found.roots["store"].state == "offline" and "0 of the 1 first-level folders" in found.roots["store"].reason
    assert found.roots["arch"].state == "reachable"
    assert [rid for rid, _ in found.runs] == ["arch"]
    # an empty root that was never scanned is simply empty: there is nothing to compare it with
    assert discovery.discover(cfg, recorded={}).roots["store"].state == "reachable"


def test_half_of_the_first_level_names_is_enough():
    assert guard.holds_enough(["a", "b"], ["a", "b", "c", "d"])
    assert not guard.holds_enough(["a"], ["a", "b", "c", "d"])
    assert guard.holds_enough([], []) and guard.holds_enough(["x"], None)


@pytest.mark.filterwarnings("ignore:.*multi-threaded.*fork.*:DeprecationWarning")     # this test asks for fork on purpose
@pytest.mark.parametrize("method", ["fork", "spawn"])
def test_a_root_that_does_not_answer_or_stalls_is_given_up_within_its_deadline(tree, method):
    roots = [{"id": "seq", "path": tree.roots["seq"]}]
    t0 = time.monotonic()
    out = guard.check_roots(roots, check_seconds=0.8, target=trees.silent_child, start_method=method)
    assert out["seq"].state == "offline" and "did not answer within 0.8 s" in out["seq"].reason and out["seq"].walk is None
    out = guard.check_roots(roots, check_seconds=10, stall_seconds=0.5, target=trees.stalling_child, start_method=method)
    assert out["seq"].state == "skipped" and "no progress for 0.5 s" in out["seq"].reason
    assert time.monotonic() - t0 < 20
    ok = guard.check_roots(roots, start_method=method)
    assert ok["seq"].state == "reachable" and len(ok["seq"].walk.runs) == 4 and ok["seq"].first_level
