"""The pass over the sequencing summary and the pass over the pore-scan table: completeness, byte-range workers,
time limit, cache, damaged files."""
import os
import time

import numpy as np
import pytest
from engine_support import problems

from flowcellar.contracts import PassOptions
from flowcellar.engine import api, passes
from flowcellar.model import partial as P
from flowcellar.sources import porescan, seqsummary, workers
from flowcellar.testing import minknow_variants as V
from flowcellar.testing.fixtures import RunSpec, make_run


@pytest.fixture(scope="module")
def big(tmp_path_factory):
    """A summary large enough to be cut into several byte ranges (30,000 reads, 7 MB)."""
    return make_run(tmp_path_factory.mktemp("big"), RunSpec(flow_cell_id="TST00060", reads=30000, hours=10.0, seed=60))


def _probe(path):
    return seqsummary.probe_task({"path": path}, lambda n: None)


def _pass(f, monkeypatch, workers_n, piece=1 << 20, **kw):
    monkeypatch.setattr(seqsummary, "MIN_BYTES_PER_WORKER", 64 << 10)
    path = f.files["minknow.sequencing_summary"]
    return seqsummary.run_pass(path, _probe(path), channels=512, ratio=0.6, jobs=workers_n, piece=piece, **kw)


def _equal(a: dict, b: dict) -> bool:
    for k in a:
        if k == "pore_type":
            if a[k] != b[k]:
                return False
        elif not np.array_equal(a[k], b[k]):
            return False
    return True


def test_a_complete_pass_reproduces_minknows_read_count_and_event_sum(big, monkeypatch):
    res = _pass(big, monkeypatch, 1)
    assert res["meta"]["state"] == "ok" and res["meta"]["complete"] is True
    assert int(res["acc"]["n"][0]) == big.expected["reads"] == 30000
    assert int(res["acc"]["ev_total"][0]) == big.expected["minknow_events"]          # equals the throughput file's last row
    assert res["meta"]["bytes_read"] == os.path.getsize(big.files["minknow.sequencing_summary"]) - _probe(big.files["minknow.sequencing_summary"])["body_start"]


@pytest.mark.parametrize("n", [2, 4, 8])
def test_results_do_not_depend_on_the_number_of_workers(big, monkeypatch, n):
    one = _pass(big, monkeypatch, 1)
    many = _pass(big, monkeypatch, n)
    assert many["meta"]["workers"] == n and many["meta"]["state"] == "ok"
    assert _equal(one["acc"], many["acc"])
    f1 = seqsummary.fragments(one["acc"], big.run_key, channels=512, ratio=0.6)
    fn = seqsummary.fragments(many["acc"], big.run_key, channels=512, ratio=0.6)
    assert f1 == fn


def test_results_do_not_depend_on_the_piece_size(big, monkeypatch):
    a = _pass(big, monkeypatch, 3, piece=1 << 20)
    b = _pass(big, monkeypatch, 3, piece=37_001)                    # pieces that end in the middle of lines
    assert _equal(a["acc"], b["acc"])


def test_every_line_belongs_to_exactly_one_range(big):
    path = big.files["minknow.sequencing_summary"]
    pr = _probe(path)
    total = 0
    for lo, hi in workers.split_ranges(pr["body_start"], pr["size"], 7):
        pieces = workers.Pieces(path, lo, hi, pr["body_start"], piece=50_000)
        total += sum(data.count(b"\n") for data in pieces)
        assert not pieces.cut_tail and not pieces.stopped
    assert total == big.expected["reads"]


def test_a_time_limited_pass_says_partial_everywhere(big, monkeypatch, tmp_path):
    monkeypatch.setattr(seqsummary, "MIN_BYTES_PER_WORKER", 64 << 10)
    monkeypatch.setattr(seqsummary, "PIECE", 200_000)
    run = api.identify(big.run_dir)
    refs = api.full_pass(run, str(tmp_path), PassOptions(jobs=2, deadline=time.monotonic() - 1.0))
    assert refs == []                                                # no time at all: nothing is read, nothing is claimed
    path = big.files["minknow.sequencing_summary"]
    res = seqsummary.run_pass(path, _probe(path), channels=512, ratio=0.6, jobs=4, piece=100_000, deadline=time.monotonic() + 0.0)
    assert res["meta"]["state"] == "partially_read" and res["meta"]["complete"] is False
    part = seqsummary.make_partial(path, _probe(path), _limited(big, monkeypatch), run.run_key, channels=512, ratio=0.6)
    assert part["state"] == "partially_read" and part["coverage"]["kind"] == "byte_ranges" and 0 < part["coverage"]["fraction"] < 1
    assert not P.is_complete(part)
    written = P.write_partial(str(tmp_path), part)
    ref = passes._ref(part, written, cached=False)
    b = api.assemble(run, [ref])
    dyn, fc = b["modules"]["dyn"], b["modules"]["fc"]
    assert dyn["state"] == "partial" and dyn["partial_reason"] == "partly_read" and dyn["coverage"]["kind"] == "byte_ranges"
    assert b["modules"]["len.est"]["state"] == "partial" and fc["state"] == "partial"
    assert "fc_channel_map" not in fc["cells"][0].get("channels", {}) and "fc_active_channels" not in fc["summary"]
    assert b["tiers"]["run_health"] == {"state": "partial", "detail": "quick", "missing_roles": []}
    inp = next(i for i in b["provenance"]["inputs"] if i["role"] == "minknow.sequencing_summary")
    assert inp["state"] == "partially_read" and inp["coverage"]["records_total"] == big.expected["reads"]
    sf = api.judge(b, [], in_library=False)
    st = sf["status"]
    assert "reads_partly_counted" in [m["code"] for m in st["modifiers"]]
    rules = {r["id"]: r for r in st["rules"]}
    assert rules["fc_active_channels"]["outcome"] == "not_evaluated"
    assert "From " in rules["dyn_signal_positive_frac"]["message"] and st["provisional"] is False
    from flowcellar.model import bundle as B
    assert problems(B.compose(b, status=st, relations=sf["relations"])) == []


def _limited(big, monkeypatch):
    """A pass that stops after its first piece in every range."""
    path = big.files["minknow.sequencing_summary"]
    pr = _probe(path)
    acc = seqsummary.new_acc(512)
    done = 0
    body = pr["size"] - pr["body_start"]
    for lo, hi in workers.split_ranges(pr["body_start"], pr["size"], 4):
        pieces = workers.Pieces(path, lo, hi, pr["body_start"], piece=150_000)
        for data in pieces:
            cols, _bad = seqsummary.parse_block(data, pr["header"])
            seqsummary.update(acc, cols, 512, 0.6)
            break
        done += pieces.done
    return {"acc": acc, "meta": {"state": "partially_read", "reads": int(acc["n"][0]), "coverage_frac": done / body, "header": pr["header"],
                                 "bad_rows": 0, "cut_tail": False, "bytes_read": done, "errors": [], "seconds": 0.1, "workers": 4}}


def test_second_full_pass_reads_nothing(fake_run, tmp_path, monkeypatch):
    run = api.identify(fake_run.run_dir)
    first = api.full_pass(run, str(tmp_path))
    assert [(r.role, r.state, r.complete, r.cached) for r in first] == [("minknow.sequencing_summary", "ok", True, False),
                                                                        ("minknow.pore_scan", "ok", True, False)]
    assert all(os.path.exists(r.path) and r.path.endswith(".fcp.json.gz") for r in first)

    def no_reading(*a, **k):
        raise AssertionError("an input was opened on the second call")
    monkeypatch.setattr(workers, "run_tasks", no_reading)
    monkeypatch.setattr(workers, "run_one", no_reading)
    second = api.full_pass(run, str(tmp_path))
    assert [(r.partial_id, r.cached, r.complete) for r in second] == [(r.partial_id, True, True) for r in first]


def test_a_changed_input_is_read_again(tmp_path):
    f = make_run(tmp_path / "run", RunSpec(flow_cell_id="TST00061", seed=61, reads=1500))
    run = api.identify(f.run_dir)
    first = api.full_pass(run, str(tmp_path / "p"))
    with open(f.files["minknow.sequencing_summary"], "a") as fh:
        fh.write("")                                                 # same bytes: only the time changes
    os.utime(f.files["minknow.sequencing_summary"], (time.time() + 5, time.time() + 5))
    again = api.full_pass(api.identify(f.run_dir), str(tmp_path / "p"))
    assert again[0].partial_id == first[0].partial_id and again[0].cached is True      # same fingerprint: the partial stands
    lines = open(f.files["minknow.sequencing_summary"]).read().split("\n")
    with open(f.files["minknow.sequencing_summary"], "w") as fh:
        fh.write("\n".join(lines[:-200]) + "\n")
    third = api.full_pass(api.identify(f.run_dir), str(tmp_path / "p"))
    assert third[0].partial_id != first[0].partial_id and third[0].cached is False


def test_partial_file_describes_itself(fake_run, full8):
    p = P.read_partial(full8.refs[0].path)
    assert p["format"] == "flowcellar.partial" and p["role"] == "minknow.sequencing_summary" and p["state"] == "ok"
    assert p["modules"] == seqsummary.MODULES and p["config_hash"] == ""
    assert p["partial_id"] == P.partial_id(p["input"]["fingerprint"], p["role"], p["modules"])
    assert p["routing"] == {"runs": {fake_run.run_key: fake_run.expected["reads"]}, "unrouted": 0, "default_run": fake_run.run_key}
    assert {f["module"] for f in p["fragments"]} == {"dyn", "len.est", "fc"} and all(f["run_key"] == fake_run.run_key for f in p["fragments"])
    dyn = next(f for f in p["fragments"] if f["module"] == "dyn")
    assert "minknow_events" in dyn["cells"][0]["scalars"] and "estimated_bases" not in dyn["cells"][0]["scalars"]      # raw units
    assert p["facts"]["reads"] == fake_run.expected["reads"] and p["facts"]["minknow_events"] == fake_run.expected["minknow_events"]
    q = P.read_partial(full8.refs[1].path)
    assert q["role"] == "minknow.pore_scan" and q["modules"] == porescan.MODULES and q["coverage"]["records"] == fake_run.expected["pore_scan_rows"]


def test_truncated_summary_is_a_partial_section_never_a_crash(tmp_path, build):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00062", seed=62))
    left = V.truncate_summary(f, 0.6)
    built = build(f)
    ref = built.refs[0]
    assert (ref.state, ref.complete) == ("truncated", False)
    b = built.bundle
    dyn = b["modules"]["dyn"]
    assert dyn["state"] == "partial" and dyn["coverage"] == {"kind": "truncated", "records": left, "records_total": f.expected["reads"],
                                                             "fraction": round(left / f.expected["reads"], 6)}
    assert sum(c["scalars"]["reads"] for c in dyn["cells"]) == left
    assert b["tiers"]["run_health"]["detail"] == "quick"
    assert next(i for i in b["provenance"]["inputs"] if i["role"] == "minknow.sequencing_summary")["state"] == "truncated"
    assert "reads_partly_counted" in [m["code"] for m in built.status["modifiers"]]
    assert problems(built.composed) == []
    # the same file again: the truncated partial is final, nothing is read twice
    assert api.full_pass(built.run, os.path.dirname(os.path.dirname(ref.path)))[0].cached is True
    # as lab history the run still counts for what the small files say (output, pores, occupancy), exactly as it
    # did in quick detail; the numbers of the damaged read table are not handed on
    facts = api.card(b, None)["facts"]
    assert b["tiers"]["run_health"]["state"] == "partial" and facts["t1_complete"] is True
    assert facts["speed_med"] is None and facts["reads"] == f.expected["reads"]
    quick = api.card(api.quick_bundle(built.run), None)["facts"]
    assert quick["t1_complete"] is True and {k: facts[k] for k in ("reads", "mux0_pores", "occ_30_90")} == \
        {k: quick[k] for k in ("reads", "mux0_pores", "occ_30_90")}


def test_summary_that_ends_at_a_line_break_but_lacks_reads_fails_the_completeness_test(tmp_path, build):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00063", seed=63))
    lines = open(f.files["minknow.sequencing_summary"]).read().split("\n")
    with open(f.files["minknow.sequencing_summary"], "w") as fh:
        fh.write("\n".join(lines[:1001]) + "\n")                     # header and 1,000 whole lines
    built = build(f)
    assert built.refs[0].state == "ok"                              # the file itself is whole
    dyn = built.bundle["modules"]["dyn"]
    assert dyn["state"] == "partial" and dyn["params"]["completeness"]["equal"] is False
    assert dyn["coverage"]["kind"] == "truncated" and dyn["coverage"]["records_total"] == f.expected["reads"]
    assert "summary_incomplete" in [n["code"] for n in dyn["notes"]]
    assert built.bundle["tiers"]["run_health"]["detail"] == "quick" and problems(built.composed) == []


def test_bad_rows_are_skipped_and_counted(tmp_path, build):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00064", seed=64))
    path = f.files["minknow.sequencing_summary"]
    lines = open(path).read().split("\n")
    lines.insert(500, "this line\thas too few columns")
    lines.insert(900, "\t".join(["-"] * 7 + ["not_a_number", "1", "x", "y", "z", "not_set", "e", "s", "signal_positive"]))
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    built = build(f)
    p = P.read_partial(built.refs[0].path)
    assert p["state"] == "ok" and p["facts"]["bad_rows"] == 2 and p["facts"]["reads"] == f.expected["reads"]
    inp = next(i for i in built.bundle["provenance"]["inputs"] if i["role"] == "minknow.sequencing_summary")
    assert inp["records"] == f.expected["reads"] and inp["records_skipped"] == 2
    assert built.bundle["modules"]["dyn"]["state"] == "ok"


def test_a_file_that_is_not_a_summary_is_an_error_section(tmp_path, build):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00065", seed=65))
    with open(f.files["minknow.sequencing_summary"], "w") as fh:
        fh.write("read_id\tsomething\n" + "abc\t1\n" * 100)
    built = build(f)
    assert built.refs[0].state == "error"
    b = built.bundle
    assert b["modules"]["dyn"]["state"] == "error" and b["modules"]["len.est"]["state"] == "error"
    assert "columns" in b["modules"]["dyn"]["notes"][0]["message"]
    assert next(i for i in b["provenance"]["inputs"] if i["role"] == "minknow.sequencing_summary")["state"] == "error"
    assert built.rules["dyn_signal_positive_frac"]["outcome"] == "not_evaluated" and built.status["pending"] == []
    assert problems(built.composed) == []
    garbage = make_run(tmp_path / "g", RunSpec(flow_cell_id="TST00066", seed=66))
    lines = open(garbage.files["minknow.sequencing_summary"]).read().split("\n")
    with open(garbage.files["minknow.sequencing_summary"], "w") as fh:
        fh.write(lines[0] + "\n" + "\n".join("garbage\x00line %d" % i for i in range(300)) + "\n")
    g = build(garbage)
    assert g.refs[0].state == "error" and g.bundle["modules"]["dyn"]["state"] == "error"


def test_summary_without_optional_columns(tmp_path, build):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00067", seed=67))
    path = f.files["minknow.sequencing_summary"]
    lines = [ln.split("\t") for ln in open(path).read().split("\n") if ln]
    keep = [i for i, name in enumerate(lines[0]) if name not in ("pore_type", "mux", "filename_fastq", "experiment_id")]
    with open(path, "w") as fh:
        fh.write("\n".join("\t".join(ln[i] for i in keep) for ln in lines) + "\n")
    built = build(f)
    assert built.refs[0].complete
    assert "run_level" not in built.bundle["modules"]["dyn"]
    assert built.bundle["modules"]["fc"]["cells"][0]["scalars"]["reads_mux_unknown"] == f.expected["reads"]
    assert problems(built.composed) == []


def test_missing_summary_means_not_assessed_not_pending(tmp_path, build):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00068", seed=68))
    os.unlink(f.files["minknow.sequencing_summary"])
    built = build(f)
    assert [r.role for r in built.refs] == ["minknow.pore_scan"]
    assert built.status["pending"] == [] and built.status["provisional"] is False
    assert built.rules["dyn_speed_est_median"]["outcome"] == "not_evaluated"
    assert "sequencing summary is not in the run folder" in built.status["detail"]
    assert "fc_well_map" in built.bundle["modules"]["fc"]["cells"][0]["channels"]                # the scan table was still read
    assert built.bundle["tiers"]["run_health"]["missing_roles"] == ["minknow.sequencing_summary"]


def test_truncated_pore_scan_table(tmp_path, build):
    f = make_run(tmp_path, RunSpec(flow_cell_id="TST00069", seed=69))
    V.truncate_summary(f, 0.7, role="minknow.pore_scan")
    built = build(f)
    assert built.refs[1].state == "truncated" and not built.refs[1].complete
    assert built.bundle["modules"]["fc"]["state"] in ("ok", "partial")
    assert problems(built.composed) == []


def test_a_worker_that_fails_or_stalls_is_reported_not_raised(tmp_path, monkeypatch):
    r = workers.run_tasks("flowcellar.sources.seqsummary:probe_task", [{"path": str(tmp_path / "no_such_file")}])
    assert r[0].state == "error" and "FileNotFoundError" in r[0].error
    r = workers.run_tasks("time:sleep_is_not_a_task", [{}])
    assert r[0].state == "error"
    t0 = time.monotonic()
    r = workers.run_tasks("flowcellar.sources.workers:_sleep_task", [{"seconds": 30}], stall_timeout_s=1.0)
    assert r[0].state == "stalled" and time.monotonic() - t0 < 10


def test_gentle_reading_uses_two_workers():
    assert seqsummary.n_workers(0, 40 << 30, True) == 2
    assert seqsummary.n_workers(0, 40 << 30, False) == min(8, os.cpu_count() or 1)
    assert seqsummary.n_workers(0, 500 << 20, False) == 1                     # one below 1 GB
    assert seqsummary.n_workers(8, 20 << 20, False) == 3                      # never more than one per 8 MB
    assert seqsummary.n_workers(4, 40 << 30, True) == 4                       # an explicit number wins


# ------------------------------------------------------------------------------------ exact N50 and median
def _events_of(run) -> np.ndarray:
    """Every read's event count, read straight from the fixture's summary (the slow, certain way)."""
    path = run.files["minknow.sequencing_summary"]
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        k = header.index("minknow_events")
        return np.array([int(line.split("\t")[k]) for line in fh if line.strip()], np.int64)


def _exact_n50(lengths: np.ndarray) -> float:
    s = np.sort(lengths)[::-1]
    cum = np.cumsum(s)
    return float(s[int(np.argmax(2 * cum >= cum[-1]))])


def test_estimated_n50_and_median_are_exact_not_read_off_the_histogram(big, tmp_path):
    """The catalogue promises an exact value from a counting array. A histogram with 40 bins per decade can sit
    several per cent away, which is enough to put a check on the other side of its limit."""
    ident = api.identify(big.run_dir)
    refs = api.full_pass(ident, str(tmp_path / "partials"), PassOptions(jobs=3))
    bundle = api.assemble(ident, refs)
    sec = bundle["modules"]["len.est"]
    ratio = sec["params"]["events_to_base_ratio"]
    ev = _events_of(big)
    assert sec["params"]["quantile_source"] == "exact"
    assert sec["summary"]["len_est_n50"] == round(_exact_n50(ev) * ratio)
    assert sec["summary"]["len_est_median"] == round(float(np.median(ev)) * ratio)
    assert not [n for n in sec.get("notes", []) if n["code"] == "quantiles_from_histogram"]
    ss = next(r for r in refs if r.role == seqsummary.ROLE)
    assert P.read_partial(ss.path)["facts"]["exact"]["reads"] == len(ev)


def test_a_result_cached_without_the_counting_array_says_its_n50_is_from_the_histogram(big, tmp_path):
    ident = api.identify(big.run_dir)
    pdir = str(tmp_path / "partials")
    refs = api.full_pass(ident, pdir, PassOptions(jobs=1))
    ss = next(r for r in refs if r.role == seqsummary.ROLE)
    old = P.read_partial(ss.path)
    old["facts"].pop("exact")                                # what an earlier build cached
    os.unlink(ss.path)
    P.write_partial(pdir, old)
    refs = api.full_pass(ident, pdir, PassOptions(jobs=1))
    assert all(r.cached for r in refs), "the cached result is still good: nothing is read again"
    sec = api.assemble(ident, refs)["modules"]["len.est"]
    assert sec["params"]["quantile_source"] == "histogram"
    assert [n["code"] for n in sec["notes"]] == ["quantiles_from_histogram"]
    exact = _exact_n50(_events_of(big)) * sec["params"]["events_to_base_ratio"]
    assert abs(sec["summary"]["len_est_n50"] - exact) / exact < 0.06
    assert problems(bundle := api.assemble(ident, refs)) == [] and bundle["modules"]["len.est"]["state"] == "ok"


def test_exact_lengths_on_small_arrays():
    counts = np.bincount(np.array([10, 10, 20, 30, 100]))
    out = seqsummary.exact_lengths(counts)
    assert out == {"reads": 5, "events_n50": 100, "events_median": 20.0, "over_cap": 0}
    assert seqsummary.exact_lengths(np.bincount(np.array([10, 20, 30, 40])))["events_median"] == 25.0
    assert seqsummary.exact_lengths(np.bincount(np.array([10, 20, 30, 40])))["events_n50"] == 30       # 40 + 30 >= 50
    assert seqsummary.exact_lengths(np.zeros(4, np.int64)) is None
