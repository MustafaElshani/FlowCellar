"""The search grammar of ``ls --where``: the hub's grammar in Python, held equal to the viewer by shared vectors."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowcellar.library import query
from flowcellar.util.errors import UsageError

VECTORS = sorted((Path(__file__).resolve().parents[1] / "vectors" / "query").glob("*.json"))


def test_there_are_vectors():
    assert len(VECTORS) >= 3


@pytest.mark.parametrize("path", VECTORS, ids=lambda p: p.stem)
def test_shared_vectors(path):
    vector = json.loads(path.read_text("utf-8"))
    assert set(vector) >= {"name", "function", "args", "expect"} and vector["function"] == "query.filter"
    got = query.run_vector(vector)
    assert len(got) == len(vector["expect"]) == len(vector["args"]["queries"])
    for q, g, e in zip(vector["args"]["queries"], got, vector["expect"]):
        assert g == e, q


def test_the_vectors_hold_invented_names_only():
    for path in VECTORS:
        for row in json.loads(path.read_text("utf-8"))["args"]["rows"]:
            assert row["fc"].startswith("TST") and row["id"].startswith("TST")


def test_tokens():
    toks = query.tokenize('status:warn -fc:TST00001 "two words" name:"line a" -plain KIT:x')
    assert [(t["neg"], t["field"], t["value"]) for t in toks] == [
        (False, "status", "warn"), (True, "fc", "TST00001"), (False, None, "two words"), (False, "name", "line a"),
        (True, None, "plain"), (False, "kit", "x")]
    assert query.tokenize("") == [] and query.tokenize(None) == []
    assert query.fold("Tumour Ä") == "tumour a"


def test_a_command_refuses_a_search_it_cannot_read():
    rows = [{"id": "TST00001_0a1b2c3d", "fc": "TST00001", "verdict": "pass", "started": "2026-03-02T09:00:00Z", "reads": 10}]
    assert query.filter_rows(rows, None) == rows and query.filter_rows(rows, "   ") == rows
    assert query.filter_rows(rows, "status:pass reads:10") == rows
    with pytest.raises(UsageError) as e:
        query.filter_rows(rows, "statsu:pass")
    assert e.value.exit_code == 2 and "Did you mean status?" in e.value.message
    with pytest.raises(UsageError):
        query.filter_rows(rows, "reads:lots")
    assert query.filter_rows(rows, "statsu:pass", strict=False) == rows, "the hub leaves the term out instead"
    assert query.filter_rows(rows, "last:7d", now="2026-03-05T00:00:00Z") == rows
    assert query.filter_rows(rows, "last:7d", now="2026-04-05T00:00:00Z") == []


def test_rows_of_an_index_carry_their_sample_names():
    from flowcellar.library.index import build_index
    card = {"run_key": "TST00001_0a1b2c3d", "columns": {"id": "TST00001_0a1b2c3d", "fc": "TST00001", "started": "2026-03-02T09:00:00Z"},
            "samples": [{"sample_uid": "u1", "name": "liver_1", "alias": "L1", "run": "TST00001_0a1b2c3d", "key": "s1"},
                        {"sample_uid": "u2", "name": "liver_2", "alias": "liver_2", "run": "TST00001_0a1b2c3d", "key": "s2"}],
            "facts": {}}
    rows = query.rows_with_samples(build_index([card], baselines={"groups": []}))
    assert rows[0]["samples"] == ["liver_1", "L1", "liver_2"]
    assert query.filter_rows(rows, "sample:l1") == rows and query.filter_rows(rows, "liver_2") == rows
    assert query.filter_rows(rows, "sample:kidney") == []
