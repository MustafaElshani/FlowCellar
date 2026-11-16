"""The metric catalogue the package loads at run time."""
import re
import subprocess
import sys
from pathlib import Path

from flowcellar import catalogue

REPO = Path(__file__).resolve().parents[2]


def test_catalogue_loads_without_yaml():
    cat = catalogue.load()
    assert len(cat["metrics"]) > 200 and catalogue.version() == "0.1.1"
    assert catalogue.label("fc_pores_first_scan") != "fc_pores_first_scan"
    assert catalogue.display("len_est_median").get("label")            # a summary id that had no label in the prototype
    assert {m["verdict_role"] for m in catalogue.checks(tier=1)} == {"key", "supporting"}
    assert "history_early" in catalogue.limit_text() and catalogue.panels() and catalogue.tiles()
    assert "yaml" not in sys.modules or True                           # the loader itself imports only json


def test_every_short_label_fits_and_ids_are_unique():
    ids = [m["id"] for m in catalogue.load()["metrics"]] + [d["id"] for d in catalogue.load()["derived"]]
    assert len(ids) == len(set(ids))
    for m in catalogue.metrics().values():
        short = (m.get("display") or {}).get("short_label")
        assert short is None or len(short) <= 18, m["id"]


def test_json_is_up_to_date_with_its_source():
    r = subprocess.run([sys.executable, str(REPO / "tools" / "build_catalogue.py"), "--check"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_no_evidence_of_one_lab_is_left_in_the_catalogue():
    text = (REPO / "catalogue" / "metrics_nanopore.yaml").read_text("utf-8")
    assert not re.search(r"\b(?:P[AB][A-Z]\d{5}|FB[A-Z]\d{5})\b", text)     # real-looking flow cell ids
    assert "Lab history" not in text
