"""``flowcellar --config FILE run RUN_DIR -o DIR``: a run without a library is judged with the limits of that file."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys

from flowcellar.testing.fixtures import RunSpec, make_run

LIMITS = """\
[profiles.cdna.rules.yield_reads_at_72h]
warn = 10000
fail = 5000

[profiles.cdna.rules.yield_estimated_bases_at_72h]
warn = 10000000
fail = 5000000

[profiles.cdna.rules.fc_occupancy_first_hour]
enabled = false
"""


def flowcellar(*args: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if not k.startswith("FLOWCELLAR_")}
    return subprocess.run([sys.executable, "-m", "flowcellar", *args], capture_output=True, text=True, timeout=300, env=env)


def _status(path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)["status"]


def test_a_run_without_a_library_is_judged_with_the_limits_of_the_file_given(tmp_path):
    fr = make_run(tmp_path / "data", RunSpec(reads=4000, hours=8, seed=5))            # a small run: far below the built-in output level
    plain = flowcellar("run", fr.run_dir, "-o", str(tmp_path / "plain"))
    assert plain.returncode == 0, plain.stdout + plain.stderr
    before = _status(tmp_path / "plain" / f"{fr.run_key}.fcb.json.gz")
    by_id = {r["id"]: r for r in before["rules"]}
    assert before["verdict"] == "fail" and by_id["yield_reads_at_6h"]["severity"] == "fail" and before["profile"]["source"] == "builtin"

    config = tmp_path / "limits.toml"
    config.write_text(LIMITS)
    own = flowcellar("--config", str(config), "run", fr.run_dir, "-o", str(tmp_path / "own"))
    assert own.returncode == 0, own.stdout + own.stderr
    after = _status(tmp_path / "own" / f"{fr.run_key}.fcb.json.gz")
    by_id = {r["id"]: r for r in after["rules"]}
    assert by_id["yield_reads_at_6h"]["severity"] == "none" and after["verdict"] == "pass", [r for r in after["rules"] if r["severity"] != "none"]
    assert after["profile"]["source"] == "library_config"
    assert "fc_occupancy_first_hour" not in by_id or by_id["fc_occupancy_first_hour"]["outcome"] != "triggered"
    assert any(m["code"] == "custom_profile" for m in after["modifiers"])               # the page says that other limits were used

    missing = flowcellar("--config", str(tmp_path / "nowhere.toml"), "run", fr.run_dir, "-o", str(tmp_path / "x"))
    assert missing.returncode == 4 and "no such file" in missing.stderr
