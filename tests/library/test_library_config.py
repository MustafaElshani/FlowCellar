"""flowcellar.toml: defaults, what is honoured, what is only warned about, and the two places it is written."""
from __future__ import annotations

import tomllib

import pytest

from flowcellar.config import loader


def test_defaults_are_generic_and_complete():
    cfg = loader.defaults()
    assert cfg["discovery"] == {"min_age_seconds": 120, "stale_after_hours": 6, "prune": [], "subtree_listing_budget": 20000,
                                "root_check_seconds": 20, "stall_timeout_seconds": 120}
    assert cfg["scan"] == {"jobs": 0, "threads": 4, "quick_first": True, "time_box_minutes": 240, "live_workers": 2}
    assert cfg["serve"] == {"host": "127.0.0.1", "port": 8737} and cfg["hub"] == {"sidecars": True, "reports": False}
    assert cfg["baselines"] == {"window": 40, "min_runs": 5, "trusted_runs": 8, "spread_guard": 0.45, "groups": []}
    assert cfg["flow_cell"] == {"restart_max_hours": 2, "restart_max_output_frac": 0.02}
    assert cfg["roots"] == [] and cfg["inputs"] == [] and cfg["link"] == {"rules": []} and cfg["profiles"] == {}
    assert cfg["samples"] == {"unclassified": ["unclassified"], "alias": {}}


def test_a_file_overrides_key_by_key_and_roots_get_ids(tmp_path):
    p = tmp_path / "flowcellar.toml"
    p.write_text('''
[library]
name = "demo"
[[roots]]
path = "/data/site_a/projects"
[[roots]]
path = "/data/site_b/projects"
project_depth = 2
live = true
[[roots]]
id = "seq"
path = "/data/minknow/"
project_pattern = "*/{project}/**"
[discovery]
prune = ["notebooks", "figures*"]
[scan]
time_box_minutes = 30
[[inputs]]
role = "align.genome_bam"
glob = "**/analysis/{sample}.bam"
exclude = ["*_tmp_*"]
[samples.alias]
"A-1" = "A1"
[profiles.cdna.rules.fc_occupancy_first_hour]
warn = 0.6
[[baselines.groups]]
match = { kit = "SQK-TEST" }
label = "test protocol"
''')
    loaded = loader.load(p)
    cfg = loaded.config
    assert loaded.warnings == [] and loaded.notices == [] and loaded.path == str(p)
    assert [r["id"] for r in cfg["roots"]] == ["site_a_projects", "site_b_projects", "seq"]
    assert cfg["roots"][1] == {"id": "site_b_projects", "path": "/data/site_b/projects", "live": True, "project_depth": 2,
                               "project_pattern": None}
    assert cfg["roots"][2]["path"] == "/data/minknow" and cfg["roots"][2]["project_pattern"] == "*/{project}/**"
    assert cfg["discovery"]["prune"] == ["notebooks", "figures*"] and cfg["discovery"]["min_age_seconds"] == 120
    assert cfg["scan"]["time_box_minutes"] == 30 and cfg["scan"]["threads"] == 4
    assert cfg["inputs"] == [{"role": "align.genome_bam", "glob": "**/analysis/{sample}.bam", "exclude": ["*_tmp_*"]}]
    assert cfg["samples"]["alias"] == {"A-1": "A1"} and cfg["samples"]["unclassified"] == ["unclassified"]
    assert cfg["profiles"]["cdna"]["rules"]["fc_occupancy_first_hour"]["warn"] == 0.6
    assert cfg["baselines"]["groups"][0]["label"] == "test protocol" and cfg["baselines"]["window"] == 40


def test_unknown_keys_warn_and_never_stop(tmp_path):
    p = tmp_path / "flowcellar.toml"
    p.write_text('[scan]\njob = 3\n[discovry]\nprune = []\n[[roots]]\npath = "/x"\nlive_mode = true\n')
    loaded = loader.load(p)
    text = "\n".join(loaded.warnings)
    assert "unknown key 'job' in [scan] (did you mean jobs?)" in text
    assert "unknown section [discovry] (did you mean discovery?)" in text
    assert "unknown key 'live_mode' in [[roots]] entry 1" in text
    assert loaded.config["scan"]["jobs"] == 0


def test_keys_of_features_that_are_not_built_are_read_and_reported_once(tmp_path):
    p = tmp_path / "flowcellar.toml"
    p.write_text('[hub]\ndefault_group = "month"\n[[hub.views]]\nname = "x"\n[[projects]]\nname = "p"\n'
                 '[[roots]]\npath = "/x"\nminknow_tmp = "/tmp/x"\nmax_read_mb_s = 50\n')
    loaded = loader.load(p)
    assert loaded.warnings == []
    assert sorted(n["key"] for n in loaded.notices) == ["hub.default_group", "hub.views", "projects", "roots.max_read_mb_s",
                                                         "roots.minknow_tmp"]
    assert {n["code"] for n in loaded.notices} == {"config_not_used"}
    assert "default_group" not in loaded.config["hub"]


@pytest.mark.parametrize("text, says", [
    ('[scan]\ntime_box_minutes = "long"\n', "[scan] time_box_minutes must be a whole number or a number"),
    ('[[roots]]\nid = "a"\n', "needs a path"),
    ('[[roots]]\npath = "/a"\nid = "x"\n[[roots]]\npath = "/b"\nid = "x"\n', "two [[roots]] entries have the id 'x'"),
    ('[[roots]]\npath = "/a"\nid = "not ok"\n', "may hold letters, digits"),
    ('[[roots]]\npath = "/a"\nproject_pattern = "*/x"\n', "must hold the capture {project}"),
    ('[[inputs]]\nrole = "x.y"\n', "needs a role and a glob or a sniff rule"),
    ('[scan\njobs = 1\n', "not valid TOML"),
    ('[discovery]\nmin_age_seconds = -1\n', "must not be negative"),
])
def test_a_value_that_cannot_be_used_is_an_error_that_names_it(tmp_path, text, says):
    p = tmp_path / "flowcellar.toml"
    p.write_text(text)
    with pytest.raises(loader.ConfigError) as e:
        loader.load(p)
    assert says in e.value.message and str(p) in e.value.message and e.value.exit_code == 6


def test_the_file_scan_writes_is_valid_and_a_root_is_appended_without_touching_the_rest(tmp_path):
    p = tmp_path / "flowcellar.toml"
    roots = [{"id": "minknow", "path": '/data/with "quotes"/minknow', "live": True, "project_depth": 1, "project_pattern": None}]
    loader.write_new(p, "My lab's \"library\"", roots)
    raw = tomllib.loads(p.read_text())
    assert raw["library"]["name"] == "My lab's \"library\"" and raw["roots"][0]["path"] == '/data/with "quotes"/minknow'
    assert raw["roots"][0]["live"] is True
    edited = p.read_text() + '\n[discovery]\nprune = ["notebooks"]   # my own line\n'
    p.write_text(edited)
    loader.append_roots(p, [{"id": "projects", "path": "/data/projects", "live": False, "project_depth": 2, "project_pattern": None}],
                        "2026-05-06T07:08:09Z")
    text = p.read_text()
    assert text.startswith(edited) and "# added by flowcellar scan on 2026-05-06" in text
    cfg = loader.load(p).config
    assert [r["id"] for r in cfg["roots"]] == ["minknow", "projects"] and cfg["roots"][1]["project_depth"] == 2
    assert cfg["discovery"]["prune"] == ["notebooks"]
    p.write_text('roots = []\n')
    with pytest.raises(loader.ConfigError):
        loader.append_roots(p, roots, "2026-05-06T07:08:09Z")
    assert p.read_text() == 'roots = []\n', "a file that cannot take the entry is left as it was"


def test_root_ids():
    assert loader.root_ids_for(["/a/minknow", "/b/projects", "/c/projects"], set()) == ["minknow", "b_projects", "c_projects"]
    assert loader.root_ids_for(["/x/projects"], {"projects"}) == ["x_projects"]
    assert loader.root_ids_for(["/x/projects"], {"projects", "x_projects"}) == ["x_projects_2"]
    assert loader.root_ids_for(["/"], set()) == ["root"]
    assert loader.root_ids_for(["/data/My Runs (2026)"], set()) == ["my_runs_2026"]


def test_the_shipped_defaults_name_no_laboratory():
    from importlib import resources
    text = resources.files("flowcellar.config").joinpath("defaults.toml").read_text("utf-8")
    cfg = tomllib.loads(text)
    assert "roots" not in cfg and "inputs" not in cfg and cfg["discovery"]["prune"] == []


def test_settings_that_have_no_effect_in_this_build_are_said_once(tmp_path):
    p = tmp_path / "flowcellar.toml"
    p.write_text('[scan]\nquick_first = false\n[project]\nfrom = "sample_sheet"\n')
    loaded = loader.load(p)
    assert sorted(n["key"] for n in loaded.notices) == ["project.from", "scan.quick_first"]
