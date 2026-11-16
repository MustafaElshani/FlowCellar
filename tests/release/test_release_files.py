"""The path from this repository to a release: metadata, licence files, the two workflows, the documentation.

These tests read the files a release is made from and fail when one of them stops saying what the release
needs. They build nothing and contact nothing. In a source distribution the workflow files are absent; the
tests about them are skipped there.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import tomllib
from importlib import metadata
from pathlib import Path

import pytest

import flowcellar

REPO = Path(__file__).resolve().parents[2]
WORKFLOWS = REPO / ".github" / "workflows"
PYPROJECT = tomllib.loads((REPO / "pyproject.toml").read_text("utf-8"))
PROJECT = PYPROJECT["project"]
# What PyPI's trusted publisher is registered for. All four must stay as they are.
PYPI_PROJECT, REPOSITORY, WORKFLOW_FILE, ENVIRONMENT = "flowcellar", "MustafaElshani/FlowCellar", "build.yml", "pypi"
# Current major versions of the actions the workflows use (checked against their release pages, November 2026).
ACTIONS = {"actions/checkout": "v7", "actions/setup-node": "v7", "actions/setup-python": "v7",
           "actions/upload-artifact": "v7", "actions/download-artifact": "v8", "pypa/gh-action-pypi-publish": "release/v1"}


def _tool(name: str):
    spec = importlib.util.spec_from_file_location(f"flowcellar_tool_{name}", REPO / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _workflow(name: str) -> dict:
    yaml = pytest.importorskip("yaml")
    path = WORKFLOWS / name
    if not path.exists():
        pytest.skip("no workflow files here (a source distribution does not carry them)")
    wf = yaml.safe_load(path.read_text("utf-8"))
    if True in wf:                       # YAML reads the key "on" as the boolean true
        wf["on"] = wf.pop(True)
    return wf


# ------------------------------------------------------------------------------------------- metadata
def test_metadata_is_fixed():
    assert PROJECT["name"] == PYPI_PROJECT
    assert PROJECT["requires-python"] == ">=3.11"
    assert PROJECT["dynamic"] == ["version", "readme"] and PYPROJECT["tool"]["hatch"]["version"]["path"] == "src/flowcellar/_version.py"
    assert PROJECT["scripts"] == {"flowcellar": "flowcellar.cli:main"}
    names = sorted(re.split(r"[<>=!~ ;\[]", d, maxsplit=1)[0].lower() for d in PROJECT["dependencies"])
    assert names == ["click", "numpy", "pyarrow"], "run-time dependencies are fixed: numpy, pyarrow, click"
    extras = PROJECT["optional-dependencies"]
    assert [d.split(">=")[0] for d in extras["bam"]] == ["pysam"]
    assert [d.split(">=")[0] for d in extras["validate"]] == ["jsonschema"]
    assert [d.split(">=")[0] for d in extras["sc"]] == ["h5py"]         # the single-cell and spatial add-on: HDF5 files only
    assert extras["all"] == ["flowcellar[bam,validate,sc]"]
    # the long description is README.md with its paths made into addresses, written by the metadata hook
    assert "readme" not in PROJECT and PYPROJECT["tool"]["hatch"]["metadata"]["hooks"]["custom"]["path"] == "hatch_build.py"
    assert (REPO / "README.md").is_file() and len(PROJECT["description"]) < 200


def test_the_licence_is_stated_once_and_its_files_exist():
    assert PROJECT["license"] == "GPL-3.0-or-later"
    assert not [c for c in PROJECT["classifiers"] if c.startswith("License ::")], "an SPDX expression replaces the licence classifier"
    for name in PROJECT["license-files"]:
        assert (REPO / name).is_file(), name
    assert PROJECT["license-files"] == ["LICENSE", "THIRD_PARTY_NOTICES.md"]
    text = (REPO / "LICENSE").read_text("utf-8")
    assert text.lstrip().startswith("GNU GENERAL PUBLIC LICENSE") and "Version 3, 29 June 2007" in text[:200]
    assert "END OF TERMS AND CONDITIONS" in text, "the whole text of the licence, as the Free Software Foundation publishes it"
    readme = (REPO / "README.md").read_text("utf-8")
    said = readme.split("## Licence and contact", 1)[1]
    assert re.search(r"GNU General Public License.{0,60}version 3.{0,40}any later version", said, re.S) and "provisional" not in readme
    assert "without any warranty" in said and "(LICENSE)" in said
    banner = (REPO / "viewer" / "build.mjs").read_text("utf-8")
    assert "GPL-3.0-or-later" in banner and "MIT licence. Includes" not in banner, "the viewer inside every page names the licence too"


def test_project_urls_point_at_the_repository():
    for label, url in PROJECT["urls"].items():
        assert url.startswith(f"https://github.com/{REPOSITORY}"), label


def test_the_installed_package_agrees_with_pyproject():
    try:
        dist = metadata.distribution("flowcellar")
    except metadata.PackageNotFoundError:
        pytest.skip("flowcellar is not installed as a distribution here")
    assert dist.version == flowcellar.__version__
    scripts = {e.name: e.value for e in dist.entry_points if e.group == "console_scripts"}
    assert scripts == {"flowcellar": "flowcellar.cli:main"}


def test_nothing_private_can_reach_a_distribution():
    sdist = PYPROJECT["tool"]["hatch"]["build"]["targets"]["sdist"]
    wheel = PYPROJECT["tool"]["hatch"]["build"]["targets"]["wheel"]
    assert wheel["packages"] == ["src/flowcellar"], "a wheel holds the package and nothing else"
    for folder in ("/design", "/prototype", "/local"):
        assert folder in sdist["exclude"], folder
        assert not [i for i in sdist["include"] if i.startswith(folder)], folder
    for target in (sdist, wheel):
        assert "src/flowcellar/viewer_dist/**" in target["artifacts"], "the built viewer is git-ignored and must be named"
    ignored = (REPO / ".gitignore").read_text("utf-8").split()
    for folder in ("design/", "prototype/", "local/"):
        assert folder in ignored, f"{folder} must stay out of git"


def test_the_build_hook_stops_a_package_without_the_viewer(tmp_path):
    spec = importlib.util.spec_from_file_location("flowcellar_hatch_build", REPO / "hatch_build.py")
    hook = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hook)
    assert PYPROJECT["tool"]["hatch"]["build"]["hooks"]["custom"]["path"] == "hatch_build.py"
    problems = hook.check_release_files(str(tmp_path))
    assert len(problems) == 1 and "the built viewer is missing" in problems[0] and "node build.mjs" in problems[0]
    dist = tmp_path / "src" / "flowcellar" / "viewer_dist"
    dist.mkdir(parents=True)
    files = {}
    for name in hook.VIEWER_FILES:
        (dist / name).write_text(name)
        files[name] = {"sha256": hashlib.sha256(name.encode()).hexdigest()}
    (dist / "manifest.json").write_text(json.dumps({"files": files}))
    assert hook.check_release_files(str(tmp_path)) == [
        "src/flowcellar/catalogue/metrics_nanopore.json is missing; write it with: python tools/build_catalogue.py",
        "src/flowcellar/sc/catalogue.json is missing; write it with: python tools/build_sc_catalogue.py"]
    (tmp_path / "src" / "flowcellar" / "catalogue").mkdir()
    (tmp_path / "src" / "flowcellar" / "catalogue" / "metrics_nanopore.json").write_text("{}")
    assert hook.check_release_files(str(tmp_path)) == ["src/flowcellar/sc/catalogue.json is missing; write it with: python tools/build_sc_catalogue.py"]
    (tmp_path / "src" / "flowcellar" / "sc").mkdir()
    (tmp_path / "src" / "flowcellar" / "sc" / "catalogue.json").write_text("{}")
    assert hook.check_release_files(str(tmp_path)) == []
    (dist / "viewer.css").write_text("changed after the build")
    assert any("viewer.css does not match" in p for p in hook.check_release_files(str(tmp_path)))


# ------------------------------------------------------------------------------------------ workflows
def _steps(job: dict) -> list[dict]:
    return job.get("steps") or []


def _uses(wf: dict) -> list[str]:
    return [s["uses"] for job in wf["jobs"].values() for s in _steps(job) if "uses" in s]


def test_the_publish_job_matches_the_trusted_publisher_exactly():
    assert (WORKFLOWS / WORKFLOW_FILE).exists() or not WORKFLOWS.exists(), f"the workflow file must be named {WORKFLOW_FILE}"
    wf = _workflow(WORKFLOW_FILE)
    assert wf["permissions"] == {"contents": "read"}, "nothing is granted to the workflow as a whole"
    publish = wf["jobs"]["publish"]
    env = publish["environment"]
    assert (env["name"] if isinstance(env, dict) else env) == ENVIRONMENT
    assert publish["permissions"] == {"id-token": "write"}, "the one permission trusted publishing needs, and no other"
    assert publish["needs"] == "build"
    cond = publish["if"]
    assert "github.event_name == 'push'" in cond and "startsWith(github.ref, 'refs/tags/v')" in cond
    assert f"github.repository == '{REPOSITORY}'" in cond
    assert " || " not in cond, "every condition must hold"
    steps = _steps(publish)
    assert [s.get("uses", "").split("@")[0] for s in steps] == ["actions/download-artifact", "pypa/gh-action-pypi-publish"]
    upload = steps[-1]
    assert upload["uses"] == "pypa/gh-action-pypi-publish@release/v1"
    assert not upload.get("with"), "no password, no user, no other index: the defaults are trusted publishing to PyPI"
    assert not any("run" in s for s in steps), "the job that can publish runs no script of this repository"
    for name, job in wf["jobs"].items():
        if name != "publish":
            assert "id-token" not in (job.get("permissions") or {}), name


def test_publishing_happens_only_for_a_version_tag():
    wf = _workflow(WORKFLOW_FILE)
    on = wf["on"]
    assert set(on) == {"push", "workflow_dispatch"}
    assert set(on["push"]) == {"tags"} and all(t.startswith("v") for t in on["push"]["tags"]), "no branch push starts a release"
    build = wf["jobs"]["build"]
    runs = "\n".join(s.get("run", "") for s in _steps(build))
    assert 'test "v${version}" = "${GITHUB_REF_NAME}"' in runs, "a tag that differs from the package version stops the build"
    for needed in ("npm ci", "node build.mjs", "python -m pytest", "python -m build", "twine check", "tools/check_dist.py", "tools/smoke_test.py"):
        assert needed in runs, needed
    assert runs.index("node build.mjs") < runs.index("python -m pytest") < runs.index("python -m build")
    assert [s["with"]["name"] for s in _steps(build) if s.get("uses", "").startswith("actions/upload-artifact")] == ["dist"]


def test_no_secret_and_no_other_way_to_publish():
    if not WORKFLOWS.exists():
        pytest.skip("no workflow files here")
    names = sorted(p.name for p in WORKFLOWS.iterdir())
    assert names == ["build.yml", "ci.yml"]
    for path in WORKFLOWS.iterdir():
        text = path.read_text("utf-8")
        code = "\n".join(line.split(" #", 1)[0] for line in text.splitlines() if not line.lstrip().startswith("#"))
        assert "secrets." not in code and "password" not in code.lower() and "api-token" not in code.lower(), path.name
        assert "twine upload" not in code and "upload.pypi.org" not in code, path.name
        if path.name != WORKFLOW_FILE:
            assert "gh-action-pypi-publish" not in code and "id-token" not in code and "environment:" not in code, \
                f"{path.name} must not be able to publish"


def test_the_actions_are_at_their_current_major_versions():
    for name in ("ci.yml", WORKFLOW_FILE):
        for use in _uses(_workflow(name)):
            action, ref = use.split("@")
            assert ACTIONS.get(action) == ref, use


def test_ci_tests_every_python_the_metadata_names_with_the_viewer_built_first():
    wf = _workflow("ci.yml")
    assert wf["permissions"] == {"contents": "read"}
    assert set(wf["on"]) == {"push", "pull_request"}
    test = wf["jobs"]["test"]
    declared = sorted(c.rsplit(" :: ", 1)[1] for c in PROJECT["classifiers"] if re.fullmatch(r"Programming Language :: Python :: 3\.\d+", c))
    assert sorted(test["strategy"]["matrix"]["python"]) == declared
    assert declared[0] == PROJECT["requires-python"].removeprefix(">=")
    order = [s.get("name") or s.get("uses") for s in _steps(test)]
    assert order.index("Build the viewer") < order.index("Tests")
    runs = "\n".join(s.get("run", "") for job in wf["jobs"].values() for s in _steps(job))
    for needed in ("ruff check", "tools/build_catalogue.py --check", "tools/build_docs.py --check", "tools/licence_check.py",
                   "python -m pytest tests", "python -m build", "tools/check_dist.py", "tools/smoke_test.py"):
        assert needed in runs, needed


# -------------------------------------------------------------------------------------- documentation
def test_the_generated_pages_are_up_to_date():
    docs = _tool("build_docs")
    page = (REPO / "docs" / "commands.md").read_text("utf-8")
    for name in docs.command_names():
        assert f"## flowcellar {name}\n" in page, name
    assert (REPO / "docs" / "metrics.md").read_text("utf-8") == docs.metrics_page(), "run python tools/build_docs.py"
    made_with = re.search(r"\(click ([0-9.]+)\)", page.splitlines()[0]).group(1)
    if made_with != metadata.version("click"):
        pytest.skip(f"the commands page was written with click {made_with}; the layout of --help differs between versions")
    assert page == docs.commands_page(), "run python tools/build_docs.py"


def test_every_option_of_every_command_is_in_the_command_reference():
    import click

    from flowcellar.cli import cli
    page = (REPO / "docs" / "commands.md").read_text("utf-8")
    ctx = click.Context(cli)
    for name in cli.list_commands(ctx):
        section = page.split(f"## flowcellar {name}\n", 1)[1].split("\n## ", 1)[0]
        for param in cli.get_command(ctx, name).params:
            for opt in getattr(param, "opts", []):
                if opt.startswith("-"):
                    assert opt in section, f"{name} {opt}"


def test_every_configuration_example_is_accepted_without_a_warning():
    from flowcellar.config import loader
    blocks = []
    for path in [REPO / "README.md", *sorted((REPO / "docs").glob("*.md"))]:
        blocks += [(path.name, b) for b in re.findall(r"```toml\n(.*?)```", path.read_text("utf-8"), re.S)]
    assert len(blocks) >= 4
    for name, block in blocks:
        loaded = loader.merge(tomllib.loads(block), source=name)
        assert loaded.warnings == [] and loaded.notices == [], (name, loaded.warnings, loaded.notices)


def test_the_links_between_the_pages_lead_somewhere():
    for path in [REPO / "README.md", *sorted((REPO / "docs").glob("*.md"))]:
        for target in re.findall(r"\]\(([^)#]+)(?:#[^)]*)?\)", path.read_text("utf-8")):
            if "://" in target:
                continue
            assert (path.parent / target).exists(), f"{path.name} links to {target}"


README = (REPO / "README.md").read_text("utf-8")
IMG = REPO / "docs" / "img"
MAX_PICTURE, MAX_PICTURES = 350_000, 3_000_000


def _outside_code(text: str) -> str:
    return "".join(re.split(r"```.*?```", text, flags=re.S))


def _picture_targets(text: str) -> list[str]:
    plain = _outside_code(text)
    return re.findall(r'\b(?:src|srcset)="([^"]+)"', plain) + re.findall(r"!\[[^\]]*\]\(([^)\s]+)\)", plain)


def test_the_readme_covers_what_a_new_user_needs():
    # what it is, how to install it, how to get from nothing to a library in the browser, what goes in and comes out
    for needed in ("## The problem", "## The idea", "## Install", "## Quick start", "## Try it without your own data", "## What you look at",
                   "## What goes in", "## What comes out", "## Works with what you have", "## Where names come from", "## What a verdict means",
                   "## Commands", "## Configuration", "## Privacy and requirements", "## Contributing", "## Licence and contact",
                   "flowcellar scan", "flowcellar serve", "flowcellar run", "flowcellar analysis", "python -m flowcellar.testing.demo",
                   "hub/index.html", "ssh -N -L 8737:127.0.0.1:8737", "Run health", "Read data", "Save as PDF", "flowcellar.toml"):
        assert needed in README, needed
    # how to get it and how to start come before the tour of the pages
    assert README.index("## Install") < README.index("## Quick start") < README.index("## What you look at")
    assert "not built yet" not in README and "CHANGELOG" not in README
    assert "github.com/MustafaElshani/FlowCellar/issues" in README.split("## Contributing", 1)[1].split("\n## ", 1)[0]
    head = README.split("\n## ", 1)[0]
    assert "flowcellar-lockup.svg" in head and "flowcellar-lockup-dark.svg" in head, "the logo, for a light and for a dark page"
    assert 'alt="FlowCellar"' in head and len(_picture_targets(head)) >= 3, "the logo and one picture of the product come first"
    # one sentence says what it is, before any heading and any picture of the product
    said = re.sub(r"<[^>]+>", " ", head.split("</p>", 1)[1].split("<p", 1)[0]).strip()
    assert said.startswith("FlowCellar is ") and len(re.findall(r"[.!?](?:\s|$)", said)) == 1 and len(said) < 320, said
    # the reader's kinds of run, and the promise that nothing is written where the data is
    for word in ("genomic DNA", "cDNA", "direct RNA", "single-cell", "spatial"):
        assert word in README.split("## The problem", 1)[1].split("\n## ", 1)[0], word
    assert "Nothing is copied from, moved in or written to the folders that are read, and no command contacts the network." in README
    # a front door, not the manual: it stays readable in a few minutes, and detail is behind links to docs/
    words = len(re.findall(r"\w+", re.sub(r"<[^>]+>|```.*?```", " ", README, flags=re.S)))
    assert 1500 < words < 3600, words
    for page in ("quick-start", "library", "configuration", "commands", "metrics", "development"):
        assert f"(docs/{page}.md)" in README, page


def test_every_picture_of_the_readme_exists_has_words_and_is_inside_the_budget():
    targets = _picture_targets(README)
    assert len(targets) >= 9
    for t in targets:
        assert "://" not in t, f"{t}: pictures are kept in the repository, by a path that works on GitHub and in a checkout"
        f = REPO / t
        assert f.is_file(), f"README shows {t}, which is not there"
        assert f.suffix in (".png", ".svg"), t
        if f.suffix == ".png":
            assert f.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", t
            assert f.stat().st_size <= MAX_PICTURE, f"{t} is {f.stat().st_size // 1000} kB; a picture is {MAX_PICTURE // 1000} kB at most"
        else:
            svg = f.read_text("utf-8")
            assert svg.lstrip().startswith("<svg") and "<script" not in svg and "http://www.w3.org/2000/svg" in svg, t
    # every picture says in words what it shows (for a reader who cannot see it, and for a page where it does not load)
    for tag in re.findall(r"<img\b[^>]*>", _outside_code(README)):
        alt = re.search(r'alt="([^"]*)"', tag)
        assert alt and len(alt.group(1)) >= 8, tag[:120]
    # a light and a dark version of what is drawn with ink: the logo and the illustration
    assert len(re.findall(r'<source media="\(prefers-color-scheme: dark\)" srcset="[^"]+-dark\.svg">', README)) == 2
    # the folder holds what the README shows and nothing else; all of it together stays small (it is in the source distribution)
    shown = {Path(t).name for t in targets if t.startswith("docs/img/")}
    kept = {f.name for f in IMG.iterdir()}
    assert kept == shown, f"docs/img and the README differ: {sorted(kept ^ shown)}"
    assert sum(f.stat().st_size for f in IMG.iterdir()) <= MAX_PICTURES


def _command_lines(text: str) -> list[str]:
    """Every line of the README that a reader would type: in a block of commands, or written as code in a sentence."""
    lines = []
    for block in re.findall(r"```([a-z]*)\n(.*?)```", text, re.S):
        if block[0] in ("", "sh", "bash", "console"):
            whole = re.sub(r"\\\n\s*", "", block[1])                # a command that goes on in the next line is one command
            lines += [ln.split("  #", 1)[0].strip() for ln in whole.splitlines()]
    lines += [c.strip() for c in re.findall(r"`([^`\n]+)`", _outside_code(text))]
    # "flowcellar 0.1.0 · bundle schema ..." is what doctor prints, not something to type
    return [ln for ln in lines if re.match(r"(?:flowcellar (?!\d)|python3? -m flowcellar\b)", ln)]


def test_every_command_line_of_the_readme_names_a_command_and_options_that_exist():
    import shlex

    import click

    from flowcellar.cli import cli
    from flowcellar.testing import demo
    ctx = click.Context(cli)
    top = {o for p in cli.params for o in getattr(p, "opts", [])}
    takes_value = {o for p in cli.params if not getattr(p, "is_flag", False) for o in getattr(p, "opts", [])}
    commands = set(cli.list_commands(ctx))
    lines = _command_lines(README)
    assert len(lines) >= 14, lines
    seen = set()
    for line in lines:
        words = shlex.split(line.replace("…", ""))
        if words[0].startswith("python"):
            assert words[1:3] == ["-m", "flowcellar.testing.demo"], line
            known = set(demo.parser()._option_string_actions)
            for w in words[3:]:
                assert not w.startswith("-") or w in known, f"{line}: no option {w}"
            continue
        rest = words[1:]
        while rest and rest[0].startswith("-"):                     # options of flowcellar itself come before the command
            assert rest[0] in top, f"{line}: flowcellar has no option {rest[0]}"
            rest = rest[2:] if rest[0] in takes_value else rest[1:]
        if not rest:
            continue
        assert rest[0] in commands, f"{line}: no command {rest[0]}"
        seen.add(rest[0])
        options = {o for p in cli.get_command(ctx, rest[0]).params for o in getattr(p, "opts", []) + getattr(p, "secondary_opts", [])}
        for w in rest[1:]:
            if w.startswith("-") and not re.fullmatch(r"-?\d+(\.\d+)?", w):
                assert w.split("=", 1)[0] in options, f"{line}: {rest[0]} has no option {w}"
    assert {"scan", "serve", "run", "analysis", "discover", "ls", "show", "init", "recompute", "doctor"} <= seen, seen


def test_on_pypi_every_path_of_the_readme_is_an_address():
    spec = importlib.util.spec_from_file_location("flowcellar_hatch_build_readme", REPO / "hatch_build.py")
    hook = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hook)
    assert hook.REPOSITORY == REPOSITORY
    # a release points at its tag, anything else at the main branch
    assert [hook.release_ref(v) for v in ("0.1.0", "1.2.3rc1", "0.1.0.dev0", "0.2.0.dev3+g1234", "")] == ["v0.1.0", "v1.2.3rc1", "main", "main", "main"]
    out = hook.absolute_readme(README, "0.1.0")
    plain = _outside_code(out)
    targets = re.findall(r'\b(?:src|srcset|href)="([^"]+)"', plain) + re.findall(r"\]\(([^)\s]+)\)", plain)
    assert len(targets) > 25
    for t in targets:
        assert t.startswith(("https://", "#")), f"{t} leads nowhere on PyPI"
    raw, page = f"https://raw.githubusercontent.com/{REPOSITORY}/v0.1.0/", f"https://github.com/{REPOSITORY}/blob/v0.1.0/"
    assert f'src="{raw}docs/img/home.png"' in out and f'src="{raw}viewer/brand/flowcellar-lockup.svg?sanitize=true"' in out
    assert f"]({page}docs/quick-start.md)" in out and f"]({page}LICENSE)" in out and "](#configuration)" in out
    assert f"/{REPOSITORY}/main/docs/img/home.png" in hook.absolute_readme(README, "0.1.0.dev0")
    # nothing else changes: the commands and the configuration example are as written
    assert re.findall(r"```.*?```", out, re.S) == re.findall(r"```.*?```", README, re.S)
    assert hook.absolute_readme("[a](./docs/x.md) ![b](docs/img/y.png) [c](https://example.org/z) <a href=\"mailto:a@example.org\">m</a>", "1.0") == \
        f"[a](https://github.com/{REPOSITORY}/blob/v1.0/docs/x.md) ![b](https://raw.githubusercontent.com/{REPOSITORY}/v1.0/docs/img/y.png) " \
        "[c](https://example.org/z) <a href=\"mailto:a@example.org\">m</a>"
    # what PyPI's own renderer makes of it: every picture and the link to the example survive
    renderer = pytest.importorskip("readme_renderer.markdown", reason="readme_renderer[md] shows what PyPI keeps")
    html = renderer.render(out)
    if html is None:
        pytest.skip("readme_renderer is installed without its Markdown extra")
    assert html.count("<img ") == len(re.findall(r"<img\b", plain)) and 'alt="FlowCellar"' in html
    assert f'src="{raw}docs/img/compare.png"' in html and "<table>" in html and "relative" not in html.split("<h2", 1)[0]


def test_the_built_wheel_carries_the_readme_with_addresses():
    wheel = os.environ.get("FLOWCELLAR_WHEEL")
    if not wheel:
        pytest.skip("set FLOWCELLAR_WHEEL to a built wheel to check its METADATA")
    import zipfile
    with zipfile.ZipFile(wheel) as z:
        meta = z.read(next(n for n in z.namelist() if n.endswith(".dist-info/METADATA"))).decode("utf-8")
    head, body = meta.split("\n\n", 1)
    assert "Description-Content-Type: text/markdown" in head
    version = re.search(r"^Version: (.+)$", head, re.M).group(1)
    ref = "main" if "dev" in version else f"v{version}"
    assert f"https://raw.githubusercontent.com/{REPOSITORY}/{ref}/docs/img/home.png" in body
    plain = _outside_code(body)
    for t in re.findall(r'\b(?:src|srcset|href)="([^"]+)"', plain) + re.findall(r"\]\(([^)\s]+)\)", plain):
        assert t.startswith(("https://", "#")), t


def test_the_changelog_has_an_entry_for_this_version():
    base = flowcellar.__version__.split(".dev")[0]
    text = (REPO / "CHANGELOG.md").read_text("utf-8")
    assert re.search(rf"^## {re.escape(base)}\b", text, re.M), base
    entry = text.split(f"## {base}", 1)[1].split("\n## ", 1)[0]
    assert "### Known limits" in entry, "an honest entry says what does not work yet"


# ----------------------------------------------------------------------------------------------- tools
def test_the_notices_name_everything_the_viewer_bundles():
    r = subprocess.run([sys.executable, str(REPO / "tools" / "licence_check.py")], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    notices = (REPO / "THIRD_PARTY_NOTICES.md").read_text("utf-8")
    for p in _tool("licence_check").bundled_packages():
        assert p["licence"] == "MIT" and f"| {p['version']} | MIT |" in notices and p["name"] in notices, p


def test_check_dist_finds_what_must_not_ship(tmp_path, monkeypatch):
    check = _tool("check_dist")
    private = tmp_path / "private_patterns.txt"
    private.write_text("# words that must not ship\nzz" + "secret\n", "utf-8")
    monkeypatch.setenv("FLOWCELLAR_PRIVATE_PATTERNS", str(private))
    files = {"flowcellar/a.py": b"x = 1\n", "flowcellar/b.txt": b"made with ZZ" + b"Secret\n",
             "flowcellar/c.json": b'{"fc": "' + b"P" + b"AX12345" + b'"}', "flowcellar/d.py": b"TST00001 0a1b2c3d\n"}
    problems = check.check_wording(files)
    assert len(problems) == 2 and "flowcellar/b.txt" in problems[0] and "flowcellar/c.json" in problems[1]
    # names of a real library: a long name anywhere, a short one only as a word of its own
    assert check.name_pattern("sample_x9").search(b"path/to/Sample_X9_sorted.bam")
    assert check.name_pattern("wt").search(b"the WT sample") and not check.name_pattern("wt").search(b"lowtide")
    assert check.NOT_A_NAME.fullmatch("barcode07") and check.NOT_A_NAME.fullmatch("no_sample") and not check.NOT_A_NAME.fullmatch("alpha")
    hits = check.search_names({"a.whl": {"x/readme.md": b"see run_alpha_7 here"}}, {"run_alpha_7": "run name", "zz": "sample"})
    assert hits == [("a.whl", "x/readme.md", "run name", "run_alpha_7")]


def test_the_network_guard_of_the_smoke_test_refuses_other_computers(tmp_path):
    smoke = _tool("smoke_test")
    (tmp_path / "sitecustomize.py").write_text(smoke.GUARD)
    log = tmp_path / "net.log"
    code = ("import socket\n"
            "for target in (('192.0.2.1', 9), ('127.0.0.1', 9)):\n"
            "    s = socket.socket(); s.settimeout(0.2); refused = False\n"
            "    try:\n"
            "        s.connect(target)\n"
            "    except OSError as e:\n"
            "        refused = 'smoke test' in str(e)\n"
            "    finally:\n"
            "        s.close()\n"
            "    print(target[0], refused)\n"
            "try:\n"
            "    socket.getaddrinfo('host.invalid', 80)\n"
            "except OSError as e:\n"
            "    print('lookup', 'smoke test' in str(e))\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60,
                       env={"PYTHONPATH": str(tmp_path), "FC_SMOKE_NETLOG": str(log), "PATH": "/usr/bin:/bin"})
    assert r.returncode == 0, r.stderr
    assert r.stdout.split() == ["192.0.2.1", "True", "127.0.0.1", "False", "lookup", "True"], r.stdout
    lines = log.read_text().splitlines()
    assert [ln.split("\t")[0] for ln in lines] == ["refused", "local", "refused"]
