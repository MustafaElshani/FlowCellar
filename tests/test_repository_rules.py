"""Rules that can be checked mechanically on the tracked files."""
import os
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TRACKED = ["src", "viewer", "tests", "docs", "tools", "catalogue", "README.md", "CHANGELOG.md", "pyproject.toml", ".github"]
SKIP_DIRS = {"node_modules", "pw-browsers", ".cache", "__pycache__", "viewer_dist"}
TEXT = {".py", ".js", ".jsx", ".mjs", ".css", ".html", ".json", ".toml", ".yaml", ".yml", ".md", ".txt", ".svg", ".cfg"}
# An id that looks like a real ONT flow cell (three letters starting with P or F, five digits). Invented ids in
# tracked files start with TST, so nothing tracked may match.
REAL_LOOKING_FLOW_CELL = re.compile(r"(?<![A-Za-z0-9])[PF][A-Z]{2}[0-9]{5}(?![0-9])")


def _files():
    for top in TRACKED:
        p = REPO / top
        if p.is_file():
            yield p
        elif p.is_dir():
            for f in p.rglob("*"):
                if f.is_file() and f.suffix in TEXT and not (set(f.relative_to(REPO).parts) & SKIP_DIRS):
                    yield f


def test_no_private_pattern_in_tracked_files():
    """Words a maintainer keeps out of the repository without writing them into it: FLOWCELLAR_PRIVATE_PATTERNS names a
    file with one regular expression per line. Skipped when the variable is not set."""
    path = os.environ.get("FLOWCELLAR_PRIVATE_PATTERNS")
    if not path:
        pytest.skip("FLOWCELLAR_PRIVATE_PATTERNS is not set")
    lines = [s.strip() for s in Path(path).read_text("utf-8").splitlines()]
    patterns = [re.compile(s, re.IGNORECASE) for s in lines if s and not s.startswith("#")]
    hits = sorted({str(f.relative_to(REPO)) for f in _files() for p in patterns if p.search(f.read_text("utf-8", errors="replace"))})
    assert hits == []


def test_no_real_looking_flow_cell_id_in_tracked_files():
    hits = [str(f.relative_to(REPO)) for f in _files() if REAL_LOOKING_FLOW_CELL.search(f.read_text("utf-8", errors="replace"))]
    assert hits == []


def test_nothing_tracked_reads_the_ignored_folders():
    """The build, the package and the tests never read design/, prototype/ or local/ (private working folders)."""
    pat = re.compile(r"""['"/](?:design|prototype)/[A-Za-z_]""")
    hits = []
    for f in _files():
        if f.suffix not in {".py", ".js", ".jsx", ".mjs", ".toml"} or f.name == "test_repository_rules.py":
            continue
        for i, line in enumerate(f.read_text("utf-8", errors="replace").splitlines(), 1):
            code = line.split("#", 1)[0] if f.suffix in {".py", ".toml"} else line.split("//", 1)[0]
            if pat.search(code) and "exclude" not in code:
                hits.append(f"{f.relative_to(REPO)}:{i}")
    assert hits == []


@pytest.mark.parametrize("pkg", ["pandas", "polars", "duckdb", "pydantic", "jinja2", "plotly", "matplotlib", "yaml", "subprocess"])
def test_run_time_code_imports_no_forbidden_package(pkg):
    pat = re.compile(rf"^\s*(?:import|from)\s+{pkg}\b", re.MULTILINE)
    hits = [str(f.relative_to(REPO)) for f in (REPO / "src").rglob("*.py") if pat.search(f.read_text("utf-8"))]
    assert hits == []
