#!/usr/bin/env python3
"""Check a built wheel and source distribution before they go anywhere.

    python tools/check_dist.py dist/*                         what every release must satisfy
    python tools/check_dist.py dist/* --list                  also print every file with its size
    python tools/check_dist.py dist/* --names-from LIBRARY    also search for the names of a real library
    python tools/check_dist.py FOLDER --names-from LIBRARY    the same search in a folder of files (a checkout)

What is checked, for each archive:

* it holds what it must: the package, the built viewer with hashes that match its manifest, the catalogue JSON,
  the bundle schema, the default configuration, the licence and the third-party notices, the console script;
  a source distribution also the viewer sources, the catalogue source, the tools, the tests and the docs;
* it holds nothing it must not: no file from design/, prototype/ or local/, no node_modules, no cache;
* its metadata says what pyproject.toml promises: name, Python 3.11 or newer, exactly numpy, pyarrow and click
  as run-time dependencies, pysam in the extra ``bam``, jsonschema in ``validate``;
* no file in it holds an id that looks like a real flow cell (invented ids in this project start with TST), or
  matches a pattern of the private list that FLOWCELLAR_PRIVATE_PATTERNS names (see ``private_patterns``);
* its long description (the README as PyPI shows it) points at every picture and page by an address, never by a
  path inside the repository, and fetches pictures from the repository only; in a source distribution the pictures
  are there, each at most 350 kB and 3 MB together.

``--names-from LIBRARY`` reads the index of a FlowCellar library made from real data (flow cell ids, run keys,
experiment, run, sample and project names, device positions, root paths) and searches every file of every
archive for each of them: a name of five characters or more anywhere in a file, a shorter one as a word of its
own, both without regard to case. A hit fails the check. What MinKNOW writes when no name was given
(``no_sample``, ``barcode01``, ``unclassified``) is not a name; a name that is an ordinary word of the product
can be set aside with ``--not-a-name WORD`` and is then reported as set aside. Parts of names (the pieces between
separators) are searched as well and reported for a person to look at; most are ordinary words, so they fail
the check only with ``--strict-parts``. Names are printed only with ``--show-names``: use it where they may be
seen. The library is only read.

Exit code 0 when everything holds, 1 otherwise. Standard library only.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tarfile
import zipfile
from pathlib import Path

VIEWER_FILES = ("viewer.core.js", "viewer.css", "template.html", "favicon.svg", "manifest.json")
PACKAGE_FILES = (
    "flowcellar/__init__.py", "flowcellar/_version.py", "flowcellar/cli/__init__.py",
    "flowcellar/catalogue/metrics_nanopore.json", "flowcellar/sc/catalogue.json", "flowcellar/config/defaults.toml",
    *(f"flowcellar/viewer_dist/{n}" for n in VIEWER_FILES),
)
SDIST_FILES = (
    "pyproject.toml", "hatch_build.py", "README.md", "LICENSE", "THIRD_PARTY_NOTICES.md", "CHANGELOG.md", "PKG-INFO",
    "catalogue/metrics_nanopore.yaml", "catalogue/metrics_singlecell_spatial.yaml", "viewer/build.mjs", "viewer/package.json", "viewer/package-lock.json",
    "viewer/src/main.jsx", "viewer/styles/tokens.css", "viewer/verify/query_vectors.mjs",
    "tools/build_catalogue.py", "tools/build_sc_catalogue.py", "tools/build_docs.py", "tools/check_dist.py", "tests/conftest.py",
    "tools/make_demo.py", "tools/make_readme_pictures.py", "tools/readme_shots.mjs", "tools/readme_illustration.py",
    "docs/index.md", "docs/commands.md", "docs/metrics.md", "docs/img/overview-light.svg", "docs/img/overview-dark.svg", "docs/img/home.png",
    "viewer/brand/flowcellar-lockup.svg", "viewer/brand/flowcellar-lockup-dark.svg",
)
# the pictures of the README travel in the source distribution: each stays small, and so does the folder
MAX_PICTURE_BYTES, MAX_PICTURES_BYTES = 350_000, 3_000_000
REPOSITORY_ADDRESSES = ("https://github.com/MustafaElshani/FlowCellar/", "https://raw.githubusercontent.com/MustafaElshani/FlowCellar/")
FORBIDDEN_PARTS = {"design", "prototype", "local", "node_modules", "pw-browsers", "__pycache__", ".venv", ".git",
                   ".github", ".pytest_cache", ".ruff_cache", ".cache", "gen", "dist"}
BASE_REQUIRES = {"numpy", "pyarrow", "click"}
EXTRAS = {"bam": {"pysam"}, "validate": {"jsonschema"}, "sc": {"h5py"}}
FORBIDDEN_TEXT = {
    "an id that looks like a real flow cell": re.compile(rb"(?<![A-Za-z0-9])[PF][A-Z]{2}[0-9]{5}(?![0-9A-Za-z])"),
}
# What a sample is called when it has no name of its own: these are words of the product, not names of a laboratory.
NOT_A_NAME = re.compile(r"unclassified|barcode\d{1,3}|no_sample(?:_id)?|none|unknown", re.IGNORECASE)


# ------------------------------------------------------------------------------------------- archives
def read_archive(path: Path) -> tuple[str, dict[str, bytes]]:
    """('wheel' | 'sdist' | 'folder', {member path without the sdist's top folder: bytes})."""
    if path.is_dir():
        return "folder", {str(f.relative_to(path)): f.read_bytes() for f in sorted(path.rglob("*")) if f.is_file()}
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as z:
            return "wheel", {n: z.read(n) for n in z.namelist() if not n.endswith("/")}
    if path.name.endswith(".tar.gz"):
        out = {}
        with tarfile.open(path) as t:
            for m in t.getmembers():
                if m.isfile():
                    out[m.name.split("/", 1)[1] if "/" in m.name else m.name] = t.extractfile(m).read()
                elif not m.isdir():
                    out[m.name + "  (not a regular file)"] = b""
        return "sdist", out
    raise SystemExit(f"{path}: neither a wheel (.whl) nor a source distribution (.tar.gz)")


def metadata_of(kind: str, files: dict[str, bytes]) -> tuple[str, dict[str, list[str]]]:
    name = next((n for n in files if n.endswith(".dist-info/METADATA")), None) if kind == "wheel" else "PKG-INFO"
    if not name or name not in files:
        return "", {}
    head = files[name].decode("utf-8").split("\n\n", 1)[0]
    fields: dict[str, list[str]] = {}
    for line in head.splitlines():
        if ": " in line and not line.startswith(" "):
            k, v = line.split(": ", 1)
            fields.setdefault(k, []).append(v)
    return name, fields


def check_description(kind: str, files: dict[str, bytes]) -> list[str]:
    """The long description (the README as PyPI shows it): every picture and every link is an address, because a
    path inside the repository leads nowhere there; in a source distribution every picture it shows is in the archive."""
    name = next((n for n in files if n.endswith(".dist-info/METADATA")), None) if kind == "wheel" else "PKG-INFO"
    if not name or name not in files:
        return []
    text = files[name].decode("utf-8")
    if "\n\n" not in text:
        return [f"{name}: no long description"]
    body = "".join(re.split(r"```.*?```", text.split("\n\n", 1)[1], flags=re.S))
    problems = []
    pictures = re.findall(r'\b(?:src|srcset)="([^"]+)"', body) + re.findall(r"!\[[^\]]*\]\(([^)\s]+)\)", body)
    links = re.findall(r'\bhref="([^"]+)"', body) + re.findall(r"(?<!!)\[[^\]]*\]\(([^)\s]+)\)", body)
    if not pictures:
        problems.append(f"{name}: the long description shows no picture")
    for t in pictures + links:
        if not t.startswith(("https://", "#")):
            problems.append(f"{name}: the long description points at {t!r}, a path that leads nowhere on PyPI")
    for t in pictures:
        if t.startswith("https://") and not t.startswith(REPOSITORY_ADDRESSES):
            problems.append(f"{name}: a picture of the long description is fetched from elsewhere: {t}")
        m = re.match(r"https://raw\.githubusercontent\.com/[^/]+/[^/]+/[^/]+/([^?#]+)", t)
        if kind == "sdist" and m and m.group(1) not in files:
            problems.append(f"{name}: the long description shows {m.group(1)}, which is not in the archive")
    return problems


def _requirement_name(req: str) -> str:
    return re.split(r"[ ;<>=!~\[(]", req.strip(), maxsplit=1)[0].lower().replace("_", "-")


def private_patterns() -> dict[str, re.Pattern]:
    """Words a maintainer keeps out of every release without writing them into the repository: one regular
    expression per line of the file that FLOWCELLAR_PRIVATE_PATTERNS names. Without the variable there are none."""
    path = os.environ.get("FLOWCELLAR_PRIVATE_PATTERNS")
    if not path:
        return {}
    lines = [s.strip() for s in Path(path).read_text("utf-8").splitlines()]
    return {f"private pattern {i}": re.compile(s.encode(), re.IGNORECASE) for i, s in enumerate(lines, 1) if s and not s.startswith("#")}


def check_wording(files: dict[str, bytes]) -> list[str]:
    problems = []
    forbidden = {**FORBIDDEN_TEXT, **private_patterns()}
    for n, data in files.items():
        for what, pat in forbidden.items():
            m = pat.search(data) or pat.search(n.encode())
            if m:
                problems.append(f"{what} in {n}: {m.group(0).decode('latin-1')!r}")
    return problems


def check_archive(path: Path, kind: str, files: dict[str, bytes]) -> list[str]:
    if kind == "folder":
        return check_wording(files)
    problems: list[str] = []
    prefix = "" if kind == "wheel" else "src/"
    for want in PACKAGE_FILES:
        if prefix + want not in files:
            problems.append(f"missing: {prefix + want}")
    if not any(re.fullmatch(re.escape(prefix) + r"flowcellar/model/schemas/bundle-[0-9.]+\.schema\.json", n) for n in files):
        problems.append("missing: the bundle schema (flowcellar/model/schemas/bundle-<version>.schema.json)")
    if kind == "sdist":
        problems += [f"missing: {want}" for want in SDIST_FILES if want not in files]
    else:
        for want in ("licenses/LICENSE", "licenses/THIRD_PARTY_NOTICES.md", "entry_points.txt", "RECORD", "WHEEL"):
            if not any(n.endswith(".dist-info/" + want) for n in files):
                problems.append(f"missing: <dist-info>/{want}")
        entry = next((files[n].decode() for n in files if n.endswith(".dist-info/entry_points.txt")), "")
        if not re.search(r"^flowcellar\s*=\s*flowcellar\.cli:main\s*$", entry, re.M):
            problems.append("the console script 'flowcellar = flowcellar.cli:main' is not declared")
        outside = sorted({n.split("/", 1)[0] for n in files} - {"flowcellar"} - {n.split("/", 1)[0] for n in files if ".dist-info/" in n})
        if outside:
            problems.append(f"the wheel holds top-level entries besides the package: {', '.join(outside)}")
    for n in files:
        parts = set(Path(n).parts[:-1]) | ({Path(n).name} if Path(n).name in FORBIDDEN_PARTS else set())
        bad = parts & FORBIDDEN_PARTS
        if bad:
            problems.append(f"must not be shipped ({', '.join(sorted(bad))}): {n}")
        if n.endswith((".pyc", ".pyo")) or "(not a regular file)" in n:
            problems.append(f"must not be shipped: {n}")
    # the viewer: every file is the one its manifest describes
    man_name = prefix + "flowcellar/viewer_dist/manifest.json"
    if man_name in files:
        try:
            listed = json.loads(files[man_name])["files"]
        except (ValueError, KeyError):
            listed = {}
            problems.append("the viewer's manifest.json cannot be read")
        for name, entry in listed.items():
            data = files.get(prefix + "flowcellar/viewer_dist/" + name)
            if data is None:
                problems.append(f"the viewer's manifest lists {name}, which is not in the archive")
            elif hashlib.sha256(data).hexdigest() != entry.get("sha256"):
                problems.append(f"the viewer file {name} does not match the hash in its manifest")
    # metadata
    meta_name, meta = metadata_of(kind, files)
    if not meta:
        problems.append("no package metadata found")
    else:
        if meta.get("Name") != ["flowcellar"]:
            problems.append(f"{meta_name}: Name is {meta.get('Name')}, not flowcellar")
        if meta.get("Requires-Python") != [">=3.11"]:
            problems.append(f"{meta_name}: Requires-Python is {meta.get('Requires-Python')}, not >=3.11")
        if meta.get("License-Expression") != ["GPL-3.0-or-later"]:
            problems.append(f"{meta_name}: License-Expression is {meta.get('License-Expression')}, not GPL-3.0-or-later")
        if sorted(meta.get("License-File", [])) != ["LICENSE", "THIRD_PARTY_NOTICES.md"]:
            problems.append(f"{meta_name}: License-File is {meta.get('License-File')}")
        base, extras = set(), {}
        for req in meta.get("Requires-Dist", []):
            m = re.search(r"extra\s*==\s*['\"]([^'\"]+)['\"]", req)
            (extras.setdefault(m.group(1), set()) if m else base).add(_requirement_name(req))
        if base != BASE_REQUIRES:
            problems.append(f"{meta_name}: run-time dependencies are {sorted(base)}, not {sorted(BASE_REQUIRES)}")
        for extra, want in EXTRAS.items():
            if extras.get(extra) != want:
                problems.append(f"{meta_name}: extra '{extra}' requires {sorted(extras.get(extra, []))}, not {sorted(want)}")
    if kind == "sdist":
        imgs = {n: len(d) for n, d in files.items() if n.startswith("docs/img/")}
        problems += [f"{n} is {size // 1000} kB; a picture of the README is {MAX_PICTURE_BYTES // 1000} kB at most" for n, size in imgs.items()
                     if size > MAX_PICTURE_BYTES]
        if sum(imgs.values()) > MAX_PICTURES_BYTES:
            problems.append(f"docs/img holds {sum(imgs.values()) / 1e6:.1f} MB; the pictures of the README are {MAX_PICTURES_BYTES // 1_000_000} MB at most")
    return problems + check_description(kind, files) + check_wording(files)


# ---------------------------------------------------------------------------------------------- names
def _column(col) -> list:
    if isinstance(col, dict) and "dict" in col:
        return [col["dict"][c] if c is not None and c >= 0 else None for c in col["codes"]]
    return list(col or [])


def library_names(library: Path, set_aside: frozenset[str] = frozenset()) -> dict[str, set[str]]:
    """{kind of name: the names} from the index of a library made from real data."""
    with gzip.open(library / "index" / "index.json.gz") as fh:
        index = json.load(fh)
    names: dict[str, set[str]] = {}

    def add(kind: str, value) -> None:
        if isinstance(value, (list, tuple)):
            for v in value:
                add(kind, v)
        elif isinstance(value, str) and value.strip():
            names.setdefault(kind, set()).add(value.strip())

    cols = (index.get("runs") or {}).get("columns") or {}
    for kind, col in (("flow cell id", "fc"), ("run key", "id"), ("experiment", "experiment"), ("run name", "run_name"),
                      ("run folder", "folder"), ("project", "project"), ("project id", "project_ids"),
                      ("device position", "position")):
        add(kind, _column(cols.get(col)))
    for row in index.get("samples") or []:
        add("sample", [row.get("name"), row.get("alias")])
        add("project", row.get("project"))
    for row in index.get("projects") or []:
        add("project", row.get("name"))
        add("project id", row.get("id"))
    for row in index.get("analyses") or []:                # single-cell and spatial analyses: their own names
        add("analysis sample", [row.get("label"), row.get("group")])
        add("analysis folder", row.get("folder"))
        add("slide", row.get("slide"))
        for link in row.get("links") or []:
            add("flow cell id", link.get("fc"))
            add("run key", link.get("run"))
    for row in index.get("roots") or []:
        path = str(row.get("path") or "")
        add("root path", path)
        add("storage mount", next((p for p in path.split("/") if p), None))
    lib = index.get("library") or {}
    if lib.get("name") != library.resolve().name:          # a library that was never named is called after its folder
        add("library name", lib.get("name"))
    add("library path", lib.get("path"))
    add("host", (index.get("progress") or {}).get("host"))
    names.get("project id", set()).discard("_none")
    for kind in list(names):                                # words of the product itself are not names of a laboratory
        names[kind] = {n for n in names[kind] if not NOT_A_NAME.fullmatch(n) and n.lower() not in set_aside}
    return names


def name_pattern(name: str) -> re.Pattern[bytes]:
    """Five characters or more: anywhere in the text. Shorter: only as a word of its own (a two-letter sample name
    is inside countless ordinary words)."""
    core = re.escape(name.encode("utf-8"))
    if len(name) >= 5:
        return re.compile(core, re.IGNORECASE)
    return re.compile(rb"(?<![A-Za-z0-9])" + core + rb"(?![A-Za-z0-9])", re.IGNORECASE)


def pieces(names: dict[str, set[str]]) -> set[str]:
    """Parts of names between separators (four characters or more, not a number, not a date)."""
    whole = {n.lower() for group in names.values() for n in group}
    out = set()
    for kind, group in names.items():
        if kind in ("root path", "library path", "storage mount", "host"):
            continue
        for n in group:
            for piece in re.split(r"[^A-Za-z0-9]+", n):
                if len(piece) >= 4 and not piece.isdigit() and piece.lower() not in whole:
                    out.add(piece)
    return out


def search_names(archives: dict[str, dict[str, bytes]], terms: dict[str, str]) -> list[tuple[str, str, str, str]]:
    """[(archive, member, kind, name)] for every name found in a file or in a file name."""
    hits = []
    compiled = [(name, kind, name_pattern(name)) for name, kind in sorted(terms.items())]
    for archive, files in archives.items():
        for member, data in files.items():
            hay = data + b"\n" + member.encode()
            low = hay.lower()
            for name, kind, pat in compiled:
                if len(name) >= 5 and name.lower().encode() not in low:
                    continue
                if pat.search(hay):
                    hits.append((archive, member, kind, name))
    return hits


# ----------------------------------------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("archives", nargs="+", type=Path)
    ap.add_argument("--list", action="store_true", help="print every file of every archive with its size")
    ap.add_argument("--names-from", type=Path, action="append", default=[], metavar="LIBRARY",
                    help="a FlowCellar library made from real data: search the archives for its names")
    ap.add_argument("--not-a-name", action="append", default=[], metavar="WORD",
                    help="a name of that library that is an ordinary word of the product: set it aside (repeatable)")
    ap.add_argument("--strict-parts", action="store_true", help="fail when a part of a name is found, too")
    ap.add_argument("--show-names", action="store_true", help="print the names that were found (they are private)")
    ap.add_argument("--unpack", type=Path, metavar="DIR", help="also write the files of each archive under DIR")
    args = ap.parse_args(argv)

    failed = False
    opened: dict[str, dict[str, bytes]] = {}
    for path in args.archives:
        kind, files = read_archive(path)
        opened[path.name] = files
        size = sum(len(b) for b in files.values())
        if kind == "folder":
            print(f"{path}: folder, {len(files)} files, {size:,} bytes")
        else:
            print(f"{path.name}: {kind}, {len(files)} files, {path.stat().st_size:,} bytes ({size:,} unpacked)")
        if args.list:
            for n in sorted(files):
                print(f"  {len(files[n]):>9,}  {n}")
        if args.unpack:
            for n, data in files.items():
                if "(not a regular file)" in n:
                    continue
                target = args.unpack / path.name / n
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
        problems = check_archive(path, kind, files)
        if kind != "folder":
            viewer = [n for n in files if "/viewer_dist/" in n or n.startswith("flowcellar/viewer_dist/")]
            print(f"  the built viewer: {len(viewer)} files, {sum(len(files[n]) for n in viewer):,} bytes")
        for p in problems:
            print(f"  PROBLEM  {p}")
        if not problems:
            print("  ok: wording" if kind == "folder" else "  ok: contents, metadata and wording")
        failed |= bool(problems)

    set_aside = frozenset(w.lower() for w in args.not_a_name)
    for library in args.names_from:
        names = library_names(library, set_aside)
        terms = {n: kind for kind, group in names.items() for n in group if len(n) >= 2}
        counts = ", ".join(f"{len(g)} {k}" for k, g in sorted(names.items()))
        print(f"names of {library}: {len(terms)} ({counts})")
        if set_aside:
            print(f"  set aside as ordinary words: {', '.join(sorted(set_aside))}")
        hits = search_names(opened, terms)
        print(f"  whole names found: {len(hits)} hits of {len({h[3] for h in hits})} names")
        for archive, member, kind, name in hits:
            print(f"  FOUND  {name if args.show_names else f'<{kind}, {len(name)} characters>'}  in {archive}: {member}")
        part_terms = {p: "part of a name" for p in pieces(names) if p.lower() not in set_aside and not NOT_A_NAME.fullmatch(p)}
        part_hits = search_names(opened, part_terms)
        found_parts = sorted({h[3] for h in part_hits}, key=str.lower)
        print(f"  parts of names found: {len(found_parts)} of {len(part_terms)} parts, in {len({(h[0], h[1]) for h in part_hits})} files"
              + (" (to look at; most are ordinary words)" if found_parts and not args.strict_parts else ""))
        if args.show_names and found_parts:
            print("    " + ", ".join(found_parts))
        failed |= bool(hits) or (args.strict_parts and bool(part_hits))
    print("FAILED" if failed else "all checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
    raise SystemExit(main())
