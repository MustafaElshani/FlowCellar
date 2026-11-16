"""Build hooks: no package without the built viewer and the catalogue JSON, and a README that works on PyPI.


The viewer is built by ``viewer/build.mjs`` into ``src/flowcellar/viewer_dist/`` and is not kept in git, so a
checkout that was never built would otherwise give a package whose reports and hub cannot be written. When the
built viewer is missing, the hook builds it if Node and ``viewer/node_modules`` are there (it never downloads
anything), and otherwise stops the build with one line that says what to run. A development install
(``pip install -e .``) is let through: it may come before the first viewer build.

The second hook writes the long description of the package. README.md points at its pictures and at the other
pages with paths inside the repository, which is right on GitHub and in a checkout and leads nowhere on PyPI. In the
metadata of a package every such path is an address: a picture at ``raw.githubusercontent.com``, a page at
``github.com/.../blob``, at the tag of the release (``v0.1.0``) or, for a development version, at ``main``.

Hatchling loads this file by its name (``[tool.hatch.build.hooks.custom]`` and ``[tool.hatch.metadata.hooks.custom]``
in pyproject.toml).
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess

try:
    from hatchling.builders.hooks.plugin.interface import BuildHookInterface
    from hatchling.metadata.plugin.interface import MetadataHookInterface
except ImportError:  # the checks below are also used by the tests, where hatchling need not be installed
    BuildHookInterface = MetadataHookInterface = object

REPOSITORY = "MustafaElshani/FlowCellar"

VIEWER_DIST = ("src", "flowcellar", "viewer_dist")
VIEWER_FILES = ("viewer.core.js", "viewer.css", "template.html", "favicon.svg")
CATALOGUE = ("src", "flowcellar", "catalogue", "metrics_nanopore.json")
SC_CATALOGUE = ("src", "flowcellar", "sc", "catalogue.json")
HOW = "build it first: cd viewer && npm ci && node build.mjs"


def check_release_files(root: str) -> list[str]:
    """What is wrong with the files a release must carry; an empty list when all is well."""
    problems: list[str] = []
    dist = os.path.join(root, *VIEWER_DIST)
    manifest_path = os.path.join(dist, "manifest.json")
    if not os.path.isfile(manifest_path):
        return [f"the built viewer is missing ({'/'.join(VIEWER_DIST)}/manifest.json); {HOW}"]
    try:
        with open(manifest_path, encoding="utf-8") as fh:
            files = json.load(fh)["files"]
    except (ValueError, KeyError, TypeError) as exc:
        return [f"{'/'.join(VIEWER_DIST)}/manifest.json cannot be read ({exc.__class__.__name__}); {HOW}"]
    for name in VIEWER_FILES:
        path = os.path.join(dist, name)
        if name not in files or not os.path.isfile(path):
            problems.append(f"the built viewer has no {name}; {HOW}")
            continue
        with open(path, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        if digest != files[name].get("sha256"):
            problems.append(f"{name} does not match the hash in the viewer's manifest.json; {HOW}")
    if not os.path.isfile(os.path.join(root, *CATALOGUE)):
        problems.append(f"{'/'.join(CATALOGUE)} is missing; write it with: python tools/build_catalogue.py")
    if not os.path.isfile(os.path.join(root, *SC_CATALOGUE)):
        problems.append(f"{'/'.join(SC_CATALOGUE)} is missing; write it with: python tools/build_sc_catalogue.py")
    return problems


def build_viewer_if_possible(root: str) -> bool:
    """Run the viewer build when the built viewer is missing and everything it needs is already on this machine."""
    if os.path.isfile(os.path.join(root, *VIEWER_DIST, "manifest.json")):
        return False
    node = shutil.which("node")
    script = os.path.join(root, "viewer", "build.mjs")
    if not node or not os.path.isfile(script) or not os.path.isdir(os.path.join(root, "viewer", "node_modules", "esbuild")):
        return False
    done = subprocess.run([node, script], cwd=os.path.join(root, "viewer"), capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError("flowcellar cannot be packaged: the viewer build failed:\n" + done.stdout + done.stderr)
    return True


class ReleaseFilesHook(BuildHookInterface):
    PLUGIN_NAME = "custom"

    def initialize(self, version: str, build_data: dict) -> None:
        if self.target_name == "wheel" and version == "editable":
            return
        if build_viewer_if_possible(self.root):
            self.app.display_info("flowcellar: built the viewer (it was missing)")
        problems = check_release_files(self.root)
        if problems:
            raise RuntimeError("flowcellar cannot be packaged: " + "; ".join(problems))


# --------------------------------------------------------------------------------------- the README on PyPI
_ABSOLUTE = re.compile(r"^(?:[a-z][a-z0-9+.-]*:|//|#)", re.I)


def release_ref(version: str) -> str:
    """Where the files of a version are: the tag ``v<version>`` of a release, the branch ``main`` of anything else."""
    return f"v{version}" if re.fullmatch(r"\d+(?:\.\d+)*(?:(?:a|b|rc)\d+)?(?:\.post\d+)?", version or "") else "main"


def absolute_readme(text: str, version: str, repository: str = REPOSITORY) -> str:
    """``text`` with every path inside the repository written as an address: pictures (``src``, ``srcset``,
    ``![..](..)``) as raw files, everything else that is linked as a page of the repository. Addresses, anchors
    inside the page and the text of code are left alone."""
    ref = release_ref(version)
    raw, page = f"https://raw.githubusercontent.com/{repository}/{ref}/", f"https://github.com/{repository}/blob/{ref}/"

    def picture(target: str) -> str:
        if _ABSOLUTE.match(target):
            return target
        path = target.lstrip("./") if target.startswith("./") else target
        return raw + path + ("?sanitize=true" if path.lower().endswith(".svg") else "")     # an SVG is served as a picture only when asked

    def link(target: str) -> str:
        if _ABSOLUTE.match(target):
            return target
        return page + (target[2:] if target.startswith("./") else target)

    out = []
    for i, part in enumerate(re.split(r"(```.*?```)", text, flags=re.S)):
        if i % 2:                                           # a block of code: as it is
            out.append(part)
            continue
        part = re.sub(r'\b(src|srcset)="([^"]+)"', lambda m: f'{m.group(1)}="{picture(m.group(2))}"', part)
        part = re.sub(r'\bhref="([^"]+)"', lambda m: f'href="{link(m.group(1))}"', part)
        part = re.sub(r"(!\[[^\]]*\]\()([^)\s]+)(\))", lambda m: m.group(1) + picture(m.group(2)) + m.group(3), part)
        part = re.sub(r"((?<!!)\[[^\]]*\]\()([^)\s]+)(\))", lambda m: m.group(1) + link(m.group(2)) + m.group(3), part)
        out.append(part)
    return "".join(out)


def package_version(root: str) -> str:
    with open(os.path.join(root, "src", "flowcellar", "_version.py"), encoding="utf-8") as fh:
        return re.search(r'__version__\s*=\s*"([^"]+)"', fh.read()).group(1)


class ReadmeHook(MetadataHookInterface):
    PLUGIN_NAME = "custom"

    def update(self, metadata: dict) -> None:
        with open(os.path.join(self.root, "README.md"), encoding="utf-8") as fh:
            text = fh.read()
        metadata["readme"] = {"content-type": "text/markdown", "text": absolute_readme(text, package_version(self.root))}
