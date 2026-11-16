# Releasing

A release is one pure-Python wheel and one source distribution, both with the built viewer and the
catalogue inside, uploaded to PyPI by `.github/workflows/build.yml` when a version tag is pushed. No token
or password is stored anywhere: PyPI trusts this one workflow of this one repository (trusted publishing).

## What must match

PyPI accepts an upload only from exactly this combination. The same four values are in `build.yml`; a test
(`tests/release/test_release_files.py`) keeps them there.

| Setting | Value |
|---|---|
| PyPI project | `flowcellar` |
| GitHub repository | `MustafaElshani/FlowCellar` |
| Workflow file | `build.yml` (in `.github/workflows/`) |
| GitHub environment | `pypi` |

## Once, before the first release

1. On PyPI, under "Publishing", add a *pending publisher* for the project name `flowcellar` with the four
   values above. The project itself is created by the first upload; until then the name is not reserved.
2. In the repository settings on GitHub, create the environment `pypi`. Add a required reviewer to it: the
   publish job then waits for a person to approve it, which is the manual gate of a release.
3. The licence is the GNU General Public License, version 3 or later (`LICENSE`, `license` in
   `pyproject.toml`, the README and the notice at the top of the viewer bundle all say so). Change it before the
   first upload or not at all: a release on PyPI cannot be taken back.
4. Make sure nothing that names a laboratory's data is in the repository (the folders `design/`,
   `prototype/` and `local/` are ignored by git for that reason and are never packaged). Words that must
   never ship can be kept in a private file, one regular expression per line: with
   `FLOWCELLAR_PRIVATE_PATTERNS=FILE` the tests and `tools/check_dist.py` search for them.

## Every release

1. Set the version in `src/flowcellar/_version.py` (for example `0.1.0`), give the entry in `CHANGELOG.md`
   its date, and regenerate the pages that quote the version: `python tools/build_docs.py`.
   In `README.md`: the line `doctor` prints (it quotes the version).
   The long description on PyPI is `README.md` with its paths written as addresses at the tag `v<version>`
   (`hatch_build.py`): the pictures it shows must be in the tagged commit, so make them again before tagging if the
   viewer changed (`python tools/make_readme_pictures.py`, see [Development](development.md)).
2. Run what the workflow will run:

   ```
   (cd viewer && npm ci && node build.mjs)
   python tools/build_catalogue.py --check && python tools/build_docs.py --check && python tools/licence_check.py
   ruff check src tests tools hatch_build.py && python -m pytest
   rm -rf dist && python -m build
   python tools/check_dist.py dist/*
   python -m venv /tmp/fc-release && /tmp/fc-release/bin/pip install "$(ls dist/*.whl)[all]"
   /tmp/fc-release/bin/python tools/smoke_test.py --installed
   FLOWCELLAR_WHEEL="$(ls dist/*.whl)" python -m pytest tests/release -k wheel     # the README inside the wheel
   ```

   Then follow the quick start of the README once, as written, in that new environment
   (`python -m flowcellar.testing.demo demo`, `flowcellar scan …`, `flowcellar serve …`).

   On a machine that holds a library made from real data, also search the two archives for its names
   (`python tools/check_dist.py dist/* --names-from LIBRARY`); the result must be zero whole names.
3. Commit, and push to `main`; wait for `ci.yml` to pass.
4. Tag the commit `v<version>` (for example `v0.1.0`) and push the tag. `build.yml` then:
   - builds the viewer with Node, runs the tests, and refuses to go on when the tag is not `v` followed by
     the version in `_version.py`;
   - builds the source distribution and the wheel, checks them (`twine check`, `tools/check_dist.py`) and
     runs `tools/smoke_test.py` on the wheel installed in an environment of its own;
   - hands both files to the job `publish`, which runs in the environment `pypi`, has the one permission
     `id-token: write`, and uploads them. This job runs only for a pushed tag that starts with `v`, and only
     in the repository named above.
5. Approve the publish job when GitHub asks (if the environment has a required reviewer), then check
   `https://pypi.org/project/flowcellar/` and install the release into an empty environment.
6. Set the version to the next development version (for example `0.2.0.dev0`).

`build.yml` can also be started by hand ("Run workflow"). It then builds and checks everything and keeps the
two files as an artifact of the run, but does not publish.

## What is in the two files

- **The wheel** (`flowcellar-<version>-py3-none-any.whl`, about 0.5 MB): the package, the built viewer
  (`flowcellar/viewer_dist/`), the catalogue JSON, the bundle schema, the default configuration, `LICENSE`
  and `THIRD_PARTY_NOTICES.md`.
- **The source distribution** (`flowcellar-<version>.tar.gz`, about 3 MB): the same package with the built viewer,
  and what is needed to build and test it again: the viewer's source and lock file, the catalogue source, the
  tools, the tests and these pages with the pictures of the README (`docs/img/`, at most 3 MB). A wheel can be built
  from it without Node.

`hatch_build.py` stops any build that would lack the built viewer or the catalogue JSON, so a release
cannot be made from a checkout in which the viewer was never built. `design/`, `prototype/` and `local/`
are in neither file; `tools/check_dist.py` fails when a file from them, the name of a storage mount or an
id that looks like a real flow cell appears.

## Not set up yet

- A conda package (bioconda) and a container image.
- Signed release notes on GitHub. (PyPI attestations are made by the publish action by default.)
- Versions taken from git tags: the version is a plain string so that it works without git history.
