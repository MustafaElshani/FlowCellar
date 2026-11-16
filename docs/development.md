# Development

## Set up

Python 3.11 or newer and Node 22 or newer (Node only builds the viewer; nothing needs it at run time).

```
python3 -m venv .venv
(cd viewer && npm ci && node build.mjs)          # builds src/flowcellar/viewer_dist/ (not kept in git)
.venv/bin/pip install -e ".[dev,bam]"
.venv/bin/python -m pytest                        # invented data only; about four minutes
.venv/bin/ruff check src tests tools hatch_build.py
```

## Where things are

| Path | What |
|---|---|
| `src/flowcellar/` | the package: `cli/` (one file per command), `engine/`, `sources/` and `modules/` (run health), `readdata/` (read data), `sc/` (single-cell and spatial analyses), `library/`, `discovery/` and `config/` (the library and the scan), `report/` and `serve/` (pages), `model/`, `aggregates/`, `catalogue/`, `util/` (shared), `testing/` (generators of invented data) |
| `src/flowcellar/viewer_dist/` | the built viewer; written only by `viewer/build.mjs`; shipped in the wheel and the source distribution |
| `viewer/` | the viewer's source: JavaScript and JSX (Preact, signals, the project's own canvas chart kit), styles, logo, build script, browser checks |
| `catalogue/metrics_nanopore.yaml` | the metric catalogue: every metric with its label, sentences and limits |
| `catalogue/metrics_singlecell_spatial.yaml` | the catalogue of the single-cell and spatial add-on |
| `tests/` | pytest; `tests/vectors/` holds cases that Python and the viewer both run |
| `tools/` | the scripts below |
| `docs/` | these pages; `commands.md` and `metrics.md` are generated |

## Tools

| Command | What it does |
|---|---|
| `node viewer/build.mjs` (or `python tools/build_viewer.py`) | builds the viewer: one script, one stylesheet, a page template and a manifest with hashes. Two builds of the same source give identical bytes. `python tools/build_viewer.py --check` says whether the built viewer is up to date |
| `python tools/build_catalogue.py` | writes the JSON the package loads from the catalogue source; `--check` fails when it is out of date. PyYAML is a development dependency only |
| `python tools/build_sc_catalogue.py` | the same for the catalogue of the single-cell and spatial add-on (`catalogue/metrics_singlecell_spatial.yaml` to `src/flowcellar/sc/catalogue.json`); `--check` fails when it is out of date |
| `python tools/build_docs.py` | writes `docs/commands.md` from the `--help` output and `docs/metrics.md` from the catalogue; `--check` fails when either is out of date |
| `python tools/make_fixtures.py OUT [--tree]` | writes fake MinKNOW run folders with invented names |
| `python tools/make_demo.py OUT [--scale 0.01]` | writes the example of the README: an invented laboratory with run folders, read data and analysis folders (the same as `python -m flowcellar.testing.demo`) |
| `python tools/make_readme_pictures.py --library DIR` | makes the pictures of the README again from a library of the example (below) |
| `python tools/licence_check.py` | the licences of what the viewer bundles and of what is installed beside FlowCellar, against an allow-list; the bundled packages must be named in `THIRD_PARTY_NOTICES.md` |
| `python tools/check_dist.py dist/*` | what a built wheel and source distribution must hold, must not hold, and must not say |
| `python tools/smoke_test.py --installed` | tries an installed FlowCellar from the outside on invented data, with the network refused |

After changing a help text, run `python tools/build_docs.py`. After changing the catalogue, run
`python tools/build_catalogue.py`, `node viewer/build.mjs` and `python tools/build_docs.py`.

## The example data and the pictures of the README

`flowcellar.testing.demo` writes the example a reader tries the product on, and every picture of the README shows
it. Nothing in it comes from a laboratory: runs, flow cell ids (`TST…`), samples and reads are drawn from a fixed
seed, so equal arguments give the same example.

- `demo.py` holds the runs (`demo_runs()`), one table of reads per run (`reads_of`) from which both the MinKNOW
  folder and the read-data files of the run are written, and the configuration of the example library
  (`config_text`). The MinKNOW folder starts from `fixtures.make_run`; BAM files are written as SAM text and
  converted by the samtools inside pysam. `demo_spatial.py` holds the Visium and Visium HD slides and their two
  pipelines.
- `--scale` sets the size of the runs as a share of full-size runs (default 0.01). Below 1 the configuration
  scales the built-in output limits by the same factor. The output limits of the engine are absolute numbers, so
  an example of small runs needs this; the block is not written at `--scale 1`.
- `tests/integration/test_example.py` writes it at a five-hundredth of full size, scans it and holds what the
  README says about it (15 runs, the verdicts, the sample on two flow cells, the two basecalling passes, the
  analyses and their flow cells, the paths the README types).

The pictures in `docs/img/` are made by a script, so that they can be made again after a change of the viewer:

```
python -m flowcellar.testing.demo /tmp/fc-example --scale 1 --jobs 22  # full size: 135 GB, 20 minutes on 22 cores
flowcellar scan /tmp/fc-example/minknow /tmp/fc-example/projects -L /tmp/fc-example/library     # 10 minutes
python tools/make_readme_pictures.py --library /tmp/fc-example/library
```

The pictures in the repository show the example at full size, so that their numbers are those of real runs and the
built-in limits judge them. `--scale 0.1` (12 GB, four minutes to write and four to scan) gives the same pages
with numbers a tenth as large and the note that limits of the library were used; it is enough to check a change.

`tools/make_readme_pictures.py` writes the illustration (`tools/readme_illustration.py`: one drawing in a light and
a dark colouring), takes the screenshots with `tools/readme_shots.mjs` (Playwright from `viewer/node_modules`,
Chromium or Firefox from `viewer/pw-browsers`), puts the three sheets of single charts together with ImageMagick,
turns the first two sheets of a run page's PDF into a picture (Chromium's own PDF, `pdftoppm`), reduces every PNG
to a palette, and ends with code 1 when a picture is over 350 kB or the folder over 3 MB. After
a change of the viewer, write the pages again first (`flowcellar hub /tmp/fc-example/library`). Take pictures from
the example only, never from a library of real runs. `tests/release/test_release_files.py` fails when the README
shows a picture that is missing or too large, or when `docs/img/` holds one the README does not show.

On PyPI the README is the long description of the package, and a path inside the repository leads nowhere there.
`hatch_build.py` therefore writes the description with every such path as an address (pictures at
`raw.githubusercontent.com`, pages at `github.com/…/blob`, at the tag of a release or at `main` for a development
version). `README.md` itself keeps relative paths, which work on GitHub and in a checkout. `tools/check_dist.py`
fails a built package whose description still holds a relative path.

## Tests

- No test reads real data. `flowcellar.testing.fixtures` writes run folders, `flowcellar.testing.readdata`
  small BAMs and tables, `flowcellar.testing.trees` a tree of three roots, `flowcellar.testing.demo` the example
  of the README; all names are invented and every flow cell id starts with `TST`.
- `tests/test_repository_rules.py` fails when a tracked file holds a storage path of a laboratory, an id
  that looks like a real flow cell, or an import of a package that must not be a run-time dependency.
- `tests/integration/` starts the real console program: a scan of three roots, the hub, a second scan that
  opens no data file, `serve`, `run`.
- Browser checks of the viewer (headless Firefox, about fifteen minutes for all of them) run on request:
  `FLOWCELLAR_BROWSER=1 .venv/bin/python -m pytest tests/viewer`, with a Playwright
  Firefox under `viewer/pw-browsers` (`cd viewer && PLAYWRIGHT_BROWSERS_PATH=pw-browsers npx playwright
  install firefox`).
- A built wheel, installed away from the source tree:
  `FLOWCELLAR_WHEEL=dist/flowcellar-<version>-py3-none-any.whl .venv/bin/python -m pytest tests/integration/test_wheel.py`.

## Rules of the code

- Run-time dependencies are `numpy`, `pyarrow` and `click`; `pysam` in the extra `bam`; `jsonschema` in the
  extra `validate`; `h5py` in the extra `sc`. Nothing else is imported at run time, no other program is started, and nothing
  contacts the network.
- The viewer never fetches: a page reads its data from script elements inside the file or from script
  files beside it, under a Content-Security-Policy with `default-src 'none'`.
- Nothing that names a laboratory's folders, samples or projects belongs in the package: such patterns
  live in a library's `flowcellar.toml`.
- Every label and every sentence a person reads comes from the catalogue.

## Continuous integration

`.github/workflows/ci.yml` runs on every push to `main` and on every pull request: lint, the generated
files (`--check` of the build tools), the licence check, the tests on every supported Python
version with the viewer built first, and a package job that builds the wheel and the source distribution,
checks them with `tools/check_dist.py` and runs `tools/smoke_test.py` on the installed wheel. It never
publishes. [Releasing](releasing.md) describes the one workflow that does.
