# The FlowCellar viewer

Everything a person sees: Home, Runs, Projects, a project, a run, Compare and the history of a flow cell. One
script and one stylesheet, built here and shipped inside the Python package; Python assembles the pages and
needs no Node.

```
viewer/
  src/            the viewer: Preact, signals, and the project's own canvas chart kit (src/charts/kit.js)
    core/         data loading, address and theme state, formats, labels, derivations from a run's data
    home/ hub/ projects/ run/ compare/ flowcell/ shell/    one folder per screen
    gen/          written by the build from the catalogue and the logo (not kept in git)
  styles/         tokens.css (the visual system, unchanged), components.css, app.css
  brand/          logo and favicon
  build.mjs       viewer/ + the package catalogue -> src/flowcellar/viewer_dist/
  verify/         browser checks (headless Firefox and Chromium through Playwright)
```

## Build

```
npm ci                    # preact, @preact/signals, esbuild, playwright -> ./node_modules only
node build.mjs            # -> ../src/flowcellar/viewer_dist/  (about 0.2 s);  --dev keeps the script readable
```

The build writes `viewer.core.js` (one classic script), `viewer.css`, `favicon.svg`, `template.html` and
`manifest.json` (size, SHA-256 and the policy hash of each file). The same sources give the same bytes: nothing in
the output depends on the time or the machine. Labels and "How to read" texts come from the catalogue inside the
package (`src/flowcellar/catalogue/metrics_nanopore.json`), which is the only home of a label; a number the
catalogue has no label for is named from its id in words (the id stays in the row's title and in its own column). The build fails if the script would end its own element, names another
origin, could ask the network for anything, or builds markup from text.

`python tools/build_viewer.py --check` says whether `viewer_dist/` is up to date.

## How a page gets its data

The page carries a boot block, `{"mode": "report" | "hub", ...}`; `src/core/data.js` reads it once.

| Page | Written by | Data |
|---|---|---|
| one-run report, `<run>.flowcellar.html` | `flowcellar.report.api.write_report` | the index (one row) and the run, base64 of gzip, inside the file |
| hub, `hub/index.html` | `flowcellar.report.api.write_hub` | none in the page: `index.js` at start-up, `data/<run key>.js` when a screen first needs a run |
| the same hub over `flowcellar serve` | `flowcellar.report.api.serve_file` | the same three kinds of file |

Nothing is fetched (no `fetch`, XHR, `import()` or worker): sibling files are loaded as scripts that call
`__FC_INDEX__` or `__FCB_SIDECAR__`, which is the one way that also works from `file://`. Every page carries a
Content-Security-Policy with `default-src 'none'`, `connect-src 'none'` and the hash of the one script; a hub adds
`'self'` for its sibling files. A run's checksum is checked before it is drawn.

Text from a run folder (run, sample and project names, messages) reaches the DOM as text nodes and attribute
values only; the viewer never builds markup from data.

## Paper

"Save as PDF" (on a run, a project, an analysis, Compare and the history of a flow cell) and the browser's own Print command
lay the page out for paper before the browser takes its picture of it. `src/core/print.js` answers `beforeprint`:
`<html>` gets the class `fc-print` and the theme "print", every screen is drawn again, every chart paints itself
again at the width it has on the sheet (three device pixels to a CSS pixel), and every heading is tied to the block
under it; `afterprint` undoes all of it. The layout is in `styles/app.css` under "Paper": it hangs on that class,
never on the width of the window, so a chart is measured at the width it is printed at. `src/shell/paper.jsx` holds
the control, the head strip and the line that says which selection is shown. No PDF library is involved and nothing
is rendered anywhere but in the browser.

## Checks

```
node verify/fixture.mjs             # the synthetic library: file://, HTTP, reports, empty library, no index (about 4 min)
node verify/fixture.mjs --quick     # one theme, one width
node verify/check.mjs --hub DIR_OR_URL --out DIR [--report FILE]     # the same checks on any hub
node verify/smoke.mjs FILE_OR_URL --hash '#/runs' --out shot.png     # one page, one screenshot
node verify/query_vectors.mjs       # the search grammar against tests/vectors/query (shared with Python)
node verify/print.mjs --fixture FILE.json --out DIR                  # "Save as PDF" in Chromium and Firefox (about 4 min)
node verify/print.mjs --out DIR --page NAME=KIND=URL                 # the same for any page of any hub or report
node verify/firefox_pdf.mjs URL OUT.pdf [a4|letter]                  # one PDF from Firefox's own print pipeline
```

Every check runs in Firefox; `FLOWCELLAR_BROWSER_ENGINE=chromium` runs it in Chromium instead (`verify/print.mjs`
always uses both).

The checks open every screen in both themes at 1,366, 1,280 and 820 px and fail on any page error, console
error, policy violation or request to another origin. At 2,560 px they hold that the application bar and the
content of every screen start and end on the same two lines (no screen is a centred column) and that the width is
used for columns; and they move the pointer across every chart of Home, a run and Compare and fail if any element
of the page changes its size or its place (a read-out has a box of fixed height: only the text in it changes). They need a Playwright Firefox on the machine (and a Playwright Chromium for the checks in Chromium):
`PLAYWRIGHT_BROWSERS_PATH` (or `FLOWCELLAR_PW_BROWSERS`) names their folder, or install them into `viewer/pw-browsers`
with `PLAYWRIGHT_BROWSERS_PATH=$PWD/pw-browsers npx playwright install firefox chromium`. The fixture is made by
`tests/viewer/viewer_make_fixture.py` from invented runs (`tests/viewer/viewer_synth.py`).

The Python side is tested by `pytest tests/viewer`; `FLOWCELLAR_BROWSER=1` adds the browser checks as one test.

Not checked: WebKit (Safari), Chrome and Edge on Windows and macOS, a screen reader, forced colours, and the print
dialog itself under a person's hand (the PDFs are made through the print pipelines of the two browsers, which fire the
same events).
