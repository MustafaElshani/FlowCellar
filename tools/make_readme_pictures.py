"""The pictures of the README, made again from a library of the example data.

    python -m flowcellar.testing.demo /tmp/fc-example --scale 0.1             # the example data (docs/development.md)
    flowcellar scan /tmp/fc-example/minknow /tmp/fc-example/projects -L /tmp/fc-example/library
    python tools/make_readme_pictures.py --library /tmp/fc-example/library     # writes docs/img/

What it writes into ``docs/img/``:

* ``overview-light.svg`` and ``overview-dark.svg``: the illustration (``tools/readme_illustration.py``);
* screenshots of the hub of that library: Home, a run page, Compare, the history of a flow cell, an analysis and the
  notice of a sample that went on two flow cells, taken by
  ``tools/readme_shots.mjs`` (Playwright from ``viewer/node_modules``, a browser from ``viewer/pw-browsers``);
* three sheets of single charts, grouped as a run page groups them and named, put together with ImageMagick;
* the first two sheets of the PDF that "Save as PDF" gives for a run page (Chromium's own PDF, turned into a picture
  with ``pdftoppm``).

Every PNG is reduced to at most 256 colours and checked against the budget (350 kB a picture, 3 MB for the folder);
the script ends with code 1 when a picture is over it. ``--only home,run`` makes some of them; ``--raw DIR`` keeps the
screenshots as taken. Needs Node, ImageMagick (``convert``, ``identify``), ``pdftoppm`` (poppler) and the example
library; nothing here runs at install time or in the tests.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import readme_illustration  # noqa: E402

MAX_PNG = 350_000
MAX_ALL = 3_000_000
INK, PAPER, MUTED = "#0a213c", "#ffffff", "#606468"
# the runs of the example the pictures show, by flow cell (flowcellar.testing.demo): a cDNA pool on one flow cell with
# every kind of read data; the pool that went on two flow cells; a genomic run; the flow cell that was used twice
RUN, PAIR, PAIR2, GENOME, REUSED, LOW = "TST00107", "TST00105", "TST00106", "TST00202", "TST00201", "TST00104"
HD = "{analysis:visium_hd_b:spaceranger}"

#: whole regions of pages: name -> (route, from, to, window width in CSS pixels, device pixels per CSS pixel, only as wide as "from")
PAGES = {
    "home": ("#/home", None, ".fc-home__block--projects", 1120, 1.6, False),
    "run": ("#/run/{run:%s}" % RUN, None, "#p-yield-output_over_time", 1120, 1.6, False),
    "compare": ("#/compare?items={run:%s},{run:%s},{run:%s},{run:%s}" % (RUN, PAIR, PAIR2, LOW), None, ".fc-difftable", 1440, 1.28, False),
    "flowcell": ("#/flowcell/%s" % REUSED, None, ".fc-tracks", 1120, 1.6, False),
    "analysis": ("#/analysis/" + HD, None, "#s-map", 1120, 1.6, False),
    # the sample that went on two flow cells: what its page says about a step that ran on the merged reads
    "merged": ("#/run/{run:%s}" % PAIR, "#s-pychopper", "#s-pychopper .fc-tile:nth-of-type(3)", 1120, 1.6, True),
}
#: CSS pixels kept around a region; less where the next thing on the page begins right below it
PADS = {"run": 2, "merged": 4}
#: the PDF of a page (what "Save as PDF" gives): name -> (route, sheets shown side by side)
PDFS = {"pdf": ("#/run/{run:%s}" % RUN, (1, 2))}
#: what is pressed on a page before its picture is taken: the run page counts by bases, so that the output track
#: shows the band of what is typical in the library
CLICKS = {"run": ['.fc-controls button:has-text("Bases")']}
#: single charts: name -> (route, selector, nth)
CHARTS = {
    "t-timeline": ("#/run/{run:%s}" % RUN, "#s-run_timeline .fc-tracks", 0),
    "t-map": ("#/run/{run:%s}" % RUN, "#p-fc-channel_map", 0),
    "t-firstscan": ("#/run/{run:%s}" % RUN, "#p-fc-first_scan", 0),
    "t-estlen": ("#/run/{run:%s}" % RUN, "#p-len-estimated_length", 0),
    "t-ended": ("#/run/{run:%s}" % RUN, "#p-dyn-why_reads_ended-bars", 0),
    "t-fate": ("#/run/{run:%s}" % RUN, "#p-read-fate", 0),
    "t-fate-samples": ("#/run/{run:%s}" % RUN, "#p-read-fate-samples", 0),
    "t-per-sample": ("#/run/{run:%s}" % RUN, "#p-yield-per_sample", 0),
    "t-len": ("#/run/{run:%s}" % RUN, "#p-len", 0),
    "t-qual": ("#/run/{run:%s}" % RUN, "#p-qual", 0),
    "t-lenqual": ("#/run/{run:%s}" % RUN, "#p-lenqual", 0),
    "t-aln-records": ("#/run/{run:%s}" % RUN, "#p-aln-records", 0),
    "t-aln-bases": ("#/run/{run:%s}" % RUN, "#p-aln-bases", 0),
    "t-aln-contigs": ("#/run/{run:%s}" % RUN, "#p-aln-contigs", 0),
    "t-aln-lengths": ("#/run/{run:%s}" % RUN, "#p-aln-lengths", 0),
    "t-ident": ("#/run/{run:%s}" % RUN, "#p-ident", 0),
    "t-calibration": ("#/run/{run:%s}" % RUN, "#p-aln-calibration", 0),
    "t-depth-cdna": ("#/run/{run:%s}" % RUN, "#p-aln-depth", 0),
    "t-depth-gdna": ("#/run/{run:%s}" % GENOME, "#p-aln-depth", 0),
    "t-time-length": ("#/run/{run:%s}" % RUN, "#p-reads-length_over_time", 0),
    "t-time-quality": ("#/run/{run:%s}" % RUN, "#p-reads-quality_over_time", 0),
    "t-time-identity": ("#/run/{run:%s}" % RUN, "#p-reads-identity_over_time", 0),
    "t-time-speed": ("#/run/{run:%s}" % RUN, "#p-reads-speed_over_time", 0),
    "t-polya": ("#/run/{run:%s}" % RUN, "#p-rna-polya", 0),
    "t-body": ("#/run/{run:%s}" % RUN, "#p-rna-body_coverage", 0),
    "t-saturation": ("#/run/{run:%s}" % RUN, "#p-rna-saturation", 0),
    "t-per-transcript": ("#/run/{run:%s}" % RUN, "#p-rna-reads_per_transcript", 0),
    "t-types": ("#/run/{run:%s}" % RUN, "#p-rna-types", 0),
    "t-tail-types": ("#/run/{run:%s}" % RUN, "#p-rna-tail", 0),
    "t-genes": ("#/run/{run:%s}" % RUN, "#p-rna-genes", 0),
    "t-topgenes": ("#/run/{run:%s}" % RUN, "#p-rna-topgenes", 0),
    "t-pych-class": ("#/run/{run:%s}" % RUN, "#p-pych-class", 0),
    "t-pych-config": ("#/run/{run:%s}" % RUN, "#p-pych-config", 0),
    "t-pych-trimmed": ("#/run/{run:%s}" % RUN, "#p-pych-trimmed", 0),
    "t-pych-tune": ("#/run/{run:%s}" % RUN, "#p-pych-tune", 0),
    "t-sp-map": ("#/analysis/" + HD, "#s-map figure", 0),
    "t-sp-units": ("#/analysis/" + HD, "#s-units figure", 0),
    "t-sp-genes": ("#/analysis/" + HD, "#s-units figure", 1),
    "t-sp-saturation": ("#/analysis/" + HD, "#s-saturation figure", 0),
}
#: the sheets: name -> rows of (group, [(chart, caption)])
SHEETS = {
    "figures-run": [
        ("Run health · from MinKNOW's own files, no basecalling needed",
         [("t-timeline", "Run timeline: five tracks on one time axis"), ("t-map", "Flow cell map"), ("t-estlen", "Estimated read length")]),
        ("", [("t-firstscan", "The first pore scan"), ("t-ended", "Why reads ended"), ("t-fate", "Read fate: where every read went"),
              ("t-per-sample", "Reads per sample")]),
        ("Read length and quality · from basecalled reads",
         [("t-len", "Read length"), ("t-qual", "Read quality"), ("t-lenqual", "Length against quality"), ("t-fate-samples", "Read fate per sample")]),
    ],
    "figures-alignment": [
        ("Alignment · from aligned BAM files",
         [("t-aln-records", "Summary of read alignment"), ("t-aln-bases", "What the aligned bases are"), ("t-aln-contigs", "Reads per chromosome")]),
        ("", [("t-aln-lengths", "Read length and aligned length"), ("t-ident", "Identity with the reference"),
              ("t-calibration", "Basecaller quality against measured quality")]),
        ("", [("t-depth-gdna", "Depth along the genome (a genomic DNA run)")]),
        ("Reads over run time · from the start time every read carries",
         [("t-time-length", "Read length over the run"), ("t-time-quality", "Read quality over the run")]),
        ("", [("t-time-identity", "Identity over the run"), ("t-time-speed", "Speed through the pore over the run")]),
    ],
    "figures-rna-spatial": [
        ("Transcripts · from an alignment to transcripts and the basecaller's poly(A) estimate",
         [("t-polya", "Poly(A) tail length"), ("t-body", "Coverage along the transcript"), ("t-saturation", "Transcripts found against reads"),
          ("t-per-transcript", "Reads per transcript")]),
        ("RNA composition · from annotated read tables",
         [("t-types", "Reads by transcript type"), ("t-tail-types", "Poly(A) tail by transcript type"), ("t-genes", "How reads spread over genes"),
          ("t-topgenes", "Genes with the most reads")]),
        ("Pychopper · from its statistics file",
         [("t-pych-class", "What pychopper made of the reads"), ("t-pych-config", "Primer arrangements found"),
          ("t-pych-trimmed", "Bases cut off each read"), ("t-pych-tune", "The cut-off pychopper chose")]),
        ("Single-cell and spatial · from Space Ranger and wf-single-cell folders",
         [("t-sp-map", "Spatial map: molecules and tissue"), ("t-sp-units", "UMIs per bin"), ("t-sp-saturation", "Saturation against reads")]),
    ],
}
SHEET_W, GAP, CAPTION_H, LABEL_H, MARGIN = 1800, 20, 46, 64, 28


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise SystemExit(f"{' '.join(cmd[:3])} ... failed:\n{r.stdout}\n{r.stderr}")
    return r


def size_of(path: str) -> tuple[int, int]:
    w, h = run(["identify", "-format", "%w %h", path]).stdout.split()
    return int(w), int(h)


def reduce_png(src: str, dest: str, *, limit: int = MAX_PNG) -> int:
    """Write ``src`` as a palette PNG: as many colours as fit under the limit, never fewer than 48."""
    for colours in (256, 160, 96, 64, 48):
        run(["convert", src, "-strip", "-dither", "None", "-colors", str(colours), "-define", "png:compression-level=9",
             "-define", "png:exclude-chunk=all", "PNG8:" + dest])
        if os.path.getsize(dest) <= limit:
            break
    return os.path.getsize(dest)


def take(hub: str, raw: str, wanted: set[str]) -> None:
    jobs: dict[tuple, dict] = {}
    for name, (route, start, end, width, scale, tight) in PAGES.items():
        if name in wanted:
            job = jobs.setdefault((route, width, scale, name in CLICKS), {"route": route, "width": width, "height": 900, "scale": scale,
                                                                          "scheme": "light", "click": CLICKS.get(name, []), "shots": []})
            job["shots"].append({"name": name, "from": start, "to": end, "tight": tight, "pad": PADS.get(name, 12)})
    for name, (route, _sheets) in PDFS.items():
        if name in wanted:
            jobs[(route, "pdf")] = {"route": route, "width": 1200, "height": 900, "scale": 1, "scheme": "light", "shots": [{"name": name, "pdf": True}]}
    need = {c for sheet, rows in SHEETS.items() if sheet in wanted for _, row in rows for c, _ in row}
    for name, (route, selector, nth) in CHARTS.items():
        if name in need:
            job = jobs.setdefault((route, 1400, 1.25), {"route": route, "width": 1400, "height": 900, "scale": 1.25, "scheme": "light", "shots": []})
            job["shots"].append({"name": name, "selector": selector, "nth": nth})
    if not jobs:
        return
    spec = os.path.join(raw, "spec.json")
    with open(spec, "w", encoding="utf-8") as fh:
        json.dump(list(jobs.values()), fh, indent=1)
    r = subprocess.run(["node", os.path.join(REPO, "tools", "readme_shots.mjs"), "--hub", hub, "--out", raw, "--spec", spec],
                       capture_output=True, text=True)
    try:
        report = json.loads(r.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        raise SystemExit(f"the screenshots could not be taken:\n{r.stdout}\n{r.stderr}") from None
    print(f"  browser: {report['browser']}; {len(report['shots'])} screenshots")
    for e in report["errors"]:
        print(f"  problem: {e}")
    if r.returncode != 0:
        raise SystemExit("some screenshots are missing (see above)")


def sheet(name: str, rows: list, raw: str, work: str) -> str:
    """Rows of charts, every row as wide as the sheet, a caption under each chart and a label above each group."""
    parts = []
    for k, (group, charts) in enumerate(rows):
        if group:
            label = os.path.join(work, f"{name}-label{k}.png")
            run(["convert", "-size", f"{SHEET_W}x{LABEL_H}", f"xc:{PAPER}", "-font", "DejaVu-Sans-Bold", "-pointsize", "30", "-fill", INK,
                 "-gravity", "SouthWest", "-annotate", "+0+14", group, label])
            parts.append(label)
        sizes = [size_of(os.path.join(raw, c + ".png")) for c, _ in charts]
        inner = SHEET_W - GAP * (len(charts) - 1)
        height = min(int(inner / sum(w / h for w, h in sizes)), 440)
        tiles = []
        for (chart, caption), (w, h) in zip(charts, sizes):
            tile = os.path.join(work, f"{name}-{chart}.png")
            width = int(round(w * height / h))
            run(["convert", os.path.join(raw, chart + ".png"), "-filter", "Lanczos", "-resize", f"{width}x{height}!", "-background", PAPER,
                 "-gravity", "North", "-extent", f"{width}x{height + CAPTION_H}", "-font", "DejaVu-Sans", "-pointsize", "23", "-fill", MUTED,
                 "-gravity", "SouthWest", "-annotate", "+2+8", caption, tile])
            tiles.append(tile)
        row = os.path.join(work, f"{name}-row{k}.png")
        cmd = ["convert"]
        for i, t in enumerate(tiles):
            cmd += [t] + (["-size", f"{GAP}x1", f"xc:{PAPER}"] if i < len(tiles) - 1 else [])
        run(cmd + ["-background", PAPER, "-gravity", "North", "+append", "-gravity", "West", "-extent", f"{SHEET_W}x{height + CAPTION_H}", row])
        parts.append(row)
    out = os.path.join(work, name + ".png")
    run(["convert", *parts, "-background", PAPER, "-gravity", "West", "-append", "-bordercolor", PAPER, "-border", str(MARGIN), out])
    return out


def sheets_picture(name: str, sheets: tuple[int, ...], raw: str, work: str) -> str:
    """Sheets of a PDF side by side, each with a thin frame, as one picture."""
    pdf = os.path.join(raw, name + ".pdf")
    pages = []
    for n in sheets:
        prefix = os.path.join(work, f"{name}-sheet{n}")
        run(["pdftoppm", "-r", "100", "-png", "-f", str(n), "-l", str(n), "-singlefile", pdf, prefix])
        framed = prefix + "-framed.png"
        run(["convert", prefix + ".png", "-bordercolor", "#cac7c2", "-border", "2", framed])
        pages.append(framed)
    out = os.path.join(work, name + ".png")
    cmd = ["convert"]
    for i, f in enumerate(pages):
        cmd += [f] + (["-size", "36x1", "xc:#f1efeb"] if i < len(pages) - 1 else [])
    run(cmd + ["-background", "#f1efeb", "-gravity", "North", "+append", "-bordercolor", "#f1efeb", "-border", "28", out])
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Make the pictures of the README from a library of the example data.")
    p.add_argument("--library", help="the library that `flowcellar scan` made from the example data (its hub/ is what is photographed)")
    p.add_argument("--out", default=os.path.join(REPO, "docs", "img"))
    p.add_argument("--only", default="", help="names of pictures to make, separated by commas (default: all)")
    p.add_argument("--raw", default=None, help="keep the screenshots as taken in this folder")
    a = p.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    for path in readme_illustration.write(a.out):
        print(f"  {os.path.relpath(path, REPO)}  {os.path.getsize(path) / 1000:.0f} kB")
    if not a.library:
        print("no --library given: only the illustration was written")
        return 0
    hub = os.path.join(a.library, "hub")
    if not os.path.isfile(os.path.join(hub, "index.html")):
        raise SystemExit(f"{hub}/index.html is missing: scan the example data first, or run `flowcellar hub {a.library}`")
    for tool in ("node", "convert", "identify", "pdftoppm"):
        if not shutil.which(tool):
            raise SystemExit(f"{tool} is needed and was not found")
    names = set(PAGES) | set(SHEETS) | set(PDFS)
    wanted = {n.strip() for n in a.only.split(",") if n.strip()} or names
    if wanted - names:
        raise SystemExit(f"no such picture: {', '.join(sorted(wanted - names))}; there are {', '.join(sorted(names))}")
    work = tempfile.mkdtemp(prefix="fc-readme-")
    raw = a.raw or os.path.join(work, "raw")
    os.makedirs(raw, exist_ok=True)
    try:
        take(hub, raw, wanted)
        over = []
        for name in sorted(wanted):
            src = (os.path.join(raw, name + ".png") if name in PAGES else sheets_picture(name, PDFS[name][1], raw, work) if name in PDFS
                   else sheet(name, SHEETS[name], raw, work))
            dest = os.path.join(a.out, name + ".png")
            n = reduce_png(src, dest)
            w, h = size_of(dest)
            print(f"  {os.path.relpath(dest, REPO)}  {w} x {h}  {n / 1000:.0f} kB")
            if n > MAX_PNG:
                over.append(name)
    finally:
        if not a.raw:
            shutil.rmtree(work, ignore_errors=True)
    total = sum(os.path.getsize(os.path.join(a.out, f)) for f in os.listdir(a.out))
    print(f"  {os.path.relpath(a.out, REPO)}: {total / 1e6:.2f} MB in {len(os.listdir(a.out))} files")
    if over or total > MAX_ALL:
        print(f"over the budget: {', '.join(over) or 'the folder as a whole'} (a picture: {MAX_PNG // 1000} kB, the folder: {MAX_ALL // 1_000_000} MB)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
