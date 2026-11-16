"""The illustration of the README: what goes in, what the library holds, what one looks at.

    python tools/readme_illustration.py            # writes docs/img/overview-light.svg and docs/img/overview-dark.svg

One drawing in two colourings, because a page on GitHub is light or dark and a picture cannot ask which: the README
shows the one that fits with a ``<picture>`` element. The text is real text (system fonts), 13 px and larger at the
drawn width of 1,000 px. ``tools/make_readme_pictures.py`` calls this before it takes the screenshots.
"""
from __future__ import annotations

import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FONT = "system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"
MONO = "ui-monospace, SFMono-Regular, Menlo, Consolas, 'DejaVu Sans Mono', monospace"
# the colours of the brand (viewer/brand, design tokens): ink and paper change places, the lamp gets lighter on dark
THEMES = {
    "light": {"ink": "#0a213c", "muted": "#606468", "panel": "#f7f6f3", "card": "#ffffff", "line": "#cac7c2", "lamp": "#2a78d6",
              "lamp_soft": "#eef5fe", "on_lamp": "#ffffff", "pass": "#1e7a4c"},
    "dark": {"ink": "#f1efeb", "muted": "#a7a49f", "panel": "#171f28", "card": "#202730", "line": "#464b51", "lamp": "#6da7ec",
             "lamp_soft": "#0a2e57", "on_lamp": "#07162a", "pass": "#6fcf97"},
}
W, H = 1000, 600

INPUTS = (
    ("folder", "MinKNOW run folders", "reports, pore scans, sequencing summary"),
    ("file", "Basecalled BAM files", "length, quality, samples, poly(A)"),
    ("file", "Aligned BAM files", "to a genome, to transcripts"),
    ("table", "pychopper statistics", "full-length cDNA, rescued, unusable"),
    ("table", "Annotated read tables", "transcript type, poly(A) tail per read"),
    ("table", "Quantification metadata", "salmon, oarfish, kallisto"),
    ("folder", "wf-single-cell output folders", "single-cell and spatial, ONT workflow"),
    ("folder", "Space Ranger output folders", "Visium and Visium HD"),
)
LIBRARY = (
    ("bundles/", "one small file per run, and one per analysis: what was measured"),
    ("status/", "each verdict with its reasons, in plain sentences"),
    ("index/", "one row per run, for the tables and the search"),
    ("hub/", "the pages: open them from a file or the built-in server"),
    ("flowcellar.toml", "your file patterns and limits: the only file you edit"),
)
PAGES = (
    ("Home", "what is there, what needs a look"),
    ("Runs · Projects · Analyses", "tables with search and filters"),
    ("A page per run", "run health and read data, every chart"),
    ("Compare", "up to eight runs side by side"),
    ("Flow cell history", "every loading of one flow cell"),
    ("One-file report · PDF", "to mail, archive or print"),
)


def _icon(kind: str, x: float, y: float, c: dict) -> str:
    """A 22 px pictogram: a folder, a file or a table."""
    s = f'<g transform="translate({x} {y})" fill="none" stroke="{c["ink"]}" stroke-width="1.5" stroke-linejoin="round" stroke-linecap="round">'
    if kind == "folder":
        s += '<path d="M1.5 5.5h6l2 2.5h11v10.5h-19z"/>'
    elif kind == "file":
        s += '<path d="M5 1.5h8l4.5 4.5v14.5h-12.500z"/><path d="M13 1.500v4.500h4.500"/><path d="M8 11h6.500M8 14.500h6.500"/>'
    else:
        s += '<path d="M2 3.500h18v15h-18z"/><path d="M2 8.500h18M2 13.500h18M8 3.500v15"/>'
    return s + "</g>"


def _text(x: float, y: float, text: str, c: dict, *, size: int = 13, weight: int = 400, fill: str = "ink", anchor: str = "start", font: str = FONT,
          spacing: str | None = None) -> str:
    extra = f' letter-spacing="{spacing}"' if spacing else ""
    return (f'<text x="{x}" y="{y}" font-family="{font}" font-size="{size}" font-weight="{weight}" fill="{c[fill]}" text-anchor="{anchor}"{extra}>'
            f"{text.replace('&', '&amp;').replace('<', '&lt;')}</text>")


def _arrow(x1: float, y: float, x2: float, c: dict) -> str:
    return (f'<path d="M{x1} {y}H{x2 - 7}" stroke="{c["lamp"]}" stroke-width="2.5" fill="none"/>'
            f'<path d="M{x2 - 9} {y - 6}L{x2} {y}L{x2 - 9} {y + 6}z" fill="{c["lamp"]}"/>')


def illustration(theme: str) -> str:
    c = THEMES[theme]
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" role="img" aria-labelledby="t d">',
           "<title id=\"t\">What FlowCellar reads, what its library holds and what you look at</title>",
           "<desc id=\"d\">On the left the files a laboratory already has: MinKNOW run folders, which are all that is needed, and seven "
           "further kinds that each add sections. In the middle the command flowcellar scan, which only reads them. Then the library, "
           "one folder with a small file per run, the verdicts, an index, the pages and the configuration. On the right the pages: Home, "
           "Runs, Projects, Analyses, a page per run, Compare, flow cell history, one-file reports and PDF.</desc>"]
    cap = dict(size=11, weight=600, fill="muted", spacing="0.06em")
    # ------------------------------------------------------------------ what you have
    ax, aw, top, row = 0, 296, 34, 62
    out.append(_text(ax + 2, 16, "WHAT YOU ALREADY HAVE", c, **cap))
    for i, (kind, title, sub) in enumerate(INPUTS):
        y = top + i * row + (20 if i else 0)
        needed = i == 0
        out.append(f'<rect x="{ax + 0.5}" y="{y + 0.5}" width="{aw - 1}" height="{row - 9}" rx="6" fill="{c["lamp_soft"] if needed else c["card"]}" '
                   f'stroke="{c["lamp"] if needed else c["line"]}" stroke-width="{1.5 if needed else 1}"/>')
        out.append(_icon(kind, ax + 13, y + 15.5, c))
        out.append(_text(ax + 46, y + 23, title, c, size=14, weight=600))
        out.append(_text(ax + 46, y + 41, sub, c, size=12, fill="muted"))
        if needed:
            out.append(f'<rect x="{aw - 88}" y="{y + 9}" width="78" height="20" rx="10" fill="{c["lamp"]}"/>')
            out.append(_text(aw - 49, y + 23.5, "all you need", c, size=11, weight=600, fill="on_lamp", anchor="middle"))
    y_opt = top + row + 5
    out.append(_text(ax + 2, y_opt + 8, "OPTIONAL · EACH KIND ADDS SECTIONS", c, **cap))
    # ------------------------------------------------------------------ scan
    mid = top + 4 * row + 20
    sx = aw + 10
    out.append(_arrow(sx, mid, sx + 132, c))
    out.append(f'<rect x="{sx - 2}" y="{mid - 44}" width="134" height="30" rx="6" fill="{c["ink"]}"/>')
    out.append(_text(sx + 65, mid - 24, "flowcellar scan", c, size=13, weight=600, fill="panel", anchor="middle", font=MONO))
    for k, line in enumerate(("reads, never writes", "nothing leaves", "the computer", "goes on where", "it stopped")):
        out.append(_text(sx + 65, mid + 26 + k * 17 + (6 if k >= 1 else 0) + (6 if k >= 3 else 0), line, c, size=12, fill="muted", anchor="middle"))
    # ------------------------------------------------------------------ the library
    lx, lw = sx + 140, 250
    ly = top + 62
    out.append(_text(lx + 2, 16, "THE LIBRARY · ONE FOLDER", c, **cap))
    out.append(f'<rect x="{lx + 0.5}" y="{top + 0.5}" width="{lw - 1}" height="{H - top - 12}" rx="8" fill="{c["panel"]}" stroke="{c["line"]}"/>')
    # the mark of the brand: a vault of nine cells, one lit (the library keeps one cell per run)
    mx, my = lx + 18, top + 16
    out.append(f'<g transform="translate({mx - 10 * 0.62} {my - 6 * 0.62}) scale(0.62)"><path fill="{c["ink"]}" d="M10.367 24V24A22 22 0 0 1 22 8.404V24ZM26 24V6.834A22 22 0 0 1 38 '
               f'6.834V24ZM42 24V8.404A22 22 0 0 1 53.633 24V24ZM10 28H22V41H10ZM26 28H38V41H26ZM10 45H22V58H10ZM26 45H38V58H26ZM42 45H54V58H42Z"/>'
               f'<path fill="{c["lamp"]}" d="M42 28H54V41H42Z"/></g>')
    out.append(_text(mx + 40, my + 14, "small, rebuildable,", c, size=13, weight=600))
    out.append(_text(mx + 40, my + 31, "never a copy of your data", c, size=13, weight=600))
    for i, (name, what) in enumerate(LIBRARY):
        y = ly + 6 + i * 96
        out.append(f'<rect x="{lx + 12.5}" y="{y + 0.5}" width="{lw - 25}" height="84" rx="6" fill="{c["card"]}" stroke="{c["line"]}"/>')
        out.append(_text(lx + 24, y + 24, name, c, size=13, weight=600, font=MONO))
        words, line, lines = what.split(), "", []
        for w in words:
            if len(line) + len(w) + 1 > 32:
                lines.append(line)
                line = w
            else:
                line = (line + " " + w).strip()
        lines.append(line)
        for k, ln in enumerate(lines[:3]):
            out.append(_text(lx + 24, y + 44 + k * 16, ln, c, size=12, fill="muted"))
    # ------------------------------------------------------------------ what you look at
    px = lx + lw + 38
    pw = W - px
    out.append(_arrow(lx + lw + 8, mid, px - 6, c))
    out.append(_text(px + 2, 16, "WHAT YOU LOOK AT · IN A BROWSER", c, **cap))
    step = (H - top - 12) / len(PAGES)
    for i, (name, what) in enumerate(PAGES):
        y = top + i * step
        out.append(f'<rect x="{px + 0.5}" y="{y + 0.5}" width="{pw - 1}" height="{step - 10}" rx="6" fill="{c["card"]}" stroke="{c["line"]}"/>')
        out.append(f'<rect x="{px + 0.5}" y="{y + 0.5}" width="5" height="{step - 10}" rx="2.5" fill="{c["lamp"]}"/>')
        out.append(_text(px + 18, y + (step - 10) / 2 - 3, name, c, size=14, weight=600))
        out.append(_text(px + 18, y + (step - 10) / 2 + 16, what, c, size=12, fill="muted"))
    out.append("</svg>")
    return "\n".join(out) + "\n"


def write(out_dir: str | None = None) -> list[str]:
    out_dir = out_dir or os.path.join(REPO, "docs", "img")
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for theme in THEMES:
        p = os.path.join(out_dir, f"overview-{theme}.svg")
        with open(p, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(illustration(theme))
        paths.append(p)
    return paths


if __name__ == "__main__":
    for path in write(sys.argv[1] if len(sys.argv) > 1 else None):
        print(path)
