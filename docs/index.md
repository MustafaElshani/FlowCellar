# FlowCellar documentation

FlowCellar is quality control for Oxford Nanopore sequencing runs. It reads the folders MinKNOW writes and,
when they exist, the basecalled and aligned reads. It keeps one small file per run in a **library**, judges
each run, and shows everything in one viewer that works offline.

| Page | What is in it |
|---|---|
| [Quick start](quick-start.md) | install, the example data, the first scan, opening the library, one run without a library, paper |
| [The library](library.md) | what a library is, roots, where names come from, projects in several places, what a scan does |
| [Configuration](configuration.md) | `flowcellar.toml` with a complete example, file patterns for read data, analyses and their flow cells, limits |
| [Commands](commands.md) | every command with its options, written from the `--help` output; exit codes |
| [Metric catalogue](metrics.md) | every check with its limits and how to read it, written from the catalogue |
| [Development](development.md) | building the viewer, tests, the tools |
| [Releasing](releasing.md) | from this repository to PyPI |

## In one minute

```
pip install "flowcellar[bam]"
flowcellar scan /data/minknow /data/projects -L ~/flowcellar-library
flowcellar serve ~/flowcellar-library --open
flowcellar run /data/minknow/demo_experiment/sample_a/20260302_0900_MD-000001_TST00001_0a1b2c3d -o qc/
```

No data at hand? `python -m flowcellar.testing.demo demo` writes an invented laboratory to try all of this on
([Quick start](quick-start.md)).

`scan` finds every run under the folders it is given, writes a verdict for each from MinKNOW's small files
within the first minutes, and then counts the reads, newest runs first. `serve` shows the library in a
browser. `run` handles one run folder without a library: it writes a bundle (`<run key>.fcb.json.gz`) and
one HTML report that opens from a file, a mail attachment or a shared drive.

Nothing is copied from, moved in or written to the folders that are read, and no command contacts the
network.

## Words used on every page

| Word | Meaning |
|---|---|
| run | one start of sequencing on one flow cell: one folder MinKNOW wrote |
| run key | flow cell id and the first eight characters of the run id, for example `TST00001_0a1b2c3d` |
| library | the folder FlowCellar owns: one bundle and one verdict per run, an index, the hub |
| root | a folder you asked FlowCellar to look under; a library has any number of them |
| project | named by the folder structure under a root; the same name under two roots is one project |
| bundle | everything measured for one run, as one compressed JSON file |
| hub | the pages of a library: Home, Runs, Projects, Analyses, a run, an analysis, a comparison, the history of a flow cell |
| analysis | one output folder of wf-single-cell or Space Ranger: a single-cell or spatial sample, with a verdict of its own |
| run health | the part of a verdict that comes from MinKNOW's own files; needs no basecalling |
| read data | the part of a verdict that comes from basecalled and aligned reads |
| quick detail, full detail | run health from the small files only, or with every read of the run counted |
