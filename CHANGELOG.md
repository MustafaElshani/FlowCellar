# Changelog

Versions of the package, the bundle schema and the metric catalogue are counted separately; each entry
names all three.

## 0.1.0 (November 2026)

The first public release. Bundle schema 0.1.1, metric catalogue 0.1.1, library layout 1.
Licence: GNU General Public License, version 3 or any later version (`GPL-3.0-or-later`).

### What it does

- **Run health** from the folder MinKNOW writes, for the file layouts of MinKNOW 5 and 6: quick detail from
  the small files in seconds, full detail from one pass over the sequencing summary and the pore-scan
  table. Assay detection, runs in progress, interrupted runs, folders with raw data only. A flow cell that
  is loaded a second time is recognised by its loadings and restarts and judged as a re-used flow cell.
- **Verdicts** in two parts, Run health and Read data, each with stored sentences that say what was
  measured and against which limit. Output is compared with the lab's own history as of the day a run
  started; a verdict is not changed by runs that come later.
- **Read data**, counted exactly and never sampled: one pass over Dorado and aligned BAMs (in parallel
  where a BAM has an index), pychopper statistics, salmon, oarfish and kallisto output, basecaller sample
  sheets. Files are attached to runs by the identifiers inside them first and by configured patterns
  second. One source of record per run, sample and question; two candidates are never added together.
  Sections: basecalled output, length, quality, length against quality, alignment, cDNA structure, read
  fate.
- **Library**: `scan` over any number of roots creates and updates it, remembers the roots, names
  projects from the folder structure and treats a project or a run found under two roots as one. Every
  pass result is cached, so a second scan opens no data file. A root that is offline is skipped and
  nothing under it changes. One scan at a time per library.
- **Viewer**: Home, Runs, Projects, a project, a run, a comparison of runs and the history of a flow
  cell. The same viewer as one offline report per run, as a static hub opened from a file, and over
  `flowcellar serve`. No page makes a request to anywhere.
- **Commands**: `scan`, `run`, `hub`, `serve`, `init`, `discover`, `ls`, `show`, `rm`,
  `recompute --status-only`, `doctor`, `validate`, `schema`.
- **Packaging**: one pure-Python wheel and one source distribution, both with the built viewer and the
  catalogue inside; Python 3.11 to 3.15; run-time dependencies numpy, pyarrow and click, with pysam in the
  extra `bam` and jsonschema in the extra `validate`.

- **Every page uses the whole width of the window**, with the same two edges on every screen, and a pointer
  over a chart changes the text of its read-out and nothing else.
- **Alignment**, from the same single pass over each BAM: a summary of read alignment as a ring with its
  table (mapped and unmapped reads, secondary and supplementary records, reads with a supplementary
  alignment); the breakdown of aligned bases (match, mismatch, insertion, deletion, soft and hard clip);
  strand balance; depth along the reference with chromosome boundaries and the mean of each chromosome; reads
  and bases per chromosome with the mitochondrial share; aligned length beside read length; the basecaller's
  quality against the quality measured from the alignment; soft-clipped share per read; accuracy by read
  length; read length against aligned length, identity against aligned length, read quality against identity.
- **Reads over run time**: output from the basecalled reads (reads and bases, passed and all), read length,
  read quality, aligned length, identity and the measured speed through the pore, as a median with the
  middle half of the reads, on the time axis of the run.
- **Transcripts**: for alignments to a transcriptome, coverage along the transcript from 5′ to 3′ by
  transcript length, transcripts and genes found against reads sequenced, reads per transcript; poly(A) tail
  length per sample from the `pt` tag, with the share of reads that have an estimate.
- **Pychopper**: a section of its own with every category of the statistics file, per sample and added up:
  classification, strand, rescue, primer hits, primer arrangements, length and quality failures, UMIs, and
  the tuning curve with the cut-off that was chosen.
- **RNA composition** from annotated read tables (role `polya.read_table`, read in one streamed pass): reads
  and poly(A) tail length by transcript type, genes detected, the share of the most abundant genes.
- **Samples**: reads and bases per sample as charts, with the reads no sample was given to.
- **Single-cell and spatial analyses**: output folders of wf-single-cell and Space Ranger (Visium, Visium HD)
  are recognised by their files, read, linked to the flow cell runs their reads came from and judged on their
  own. New command `flowcellar analysis`; an Analyses tab and a page per analysis with a spatial map, read
  fate, saturation, distributions per spot, bin or cell, and the two pipelines side by side; `ls --analyses`,
  `show an_…`. The extra `sc` adds h5py for Space Ranger's HDF5 files.
- **Scan**: after an upgrade that changes what a pass counts, the files an earlier version had read are read
  again first (newest runs first), then files never read; no file twice in one scan. `scan --read
  bam|tables|analyses` chooses which kind of large file a scan reads.
- **Index**: column `t2_full_length` (full-length share of a cDNA run), shown in the peek panel of Runs and in
  Compare; member `analyses`.
- **Wide windows**: the opening charts of Alignment lie four side by side; the spatial map of an analysis shows
  a second picture of the same slide beside the first (another colouring, the tissue image where there is one),
  and the two share the cursor.
- **`flowcellar run`**: the genome and the transcriptome alignment of one sample are taken for one sample
  whatever their file names end in (`_genome`, `_gnome`, `_transcriptome`, `_txome`, …).
- Fixed: a roll-up over cells that hold curves no longer fails on the third cell; merged count tables keep
  the name of their weight and the count that was folded into "other"; a share that is not zero is no longer
  written "0.0 %"; a file read by an earlier version of a pass no longer counts as read for ever; a folder
  into which the reports of several runs were copied is no longer taken for a run folder.

- **Save as PDF.** The page of a run, of a single-cell or spatial analysis, a comparison of runs and the
  history of a flow cell have a control "Save as PDF", in the hub and in the one-file reports. It opens the
  print dialog of the browser with a layout made for paper (A4 or Letter, upright): the first sheet says what
  the sheets are, when the numbers were computed and by which version, then the verdict with both parts, the
  key numbers and the selection that is shown; every section follows in the order of the page, without
  navigation, controls and links, in the light theme whatever the screen shows. No chart, row of key numbers
  or table row is cut where a sheet ends, no heading stands alone at the foot of one, and a table that is
  wider than the sheet (the differences of eight compared runs) is scaled down until it fits. A PDF holds the
  page as it is shown (open groups are in it, closed ones are counted as left out); a run page also offers
  "PDF with every chart". The browser's own Print command gives the same sheets. Nothing is rendered outside
  the browser and no dependency was added.
- **One analysis of merged reads, said as that.** A pychopper table or an annotated read table of a sample
  that was sequenced on two flow cells is described as what it is: one analysis of reads that were merged
  before the tool ran, shown on the page of each flow cell. A single-cell or spatial analysis made from several
  flow cells says the same in its header, under its verdict and on the run page. The notes a bundle stores and
  the sources of read fate use the same words.
- **A project as one PDF.** The page of a project has "Save as PDF" too. On paper it is a report of the project:
  its numbers, every run of it oldest first with its verdict, run time, output, reads, pores, occupancy and read
  N50, the reason of each run that did not pass, and the single-cell and spatial analyses made from its runs. The
  folders the project lies in are left out of the PDF.
- A change of theme now reaches the colours of charts and legends that are already drawn (they kept the
  colours of the theme the page was opened in).
- The browser checks run in Chromium as well (`FLOWCELLAR_BROWSER_ENGINE=chromium`); `viewer/verify/print.mjs`
  checks the PDFs of both browsers.
- **An example to try it on.** `python -m flowcellar.testing.demo DIR` writes an invented laboratory (15 runs of 13
  flow cells in four projects, with BAM files, pychopper statistics, annotated read tables, quantification
  metadata and the analysis folders of a Visium and a Visium HD slide) and the `flowcellar.toml` of a library for
  it. The README and its pictures are made from it; `tools/make_readme_pictures.py` makes the pictures again.
- **Limits of read data checks can be set per kind of run**, like those of run health:
  `[profiles.<assay group>.rules.<check id>]` with `warn`, `fail` or `enabled = false`.
  `flowcellar --config FILE run RUN_FOLDER -o DIR` judges a run without a library with the limits of that file.
- An analysis folder is no longer left out because a folder *above* the scanned folder is named `work` (data
  under `/work/...`): scratch folders are told by the names below the scanned folder only.
- "Typical here is 727 Mb to 1.1 Gb": the lower end of a range keeps its unit when the two ends differ in unit.
- On PyPI the long description is the README with every path written as an address (`hatch_build.py`).

### Known limits

What is measured

- The alignment, transcript, pychopper and RNA composition charts are shown, not judged: no check and no limit
  uses their numbers yet.
- Charts over run time are binned on the clock (ten minutes for output, one hour for the medians), because a
  pass over a file does not know when the run started. A step can begin up to that long before a whole step
  of the run timeline; every point is drawn at its true time and the page says so. A rate per hour of the first
  and the last step is worked out over the part of the step the run covers.
- A pychopper statistics table or an annotated read table that covers a sample on two flow cells cannot be
  split: its counts are shown in full on both runs, and each page says so.
- A pychopper table that lacks a category stores zero for the numbers read from it; the section names them.
- Depth along the reference counts positions under a match, mismatch or deletion of primary alignments;
  introns are not covered, secondary and supplementary records are left out.
- The curve of transcripts found against reads sequenced is the exact expectation from the read count of each
  transcript, not a replay of the run in time. Reads per transcript are kept only while samples times
  reference sequences stay under 64 million.
- Whether a BAM is aligned to a genome or to a transcriptome is told from its header; an unusual reference
  can be misjudged (the evidence is stored with the section).
- The measured speed reads a little low on reads that were trimmed after basecalling: the duration covers the
  whole read.
- For a choice of several samples, medians are read from stored histograms and can be one bin off; for all
  samples and for one sample they are exact. Gene numbers exist per sample only.
- In an annotated read table a read that pychopper cut into pieces counts once per piece.
- From 13 samples on, two-dimensional charts, depth and curves over time are kept once for the run, not per
  sample.
- Analyses: a link to a run that rests on a path or a sample name says so and is of low or medium confidence;
  for wf-single-cell 3.x the molecules of the gene matrix and of the per-read table differ, so the saturation
  tile and the saturation curve end a few points apart (the page says so); limits are provisional and not
  tuned on any history.

- The estimated read N50 and median (before basecalling) are exact: the pass over the sequencing summary counts
  reads per number of events. A result cached by a build that kept no such count is read off a histogram with
  40 bins per factor of ten instead (by rank inside a bin: typically within a fraction of a per cent of the
  exact value, at worst the width of one bin, about 6 %), and the run page says so ("approximate");
  `flowcellar scan --force --tier 1` counts such runs again.
- Read N50, median read length and median read quality of the read data, and the speed median, are read off
  histograms by rank inside a bin (read lengths: 40 bins per factor of ten; speed: 5 bases per second). They
  are close to the exact value (at worst one bin away), not exact.
- `scan` does not pick up the BAM files MinKNOW writes inside a run folder when it basecalls during the
  run. `flowcellar run` on that folder does.
- The time box stops new passes from starting. A pass that is running when the box ends is finished, so a
  scan can run a few minutes longer than its box.
- Home counts the read-data files that wait per run, the last line of `scan` counts them per file; the two
  numbers can differ.

What is shown

- Of the 238 metrics in the catalogue, 107 have a label (and all 70 derived numbers, and all 57 metrics of
  the single-cell and spatial catalogue); the others are shown under a name made from their id.
- The viewer is checked in Firefox and in Chromium (Playwright builds on Linux): not in Safari, not in Chrome
  or Edge on Windows or macOS, not with a screen reader.
- A PDF comes from the print dialog of the browser: there is no `flowcellar` command that writes one (a
  headless Chrome can print a report file from the command line). Home, Runs and Projects have no layout for
  paper. Firefox writes charts into a PDF as pictures (three pixels to one, so they stay sharp); Chromium
  writes them as drawings. Sheet numbers are the browser's.
- A report file names the licence of the viewer inside it and where its source is, and not the full licence
  texts (FlowCellar's and Preact's); the package holds them, in `LICENSE` and `THIRD_PARTY_NOTICES.md`.

Where it runs

- Tested on Linux only. macOS is untried; on Windows pysam cannot be installed, so read data is not
  available there, and the rest is untried.
- The two GitHub workflows have been checked step by step on a developer's machine but have never run on
  GitHub.
- Used as a Python library, the read-data pass cannot start its worker processes from a program that is
  given on standard input. The command line is not affected.
- With a temporary folder (`TMPDIR`) whose path is longer than about 75 characters, a pass that runs in
  worker processes fails with "AF_UNIX path too long" (seen on the pass that counts the reads: Python keeps
  the socket of its worker server in that folder). The scan names each run it could not read and ends with
  exit code 5; with a shorter `TMPDIR` it runs.

### Not built yet

Adaptive sampling; duplex; modified bases; what is particular to direct RNA beyond run health and poly(A)
tails; barcode balance; limits, hub columns and search keys for the numbers of those charts; analyses in
the search and in Compare; notes and tags on runs; saved views; editing limits in the viewer; a per-sample page; a project report file; a forecast for a run in
progress; export to tables and to MultiQC; import of bundles; computing parts of a run elsewhere and
merging them; migration between schema versions; a conda package.

### Schema

Bundle schema 0.1.1. What changed against the draft 0.1.0 is listed in the schema file itself
(`flowcellar schema`, member `x-changes`). Readers accept any 0.1.x; no earlier version was released, so
there is nothing to migrate.
