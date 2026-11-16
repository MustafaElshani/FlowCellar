# The library

A library is one folder that FlowCellar owns. It holds, for every run it has seen, a small **bundle** (the
measured numbers, some tens to a few hundred kB compressed), a **status file** (the verdict and its
reasons), and an **index** that the viewer loads first. It never holds a copy of your data, and nothing is
ever copied from, moved in or written to the folders it was built from.

```
flowcellar scan /data/minknow /data/projects -L ~/flowcellar-library    # creates the library; two roots
flowcellar scan /archive/2025 -L ~/flowcellar-library                   # later: adds a third root, scans all
flowcellar scan -L ~/flowcellar-library                                 # later: looks for what is new
flowcellar run /scratch/one_run_folder -L ~/flowcellar-library          # adds one run from anywhere
```

What is inside:

```
<library>/
  flowcellar.toml      the configuration; the only file a person edits
  LIBRARY              marker: this folder is a FlowCellar library
  bundles/             <run key>.fcb.json.gz, one per run: what was measured (an_….fcb.json.gz: one per analysis)
  status/              <run key>.status.json: verdict, checks, sentences
  cards/  state/       one row of the index per run; bookkeeping of the scans
  partials/  probes/   caches of pass results and of file identities; safe to delete
  index/               index.json.gz: everything Home, Runs and Projects need
  hub/                 index.html, index.js, data/<run key>.js: the pages
  logs/                one file per scan: every decision and every error
  trash/               what `flowcellar rm` moved away
```

## Where its content comes from: roots

A **root** is a folder you ask FlowCellar to look under: the output folder of the sequencer, a project share,
an archive. A library can have any number of roots; `scan` remembers each one as a `[[roots]]` entry and Home
shows, per root, how many runs were found, when it was last scanned and whether it is reachable.

- **A run is recognised by the files in its folder, never by its name**: a final summary, a report, a
  throughput file, a sequencing summary. A folder with raw data only (POD5 or FAST5 and no report) becomes an
  *Incomplete* entry so that the run is not forgotten. A folder that is merely named like a run but holds only
  basecaller output is not a run.
- **A run is known by its run key**: the flow cell id plus the first eight characters of the protocol run id,
  for example `TST00001_0a1b2c3d`. The same run found under two roots is one run with two paths. A run folder
  that was moved is found again at its new place and nothing is computed again.
- **A root that is offline or does not answer is skipped.** Nothing known under it changes, Home shows one
  notice, and the scan ends with exit code 5. The next scan picks it up again.
- **A folder that cannot be read is named.** A folder under a root that you have no permission to list is
  skipped with a warning and a notice on Home, and the scan ends with exit code 5; runs the library already
  holds from there stay as they were (they are not taken for gone).
- **A large file that cannot be opened is tried again.** A run whose sequencing summary could not be opened is
  not marked as counted: the scan warns, ends with exit code 5, shows the run as waiting, and the next scan
  tries again by itself.
- **A run folder that disappeared keeps its bundle.** Its row says the source is missing. Before saying so the
  scan looks at the folder itself: a run that was added by hand (`flowcellar run -L`) from a place a scan does
  not walk, for instance under a folder named `work`, is still there and stays as it was, with a notice.
- **`live = true`** marks the folder a sequencer is writing to: it is read with two workers and at low
  priority, and a run in progress appears as *Running*.
- A file or run folder that changed in the last two minutes is left for the next scan, so that a folder that
  is still being copied is not read half-way (not applied to a `live` root).

## Projects

Each root says which folder level names the project:

```toml
[[roots]]
id = "minknow"
path = "/data/minknow"
project_depth = 1          # /data/minknow/<project>/...

[[roots]]
id = "projects"
path = "/data/projects"
project_depth = 2          # /data/projects/<group>/<project>/...
```

`project_depth = 1` (the default) takes the first folder under the root. `project_pattern =
"*/{project}/**"` does the same with a pattern. **A project found under two roots is one project**: its page
and its row on Home list both roots. A run belongs to the project of its run folder and to the projects of the
read-data files that were counted for it. A run directly under a root has no project.

Changing `project_depth` and scanning again regroups the runs; nothing is read or computed again.

The whole file, with every section, is described in [Configuration](configuration.md).

## Where names come from

FlowCellar asks for no table of names. It takes them from the files, in a fixed order, and says where each one
came from (the page of a run lists its sources under "Data sources"; `flowcellar discover -L DIR --explain` prints
the reason for every file).

**The run.** Flow cell id, run id, start time, kit, device, experiment and sample id come from MinKNOW's report
in the run folder, then from its final summary, then from its sample sheet, then from the names of its files, in
this order; the name of the folder is the last resort (it is all a folder with raw data only has). When two of
these name different flow cells the page carries a warning and the first in this order is used.

**The samples of a run.**

1. A sample sheet of the basecaller: a CSV with the columns `barcode`, `alias` and `flow_cell_id` or
   `position_id`, in any order, as Dorado takes it. It is recognised by its columns, whatever it is called, in a
   folder that holds a run folder or a basecaller output folder (`bam_pass/`, or a folder one of your patterns
   starts at); `flowcellar run` takes one with `--sample-sheet`. The rows of that flow cell name its samples: the
   `alias`, else the `sample_id`, else the barcode. The sample sheet MinKNOW keeps *inside* a run folder names the
   run, not its barcodes.
2. Without a sheet: the place of a file. `bam_pass/<name>/…` and `{sample}` in an `[[inputs]]` pattern give the
   name; a file of your pipeline is named by the pattern that found it.
3. Else the one label the reads of a file carry (the barcode or alias the basecaller wrote on them).

What is inside a file wins over its name: a file that is named after one sample and whose reads carry the alias of
another sample of the sheet is shown under the sample of the sheet, with a note. A name that is in the data and not
in the sheet is a sample of its own, with a note; when it holds under 0.1 % of the reads its counts are kept and
its distributions are not. Two sheets that name one barcode differently give a note, and the first is used.
`[samples.alias]` in `flowcellar.toml` renames a sample; `[samples] unclassified` lists the names that mean "no
barcode found".

**A sample on several flow cells.** Each flow cell keeps its own page, with its own run health and its own reads.
Steps that ran on the merged reads of the sample (pychopper statistics, annotated read tables) are one analysis,
shown on the page of each of those flow cells and said there to be one analysis of reads merged before that step.

**Several basecalling passes of one sample** (a second pass with a modified-base model, for instance): one is
counted, the one named by `[variants] primary`; the others are listed with the reason. See "Read data: one source
of record" below.

**Projects** come from the folder structure (next section). **Analyses** are linked to the flow cells their reads
came from by the identifiers inside their files: see [Configuration](configuration.md).

## What a scan does, in order

1. takes the library's lease (one scan at a time; a lease left by a crashed scan is taken over with a notice,
   also when its process number has since been given to another program: a running scan touches the lease
   every minute, and a lease nobody touched for three minutes belongs to nobody). The time box starts when the
   lease is held, not while waiting for another scan;
2. checks every root and lists its folders in a separate process with a deadline;
3. for every new run: a quick bundle and a verdict from the small files; then the hub is written;
4. links read-data files to runs, by the identifiers inside the files first and by configured patterns second;
5. full passes over sequencing summaries and pore-scan tables, newest runs first, inside the time box;
6. read-data passes inside what is left of the time box: first the files an earlier version of FlowCellar had
   read, when this version counts more from each file, then the files never read; newest runs first each time,
   and no file twice in one scan;
7. single-cell and spatial analyses: each is read from its small files, linked to its runs and judged; then the
   per-read tables of wf-single-cell are read, one analysis at a time, while the time box lasts;
8. writes the index and the hub and releases the lease.

A failed pass is not tried again at every scan: `scan --retry-failed` does that.

`scan --read bam`, `--read tables` and `--read analyses` (repeatable) choose which kind of large read-data file
a scan reads: BAM files take hours in a large archive, annotated read tables and the per-read tables of
analyses take minutes. What a scan leaves out waits for a scan that reads it. Until then the page of the run
says how many files wait and gives the command that reads them (`flowcellar scan -L DIR --only RUN`).

After an upgrade of FlowCellar a scan may link the files again or put the pages of the runs together again
(when the new version attaches files or builds a section differently). Neither opens a data file: both work
from what the library already holds, and take seconds.

## Read data: one source of record

A run can have several files that describe the same reads (a Dorado BAM, an aligned BAM made from it, a second
basecalling). FlowCellar picks **one source of record per run, sample and question** and never adds two
candidates together. When two candidates cannot be told apart, neither is counted, both are named in a notice,
and the library waits for a rule in `flowcellar.toml`. `flowcellar discover --explain -L DIR` prints, for every
file, which run it was linked to and on what evidence.

## Deleting and rebuilding

Everything in a library is derived from the roots and from `flowcellar.toml`. Deleting `partials/` and
`probes/` only costs time: the next scan reads the inputs again. Deleting the whole library and scanning
again builds it anew: runs are judged in the order they were started, each against the runs before it.
What a new library does not know: runs that were added one by one with `run -L` from outside the roots, and
removals made with `rm`. Keep a copy of `flowcellar.toml` when it holds patterns and limits you wrote.

## Housekeeping

```
flowcellar ls -L DIR --where "status:warn kit:SQK-PCB114-24"      the Runs table as text, same search as the hub
flowcellar ls -L DIR --analyses                                   the single-cell and spatial analyses
flowcellar show RUN_KEY -L DIR --status                           one run: every check and its sentence
flowcellar rm RUN_KEY -L DIR                                      move a run to the library's trash (never the data)
flowcellar recompute --status-only --all -L DIR                   judge again with the current limits; reads no input
flowcellar doctor -L DIR                                          the installation, the library, its roots
flowcellar validate DIR                                           every bundle against the schema
flowcellar hub DIR                                                write the pages again
```
