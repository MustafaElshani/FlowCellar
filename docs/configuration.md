# Configuration: `flowcellar.toml`

The file lies in the library folder. `scan` creates it and only ever *adds* a `[[roots]]` entry to it;
everything else is yours to edit. `flowcellar init DIR --from-scan ROOT` (the option can be given several times) drafts one with every section and
with proposals for file patterns, written as comments for you to check. A key that is left out keeps its
default; an unknown key is a warning, never an error.

The package knows the layouts MinKNOW and Dorado write themselves. Everything that is particular to a
laboratory (where its pipeline puts BAMs and tables) belongs in this file. See what the patterns find, before
anything is computed, with `flowcellar discover -L DIR --explain`.

## A complete example

```toml
[library]
name = "our nanopore runs"

[[roots]]
id = "minknow"
path = "/data/minknow"
live = true                       # the sequencer writes here: read gently
project_depth = 1

[[roots]]
id = "projects"
path = "/data/projects"
project_depth = 2

[discovery]
min_age_seconds = 120             # a file changed more recently is left for the next scan
stale_after_hours = 6             # a run without a final summary counts as running while files still change
prune = ["notebooks", "figures"]  # folder names (shell patterns) never listed, added to the built-in list
subtree_listing_budget = 20000    # folders listed under one first-level folder before it is given up
root_check_seconds = 20
stall_timeout_seconds = 120

# ---- read data: where the pipeline leaves its files, relative to each root -------------------------------
[[inputs]]
role = "dorado.bam"
glob = "**/*_basecalled/**/bam_pass/{sample}/*.bam"

[[inputs]]
role = "align.genome_bam"
glob = "**/analysis/alignments/{sample}_genome_sorted.bam"
exclude = ["*_draft_*"]

[[inputs]]
role = "pychopper.stats"
glob = "**/analysis/pychopper_reports/{sample}_pychopper_stats.tsv"

[[inputs]]
role = "polya.read_table"
glob = "**/analysis/polya/{sample}_readid_annotated.csv"

[[inputs]]
role = "quant.salmon"
glob = "**/analysis/counts/salmon/{sample}_salmon/aux_info/meta_info.json"

[[inputs]]
role = "basecaller.sample_sheet"
sniff = "csv_columns:barcode,alias,flow_cell_id|position_id"

[[link.rules]]                    # tables without identifiers inside: linked through the sample name
name = "experiment-folder"
roles = ["pychopper.stats", "polya.read_table", "quant.salmon"]
via = "sample"
within = "experiment"

[variants]
primary = "main"                  # when a sample has several basecallings, this one is counted

[samples]
unclassified = ["unclassified"]
[samples.alias]
"name as recorded" = "name to show"

[scan]
jobs = 0                          # worker processes per large file; 0 chooses
threads = 4                       # decompression threads for a BAM without an index
time_box_minutes = 240
live_workers = 2

[hub]
sidecars = true                   # one data file per run beside the hub page, so that it opens without a server
reports = false                   # also write one offline report per run into reports/

[serve]
host = "127.0.0.1"
port = 8737
```

## `[[roots]]`

| Key | Meaning |
|---|---|
| `id` | a short name shown on Home and stored with every run; do not change it later |
| `path` | the folder to look under |
| `live` | `true` for the folder a sequencer is writing to: it is read with few workers and at low priority |
| `project_depth` | which folder level under the root names the project: 1 = the first (the default), 0 = none |
| `project_pattern` | instead of a depth: a pattern with the capture `{project}`, for example `"*/{project}/**"` |

A folder name that an `[[inputs]]` pattern spells out is never taken for a project. With the pattern
`**/analysis/alignments/{sample}.bam`, a file at `<root>/alpha_project/analysis/alignments/s1.bam` under a root
with `project_depth = 2` has `analysis` at project depth: that names no project, and the run keeps the project
of its own folder. More names can be listed:

```toml
[project]
exclude = ["scratch", "old"]      # folder names that are never a project
```

## `[[inputs]]`

| Key | Meaning |
|---|---|
| `role` | what the file is: `dorado.bam`, `align.genome_bam`, `align.transcriptome_bam`, `pychopper.stats`, `polya.read_table` (an annotated read table: one row per read with its transcript, transcript type, gene and poly(A) tail), `quant.salmon`, `quant.oarfish`, `quant.kallisto`, `basecaller.sample_sheet`, `dorado.summary`; for analyses `analysis.wf_single_cell` and `analysis.spaceranger` (below) |
| `glob` | a path pattern relative to each root. `**` crosses folders; `{sample}` captures the sample name from one path component; `{variant:a\|b}` captures one of the listed words |
| `exclude` | file name patterns to leave out |
| `variant` | the basecalling variant these files belong to |
| `sniff` | recognise small tables by their content instead of their path, for example a basecaller sample sheet by its column names |

A BAM is attached to a run by what is inside it: the read group of its header, the `RG` and `fn` tags of its
reads. A file that holds reads of two runs is split read by read. Patterns only say where to look and which
sample a file belongs to. Small tables have no identifiers; a `[[link.rules]]` entry says how they follow
their sample (`via = "sample"`, inside the same experiment folder).

Three things follow from this:

- **A sample on two flow cells in one file.** When the sample sheet of each flow cell names the sample its own
  way and the reads carry that name, the reads of each flow cell are shown under the sample of its own sheet,
  not under a second sample named after the file.
- **A sample sequenced on several flow cells.** Each flow cell keeps its own page, with its own run health and its
  own reads (an aligned BAM of the merged sample is split read by read, by the tags its reads carry). A step that
  ran on the merged reads of the sample (pychopper statistics, an annotated read table, a quantification) is one
  analysis: it is shown on the page of each of those flow cells, and the page says that it is one analysis of reads
  merged before that step. A check that needs the count of one flow cell alone from such a table (the usable share
  of read fate, when the quantification covers several flow cells) reads "Not assessed" with that reason.
- **A file without a variant in its name is the primary variant.** With `[variants] primary = "main"`, the
  alignment `sample_genome_sorted.bam` belongs to the basecalling `main/sample.bam`, and
  `sample_other_genome_sorted.bam` (found by a pattern with `variant = "other"`) to `other/sample.bam`. The
  primary one is counted; a page that shows another version says which.

## Single-cell and spatial analyses

An analysis folder is recognised by one file inside it: `wf-single-cell-report.html` in the output folder of
wf-single-cell, `metrics_summary.csv` in the `outs/` folder of Space Ranger. The package looks for both file
names by itself in every folder a scan lists. A scan does not list below a folder that one of your own patterns
starts at (an `analysis` folder, say): what lies there is found through your patterns. So in an archive whose
analyses lie under such folders, and in any deep archive, say where they are, and keep the pipeline's work
folders out:

```toml
[discovery]
prune = ["work", "SPATIAL_RNA_COUNTER_CS"]

[[inputs]]
role = "analysis.wf_single_cell"
glob = "**/analysis/*/wf-single-cell-report.html"

[[inputs]]
role = "analysis.spaceranger"
glob = "**/analysis/*/outs/metrics_summary.csv"
```

Write one pattern per depth instead of a `**` in the middle of a pattern: `**/analysis/**/outs/…` makes the
walk descend below every `analysis` folder of the root. `flowcellar discover -L DIR` says how many analysis
folders the patterns find.

A single-cell or spatial sample is usually sequenced on several flow cells. An analysis is one analysis of their
merged reads, and it is linked to every run its reads came from:

- by what is inside its files: the names of the raw-data files on the reads of the workflow's tagged BAM
  (`<flow cell>_<run id>_…pod5`) name the runs and give the share of the reads from each; a read group that kept the
  flow cell and the run id names the runs without a share. A Space Ranger folder holds no such names and takes the
  links of the wf-single-cell analysis made from the same reads;
- else by the place and the name: a run folder that lies beside the analysis and carries its sample name, or the
  sample name alone when exactly one run of the library carries it;
- else by a person: `flowcellar analysis FOLDER --run FLOWCELL -L DIR`.

The page of the analysis names those flow cells with their shares and verdicts; the page of each flow cell lists
the analyses it fed. A flow cell that the reads name and the library does not hold is shown as "run not in
library" and becomes a link once that run is scanned. The verdict of an analysis never changes the verdict of a run.

## Limits and comparison groups

```toml
[profiles.cdna.rules.fc_occupancy_first_hour]
warn = 0.65                       # override one limit for one assay group
# fail = 0.4      role = "supporting"      enabled = false

[baselines]
window = 40                       # how many earlier runs form the lab history
min_runs = 5                      # history is used from this many runs
trusted_runs = 8
spread_guard = 0.45               # a history that scatters more than this is not used for verdicts

[[baselines.groups]]              # give runs a comparison group of their own
match = { kit = "SQK-PCB114-24", started_after = "2026-01-01" }
label = "2026 protocol"

[flow_cell]
restart_max_hours = 2             # a start this short on a flow cell that is loaded again counts as a restart
restart_max_output_frac = 0.02
```

After changing limits: `flowcellar recompute --status-only --all --rebaseline -L DIR` judges every run again
without reading any input.

The `[profiles.<assay group>.rules.<check id>]` tables change the limits of a check for one kind of run (`gdna`,
`cdna`, `drna`): `warn` and `fail` replace the built-in levels, and `enabled = false` switches the check off, which
the page shows as "Not assessed". This holds for the checks of run health and for those of read data; the check ids
and their defaults are in the [metric catalogue](metrics.md). A page judged with limits of the library says so
under its verdict. `[profiles.<assay group>]` also takes `outcome_unit = "reads"` or `"bases"`: the unit in which
the output of a run is judged.

```toml
[profiles.gdna.rules.len_n50]
warn = 15000                      # a read data check: read N50 of genomic runs
fail = 8000

[profiles.cdna.rules.yield_basecalled_bases]
enabled = false                   # do not judge basecalled output of cDNA runs
```

The output of a run is judged at the longest of 6, 12, 24, 48 and 72 hours that it reached. The levels given for
`yield_reads_at_72h` and `yield_estimated_bases_at_72h` are scaled to the shorter horizons by the shape of a typical
run, so one pair of numbers is enough.

A single run without a library is judged with the built-in limits; `flowcellar --config FILE run RUN_FOLDER -o DIR`
judges it with the `[profiles]` and `[flow_cell]` tables of that file.
