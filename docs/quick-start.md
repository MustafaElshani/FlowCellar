# Quick start

## 1. Install

Python 3.11 or newer, on Linux.

```
python3 -m venv ~/flowcellar-env
~/flowcellar-env/bin/pip install "flowcellar[all]"
~/flowcellar-env/bin/flowcellar doctor
```

`doctor` prints the versions and says whether the viewer is inside the package. The extra `bam` adds
pysam, which reads BAM files (it has no Windows build; run health works without it); `sc` adds h5py, which
reads the HDF5 files of a Space Ranger analysis; `validate` adds jsonschema, which `flowcellar validate`
needs; `all` is the three together.

To work on FlowCellar itself, use a checkout as [Development](development.md) describes.

## Without your own data: the example

```
python -m flowcellar.testing.demo demo
flowcellar scan demo/minknow demo/projects -L demo/library
flowcellar serve demo/library --open
```

The first command writes an invented laboratory: 15 runs of 13 flow cells in four projects, with BAM files,
pychopper statistics, annotated read tables, quantification metadata and the analysis folders of a Visium and a
Visium HD slide. Everything on this page can be tried on it, with `demo/minknow` and `demo/projects` in the place
of `/data/minknow` and `/data/projects`, and `demo/library` in the place of `~/flowcellar-library`.

- It takes about a minute and 1.5 GB; the scan takes about a minute more. BAM files need the extra `bam` and the
  analysis folders the extra `sc`; without one of them those files are left out and the command says so.
- The runs hold one hundredth of the reads of full-size runs. `demo/library/flowcellar.toml`, which the command
  writes, lowers the output limits by the same factor (a `[profiles]` table per kind of run), so that the runs are
  judged as full-size runs would be. `--scale 0.1` writes runs of a tenth of full size: 12 GB, some minutes.
- The same file holds the `[[inputs]]` patterns of the example: a working model of what to write for your own
  pipeline.
- Everything is invented and derived from a fixed seed: equal arguments give the same runs, names and reads.

## 2. The first scan

```
flowcellar scan /data/minknow /data/projects -L ~/flowcellar-library
```

The library folder is created when it does not exist, and the folders are remembered in
`~/flowcellar-library/flowcellar.toml`. The command prints, in this order:

1. `created library …` and what it found: runs with MinKNOW reports, folders with raw data only, and how
   many gigabytes of read tables there are to count;
2. one line per run with its verdict from the small files, newest first; this takes seconds per run;
3. `hub ready: …/hub/index.html`: from this moment the library can be opened;
4. `counting reads: …`: one pass over each sequencing summary and pore-scan table, newest runs first;
5. read data, when file patterns are configured (see [Configuration](configuration.md));
6. one line per single-cell or spatial analysis that was found, with its verdict;
7. `done: 6 runs · 5 Pass, 1 Warning, 0 Fail, 0 Incomplete · …`.

Useful to know:

- **`--quick`** stops after step 3: verdicts from the small files only. A later scan without it counts the
  reads.
- **`--time-box MIN`** (240 by default) stops *starting* new passes after that many minutes; a pass that
  is running finishes. The next scan goes on where this one stopped.
- **`--read bam`, `--read tables`, `--read analyses`** read only that kind of large read-data file in this
  scan; the others wait.
- **A second scan is cheap.** Every pass result is cached inside the library, keyed by a fingerprint of
  the input file. A scan with nothing new lists folders, opens no data file and takes seconds.
- **Ctrl-C is safe**, and so is a plain `kill`. Whatever was finished is kept, and the pages say that no scan
  is running any more.
- **A file that changed in the last two minutes is left for the next scan**, so a folder that is still
  being copied is not read half-way.
- `flowcellar discover /data/projects` lists what a scan would find and computes nothing;
  `flowcellar scan … --dry-run` says what a scan would do and writes nothing.

## 3. Open the library

**On the same computer**

```
flowcellar serve ~/flowcellar-library --open
```

or open `~/flowcellar-library/hub/index.html` in a browser. Both show the same pages: Home (what is in the
library, what needs attention, where it came from), Runs, Projects, one page per run, a comparison of runs
and the history of a flow cell.

**On a server you reach over SSH**

```
# on the server
flowcellar serve ~/flowcellar-library             # serves on 127.0.0.1:8737 of the server

# on your computer
ssh -N -L 8737:127.0.0.1:8737 user@server.example.org
# then open http://127.0.0.1:8737/ in your browser
```

If 8737 is taken on your computer, forward another local port (`-L 9000:127.0.0.1:8737`) and open
`http://127.0.0.1:9000/`. The server answers requests addressed to `127.0.0.1` or `localhost`, whatever the
port; any other name gets "403" and one line that names the address to use. An editor with remote support
does the same with its "forward a port" command.

Times on every page are shown on the clock of the computer the page is read on, and say so: the start of a
run reads "17 Apr 2025, 13:45 UTC+1", and every other time carries the same moment in UTC as its tooltip.
For a lab that sits beside its sequencer this is the clock in the name of the run folder.

The server is read-only and has no password. It listens on the server's own loopback address, so only
people who can log in to the server can forward the port. `--allow-lan` together with `--host` opens it to
the network; do that only when everyone on that network may read the library.

**Without a server at all**

Copy the folder `hub/` to your computer (`rsync -a user@server.example.org:flowcellar-library/hub/ hub/`)
and open `hub/index.html`. The page reads `index.js` and `data/<run key>.js` beside it; it makes no
request to anywhere else.

## 4. One run, no library

```
flowcellar run /data/minknow/demo_experiment/sample_a/20260302_0900_MD-000001_TST00001_0a1b2c3d -o qc/
```

writes `qc/TST00001_0a1b2c3d.fcb.json.gz` (the bundle: every number, machine-readable) and
`qc/TST00001_0a1b2c3d.flowcellar.html` (one self-contained page). A quick report appears first and is
replaced when the reads are counted; `--quick` stops after the quick report.

- Read data is given as further paths, or with `--input ROLE=PATH` and `--sample-sheet CSV`:
  `flowcellar run RUN_FOLDER basecalled/ -o qc/`.
- `--fail-on warn|fail|incomplete` makes the command exit with code 3 when the verdict is at or above
  that level, for pipelines.
- `--reused` or `--flow-cell-use 2` says that the flow cell was loaded before; without a library
  FlowCellar cannot know.
- `flowcellar --config FILE run …` judges the run with the limits of that `flowcellar.toml` (its `[profiles]` and
  `[flow_cell]` tables) and names its read data by its `[variants]` and `[samples]` tables; without it the
  built-in limits are used.
- `-L LIBRARY` also adds the run to a library.

## 5. One single-cell or spatial analysis

```
flowcellar analysis /data/projects/alpha_project/analysis/spaceranger/outs -o qc/
flowcellar analysis /data/projects/alpha_project/analysis/wf-single-cell -o qc/ --run TST00001
```

writes `qc/an_….fcb.json.gz` and `qc/an_….flowcellar.html`: the key numbers of the analysis, a spatial map, what
became of the reads, saturation, and the pipeline's own alarms in its own words. `--run` names the flow cell
run the reads came from when the files do not say. With `-L LIBRARY` the analysis is added to a library; a scan
finds analysis folders under its roots by itself.

## Paper

A run page, a project page, an analysis page, Compare and the history of a flow cell have a button **Save as PDF**. It opens the
print dialog of the browser with a layout made for paper (A4 or Letter, upright): choose the PDF printer of that
dialog ("Save as PDF" in Chrome and Edge, "Save to PDF" in Firefox). The PDF holds the page as it is shown: the
header, the verdict with its reasons, the key numbers, a line that says which samples, reads and time span were
chosen, and every section with its charts; a group that is closed on screen is left out, and the line says how
many were. On a run page **PDF with every chart** also opens the groups of further charts. Text in the PDF is
text; nothing is cut across two sheets. It was checked in Firefox and in Chromium on Linux; other browsers
were not tried.

FlowCellar itself writes no PDF file and starts no browser. For many runs at once, write the one-file reports
(`flowcellar run`, or `[hub] reports = true`) and let a Chrome, Chromium or Edge without a window print each:

```
chrome --headless --no-pdf-header-footer --virtual-time-budget=10000 --print-to-pdf=OUT.pdf REPORT.flowcellar.html
```

`chrome` stands for the program on your machine (`google-chrome`, `chromium`, `msedge`). `--virtual-time-budget`
makes the browser wait until the charts are drawn: without it the PDF holds only the text summary of the report.
The report is printed as it opens (all samples, the whole run, groups closed), on the paper size of the system.
This was tried with one Chromium on Linux.

## 6. Add more later

```
flowcellar scan -L ~/flowcellar-library                      # look for what is new under the remembered folders
flowcellar scan /archive/2025 -L ~/flowcellar-library        # add a folder and scan all of them
flowcellar run /scratch/one_run_folder -L ~/flowcellar-library       # add one run
```

Run `flowcellar scan -L ~/flowcellar-library` from a timer (cron, a systemd timer) to keep a library up to
date. Only one scan works on a library at a time: a second one waits for the first when it runs unattended,
and exits with code 5 when it was started on a terminal (`--wait` and `--no-wait` choose).

## Keep it running

On a computer with systemd, three small files open the library after every restart and keep it up to date,
with nothing to start by hand. Write them to `~/.config/systemd/user/` (`which flowcellar` prints the full
path of the command, which the files need):

```
# flowcellar.service: the pages, at http://127.0.0.1:8737/
[Unit]
Description=FlowCellar library viewer
[Service]
ExecStart=%h/flowcellar-env/bin/flowcellar serve %h/flowcellar-library
Restart=on-failure
[Install]
WantedBy=default.target
```

```
# flowcellar-scan.service: one scan of the folders the library knows
[Unit]
Description=FlowCellar scan of the library
[Service]
Type=oneshot
ExecStart=%h/flowcellar-env/bin/flowcellar scan -L %h/flowcellar-library
```

```
# flowcellar-scan.timer: when that scan runs
[Unit]
Description=FlowCellar scan, once a week
[Timer]
OnCalendar=Sun 01:00
Persistent=true
[Install]
WantedBy=timers.target
```

```
systemctl --user daemon-reload
systemctl --user enable --now flowcellar.service flowcellar-scan.timer
loginctl enable-linger "$USER"                 # start them when the computer starts, not at your first login
```

Bookmark `http://127.0.0.1:8737/`. The Home page is titled "FlowCellar - " and the name of the library
(`[library] name` in `flowcellar.toml`), so that is what the bookmark is called. `OnCalendar` takes any
schedule (`daily`, `hourly`, `Mon..Fri 19:00`). A scan reads large files for four hours at most and goes on at
the next scan (`--time-box`); with `--quick` it only looks for new runs. With cron instead of a timer:
`0 1 * * 0  $HOME/flowcellar-env/bin/flowcellar scan -L $HOME/flowcellar-library`.

## 7. When something looks wrong

```
flowcellar doctor -L ~/flowcellar-library         # the installation, the library, whether each root answers
flowcellar show TST00001_0a1b2c3d -L ~/flowcellar-library --status      # every check of one run with its sentence
flowcellar discover -L ~/flowcellar-library --explain                   # which file was linked to which run, and why
```

Each scan writes every decision and every error to one file under `~/flowcellar-library/logs/`.
