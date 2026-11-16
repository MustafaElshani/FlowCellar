"""The single-cell and spatial add-on: quality control of wf-single-cell and Space Ranger outputs.

An *analysis* is one output folder of one pipeline for one sample: a wf-single-cell sample folder, or a Space
Ranger ``outs/``. It gets a bundle of its own (``kind: "analysis"``), because one analysis pools the reads of
several flow cells; the bundle names the runs it came from with the evidence for each link. An analysis has its
own verdict and never changes the verdict of a flow cell run.

Nothing of this package is imported until an analysis is looked for or shown: ``import flowcellar``, the run
commands and a run report do not load it. Importing this package itself loads nothing heavy either: ``h5py``
(the ``sc`` extra) is imported inside the functions that open an HDF5 file, ``pysam`` (the ``bam`` extra) inside
the function that tallies the read tags of a BAM. Without ``h5py`` the files that are tables or JSON are still
read and the bundle says what is missing.

    flowcellar.sc.detect       which folders are analyses (by their files, never by their names)
    flowcellar.sc.collect      one analysis folder -> the measured bundle
    flowcellar.sc.link         the evidence ladder from an analysis to the flow cell runs it came from
    flowcellar.sc.judge        the limits and the verdict of an analysis
    flowcellar.sc.store        analyses inside a library, the index table, what a scan does
    flowcellar.sc.report       the one-file report and the hub's data file of an analysis
"""
from __future__ import annotations

import re

MODULE_ID = "sc"
#: version of the section this package writes; a change makes every cached collection stale
MODULE_VERSION = "0.1.0"
SECTION_SCHEMA_VERSION = 1

PIPELINE_SPACERANGER = "spaceranger"
PIPELINE_WF = "wf-single-cell"
#: how a pipeline is spelled inside a bundle id (letters only: the id is a file name and an address)
PIPELINE_CODE = {PIPELINE_SPACERANGER: "spaceranger", PIPELINE_WF: "wfsinglecell"}
PIPELINE_WORD = {PIPELINE_SPACERANGER: "Space Ranger", PIPELINE_WF: "wf-single-cell"}

#: the id of an analysis bundle as it may appear in a file name, an element id and an address
ANALYSIS_ID_RE = re.compile(r"^an_[a-z]{2,16}_[0-9a-f]{12}$")

#: the roles under which the discovery walk reports the anchor file of an analysis
ROLE_SPACERANGER = "analysis.spaceranger"
ROLE_WF = "analysis.wf_single_cell"
ROLE_PREFIX = "analysis."

ASSAY_WORD = {"visium": "Visium", "visium_hd": "Visium HD", "singlecell": "Single cell", None: "Assay not known"}

#: size budget of one analysis section as JSON text (bytes): the default aim and the cap a bundle may not pass
SECTION_BUDGET_BYTES = 150_000
SECTION_CAP_BYTES = 400_000


def bundle_id(pipeline: str, analysis_hash: str) -> str:
    """``an_<pipeline>_<first 12 hex of the content hash>``."""
    return f"an_{PIPELINE_CODE[pipeline]}_{analysis_hash[:12]}"


def is_analysis_id(text) -> bool:
    return isinstance(text, str) and bool(ANALYSIS_ID_RE.match(text))
