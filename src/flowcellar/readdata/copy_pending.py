"""Display copy of the read-data checks.

The label of each check and the words its sentence starts with come from the catalogue
(``catalogue/metrics_nanopore.yaml``: ``display.label`` and ``display.message``). What is still kept here, in one
place so that it can be found and moved: the sentence forms for several samples, the reasons a check is not
assessed, and the words for populations, assays and stages. Plain words only: no internal term appears in them.
"""
from __future__ import annotations

from flowcellar import catalogue as _catalogue

#: metric id -> how its value is written
CHECK_FORMAT: dict[str, str] = {
    "yield_pass_frac": "pct",
    "yield_basecalled_bases": "bases",
    "len_n50": "length",
    "len_median": "length",
    "len_short_read_frac": "pct",
    "qual_median": "q",
    "aln_mapped_frac": "pct",
    "aln_identity_median": "pct2",
    "aln_supplementary_frac": "pct",
    "cdna_full_length_frac": "pct",
    "cdna_rescued_frac": "pct2",
    "funnel_usable_frac": "pct",
    "funnel_basecall_loss_frac": "pct",
}


def _subject(metric_id: str) -> str:
    """The words a check's sentence starts with. They live in the catalogue: the part of ``display.message`` before
    the value, else the label."""
    d = _catalogue.display(metric_id)
    message = d.get("message") or ""
    tail = ": {value}."
    return message[:-len(tail)] if message.endswith(tail) else d.get("label") or metric_id.replace("_", " ")


#: metric id -> label (from the catalogue) and how the value is written (``fmt``)
CHECK_COPY: dict[str, dict] = {k: {"label": _subject(k), "fmt": fmt} for k, fmt in CHECK_FORMAT.items()}

SENTENCES = {
    "one": "{label}: {value}.",
    "all_samples": "{label}: {value} over all samples; {extreme} sample {worst}.",
    "one_sample_off": "{label}: {worst} in sample {names}.",
    "samples_off": "{label} is {side} the limit in {k} of {n} samples ({names}); the {extreme} is {worst}.",
    "usable": "Reads that ended up {stage}: {value} of the reads sequenced.",
}

#: why a check is not assessed
REASONS = {
    "population": "it is defined for {wanted}, and the reads on this machine are {have}",
    "no_limit": "no limit is defined for this kind of run",
    "not_counted_yet": "not every read of the linked files has been counted yet",
    "weak_link": "the files behind it are linked by a rule, not by what they contain, and the link is not confirmed",
    "mixed_sources": "the samples were read from different kinds of file",
    "two_sources": "two candidate sources describe the same reads and no rule chooses between them",
    "mapped_only": "the alignment file holds mapped reads only, so the share is 100 % by construction",
    "both_ends": "the reads given to the cDNA step were already limited to reads with a barcode at both ends, so nearly all are full-length by construction ({value} here)",
    "stopped_early": "the run did not reach its planned length; its output is judged with the run's health",
    "no_value": "the number is not known",
}

POPULATION_WORDS = {
    "basecalled_all": "every read the basecaller wrote",
    "basecalled_pass": "reads that passed the quality filter",
    "classified": "reads given to a sample",
    "full_length_trimmed": "full-length reads, trimmed",
    "aligned_primary": "aligned reads only",
}

ASSAY_WORDS = {"cdna": "cDNA", "drna": "direct RNA", "gdna": "genomic DNA"}
STAGE_WORDS = {"quantified": "counted towards a transcript", "genome_mapped": "mapped to the genome"}
