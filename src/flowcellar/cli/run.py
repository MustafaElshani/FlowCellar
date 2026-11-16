"""flowcellar run: one run folder -> a bundle and one offline HTML report."""
from __future__ import annotations

import logging
import os
import sys
import time

import click

from flowcellar.cli._options import library_option, resolve_library
from flowcellar.util import log as _log
from flowcellar.util.errors import ExitCode, FlowCellarError, InputError, NotBuiltError, UsageError

_DEVICE_WORDS = {"P2_SOLO": "PromethION P2 Solo", "P2I": "PromethION P2i", "PROMETHION": "PromethION", "P24": "PromethION 24",
                 "P48": "PromethION 48", "MINION": "MinION", "MINION_MK1B": "MinION Mk1B", "MINION_MK1C": "MinION Mk1C",
                 "MINION_MK1D": "MinION Mk1D", "GRIDION": "GridION", "FLONGLE": "Flongle"}
_ASSAY_WORDS = {"cdna": "cDNA", "drna": "direct RNA", "gdna": "genomic DNA", "unknown": "assay not known"}
_LEVELS = {"incomplete": ("incomplete", "running", "unknown", "warn", "fail"), "warn": ("warn", "fail"), "fail": ("fail",)}


def _secs(s: float) -> str:
    s = int(round(s))
    return f"{s} s" if s < 60 else f"{s // 60} min {s % 60} s" if s < 3600 else f"{s // 3600} h {s % 3600 // 60} min"


def _size(n: int) -> str:
    return f"{n / 1e9:.1f} GB" if n >= 1e9 else f"{n / 1e6:.0f} MB" if n >= 1e6 else f"{max(n // 1000, 1)} kB"


def _say(event: str, label: str, text: str, **fields) -> None:
    _log.event(f"run.{event}", f"  {label:<14} {text}" if label else text, **fields)


def _health_line(status: dict) -> str:
    """``Pass so far · 9 checks done, 6 wait for the read count`` or ``Pass · 15 checks · reads at 72 h: 76.4 M``."""
    from flowcellar.engine.card import outcome_of
    from flowcellar.model.status import VERDICT_WORD
    from flowcellar.util.text import fmt_bases, fmt_count
    health = (status.get("by_tier") or {}).get("run_health", status.get("verdict"))
    word = VERDICT_WORD.get(health, health)
    done = sum(1 for r in status.get("rules", []) if r.get("outcome") != "not_evaluated" and r.get("tier") == "run_health")
    pending = len(status.get("pending") or [])
    if status.get("provisional"):
        return f"{word} so far · {done} checks done, {pending} wait for the read count"
    parts = [word]
    if done:
        parts.append(f"{done} checks")
    out = outcome_of(status)
    if out.get("unit") == "reads":
        parts.append(f"reads at {out['horizon_h']} h: {fmt_count(out['value'])}")
    elif out.get("metric", "").startswith("yield_estimated_bases"):
        parts.append(f"estimated bases at {out['horizon_h']} h: {fmt_bases(out['value'])}")
    elif out.get("metric", "").startswith("yield_per_pore"):
        parts.append(f"output per starting pore at {out['horizon_h']} h: {out['value']:.1f} Mb")
    if health in ("incomplete", "warn", "fail", "running") and status.get("headline"):
        parts.append(status["headline"])
    return " · ".join(parts)


def _data_line(status: dict) -> str:
    from flowcellar.model.status import VERDICT_WORD
    data = (status.get("by_tier") or {}).get("read_level", "not_available")
    rules = [r for r in status.get("rules", []) if r.get("tier") == "read_level" and r.get("role") != "info"]
    n = sum(1 for r in rules if r.get("outcome") != "not_evaluated")
    so_far = any(m.get("code") == "read_data_so_far" for m in status.get("modifiers") or [])
    if so_far and data in ("pass", "warn", "fail"):
        return f"{VERDICT_WORD.get(data, data)} so far · {n} of {len(rules)} checks assessed; files of this run still wait to be read"
    return VERDICT_WORD.get(data, data) + (f" · {n} checks" if n else "")


def _first_line(run, bundle: dict) -> str:
    ident = bundle.get("identity") or {}
    fc, dev, kit, assay = ident.get("flow_cell") or {}, ident.get("device") or {}, ident.get("kit") or {}, ident.get("assay") or {}
    hours = (((bundle.get("modules") or {}).get("run") or {}).get("summary") or {}).get("run_duration_h")
    parts = [_DEVICE_WORDS.get(dev.get("type") or "", dev.get("type") or "device not known"), fc.get("product_code") or "flow cell not known",
             kit.get("sequencing_kit") or "kit not known", _ASSAY_WORDS.get(assay.get("group") or "unknown", "assay not known")]
    if run.kind == "run_raw_only":
        parts = ["raw data only"]
    elif hours is not None:
        parts.append(f"{hours:.0f} h" if hours >= 10 else f"{hours:.1f} h" if hours >= 1 else f"{round(hours * 60)} min")
    return f"{run.run_key}  " + " · ".join(parts)


class _Reading:
    """Progress of a pass: one line per quarter read, or one updated line on a terminal."""

    def __init__(self):
        self.tty = _log.is_tty() and _log.log.isEnabledFor(logging.INFO) and not any(
            type(h.formatter).__name__ == "_JsonFormatter" for h in _log.log.handlers)
        self.step: dict[str, int] = {}
        self.open_line = False

    def __call__(self, ev) -> None:
        name = os.path.basename(ev.path) if ev.path else ""
        if ev.what == "reading" and ev.fraction is not None:
            rate = f"  {ev.mb_per_s:.0f} MB/s" if ev.mb_per_s else ""
            text = f"  {'reading':<14} {name}   {ev.fraction * 100:.0f} %{rate}"
            if self.tty:
                sys.stderr.write("\r" + text + "   ")
                sys.stderr.flush()
                self.open_line = True
                return
            q = int(ev.fraction * 4)
            if q > self.step.get(name, -1):
                self.step[name] = q
                _log.event("run.reading", text, path=ev.path, fraction=round(ev.fraction, 3), mb_per_s=ev.mb_per_s)
        elif ev.what in ("error", "skipped"):
            self.close()
            _log.event("run." + ev.what, f"  {'not read':<14} {name}   {ev.message or ''}", level=logging.WARNING, path=ev.path)
        elif ev.what == "done":
            self.close()

    def close(self) -> None:
        if self.open_line:
            sys.stderr.write("\r" + " " * 100 + "\r")
            sys.stderr.flush()
            self.open_line = False


def _parse_inputs(inputs: tuple[str, ...]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for item in inputs:
        role, sep, path = item.partition("=")
        if not sep or not role or not path:
            raise UsageError(f"--input takes ROLE=PATH, got {item!r}", hint="for example --input dorado.bam=bam_pass/")
        if not os.path.exists(path):
            raise InputError(f"--input {role}={path}: no such file or folder")
        out.setdefault(role, []).append(os.path.abspath(path))
    return out


@click.command("run", short_help="One run folder: a bundle and one offline HTML report.")
@click.argument("run_dir", type=click.Path(file_okay=False))
@click.argument("tier2_paths", nargs=-1, type=click.Path())
@click.option("-o", "--outdir", type=click.Path(file_okay=False), help="Where the bundle and the report are written.")
@click.option("--prefix", metavar="TEXT", help="File name prefix. [default: the run key]")
@click.option("--input", "inputs", multiple=True, metavar="ROLE=PATH",
              help="A read-data file or folder with its role, for example dorado.bam=bam_pass/ (repeatable).")
@click.option("--sample-sheet", type=click.Path(dir_okay=False), help="Basecaller sample sheet (barcode, alias).")
@click.option("--quick", is_flag=True, help="Stop after the quick report: do not count the reads.")
@click.option("--assay", type=click.Choice(["gdna", "amplicon", "cdna_pcr", "cdna_direct", "direct_rna"]),
              help="State the assay instead of detecting it.")
@click.option("--flow-cell-use", type=click.IntRange(min=1), metavar="N", help="Which loading of the flow cell this is (1 = first use, stated).")
@click.option("--reused", is_flag=True, help="The flow cell was used before (the same as --flow-cell-use 2).")
@click.option("--flow-cell-check-pores", type=click.IntRange(min=0), metavar="N", help="Pore count of MinKNOW's flow cell check before loading.")
@click.option("--baseline-group", metavar="LABEL", help="Comparison group label for this run.")
@click.option("--no-report", is_flag=True, help="Write the bundle only.")
@click.option("--no-bundle", is_flag=True, help="Write the report only.")
@click.option("--deterministic", is_flag=True, help="Pin times and leave timings out: identical inputs give identical bytes.")
@click.option("--fail-on", type=click.Choice(["warn", "fail", "incomplete"]), help="Exit with code 3 when the verdict is at or above this level.")
@library_option()
@click.pass_context
def cmd(ctx, run_dir, tier2_paths, outdir, prefix, inputs, sample_sheet, quick, assay, flow_cell_use, reused, flow_cell_check_pores,
        baseline_group, no_report, no_bundle, deterministic, fail_on, library):
    """Read one MinKNOW run folder RUN_DIR (and read-data files under TIER2_PATHS), write
    <prefix>.fcb.json.gz and <prefix>.flowcellar.html into --outdir, and with -L add the run to a library."""
    from flowcellar.contracts import Declared, PassOptions, RunContext
    from flowcellar.engine import api as engine
    from flowcellar.model import bundle as B
    from flowcellar.model.runkey import is_safe_id

    lib_path = resolve_library(ctx, library, required=False)
    if not os.path.isdir(run_dir):
        raise InputError(f"{run_dir}: no such folder", hint="RUN_DIR is the folder MinKNOW wrote for one run")
    for p in tier2_paths:
        if not os.path.exists(p):
            raise InputError(f"{p}: no such file or folder")
    if sample_sheet and not os.path.isfile(sample_sheet):
        raise InputError(f"--sample-sheet {sample_sheet}: no such file")
    if not outdir and not lib_path:
        raise UsageError("nowhere to write: give -o DIR, or -L DIR to add the run to a library")
    if reused and flow_cell_use not in (None, 2):
        raise UsageError("--reused means --flow-cell-use 2; give one of the two")
    if no_report and no_bundle and not lib_path:
        raise UsageError("--no-report and --no-bundle together leave nothing to write")
    if prefix and not is_safe_id(prefix):
        raise UsageError(f"--prefix {prefix!r} cannot be used in a file name", hint="letters, digits, dot, dash and underscore")
    explicit = _parse_inputs(inputs)
    declared = Declared(flow_cell_use=2 if reused else flow_cell_use, flow_cell_check_pores=flow_cell_check_pores,
                        baseline_group=baseline_group)
    stated = any(v is not None for v in (declared.flow_cell_use, declared.flow_cell_check_pores, declared.baseline_group))
    obj = ctx.obj or {}
    t_start = time.monotonic()
    problems = 0

    run = engine.identify(run_dir)                                   # InputError -> exit code 4
    prefix = prefix or run.run_key
    rctx = RunContext(assay_override=assay, deterministic=deterministic)
    bundle_path = os.path.join(outdir, B.bundle_filename(prefix)) if outdir else None
    report_path = os.path.join(outdir, B.report_filename(prefix)) if outdir else None
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    def compose(measured: dict, history=(), in_library: bool = False, **kw) -> tuple[dict, dict]:
        sf = engine.judge(measured, history, declared=declared, in_library=in_library, **kw)
        if deterministic and not os.environ.get("SOURCE_DATE_EPOCH", "").isdigit():
            sf["status"]["evaluated"] = measured["provenance"]["created"]
            sf["status"]["judged"]["at"] = measured["provenance"]["created"]
        return B.compose(measured, status=sf["status"], relations=sf["relations"]), sf

    def write_report(composed: dict, label: str, t0: float) -> None:
        nonlocal problems
        if no_report or not report_path:
            return
        try:
            from flowcellar.report import api as report
            report.write_report(composed, report_path)
            _say("report", label, f"{report_path}      ({_secs(time.monotonic() - t0)})", path=report_path)
        except NotBuiltError as e:
            problems += 1
            _log.event("run.report_missing", f"  {'no report':<14} {e.message}", level=logging.WARNING)
        except FlowCellarError as e:
            problems += 1
            _log.event("run.report_missing", f"  {'no report':<14} {e.message}" + (f" ({e.hint})" if e.hint else ""), level=logging.WARNING)

    # ---------------------------------------------------------------- into a library
    if lib_path:
        from flowcellar.library import api as libapi
        lib = libapi.open_library(lib_path, create=True)
        key = libapi.add_run(lib, run.run_dir, quick=quick)
        composed = lib.composed(key)
        status = composed.get("status") or {}
        _log.event("run.start", _first_line(run, composed), run_key=key)
        _say("library", "library", f"{lib.path}   {key} added ({_secs(time.monotonic() - t_start)})", library=lib.path)
        if stated:
            _log.event("run.declared", f"  {'note':<14} what was stated about the flow cell applies to the files written with -o; "
                       "a library does not store it in this version", level=logging.WARNING)
        if outdir:
            if stated or assay:
                history = [c for c in (lib.card(k) for k in lib.run_keys() if k != key) if c]
                cfg = lib.config or {}
                composed, sf = compose(lib.read_measured(key), history, in_library=True, profile=cfg.get("profiles"),
                                       baseline_groups=(cfg.get("baselines") or {}).get("groups") or (), flow_cell=cfg.get("flow_cell"))
                status = sf["status"]
            if not no_bundle:
                B.write_bundle(bundle_path, composed, with_hash=False)
                _say("bundle", "bundle", bundle_path, path=bundle_path)
            write_report(composed, "report", t_start)
        _say("health", "Run health", _health_line(status))
        _say("data", "Read data", _data_line(status))
        return _finish(ctx, status, fail_on, problems)

    # ---------------------------------------------------------------- without a library
    # --config FILE: the limits and the flow cell settings of that file judge the run, and its [variants] and
    # [samples] tables name its read data, as they would in the library the file belongs to
    own: dict = {}
    cfg = None
    if obj.get("config"):
        if not os.path.isfile(obj["config"]):
            raise InputError(f"--config {obj['config']}: no such file")
        from flowcellar.config import loader
        cfg = loader.load(obj["config"]).config
        own = {"profile": cfg.get("profiles") or None, "flow_cell": cfg.get("flow_cell")}
    measured = engine.quick_bundle(run, rctx)
    composed, sf = compose(measured, **own)
    status = sf["status"]
    _log.event("run.start", _first_line(run, measured), run_key=run.run_key)
    t0 = time.monotonic()
    large = sum(run.sizes.get(r, 0) for r in ("minknow.sequencing_summary", "minknow.pore_scan"))
    wants_read_data = bool(tier2_paths or explicit or sample_sheet)
    final = (quick or run.kind != "run" or run.run_state == "running" or not large) and not wants_read_data
    if final and not no_bundle:
        B.write_bundle(bundle_path, composed, with_hash=False)
    write_report(composed, "report" if final else "quick report", t_start)
    line = _health_line(status)
    if status.get("provisional") and not quick and large:
        line += f" ({_size(large)} to read)"
    elif status.get("provisional") and quick:
        line += " (--quick: the reads are not counted)"
    _say("health", "Run health", line)
    if final:
        if not no_bundle:
            _say("bundle", "bundle", bundle_path, path=bundle_path)
        _say("data", "Read data", _data_line(status))
        return _finish(ctx, status, fail_on, problems)

    partials_dir = os.path.join(outdir, ".flowcellar-partials")
    reading = _Reading()
    opts = PassOptions(jobs=obj.get("jobs", 0), threads=obj.get("threads", 4), progress=reading)
    refs = [] if quick else engine.full_pass(run, partials_dir, opts)
    reading.close()
    problems += sum(1 for r in refs if r.state == "error")
    # ---- read data, through READ DATA's contract
    read_data = None
    try:
        from flowcellar.readdata import api as readdata
        cands = readdata.candidates_for_run(run, paths=[os.path.abspath(p) for p in tier2_paths], explicit=explicit or None,
                                            sample_sheet=os.path.abspath(sample_sheet) if sample_sheet else None)
        if cands:
            linked = readdata.link(cands, [run], config=cfg, probes_dir=os.path.join(partials_dir, "probes"))
            for inp in linked:
                if run.run_key in (inp.run_keys or []) and not inp.not_counted_reason:
                    ref = readdata.read_pass(inp, partials_dir, [run.run_key], opts)
                    reading.close()
                    if ref is not None and ref.state == "error":
                        problems += 1
            tier1 = engine.assemble(run, refs, None, rctx)
            read_data = readdata.read_data(run, linked, partials_dir, tier1=tier1, config=cfg)
        elif wants_read_data:
            _log.event("run.no_read_data", f"  {'Read data':<14} no read-data file was recognised in what was given", level=logging.WARNING)
            problems += 1
    except FlowCellarError as e:
        problems += 1
        _log.event("run.read_data_error", f"  {'Read data':<14} not read: {e.message}" + (f" ({e.hint})" if e.hint else ""), level=logging.WARNING)
    measured = engine.assemble(run, refs, read_data, rctx)
    composed, sf = compose(measured, **own)
    status = sf["status"]
    if not no_bundle:
        B.write_bundle(bundle_path, composed, with_hash=False)
        _say("bundle", "bundle", bundle_path, path=bundle_path)
    write_report(composed, "full report" if not quick else "report", t0)
    _say("health", "Run health", _health_line(status))
    _say("data", "Read data", _data_line(status))
    return _finish(ctx, status, fail_on, problems)


class GateTripped(FlowCellarError):
    """``--fail-on``: the verdict is at or above the given level (exit code 3)."""
    exit_code = ExitCode.GATE


class FinishedWithProblems(FlowCellarError):
    """The files were written, but an input or the report failed (exit code 5)."""
    exit_code = ExitCode.PARTIAL


def _finish(ctx, status: dict, fail_on: str | None, problems: int) -> None:
    from flowcellar.model.status import VERDICT_WORD
    if fail_on and status.get("verdict") in _LEVELS[fail_on]:
        raise GateTripped(f"the verdict is {VERDICT_WORD.get(status.get('verdict'), status.get('verdict'))} (--fail-on {fail_on})")
    if problems:
        raise FinishedWithProblems("finished, but " + ("one step" if problems == 1 else f"{problems} steps") + " failed or was left out (see the lines above)")
