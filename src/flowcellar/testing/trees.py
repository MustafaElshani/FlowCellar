"""Trees of fake run folders under several roots, for discovery and library tests.

Everything is written with ``flowcellar.testing.fixtures.make_run``: invented names only. The trees hold the
cases a library must survive: one project under two roots, the same run under two roots, a moved
folder, a raw-only folder, basecaller output that is named like a run, a folder pruned by name, a subtree that
exhausts the listing budget, a root that became an empty mount point or disappeared.

    tree = make_three_roots(tmp_path)
    tree.roots            {"seq": ".../root_seq", "store": ".../root_store", "arch": ".../root_arch"}
    tree.runs["a1"]       FixtureRun (run_dir, run_key, expected)
"""
from __future__ import annotations

import os
import shutil
import time
from dataclasses import dataclass, field

from flowcellar.testing.fixtures import FixtureRun, RunSpec, make_run


@dataclass
class Tree:
    base: str
    roots: dict[str, str] = field(default_factory=dict)          # short name -> folder
    runs: dict[str, FixtureRun] = field(default_factory=dict)    # label -> run with MinKNOW reports
    raw_only: dict[str, str] = field(default_factory=dict)       # run key -> folder
    lookalikes: list[str] = field(default_factory=list)
    hidden: list[str] = field(default_factory=list)              # run folders the walk must NOT find

    def root_list(self) -> list[str]:
        return list(self.roots.values())


def age_tree(path: str | os.PathLike, seconds: float = 3600.0) -> None:
    """Set every file and folder under ``path`` to ``seconds`` ago: fixtures are written a moment before a test
    looks at them, and a scan leaves files that changed in the last two minutes for the next scan."""
    t = time.time() - seconds
    for d, _dirs, files in os.walk(path):
        for n in files:
            os.utime(os.path.join(d, n), (t, t))
        os.utime(d, (t, t))


def make_raw_only(parent: str | os.PathLike, flow_cell_id: str = "TST00009", run8: str = "0f0e0d0c",
                  stamp: str = "20260110_1200", position: str = "MD-000001") -> tuple[str, str]:
    """A folder with raw data only (no MinKNOW report of any kind). Returns (run key, folder)."""
    d = os.path.join(os.fspath(parent), f"{stamp}_{position}_{flow_cell_id}_{run8}")
    os.makedirs(os.path.join(d, "pod5"), exist_ok=True)
    os.makedirs(os.path.join(d, "other_reports"), exist_ok=True)
    with open(os.path.join(d, "pod5", f"{flow_cell_id}_{run8}_aaaaaaaa_0.pod5"), "w") as fh:
        fh.write("")
    return f"{flow_cell_id}_{run8}", d


def make_lookalike(parent: str | os.PathLike, flow_cell_id: str = "TST00001", run8: str = "0a0b0c0d",
                   stamp: str = "20260302_0900") -> str:
    """Basecaller output in a folder that is named like a run folder: it must never become a run."""
    d = os.path.join(os.fspath(parent), f"{stamp}_0_{flow_cell_id}_{run8}")
    os.makedirs(os.path.join(d, "bam_pass", "barcode01"), exist_ok=True)
    with open(os.path.join(d, "bam_pass", "barcode01", f"{flow_cell_id}_pass_barcode01_{run8}_aaaaaaaa_0.bam"), "w") as fh:
        fh.write("not a real BAM\n")
    return d


def make_deep_tree(parent: str | os.PathLike, folders: int = 60, name: str = "many_folders") -> str:
    """``folders`` nested and sibling folders with no run in them: more listings than a small budget allows."""
    top = os.path.join(os.fspath(parent), name)
    for i in range(folders):
        os.makedirs(os.path.join(top, f"d{i % 6}", f"e{i // 6}", "f"), exist_ok=True)
    return top


def move_run(run: FixtureRun, new_parent: str | os.PathLike) -> str:
    """Move a run folder (file times are kept, as ``mv`` and ``rsync -a`` keep them). Returns the new folder."""
    os.makedirs(new_parent, exist_ok=True)
    dest = os.path.join(os.fspath(new_parent), os.path.basename(run.run_dir))
    shutil.move(run.run_dir, dest)
    return dest


def mirror_run(run: FixtureRun, new_parent: str | os.PathLike, keep_times: bool = True) -> str:
    """Copy a run folder to a second place: the same run under two roots. Returns the copy."""
    os.makedirs(new_parent, exist_ok=True)
    dest = os.path.join(os.fspath(new_parent), os.path.basename(run.run_dir))
    shutil.copytree(run.run_dir, dest, copy_function=shutil.copy2 if keep_times else shutil.copy)
    if keep_times:
        for d, _dirs, _files in os.walk(run.run_dir):
            rel = os.path.relpath(d, run.run_dir)
            st = os.stat(d)
            os.utime(os.path.join(dest, rel) if rel != "." else dest, (st.st_atime, st.st_mtime))
    return dest


def take_offline(root: str | os.PathLike, how: str = "empty") -> str:
    """Make a root look the way a lost mount looks. ``empty``: an empty folder in its place (an automount that
    answers with nothing); ``missing``: no folder at all. Returns where the content went (for ``bring_back``)."""
    root = os.fspath(root)
    away = root.rstrip("/") + ".away"
    os.rename(root, away)
    if how == "empty":
        os.makedirs(root)
    return away


def bring_back(root: str | os.PathLike) -> None:
    root = os.fspath(root)
    away = root.rstrip("/") + ".away"
    if os.path.isdir(root):
        os.rmdir(root)
    os.rename(away, root)


def make_three_roots(out_dir: str | os.PathLike, *, reads: int = 400, hours: float = 2.0, aged: bool = True) -> Tree:
    """Three roots, five runs with reports, two projects of which one lies under two roots.

    ``root_seq`` (project = first folder)::

        alpha_project/pool_a/<a1: TST00001, 2026-03-02>
        alpha_project/pool_a_reload/<a2: TST00001 again, 2026-03-03>      the flow cell loaded a second time
        beta_project/pool_b/<b1: TST00002, 2026-03-05>
        beta_project/pool_d/<b2: TST00004, 2026-02-01>                    the oldest run
        work/hidden_project/pool_x/<a run the walk must not find: "work" is pruned by name>

    ``root_store`` (project = SECOND folder)::

        group_one/alpha_project/pool_c/<a3: TST00003, 2026-04-10>         alpha_project under a second root
        group_one/alpha_project/called_basecalled/<a folder named like a run that holds only bam_pass>

    ``root_arch`` (project = first folder)::

        <raw-only folder directly under the root: it has no project>

    So: project ``alpha_project`` lies under two roots, ``beta_project`` under one, and one run has no project.
    """
    base = os.fspath(out_dir)
    t = Tree(base=base, roots={"seq": os.path.join(base, "root_seq"), "store": os.path.join(base, "root_store"),
                               "arch": os.path.join(base, "root_arch")})
    seq, store, arch = t.roots["seq"], t.roots["store"], t.roots["arch"]

    def run(root: str, label: str, **kw) -> FixtureRun:
        r = make_run(root, RunSpec(reads=reads, hours=hours, **kw))
        t.runs[label] = r
        return r

    run(seq, "a1", flow_cell_id="TST00001", experiment="alpha_project", sample="pool_a", start="2026-03-02T09:00:00Z", seed=11)
    run(seq, "a2", flow_cell_id="TST00001", experiment="alpha_project", sample="pool_a_reload", start="2026-03-03T10:00:00Z",
        seed=12, first_scan_pore_frac=0.35)
    run(seq, "b1", flow_cell_id="TST00002", experiment="beta_project", sample="pool_b", start="2026-03-05T08:30:00Z", seed=13,
        kit="SQK-NBD114-24")
    hidden = make_run(os.path.join(seq, "work"), RunSpec(reads=50, hours=1.0, flow_cell_id="TST00008", experiment="hidden_project",
                                                         sample="pool_x", start="2026-01-05T08:00:00Z", seed=18))
    t.hidden.append(hidden.run_dir)
    run(os.path.join(store, "group_one"), "a3", flow_cell_id="TST00003", experiment="alpha_project", sample="pool_c",
        start="2026-04-10T11:00:00Z", seed=14)
    t.lookalikes.append(make_lookalike(os.path.join(store, "group_one", "alpha_project", "called_basecalled")))
    run(seq, "b2", flow_cell_id="TST00004", experiment="beta_project", sample="pool_d", start="2026-02-01T07:00:00Z", seed=15)
    key, d = make_raw_only(arch)
    t.raw_only[key] = d
    if aged:
        age_tree(base)
    return t


def stalling_child(conn, root_id: str, path: str, rules, settings, do_walk: bool, *, answer: bool = True, seconds: float = 30.0) -> None:
    """Stands in for the walking child process in tests of the root check: answers the first-level listing (or
    not) and then hangs, the way a process blocked on a dead network mount does."""
    if answer:
        conn.send(("root", sorted(n for n in os.listdir(path) if not n.startswith("."))))
    time.sleep(seconds)


def silent_child(conn, root_id: str, path: str, rules, settings, do_walk: bool) -> None:
    """A child that never answers at all (the root check must give up on it)."""
    stalling_child(conn, root_id, path, rules, settings, do_walk, answer=False)
