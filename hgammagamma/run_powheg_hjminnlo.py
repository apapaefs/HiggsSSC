#!/usr/bin/env python3
"""Run the HJMiNNLO POWHEG signal sample in a pipeline-friendly way."""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import os
import platform
import re
import shlex
import subprocess
import time
from pathlib import Path
from typing import Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_POWHEG_DIR = REPO_ROOT / "POWHEG-BOX-V2" / "HJ" / "HJMiNNLO"
DEFAULT_CARD = SCRIPT_DIR / "HOAnalysis" / "powheg-hjminnlo-ssc40-nnpdf40nnloqed.input"
DEFAULT_HERWIG_ENV = Path.home() / "Projects/Herwig/Herwig-REAL-stable-gcc-full/bin/activate"
DEFAULT_LINUX_HERWIG_MODULE = "herwig/stable"
DEFAULT_DARWIN_HERWIG_MODULE = "herwig/730"
DEFAULT_JOBS = min(os.cpu_count() or 1, 4)
MERGED_LHE_NAME = "powheg-hjminnlo-merged.lhe"
LHE_MANIFEST_NAME = "powheg-lhe-files.txt"
SEED_LHE_MANIFEST_NAME = "powheg-lhe-seed-files.txt"
SUMMARY_NAME = "powheg-run-summary.txt"
LHE_EVENT_START_RE = re.compile(r"^\s*<event\b", re.IGNORECASE)
LHE_CLOSE_RE = re.compile(r"^\s*</LesHouchesEvents\s*>", re.IGNORECASE)
LHE_SEED_RE = re.compile(r"^pwgevents-(?P<seed>\d+)\.lhe$")
RUN_LOG_RE = re.compile(r"^run-(?P<label>.+)-(?P<seed>\d+)\.log$")
ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
STAGE_LABELS = ("st1-xg1", "st1-xg2", "st2", "st3", "st4")
STAGE_DONE_MARKERS = {
    "st1-xg1": ("Importance sampling x grids generated and stored",),
    "st1-xg2": ("Importance sampling x grids generated and stored",),
    "st2": ("negative weight fraction:", "total (btilde+remnants) cross section"),
    "st3": ("Normalization of upper bounding function for radiation computed and stored",),
}


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return parsed


def default_cxx() -> str:
    if platform.system() == "Darwin" and Path("/opt/homebrew/bin/g++-15").exists():
        return "/opt/homebrew/bin/g++-15"
    return os.environ.get("CXX", "g++")


def default_cc() -> str:
    if platform.system() == "Darwin" and Path("/opt/homebrew/bin/gcc-15").exists():
        return "/opt/homebrew/bin/gcc-15"
    return os.environ.get("CC", "gcc")


def maybe_path(value: str | None) -> Path | None:
    if not value:
        return None
    return Path(value).expanduser()


def clean_optional_text(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    return value or None


def env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.lower() in {"1", "true", "yes", "on"}


def normalize_herwig_env(path: Path | None) -> Path | None:
    if path is None:
        return None
    expanded = path.expanduser()
    if expanded.is_dir():
        activate = expanded / "bin" / "activate"
        if activate.exists():
            return activate
    return expanded


def q(value: object) -> str:
    return shlex.quote(str(value))


def clean_shell_env(env: dict[str, str] | None = None) -> dict[str, str]:
    shell_env = dict(os.environ if env is None else env)
    for key in list(shell_env):
        if key.startswith("BASH_FUNC_") and key.endswith("%%"):
            shell_env.pop(key, None)
    return shell_env


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nevents",
        type=positive_int,
        required=True,
        help="total requested POWHEG LHE events across all seed jobs",
    )
    parser.add_argument(
        "--jobs",
        type=positive_int,
        default=positive_int(os.environ.get("POWHEG_JOBS", str(DEFAULT_JOBS))),
        help="number of POWHEG seed jobs to run in parallel",
    )
    parser.add_argument(
        "--ebeam",
        type=positive_float,
        default=positive_float(os.environ.get("POWHEG_EBEAM", "20000")),
        help="beam energy in GeV; default 20000 gives pp collisions at 40 TeV",
    )
    parser.add_argument("--run-dir", type=Path, default=None, help="run directory to create")
    parser.add_argument("--powheg-dir", type=Path, default=maybe_path(os.environ.get("POWHEG_DIR")) or DEFAULT_POWHEG_DIR)
    parser.add_argument("--card", type=Path, default=maybe_path(os.environ.get("POWHEG_CARD")) or DEFAULT_CARD)
    parser.add_argument("--herwig-env", type=Path, default=maybe_path(os.environ.get("HERWIG_ENV")))
    parser.add_argument("--no-herwig-env", action="store_true", help="do not source a Herwig/LHAPDF environment")
    parser.add_argument(
        "--herwig-module",
        default=os.environ.get("HERWIG_MODULE"),
        help=(
            "environment module to load before POWHEG, e.g. herwig/730 on the laptop "
            "or herwig/stable on timur"
        ),
    )
    parser.add_argument(
        "--no-herwig-module",
        action="store_true",
        default=env_bool("NO_HERWIG_MODULE", False),
        help="disable the automatic platform Herwig module default",
    )
    parser.add_argument("--skip-build", action="store_true", help="require an existing pwhg_main executable")
    parser.add_argument(
        "--status",
        action="store_true",
        help="show progress for the selected run directory and exit without running POWHEG",
    )
    parser.add_argument(
        "--watch-status",
        type=nonnegative_float,
        default=None,
        metavar="SECONDS",
        help="repeat --status every SECONDS until the run is complete; use Ctrl-C to stop",
    )
    parser.add_argument(
        "--merged-lhe-name",
        default=os.environ.get("POWHEG_MERGED_LHE", MERGED_LHE_NAME),
        help="name of the merged LHE file written in the run directory",
    )
    parser.add_argument(
        "--no-merge-lhe",
        action="store_true",
        default=env_bool("POWHEG_NO_MERGE_LHE", False),
        help="leave per-seed LHE files unmerged and list them in the main manifest",
    )
    parser.add_argument("--cxx", default=os.environ.get("POWHEG_CXX", default_cxx()))
    parser.add_argument("--cc", default=os.environ.get("POWHEG_CC", default_cc()))
    parser.add_argument("--stdclib", default=os.environ.get("POWHEG_STDCLIB", "-lstdc++"))
    parser.add_argument(
        "--resume",
        action="store_true",
        help="allow an existing run directory; useful after a failed integration stage",
    )
    parser.add_argument("--dry-run", action="store_true", help="print the planned stages without running POWHEG")
    return parser


def run_dir_for(powheg_dir: Path, nevents: int) -> Path:
    return powheg_dir / f"run-ssc40-hjminnlo-nnpdf40nnloqed-{nevents}ev"


def module_setup_snippet() -> str:
    return (
        "if ! type module >/dev/null 2>&1; then "
        "for module_init in "
        "/opt/homebrew/opt/modules/init/bash "
        "/etc/profile.d/modules.sh "
        "/usr/share/Modules/init/bash "
        "/usr/share/lmod/lmod/init/bash; do "
        "[ -r \"$module_init\" ] && source \"$module_init\" && break; "
        "done; "
        "fi"
    )


def runtime_prefix(herwig_env: Path | None, herwig_module: str | None) -> list[str]:
    commands = ["set -e"]
    if herwig_env:
        commands.append(f"source {q(herwig_env)}")
    if herwig_module:
        commands.append(module_setup_snippet())
        commands.append(f"module load {q(herwig_module)}")
    return commands


def shell_command(command: str, herwig_env: Path | None, herwig_module: str | None) -> str:
    return "; ".join([*runtime_prefix(herwig_env, herwig_module), command])


def run_shell(command: str, cwd: Path, herwig_env: Path | None, herwig_module: str | None, dry_run: bool) -> None:
    full = shell_command(command, herwig_env, herwig_module)
    print(f"+ cd {q(cwd)} && {full}")
    if dry_run:
        return
    subprocess.run(["bash", "-lc", full], cwd=cwd, env=clean_shell_env(), check=True)


def tail_text(path: Path, nlines: int = 40) -> str:
    if not path.exists():
        return "<log file was not created>"
    lines = path.read_text(errors="replace").splitlines()
    return "\n".join(lines[-nlines:]) if lines else "<log file is empty>"


def count_lhe_events(path: Path) -> int:
    count = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            if LHE_EVENT_START_RE.match(line):
                count += 1
    return count


def seed_entries_in_file(run_dir: Path) -> int | None:
    seed_file = run_dir / "pwgseeds.dat"
    if not seed_file.exists():
        return None
    seeds = [line for line in seed_file.read_text(errors="replace").splitlines() if line.strip()]
    return len(seeds) or None


def log_index(run_dir: Path) -> dict[tuple[str, int], Path]:
    logs: dict[tuple[str, int], Path] = {}
    for path in run_dir.glob("run-*.log"):
        match = RUN_LOG_RE.match(path.name)
        if not match:
            continue
        label = match.group("label")
        if label not in STAGE_LABELS:
            continue
        logs[(label, int(match.group("seed")))] = path
    return logs


def seed_lhe_path(run_dir: Path, seed: int) -> Path:
    return run_dir / f"pwgevents-{seed:04d}.lhe"


def seed_from_lhe_path(path: Path) -> int | None:
    match = LHE_SEED_RE.match(path.name)
    return int(match.group("seed")) if match else None


def log_has_failure(path: Path) -> bool:
    if not path.exists():
        return False
    text = path.read_text(errors="replace")
    failure_markers = (
        "Abort",
        "Cannot open",
        "Fatal Error",
        "Segmentation fault",
        "Traceback",
        "error termination",
        "returned non-zero",
    )
    return any(marker in text for marker in failure_markers)


def log_is_complete(label: str, seed: int, path: Path, run_dir: Path) -> bool:
    if label == "st4":
        return seed_lhe_path(run_dir, seed).exists()
    markers = STAGE_DONE_MARKERS.get(label, ())
    if not markers or not path.exists():
        return False
    text = path.read_text(errors="replace")
    return all(marker in text for marker in markers)


def age_text(path: Path) -> str:
    age = max(0, time.time() - path.stat().st_mtime)
    if age < 60:
        return f"{age:.0f}s ago"
    if age < 3600:
        return f"{age / 60:.1f}m ago"
    return f"{age / 3600:.1f}h ago"


def last_meaningful_line(path: Path) -> str:
    lines = tail_text(path, 80).splitlines()
    for line in reversed(lines):
        line = ANSI_RE.sub("", line).strip()
        if line:
            return line
    return "<empty log>"


def timing_labels(run_dir: Path) -> list[str]:
    timings = run_dir / "Timings.txt"
    if not timings.exists():
        return []
    labels: list[str] = []
    for line in timings.read_text(errors="replace").splitlines():
        if line.strip():
            labels.append(line.split(maxsplit=1)[0])
    return labels


def print_run_status(
    run_dir: Path,
    requested_events: int,
    groups: list[tuple[int, list[int]]],
    merged_lhe_name: str,
) -> bool:
    print(f"POWHEG run directory: {run_dir}")
    if not run_dir.exists():
        print("Status: run directory has not been created yet.")
        return False

    expected_seeds = active_seed_count(groups)
    event_targets = seed_event_targets(groups)
    seed_file_entries = seed_entries_in_file(run_dir)
    labels = timing_labels(run_dir)
    print(f"Requested events: {requested_events}")
    print(f"Seed jobs: {expected_seeds}")
    if seed_file_entries is not None and seed_file_entries != expected_seeds:
        print(f"Seed file entries: {seed_file_entries} (using --jobs-derived expected seed count)")
    print(f"Timing checkpoints: {', '.join(labels) if labels else '<none yet>'}")

    logs = log_index(run_dir)
    seed_lhe_files = sorted(run_dir.glob("pwgevents-*.lhe"))
    seed_lhe_counts: dict[int, int] = {}
    for path in seed_lhe_files:
        seed = seed_from_lhe_path(path)
        if seed is not None:
            seed_lhe_counts[seed] = count_lhe_events(path)

    print("\nStage progress:")
    for label in STAGE_LABELS:
        stage_logs = {seed: path for (stage, seed), path in logs.items() if stage == label}
        failed = sum(log_has_failure(path) for path in stage_logs.values())
        if label == "st4":
            complete = sum(seed_lhe_counts.get(seed, 0) >= target for seed, target in event_targets.items())
        else:
            complete = sum(log_is_complete(label, seed, path, run_dir) for seed, path in stage_logs.items())
        missing = max(0, expected_seeds - len(stage_logs))
        running = max(0, len(stage_logs) - complete - failed)
        print(
            f"  {label:7} complete {complete:>2}/{expected_seeds:<2}  "
            f"running/unknown {running:>2}  missing {missing:>2}  failed {failed:>2}"
        )

    seed_lhe_events = sum(seed_lhe_counts.values())
    complete_seed_lhes = sum(seed_lhe_counts.get(seed, 0) >= target for seed, target in event_targets.items())
    print("\nLHE output:")
    print(f"  seed LHE files: {len(seed_lhe_files)}/{expected_seeds} ({complete_seed_lhes} complete)")
    print(f"  seed LHE events: {seed_lhe_events}/{requested_events}")
    if seed_lhe_counts and expected_seeds <= 16:
        per_seed = ", ".join(
            f"{seed}:{seed_lhe_counts.get(seed, 0)}/{event_targets[seed]}" for seed in sorted(event_targets)
        )
        print(f"  per-seed events: {per_seed}")
    merged_lhe = resolve_merged_lhe(run_dir, merged_lhe_name)
    if merged_lhe.exists():
        print(f"  merged LHE: {merged_lhe} ({count_lhe_events(merged_lhe)} events)")
    else:
        print(f"  merged LHE: {merged_lhe} <not written yet>")
    manifest = run_dir / LHE_MANIFEST_NAME
    summary = run_dir / SUMMARY_NAME
    print(f"  manifest: {manifest if manifest.exists() else '<not written yet>'}")
    print(f"  summary: {summary if summary.exists() else '<not written yet>'}")

    latest_lhes = sorted(seed_lhe_files, key=lambda path: path.stat().st_mtime, reverse=True)[:4]
    if latest_lhes:
        print("\nLatest LHE files:")
        for path in latest_lhes:
            seed = seed_from_lhe_path(path)
            count = seed_lhe_counts.get(seed or -1, 0)
            target = event_targets.get(seed or -1, "?")
            print(f"  {path.name} ({age_text(path)}): {count}/{target} events")

    latest_logs = sorted(logs.values(), key=lambda path: path.stat().st_mtime, reverse=True)[:4]
    if latest_logs:
        print("\nLatest logs:")
        for path in latest_logs:
            print(f"  {path.name} ({age_text(path)}): {last_meaningful_line(path)}")

    return "end" in labels and manifest.exists()


def ensure_build(
    powheg_dir: Path,
    args: argparse.Namespace,
    herwig_env: Path | None,
    herwig_module: str | None,
) -> Path:
    exe = powheg_dir / "pwhg_main"
    if exe.exists() and not args.skip_build:
        return exe
    if exe.exists() and args.skip_build:
        return exe
    if args.skip_build:
        raise SystemExit(f"missing POWHEG executable: {exe}")

    make_cmd = " ".join(
        [
            "make",
            "pwhg_main",
            f"CXX={q(args.cxx)}",
            f"CC={q(args.cc)}",
            f"STDCLIB={q(args.stdclib)}",
        ]
    )
    run_shell(make_cmd, powheg_dir, herwig_env, herwig_module, args.dry_run)
    return exe


def set_powheg_value(text: str, key: str, value: object) -> str:
    key_re = re.escape(key)
    pattern = re.compile(rf"^(\s*{key_re}\s+)(\S+)(.*)$")
    out: list[str] = []
    replaced = False
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith(("!", "#")):
            out.append(line)
            continue
        match = pattern.match(line)
        if match:
            prefix, _old, tail = match.groups()
            out.append(f"{prefix}{value}{tail}")
            replaced = True
        else:
            out.append(line)
    if not replaced:
        out.append(f"{key} {value}")
    return "\n".join(out) + "\n"


def patched_input(base_text: str, updates: dict[str, object]) -> str:
    text = base_text
    for key, value in updates.items():
        text = set_powheg_value(text, key, value)
    return text


def event_groups(nevents: int, jobs: int) -> list[tuple[int, list[int]]]:
    active_jobs = min(nevents, jobs)
    base = nevents // active_jobs
    remainder = nevents % active_jobs
    groups: list[tuple[int, list[int]]] = []
    if remainder:
        groups.append((base + 1, list(range(1, remainder + 1))))
    rest = list(range(remainder + 1, active_jobs + 1))
    if rest:
        groups.append((base, rest))
    return groups


def active_seed_count(groups: Iterable[tuple[int, list[int]]]) -> int:
    return sum(len(seeds) for _events, seeds in groups)


def seed_event_targets(groups: Iterable[tuple[int, list[int]]]) -> dict[int, int]:
    return {seed: events for events, seeds in groups for seed in seeds}


def copy_or_make_seeds(powheg_dir: Path, run_dir: Path, nseeds: int, dry_run: bool) -> None:
    source = powheg_dir / "suggested_run" / "pwgseeds.dat"
    target = run_dir / "pwgseeds.dat"
    if dry_run:
        print(f"+ write {target} with {nseeds} seed entries")
        return

    seeds: list[str] = []
    if source.exists():
        for line in source.read_text().splitlines():
            line = line.strip()
            if line:
                seeds.append(line)
            if len(seeds) == nseeds:
                break
    if len(seeds) < nseeds:
        start = 31122002
        seeds.extend(str(start + 104729 * idx) for idx in range(len(seeds), nseeds))
    target.write_text("\n".join(seeds[:nseeds]) + "\n")


def write_input(run_dir: Path, base_text: str, updates: dict[str, object], dry_run: bool) -> None:
    path = run_dir / "powheg.input"
    if dry_run:
        update_text = ", ".join(f"{key}={value}" for key, value in updates.items())
        print(f"+ write {path} ({update_text})")
        return
    path.write_text(patched_input(base_text, updates))


def append_timing(run_dir: Path, label: str, dry_run: bool) -> None:
    stamp = dt.datetime.now().astimezone().strftime("%a %b %d %H:%M:%S %Z %Y")
    if dry_run:
        print(f"+ timing {label} {stamp}")
        return
    with (run_dir / "Timings.txt").open("a") as handle:
        handle.write(f"{label} {stamp}\n")


def run_seed(
    exe: Path,
    run_dir: Path,
    seed_index: int,
    log_name: str,
    herwig_env: Path | None,
    herwig_module: str | None,
    dry_run: bool,
) -> None:
    log_path = run_dir / log_name
    command = shell_command(q(exe), herwig_env, herwig_module)
    if dry_run:
        print(f"+ cd {q(run_dir)} && printf '%s\\n' {seed_index} | bash -lc {q(command)} > {log_path} 2>&1")
        return

    with log_path.open("w") as log:
        try:
            subprocess.run(
                ["bash", "-lc", command],
                cwd=run_dir,
                env=clean_shell_env(),
                input=f"{seed_index}\n",
                text=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"{log_name} failed with exit code {exc.returncode}. "
                f"See {log_path}.\nLast log lines:\n{tail_text(log_path)}"
            ) from exc


def run_seed_group(
    exe: Path,
    run_dir: Path,
    seeds: list[int],
    label: str,
    herwig_env: Path | None,
    herwig_module: str | None,
    dry_run: bool,
) -> None:
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(seeds)) as pool:
        futures = {}
        for seed in seeds:
            log_name = f"run-{label}-{seed}.log"
            future = pool.submit(run_seed, exe, run_dir, seed, log_name, herwig_env, herwig_module, dry_run)
            futures[future] = (seed, log_name)
        completed = 0
        for future in concurrent.futures.as_completed(futures):
            seed, log_name = futures[future]
            future.result()
            completed += 1
            if not dry_run:
                print(f"[{label}] completed seed {seed} ({completed}/{len(seeds)}); log {run_dir / log_name}", flush=True)


def stage_updates(
    stage: int,
    nseeds: int,
    numevts: int,
    *,
    xgriditeration: int | None = None,
    use_old_grid: int,
    use_old_ubound: int,
) -> dict[str, object]:
    updates: dict[str, object] = {
        "numevts": numevts,
        "manyseeds": 1,
        "maxseeds": nseeds,
        "parallelstage": stage,
        "use-old-grid": use_old_grid,
        "use-old-ubound": use_old_ubound,
    }
    if xgriditeration is not None:
        updates["xgriditeration"] = xgriditeration
    return updates


def resolve_merged_lhe(run_dir: Path, merged_lhe_name: str) -> Path:
    path = Path(merged_lhe_name).expanduser()
    return path if path.is_absolute() else run_dir / path


def merge_lhe_files(lhe_files: list[Path], output: Path, dry_run: bool) -> int:
    if dry_run:
        sources = " ".join(q(path) for path in lhe_files) or "<stage-4 LHE files>"
        print(f"+ merge {sources} -> {output}")
        return 0
    if not lhe_files:
        raise SystemExit("no per-seed LHE files found to merge")

    event_count = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as out:
        wrote_header = False
        for source in lhe_files:
            in_events = False
            saw_event = False
            with source.open(errors="replace") as handle:
                for line in handle:
                    if LHE_CLOSE_RE.match(line):
                        break
                    if not in_events:
                        if LHE_EVENT_START_RE.match(line):
                            in_events = True
                            saw_event = True
                            event_count += 1
                            out.write(line)
                        elif not wrote_header:
                            out.write(line)
                        continue
                    if LHE_EVENT_START_RE.match(line):
                        event_count += 1
                    out.write(line)
            if not saw_event:
                raise SystemExit(f"no <event> blocks found in {source}")
            wrote_header = True
        out.write("</LesHouchesEvents>\n")
    return event_count


def write_manifest(
    run_dir: Path,
    requested_events: int,
    groups: list[tuple[int, list[int]]],
    ebeam: float,
    merged_lhe_name: str,
    merge_lhe: bool,
    dry_run: bool,
) -> None:
    manifest = run_dir / LHE_MANIFEST_NAME
    seed_manifest = run_dir / SEED_LHE_MANIFEST_NAME
    summary = run_dir / SUMMARY_NAME
    merged_lhe = resolve_merged_lhe(run_dir, merged_lhe_name)
    if dry_run:
        if merge_lhe:
            merge_lhe_files([], merged_lhe, dry_run)
        print(f"+ write {seed_manifest}")
        print(f"+ write {manifest}")
        print(f"+ write {summary}")
        return

    lhe_files = sorted(path for path in run_dir.glob("pwgevents-*.lhe") if path.resolve() != merged_lhe.resolve())
    if len(lhe_files) != active_seed_count(groups):
        raise SystemExit(f"expected {active_seed_count(groups)} LHE files, found {len(lhe_files)}; check run-st4-*.log")

    seed_manifest.write_text("".join(f"{path}\n" for path in lhe_files))
    merged_events: int | None = None
    if merge_lhe:
        merged_events = merge_lhe_files(lhe_files, merged_lhe, dry_run)
        manifest.write_text(f"{merged_lhe}\n")
    else:
        manifest.write_text(seed_manifest.read_text())

    distribution = ", ".join(f"{len(seeds)} seed(s) x {events} events" for events, seeds in groups)
    summary_lines = [
        f"requested_events: {requested_events}",
        f"seed_jobs: {active_seed_count(groups)}",
        f"event_distribution: {distribution}",
        f"ebeam_gev: {ebeam:g}",
        f"sqrt_s_gev: {2 * ebeam:g}",
        f"seed_lhe_manifest: {seed_manifest}",
    ]
    if merge_lhe:
        summary_lines.extend(
            [
                f"merged_lhe: {merged_lhe}",
                f"merged_events: {merged_events}",
            ]
        )
    summary_lines.extend([f"lhe_manifest: {manifest}", ""])
    summary.write_text("\n".join(summary_lines))


def main() -> None:
    args = build_parser().parse_args()
    powheg_dir = args.powheg_dir.expanduser().resolve()
    card = args.card.expanduser().resolve()
    groups = event_groups(args.nevents, args.jobs)
    nseeds = active_seed_count(groups)
    max_events_per_seed = max(events for events, _seeds in groups)
    run_dir = (args.run_dir.expanduser() if args.run_dir else run_dir_for(powheg_dir, args.nevents)).resolve()

    if args.watch_status is not None:
        interval = args.watch_status
        while True:
            done = print_run_status(run_dir, args.nevents, groups, args.merged_lhe_name)
            if done or interval == 0:
                return
            print(f"\nRefreshing in {interval:g} seconds. Press Ctrl-C to stop.\n")
            time.sleep(interval)
    if args.status:
        print_run_status(run_dir, args.nevents, groups, args.merged_lhe_name)
        return

    raw_herwig_env = None if args.no_herwig_env or args.herwig_env is None else args.herwig_env
    herwig_env = normalize_herwig_env(raw_herwig_env)
    herwig_module = None if args.no_herwig_module else clean_optional_text(args.herwig_module)
    if (
        herwig_env is None
        and herwig_module is None
        and not args.no_herwig_module
    ):
        system = platform.system()
        if system == "Linux":
            herwig_module = DEFAULT_LINUX_HERWIG_MODULE
        elif system == "Darwin":
            herwig_module = DEFAULT_DARWIN_HERWIG_MODULE
    if herwig_env is None and herwig_module is None and DEFAULT_HERWIG_ENV.exists() and not args.no_herwig_env:
        herwig_env = DEFAULT_HERWIG_ENV
    if herwig_env and not herwig_env.exists():
        raise SystemExit(
            f"missing Herwig/LHAPDF activation script: {herwig_env}. "
            "Pass --herwig-env /path/to/bin/activate or a prefix containing bin/activate."
        )
    if herwig_env and herwig_env.is_dir():
        raise SystemExit(
            f"Herwig/LHAPDF environment points to a directory without bin/activate: {herwig_env}. "
            "Pass --herwig-env /path/to/bin/activate or --no-herwig-env."
        )
    if not powheg_dir.exists():
        raise SystemExit(f"missing POWHEG HJMiNNLO directory: {powheg_dir}")
    if not card.exists():
        raise SystemExit(f"missing POWHEG input card: {card}")

    if run_dir.exists() and not args.resume:
        raise SystemExit(f"run directory already exists: {run_dir}; pass --run-dir or --resume")
    if args.dry_run:
        print(f"+ mkdir -p {run_dir}")
    else:
        run_dir.mkdir(parents=True, exist_ok=True)

    exe = ensure_build(powheg_dir, args, herwig_env, herwig_module)
    base_text = patched_input(
        card.read_text(),
        {
            "numevts": max_events_per_seed,
            "ebeam1": f"{args.ebeam:g}d0",
            "ebeam2": f"{args.ebeam:g}d0",
            "manyseeds": 1,
            "maxseeds": nseeds,
        },
    )
    if args.dry_run:
        print(f"+ write {run_dir / 'powheg.input-save'}")
    else:
        (run_dir / "powheg.input-save").write_text(base_text)
        (run_dir / "Timings.txt").write_text("")
    copy_or_make_seeds(powheg_dir, run_dir, nseeds, args.dry_run)

    print(f"Requested {args.nevents} LHE events over {nseeds} POWHEG seed job(s).")
    print(f"Beam setup: pp at sqrt(s) = {2 * args.ebeam:g} GeV ({args.ebeam:g} GeV per beam).")
    for events, seeds in groups:
        print(f"  seeds {seeds[0]}..{seeds[-1]}: {events} events each")

    for igrid in (1, 2):
        label = f"st1-xg{igrid}"
        append_timing(run_dir, label, args.dry_run)
        write_input(
            run_dir,
            base_text,
            stage_updates(1, nseeds, max_events_per_seed, xgriditeration=igrid, use_old_grid=0, use_old_ubound=0),
            args.dry_run,
        )
        run_seed_group(exe, run_dir, list(range(1, nseeds + 1)), label, herwig_env, herwig_module, args.dry_run)

    append_timing(run_dir, "st2", args.dry_run)
    write_input(run_dir, base_text, stage_updates(2, nseeds, max_events_per_seed, use_old_grid=0, use_old_ubound=0), args.dry_run)
    run_seed_group(exe, run_dir, list(range(1, nseeds + 1)), "st2", herwig_env, herwig_module, args.dry_run)

    append_timing(run_dir, "st3", args.dry_run)
    write_input(run_dir, base_text, stage_updates(3, nseeds, max_events_per_seed, use_old_grid=1, use_old_ubound=0), args.dry_run)
    run_seed_group(exe, run_dir, list(range(1, nseeds + 1)), "st3", herwig_env, herwig_module, args.dry_run)

    append_timing(run_dir, "st4", args.dry_run)
    for events, seeds in groups:
        write_input(run_dir, base_text, stage_updates(4, nseeds, events, use_old_grid=1, use_old_ubound=1), args.dry_run)
        run_seed_group(exe, run_dir, seeds, "st4", herwig_env, herwig_module, args.dry_run)
    append_timing(run_dir, "end", args.dry_run)

    write_manifest(
        run_dir,
        args.nevents,
        groups,
        args.ebeam,
        args.merged_lhe_name,
        not args.no_merge_lhe,
        args.dry_run,
    )
    print(f"POWHEG run directory: {run_dir}")
    if not args.no_merge_lhe:
        print(f"Merged LHE: {resolve_merged_lhe(run_dir, args.merged_lhe_name)}")
    print(f"LHE manifest: {run_dir / LHE_MANIFEST_NAME}")


if __name__ == "__main__":
    main()
