#!/usr/bin/env python3
"""Attach shared four-lepton events and run analysis without event generation."""
from __future__ import annotations

import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

from hfourlepton.shared_events import DEFAULT_DATASET, digest, load_dataset, name, safe_path, save

REPO = Path(__file__).resolve().parent
DEFAULT_WORKSPACE = REPO / "hfourlepton/student_runs/ssc40"


def attach(dataset, workspace):
    dataset, workspace = Path(dataset).resolve(), Path(workspace).resolve()
    if workspace.is_relative_to(dataset) or dataset.is_relative_to(workspace):
        raise ValueError("student workspace and shared dataset must be separate directories")
    index = load_dataset(dataset)
    binding = dict(schema_version=1, dataset=str(dataset), dataset_sha256=digest(dataset / "dataset.json"))
    record = workspace / "shared_events.json"
    if record.exists():
        if json.loads(record.read_text()) != binding:
            raise ValueError("workspace is attached to another dataset; choose a new --workspace")
        print(f"Already attached: {workspace}")
        return
    if workspace.exists() and any(workspace.iterdir()):
        raise ValueError("workspace is nonempty; choose a new --workspace")
    workspace.mkdir(parents=True, exist_ok=True)
    # Only the small analysis cards are copied. Events stay in the read-only share.
    for relative in index["cards"].values():
        target = safe_path(workspace, relative)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(safe_path(dataset, relative), target)
    (workspace / "events").symlink_to(dataset, target_is_directory=True)
    save(record, binding)
    print(f"Attached {len(index['samples'])} samples from {dataset}\nWorkspace: {workspace}\n"
          f"Editable cards: {workspace / 'cards'}\nNo events were generated or copied into your workspace.")


def bound(workspace, *, full=False):
    workspace = Path(workspace).resolve()
    path = workspace / "shared_events.json"
    if not path.is_file():
        raise ValueError("workspace is not attached; run student_events.py attach first")
    binding = json.loads(path.read_text())
    dataset = Path(binding["dataset"]).resolve()
    if digest(dataset / "dataset.json") != binding["dataset_sha256"]:
        raise ValueError("the attached dataset release changed; use a new workspace")
    return dataset, load_dataset(dataset, full=full)


def environment(workspace):
    env = os.environ.copy()
    env.update(PYTHONNOUSERSITE="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               MPLCONFIGDIR=str(Path(workspace) / ".matplotlib"))
    return env


def setup(workspace):
    bound(workspace)
    python = Path(workspace) / ".venv/bin/python"
    env = environment(workspace)
    if not python.is_file():
        subprocess.run([sys.executable, "-m", "venv", str(Path(workspace) / ".venv")], check=True, env=env)
    subprocess.run([str(python), "-m", "pip", "install", "-r", str(REPO / "requirements-student.txt")],
                   check=True, env=env)
    subprocess.run([str(python), "-c", "import numpy,matplotlib,yaml,uproot,awkward,xgboost,sklearn; "
                    "print('Analysis dependencies ready')"], check=True, env=env)
    print("Analysis environment ready. MG5, POWHEG, ihixs and shower setup are not required.")


def command(args, dataset, index, workspace, python, output):
    profile = args.profile
    manifests = [safe_path(dataset, index["manifests"][p]) for p in
                 (("perfect", "ssc") if profile == "both" else (profile,))]
    if args.action == "report":
        return [str(python), str(REPO / "hfourlepton/make_fourlepton_report.py"),
                *map(str, manifests), "--output-dir", str(output), "--observables", "--strict",
                "--luminosity-fb", str(args.luminosity_fb)]
    if profile == "both":
        raise ValueError("choose --profile ssc or --profile perfect for cuts/XGBoost")
    card = args.config.resolve() if args.config else workspace / index["cards"][f"{args.action}_{profile}"]
    if not card.is_file():
        raise ValueError(f"missing analysis card: {card}")
    return [str(python), str(REPO / "hfourlepton/analyze_fourlepton.py"), args.action,
            "--manifest", str(manifests[0]), "--config", str(card), "--output-dir", str(output),
            "--luminosity-fb", str(args.luminosity_fb)]


def analyze(args):
    workspace = args.workspace.resolve()
    dataset, index = bound(workspace)
    if not math.isfinite(args.luminosity_fb) or args.luminosity_fb <= 0:
        raise ValueError("luminosity must be positive and finite")
    label = name(args.name)
    output = safe_path(workspace, Path("results") / label / f"{args.action}_{args.profile}")
    if output.resolve().is_relative_to(dataset):
        raise ValueError("analysis output must not be inside the shared release")
    if ((output.exists() and any(output.iterdir())) or output.with_suffix(".run.json").exists()
            or output.with_suffix(".log").exists()):
        raise ValueError(f"output is nonempty: {output}; choose a new --name")
    python = workspace / ".venv/bin/python"
    if not python.is_file():
        raise ValueError("analysis environment is missing; run student_events.py setup first")
    cmd = command(args, dataset, index, workspace, python, output)
    env = environment(workspace)
    # Store commands and identities outside the output: analysis CLIs require an
    # empty output directory and refuse to overwrite completed results.
    output.parent.mkdir(parents=True, exist_ok=True)
    log = output.with_suffix(".log")
    receipt = dict(command=cmd, dataset=str(dataset), dataset_sha256=digest(dataset / "dataset.json"),
                   code_commit=subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                       capture_output=True, text=True).stdout.strip(),
                   script_sha256=digest(cmd[1]), status="running")
    receipt["analysis_sources"] = {str(p.relative_to(REPO)): digest(p) for p in
        [*sorted((REPO / "hfourlepton").glob("*.py")), REPO / "ho_xgboost_core.py",
         REPO / "hgammagamma/ho_shower_completion.py", REPO / "hgammagamma/ho_signal_normalization.py"]}
    if "--config" in cmd:
        card = Path(cmd[cmd.index("--config") + 1])
        receipt["card_sha256"] = digest(card)
        receipt["card"] = json.loads(card.read_text()) if card.suffix == ".json" else card.read_text()
    receipt_path = output.with_suffix(".run.json")
    save(receipt_path, receipt)
    print(f"Running {args.action} ({args.profile}); log: {log}", flush=True)
    with log.open("w") as stream:
        stream.write("COMMAND: " + json.dumps(cmd) + "\n")
        stream.flush()
        result = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, env=env, cwd=REPO)
    receipt.update(status="complete" if result.returncode == 0 else "failed", exit_code=result.returncode)
    save(receipt_path, receipt)
    if result.returncode:
        raise ValueError(f"analysis failed (exit {result.returncode}); inspect {log}")
    print(f"Finished: {output / 'index.html'}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("attach", "setup", "check", "report", "cuts", "xgboost"))
    parser.add_argument("--workspace", type=Path, default=DEFAULT_WORKSPACE)
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET, help="Shared release to attach.")
    parser.add_argument("--profile", choices=("ssc", "perfect", "both"), default=None)
    parser.add_argument("--name", default="baseline", help="New result label; existing outputs are never overwritten.")
    parser.add_argument("--config", type=Path, help="Private analysis card for cuts/XGBoost.")
    parser.add_argument("--luminosity-fb", type=float, default=100.)
    parser.add_argument("--full", action="store_true", help="With check, also hash every large event file.")
    args = parser.parse_args(argv)
    args.profile = args.profile or ("both" if args.action == "report" else "ssc")
    if args.action == "attach":
        attach(args.dataset, args.workspace)
    elif args.action == "setup":
        setup(args.workspace)
    elif args.action == "check":
        dataset, index = bound(args.workspace, full=args.full)
        print(f"Ready: {len(index['samples'])} samples; {dataset}\n"
              f"Checks: {'all SHA-256 hashes' if args.full else 'all file sizes/access and metadata SHA-256 hashes'}")
    else:
        # Keep a student's concurrent commands from racing on result paths.
        bound(args.workspace)
        with (args.workspace / ".analysis.lock").open("a") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise ValueError("another analysis is using this workspace") from error
            analyze(args)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        sys.exit(f"Error: {error}")
