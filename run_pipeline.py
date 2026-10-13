#!/usr/bin/env python3
"""Set up and run the student HiggsSSC workflows on Timur.

Run commands from the repository root. A run includes event generation,
showering, detector analysis and the HTML report; cuts and ML are optional.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import math
from pathlib import Path
import re
import subprocess

from hgammagamma import run_gammagamma_ho_pipeline as shared

REPO = Path(__file__).resolve().parent
WORKFLOWS = ("gammagamma-lo", "gammagamma-ho", "fourlepton-lo")
ACTIONS = ("setup", "check", "run", "status", "report", "cuts", "xgboost")


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=ACTIONS)
    p.add_argument("--workflow", required=True, choices=WORKFLOWS)
    p.add_argument("--work-dir", type=Path, default=REPO / "inputs/pipeline")
    p.add_argument("--run-tag", help="name of this campaign; required for run/status/report/cuts/xgboost")
    p.add_argument("--nevents", type=int, default=30)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--luminosity-fb", type=float, default=100.)
    p.add_argument("--herwig-module", default=None)
    p.add_argument("--herwig-env", type=Path)
    p.add_argument("--cc", default=None)
    p.add_argument("--cxx", default=None)
    p.add_argument("--fc", default=None)
    p.add_argument("--ninja-dir", type=Path)
    p.add_argument("--with-ml", action="store_true", help="include optional XGBoost dependencies during setup/check")
    p.add_argument("--signal-lhe", type=Path, help="HO only: reuse a compatible signal LHE")
    p.add_argument("--signal-normalization", type=Path, help="HO only: reuse a validated ihixs record")
    p.add_argument("--config", type=Path, help="optional cuts/XGBoost analysis card")
    p.add_argument("--dry-run", action="store_true", help="print the plan without installing or running programs")
    args = p.parse_args(argv)
    if args.action not in ("setup", "check") and not args.run_tag:
        p.error("--run-tag is required for this action")
    if args.run_tag and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.run_tag):
        p.error("run tag must contain only letters, digits, underscores and hyphens")
    if args.nevents < 1 or args.jobs < 1 or not math.isfinite(args.luminosity_fb) or args.luminosity_fb <= 0:
        p.error("nevents, jobs and luminosity must be positive and finite")
    if args.herwig_env and args.herwig_module:
        p.error("select either --herwig-env or --herwig-module")
    if args.workflow != "gammagamma-ho" and (args.signal_lhe or args.signal_normalization):
        p.error("signal input overrides apply only to gammagamma-ho")
    if args.workflow == "fourlepton-lo" and args.action in ("cuts", "xgboost"):
        p.error("cuts and xgboost currently support the diphoton workflows")
    if args.config and args.action not in ("cuts", "xgboost"):
        p.error("--config applies to cuts/xgboost")
    for key in ("work_dir", "herwig_env", "ninja_dir", "signal_lhe", "signal_normalization", "config"):
        value = getattr(args, key)
        if value is not None:
            setattr(args, key, value.expanduser().resolve())
    for path in (REPO, args.work_dir, args.signal_lhe, args.signal_normalization):
        if path and re.search(r"\s|[;#\"']", str(path)):
            p.error("generator paths cannot contain whitespace or card metacharacters")
    return args


def locations(args):
    if args.workflow == "fourlepton-lo":
        directory = REPO / "hfourlepton/runs" / args.run_tag
        return directory, directory / "report"
    order = "HOAnalysis" if args.workflow == "gammagamma-ho" else "LOAnalysis"
    base = REPO / "hgammagamma" / order
    return base / "runs" / args.run_tag, base / "plots" / args.run_tag


def state_path(args):
    return args.work_dir / "campaigns" / args.workflow / (args.run_tag + ".json")


def read_state(args):
    path = state_path(args)
    if not path.is_file():
        raise ValueError(f"no managed campaign {args.run_tag!r} for {args.workflow}; check --work-dir and --run-tag")
    return json.loads(path.read_text())


def preflight_inputs(args):
    for path in (args.signal_lhe, args.signal_normalization, args.config):
        if path and (not path.is_file() or not path.stat().st_size):
            raise ValueError(f"supplied input is missing or empty: {path}")
    if args.signal_normalization:
        shared.norm.load_record(args.signal_normalization)


def identity(args, config):
    result = {"workflow": args.workflow, "run_tag": args.run_tag, "environment": config,
              "nevents": args.nevents, "jobs": args.jobs}
    for key in ("signal_lhe", "signal_normalization"):
        path = getattr(args, key)
        result[key] = {"path": str(path), "sha256": shared.digest(path)} if path else None
    if args.workflow == "gammagamma-ho":
        result["powheg_card_sha256"] = shared.digest(shared.HO / "powheg-hjminnlo-ssc40-nnpdf40nnloqed.input")
    return result


def prepare_run(args, config):
    expected = identity(args, config)
    path = state_path(args)
    if path.exists():
        state = read_state(args)
        if state["identity"] != expected:
            raise ValueError("campaign configuration changed; repeat the original options or choose a new --run-tag")
        return state
    directory, report = locations(args)
    if directory.exists() and any(directory.iterdir()):
        raise ValueError(f"unmanaged nonempty campaign: {directory}; use its original runner or a new tag")
    state = {"schema_version": 1, "identity": expected, "directory": str(directory),
             "report": str(report), "created_at": now(), "state": "prepared", "stage": "preflight"}
    shared.save(path, state)
    return state


@contextmanager
def installation_lock(work):
    work.mkdir(parents=True, exist_ok=True)
    with (work / ".pipeline.lock").open("a") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("another command is using this installation; use status to inspect it") from exc
        yield


def progress(args, **values):
    if not args.run_tag or not state_path(args).is_file():
        return
    state = read_state(args)
    state.update(values, updated_at=now())
    shared.save(state_path(args), state)


class StudentPipeline(shared.Pipeline):
    def __init__(self, args, config):
        super().__init__(args, config)
        # These legacy shell defaults must not change a saved student campaign.
        from hgammagamma import run_gammagamma_campaign as lo
        from hfourlepton import run_four_lepton_campaign as four
        for module in (lo, four):
            source = Path(module.__file__).read_text()
            for key in re.findall(r'(?:env_(?:text|int|float|bool)|os\.environ\.get)\(\s*"([A-Z][A-Z0-9_]*)"', source):
                if key != "PATH":
                    self.env.pop(key, None)
        self.env["PATH"] = str(self.python.parent) + ":" + self.env["PATH"]
        self.env["PYTHONNOUSERSITE"] = "1"

    def call(self, label, argv, cwd=None, extra_env=None, mg5=False):
        progress(self.args, state="running", stage=label, log=str(self.logs / (label + ".log")), error=None)
        try:
            result = super().call(label, argv, cwd=cwd, extra_env=extra_env, mg5=mg5)
        except (Exception, KeyboardInterrupt) as exc:
            progress(self.args, state="failed", error=str(exc) or "command interrupted")
            raise
        progress(self.args, last_completed_stage=label)
        return result


def campaign_command(pipeline, workflow, profile=None):
    args = pipeline.args
    script = "hfourlepton/run_four_lepton_campaign.py" if workflow == "fourlepton-lo" else "hgammagamma/run_gammagamma_campaign.py"
    command = [pipeline.python, REPO / script, "--run-tag", args.run_tag,
               "--nevents", args.nevents, "--nb-core", args.jobs,
               "--mg5-dir", pipeline.mg5, "--mg5-python", pipeline.python,
               "--no-herwig-module", "--no-herwig-env", "--ebeam", 20000]
    if workflow == "fourlepton-lo":
        command.extend(["--sample-set", "irreducible", "--detector-response", profile])
    else:
        command.extend(["--stage", "all", "--resume", "--output-dir", locations(args)[0],
                        "--run-samples", ",".join(shared.SAMPLES), "--detector-response", "ssc"])
    return command


def physics_artifacts(args):
    """Check saved analysis provenance without invoking an analyzer or generator.

    Successful student runs retain per-campaign hashes in their saved state.
    Older HO manifests record summary values instead of output hashes. Retain
    that contract, and use existing pilot hashes when they cover this campaign.
    """
    directory, _ = locations(args)
    paths = []
    expected_hashes = {}
    ho_summaries = []

    def expect_hash(path, value):
        if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
            raise ValueError(f"completed analysis lacks a recorded SHA-256: {path}")
        key = str(Path(path).resolve())
        if key in expected_hashes and expected_hashes[key] != value:
            raise ValueError(f"conflicting completed analysis provenance: {path}")
        expected_hashes[key] = value

    if args.workflow == "fourlepton-lo":
        for profile in ("perfect", "ssc"):
            path = directory / f"manifest_{profile}.json"
            manifest = json.loads(path.read_text())
            if len(manifest["samples"]) != 9 or any(s["status"] != "complete" for s in manifest["samples"]):
                raise ValueError(f"incomplete four-lepton {profile} campaign")
            paths.append(path)
            for sample in manifest["samples"]:
                paths.extend(Path(sample[key]) for key in ("root_file", "analysis_summary"))
                expect_hash(sample["root_file"], sample.get("output_root_sha256"))
                expect_hash(sample["analysis_summary"], sample.get("analysis_summary_sha256"))
    else:
        for sample in shared.SAMPLES:
            category = "Signal" if sample.startswith("signal_") else "Backgrounds"
            folder = directory / category / "events" / sample
            path = folder / "campaign.json"
            manifest = json.loads(path.read_text())
            if not {"generate", "shower", "analyze"}.issubset(manifest.get("completed", [])):
                raise ValueError(f"analysis is incomplete for {sample}")
            paths.append(path)
            outputs = [folder / f"{sample}_hwsim_roots-{args.run_tag}{suffix}"
                       for suffix in (".dat", ".top", "_var.root")]
            paths.extend(outputs)
            if args.workflow == "gammagamma-ho":
                paths.append(folder / f"normalization-{args.run_tag}.json")
                ho_summaries.append((folder, outputs[0]))
            else:
                recorded = {str(Path(entry["path"]).resolve()): entry.get("sha256")
                            for entry in manifest.get("artifacts", {}).get("analyze", [])}
                for output in outputs:
                    expect_hash(output, recorded.get(str(output.resolve())))
    for path in paths:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"missing or empty analysis product: {path}")
    if ho_summaries:
        from hgammagamma.make_gammagamma_report import validate_ho_analysis_summary
        for folder, summary in ho_summaries:
            validate_ho_analysis_summary(folder, summary)
    current = {str(path.resolve()): shared.digest(path) for path in paths}
    if args.workflow == "gammagamma-ho":
        shared.validate_analysis_snapshot(directory, current)
    if state_path(args).is_file():
        recorded = read_state(args).get("analysis_artifacts")
        if recorded is not None:
            if not isinstance(recorded, dict) or {str(Path(name).resolve()) for name in recorded} != set(current):
                raise ValueError("completed campaign analysis-artifact inventory differs from saved evidence")
            for name, value in recorded.items():
                expect_hash(name, value)
    snapshots = [pilot_path(args)]
    if args.workflow == "gammagamma-ho":
        snapshots.append(args.work_dir / "pilot-complete.json")
    for snapshot in snapshots:
        if snapshot.is_file():
            recorded = json.loads(snapshot.read_text()).get("artifacts", {})
            for name, value in recorded.items():
                # Reporting a production tag does not depend on another tag's pilot.
                if str(Path(name).resolve()) in current:
                    expect_hash(name, value)
    for name, expected in expected_hashes.items():
        if current[name] != expected:
            raise ValueError(f"completed analysis product changed: {name}; preserve it and use a new run tag")
    return {str(path): current[str(path.resolve())] for path in paths}


def pilot_path(args):
    return args.work_dir / "pilots" / (args.workflow + ".json")


def verify_pilot(args, config):
    path = pilot_path(args)
    if not path.is_file():
        raise ValueError(f"first complete a 30-event {args.workflow} pilot with a separate run tag")
    pilot = json.loads(path.read_text())
    if pilot["environment"] != config:
        raise ValueError("pilot used a different environment; complete a new pilot")
    for name, expected in pilot["artifacts"].items():
        if not Path(name).is_file() or shared.digest(name) != expected:
            raise ValueError(f"pilot physics artifact changed: {name}")


def report(pipeline):
    args = pipeline.args
    physics_artifacts(args)  # Never generate to repair a missing analysis.
    directory, output = locations(args)
    legacy = args.work_dir / "pilot-complete.json"
    legacy_pilot = None
    if args.workflow == "gammagamma-ho" and legacy.exists():
        candidate = json.loads(legacy.read_text())
        if str(output / "index.html") in candidate.get("artifacts", {}):
            shared.verify_pilot(args.work_dir)
            legacy_pilot = candidate
    if args.workflow == "fourlepton-lo":
        command = [pipeline.python, REPO / "hfourlepton/make_fourlepton_report.py",
                   directory / "manifest_perfect.json", directory / "manifest_ssc.json",
                   "--output-dir", output, "--luminosity-fb", args.luminosity_fb, "--strict"]
    else:
        command = [pipeline.python, REPO / "hgammagamma/make_gammagamma_report.py",
                   "--analysis-root", directory, "--run-tag", args.run_tag, "--output-dir", output,
                   "--no-density", "--normalization", "event_xsec", "--luminosity-fb", args.luminosity_fb]
    pipeline.call(args.run_tag + "-report", command)
    if not (output / "index.html").is_file():
        raise ValueError(f"report did not produce {output / 'index.html'}")
    if legacy_pilot is not None:
        legacy_pilot["artifacts"][str(output / "index.html")] = shared.digest(output / "index.html")
        shared.save(legacy, legacy_pilot)
    print(f"Report: {output / 'index.html'}")


def downstream_analysis(pipeline):
    args = pipeline.args
    physics_artifacts(args)
    if args.action == "xgboost":
        try:
            pipeline.call("check-python-ml", [pipeline.python, "-c", "import xgboost,sklearn,tqdm"])
        except RuntimeError as exc:
            raise ValueError(f"ML dependencies are unavailable; run setup --workflow {args.workflow} --with-ml") from exc
    prefix = "ho_" if args.workflow == "gammagamma-ho" else ""
    basename = "baseline_cuts.yaml" if args.action == "cuts" else "xgboost_baseline.yaml"
    card = args.config or REPO / "hgammagamma/analysis_cards" / (prefix + basename)
    directory, _ = locations(args)
    output = directory.parent.parent / "analyses" / args.run_tag / args.action
    pipeline.call(args.run_tag + "-" + args.action, [pipeline.python, REPO / "analyze_lo_varfiles.py",
                  args.action, "--config", card, "--run-tag", args.run_tag,
                  "--analysis-root", directory, "--output-dir", output, "--luminosity-fb", args.luminosity_fb])


def validate_saved_run_artifacts(args):
    """Validate completed evidence before a resumed run can rewrite reports."""
    saved = read_state(args) if state_path(args).is_file() else None
    completed = saved is not None and (
        saved.get("state") == "complete" or saved.get("analysis_artifacts") is not None
    )
    directory = locations(args)[0].resolve()
    if args.workflow == "gammagamma-ho" and shared.analysis_snapshot_path(directory).is_file():
        completed = True
    snapshots = [pilot_path(args)]
    if args.workflow == "gammagamma-ho":
        snapshots.append(args.work_dir / "pilot-complete.json")
    for snapshot in snapshots:
        if snapshot.is_file():
            artifacts = json.loads(snapshot.read_text()).get("artifacts", {})
            if any(Path(name).resolve().is_relative_to(directory) for name in artifacts):
                completed = True
    if completed:
        physics_artifacts(args)


def run(pipeline, *, saved_artifacts_checked=False):
    args = pipeline.args
    if not saved_artifacts_checked:
        validate_saved_run_artifacts(args)
    if args.workflow == "gammagamma-ho":
        pipeline.run()
    else:
        pipeline.check()
        if args.workflow == "fourlepton-lo":
            for profile in ("perfect", "ssc"):
                pipeline.call(args.run_tag + "-" + profile, campaign_command(pipeline, args.workflow, profile))
        else:
            pipeline.call(args.run_tag + "-campaign", campaign_command(pipeline, args.workflow))
        report(pipeline)
    artifacts = physics_artifacts(args)
    if args.nevents <= 30:
        shared.save(pilot_path(args), {"environment": pipeline.config, "artifacts": artifacts})
    progress(args, state="complete", stage="complete", error=None, analysis_artifacts=artifacts)
    print(f"Pipeline complete: {locations(args)[1] / 'index.html'}")


def status(args):
    state = read_state(args)
    print(f"{args.workflow}: {args.run_tag}\nState: {state['state']}\nLast recorded stage: {state['stage']}")
    if state.get("updated_at"):
        print("Updated:", state["updated_at"])
    log = Path(state["log"]) if state.get("log") else None
    if log:
        print("Log:", log)
        if log.exists():
            print("Last log activity:", datetime.fromtimestamp(log.stat().st_mtime, timezone.utc).isoformat(timespec="seconds"))
    directory, output = locations(args)
    for path in sorted(directory.glob("*/events/*/campaign.json")):
        sample = json.loads(path.read_text())
        print(f"  {sample['sample']}: {', '.join(sample.get('completed', [])) or 'prepared'}")
    for path in sorted(directory.glob("manifest_*.json")):
        samples = json.loads(path.read_text()).get("samples", [])
        print(f"  {path.stem}: {sum(s.get('status') == 'complete' for s in samples)}/{len(samples)} samples complete")
    print("Report:", output / "index.html" if (output / "index.html").exists() else "not created")
    if state.get("error"):
        print("Error:", state["error"])
    if state["state"] == "complete":
        print("Next: inspect the report and event/weight checks before choosing a larger run tag.")
    elif state["state"] == "running":
        print("Next: inspect the log. Integrations may be quiet; this is the last recorded state, not a liveness check.")
    else:
        print("Next: inspect the error/log, then repeat the original run command to resume valid stages.")
    return 0


def main(argv=None):
    args = parse_args(argv)
    if args.dry_run:
        stages = {"gammagamma-lo": ["MG5 signal and three backgrounds", "Herwig showers", "SSC analysis", "report"],
                  "gammagamma-ho": ["ihixs build/benchmark/109-point calculation", "POWHEG integration and signal", "MG5 backgrounds", "Herwig showers", "SSC analysis", "report"],
                  "fourlepton-lo": ["seven hard-process integrations / nine samples", "MadSpin signal decays", "perfect and SSC showers/analysis", "strict comparison report"]}
        print(json.dumps({"action": args.action, "workflow": args.workflow, "work_dir": str(args.work_dir),
                          "run_tag": args.run_tag, "nevents": args.nevents, "jobs": args.jobs,
                          "stages": stages[args.workflow] if args.action == "run" else [args.action],
                          "note": "Small samples still require generator integration; HO includes full normalization."}, indent=2))
        return 0
    if args.action == "status":
        return status(args)
    preflight_inputs(args)
    config = shared.configuration(args)
    if args.action != "setup" and not (args.work_dir / "environment.json").is_file():
        raise ValueError("run setup for this workflow first")
    if args.action == "setup" and args.work_dir.exists() and not (args.work_dir / "environment.json").exists() and any(args.work_dir.iterdir()):
        raise ValueError("unmanaged nonempty installation; choose a new --work-dir")
    with installation_lock(args.work_dir):
        shared.frozen(args.work_dir / "environment.json", config)
        if args.action == "run":
            if args.nevents > 30:
                verify_pilot(args, config)
            prepare_run(args, config)
            validate_saved_run_artifacts(args)
        elif args.action in ("report", "cuts", "xgboost"):
            saved = read_state(args)
            if saved["identity"]["environment"] != config:
                raise ValueError("campaign belongs to a different environment")
        try:
            pipeline = StudentPipeline(args, config)
            if args.action in ("setup", "check"):
                getattr(pipeline, args.action)()
            elif args.action == "run":
                run(pipeline, saved_artifacts_checked=True)
            else:
                previous = read_state(args)
                if args.action == "report":
                    report(pipeline)
                else:
                    downstream_analysis(pipeline)
                progress(args, state=previous["state"], stage=previous["stage"], error=None)
        except (Exception, KeyboardInterrupt) as exc:
            progress(args, state="failed", error=str(exc) or "command interrupted")
            raise
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"Error: {exc}")
