#!/usr/bin/env python3
"""Audit and restart error-limited qq showers, using completed hard events only.

The recovery is deliberately separate from the physics adapters: changing a
shower error budget must not invalidate or regenerate completed MG5 production.
Preparation archives failed showers. Running applies the recorded card change
in both the coordinator and its isolated workers, then analyzes perfect/SSC.
"""
import argparse
from contextlib import contextmanager
from dataclasses import replace
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import sys
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hfourlepton import higher_order as ho, minimal_runtime as runtime
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton.repair_loop_production import completed_production

SAMPLES = {"bkg_qq_4e", "bkg_qq_4mu", "bkg_qq_2e2mu"}
RECEIPT = "shower_error_recovery.json"
LIMIT = r"(?m)^set EventGenerator:MaxErrors ([0-9]+)$"


def write_json(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=runtime._encode) + "\n")
    temporary.replace(path)


def change_budget(card, budget):
    limits = re.findall(LIMIT, card)
    if len(limits) != 1 or not 0 < int(limits[0]) < budget:
        raise ValueError("recovery requires exactly one positive, increasing MaxErrors setting")
    return re.sub(LIMIT, f"set EventGenerator:MaxErrors {budget}", card)


def configuration(run):
    jobs = sorted((run / "workers/perfect").glob("*.job.json"))
    if not jobs:
        raise ValueError("no original perfect-profile worker configuration")
    configs = [json.loads(p.read_text(), object_hook=runtime._decode)["config"] for p in jobs]
    cfg = replace(four.Config(**configs[0]), manifest_override=None, postprocess_only=False)
    for value in configs[1:]:
        if replace(four.Config(**value), manifest_override=None, postprocess_only=False) != cfg:
            raise ValueError("original worker configurations disagree")
    if (cfg.run_dir.resolve() != run or cfg.force or cfg.dry_run or cfg.order != "ho"
            or cfg.preset != "minimal" or cfg.detector_response != "perfect"
            or cfg.selection_profile != "ssc"):
        raise ValueError("recovery requires the original minimal HO/SSC perfect-profile campaign")
    if {s.name for s in four.selected_samples(cfg)} != SAMPLES:
        raise ValueError("this recovery supports exactly the three qq continuum samples")
    return cfg


def check_identity(cfg, manifest):
    versions = {key: manifest["software"][key + "_version"]
                for key in ("python", "mg5", "herwig", "root")}
    if (four.campaign_tracked_inputs(cfg) != manifest["tracked_inputs"]
            or four.configuration_fingerprint(cfg, four.selected_samples(cfg),
                software_versions=versions) != manifest["configuration_fingerprint"]):
        raise ValueError("campaign configuration or tracked physics software changed")


def require_idle(run):
    """Refuse to move files while a local shower/analysis is using this run."""
    for directory in Path("/proc").glob("[0-9]*"):
        if int(directory.name) == os.getpid():
            continue
        try:
            if directory.stat().st_uid != os.getuid():
                continue
            cwd = (directory / "cwd").resolve(strict=True)
            if cwd == run or run in cwd.parents:
                raise ValueError(f"campaign is in use by PID {directory.name}: {cwd}")
        except (OSError, RuntimeError):
            continue


def audit(run, budget):
    require_idle(run)
    cfg = configuration(run)
    manifest = json.loads(cfg.manifest_path.read_text())
    check_identity(cfg, manifest)
    if (run / "manifest_ssc.json").exists():
        raise ValueError("SSC analysis already exists; recovery scope differs")
    samples = four.selected_samples(cfg)
    records = {}
    for sample in samples:
        print(f"Auditing {sample.name}: completed production and failed shower", flush=True)
        entry = four.manifest_sample(manifest, sample.name)
        journal = run / "workers/perfect" / (sample.name + ".json")
        private = json.loads(journal.read_text())
        if (entry["status"] != "failed" or four.expected_analysis_output(sample, cfg).exists()
                or private["configuration_fingerprint"] != manifest["configuration_fingerprint"]
                or four.manifest_sample(private, sample.name) != entry):
            raise ValueError(f"not an unanalyzed failed shower with matching journal: {sample.name}")
        production = completed_production(sample, samples, cfg)
        directory = four.shower_work_dir(sample, cfg) / "herwig"
        if (directory / "matched_completion.json").exists():
            raise ValueError(f"completed shower cannot be replaced: {directory}")
        card = (directory / (sample.name + ".in")).read_text()
        source = Path(ho.completion._audit_source_card(directory / (sample.name + ".in")))
        expected = ho.herwig_card(sample, source, four.seed_for(sample.name + ":herwig", cfg),
                                  cfg.events_for(sample), cfg)
        if card != expected:
            raise ValueError(f"saved shower card differs: {sample.name}")
        updated = change_budget(card, budget)
        parsed = ho.completion._parse_logs(directory, sample.name, "error_limit")
        weights = ho.completion.source_weight_summary(source, parsed["attempted_events"])
        if weights["lhe_events"] != cfg.events_for(sample):
            raise ValueError(f"unexpected hard-event count: {sample.name}")
        old_budget = ho.ho.validate_error_limit_stop(parsed, {
            "configuration": {"herwig_card": card}, "nevents_requested": cfg.events_for(sample),
            "lhe": {"events": weights["lhe_events"]}})
        metadata = json.loads((directory / "root-metadata.log").read_text().splitlines()[-1])
        if ho.completion._root_summary(metadata)["entries"] != parsed["saved_events"]:
            raise ValueError(f"failed ROOT entry count disagrees: {sample.name}")
        roots = {str(p.resolve()) for p in (directory / "events").rglob("*.root")}
        if roots != set(parsed["root_writers"]) or roots != {v["path"] for v in metadata}:
            raise ValueError(f"failed ROOT inventory disagrees: {sample.name}")
        for value in metadata:
            stat = Path(value["path"]).stat()
            if (stat.st_size, stat.st_mtime_ns) != (value["size"], value["mtime_ns"]):
                raise ValueError(f"failed ROOT changed: {value['path']}")
        artifacts = []
        for path in sorted(directory.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"cannot archive shared shower artifacts: {path}")
            if path.is_file():
                artifacts.append(four.artifact_stage_record(path))
        records[sample.name] = dict(production=production, source=four.artifact_stage_record(source),
            previous_max_errors=old_budget, max_errors=budget, old_card=card, card=updated,
            source_weights=weights, failed_shower=parsed, artifacts=artifacts)
    return dict(schema_version=1, run_dir=str(run), config=vars(cfg),
        configuration_fingerprint=manifest["configuration_fingerprint"], samples=records,
        manifest=four.artifact_stage_record(cfg.manifest_path),
        implementation=four.artifact_stage_record(Path(__file__).resolve()))


def prepare(run, budget):
    receipt = run / RECEIPT
    if receipt.exists():
        raise ValueError(f"recovery already prepared: {receipt}; use run")
    record = audit(run, budget)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = run / "shower-recovery" / stamp
    archive.mkdir(parents=True, exist_ok=False)
    moved = []
    try:
        shutil.copy2(run / "manifest_perfect.json", archive / "manifest_perfect.json")
        shutil.copytree(run / "workers", archive / "workers")
        for name in sorted(record["samples"]):
            source = run / "samples/shared" / name / "herwig"
            target = archive / name
            source.rename(target)
            moved.append((source, target))
        record.update(archive=str(archive), prepared_at_utc=stamp)
        write_json(receipt, record)
    except BaseException:
        for source, target in reversed(moved):
            target.rename(source)
        raise
    print(f"Prepared recovery: {receipt}\nFailed showers preserved in {archive}", flush=True)
    return record


def load_receipt(run):
    path = run / RECEIPT
    receipt = json.loads(path.read_text(), object_hook=runtime._decode)
    if receipt["run_dir"] != str(run) or set(receipt["samples"]) != SAMPLES:
        raise ValueError("recovery receipt has a different campaign or sample inventory")
    if receipt["implementation"] != four.artifact_stage_record(Path(__file__).resolve()):
        raise ValueError("recovery implementation changed since preparation")
    cfg = four.Config(**receipt["config"])
    archived = Path(receipt["archive"]) / "manifest_perfect.json"
    if (four.sha256_file(archived) != receipt["manifest"]["sha256"]
            or archived.stat().st_size != receipt["manifest"]["size"]):
        raise ValueError("archived original manifest changed")
    original = json.loads(archived.read_text())
    check_identity(cfg, original)
    return receipt, cfg


@contextmanager
def recovery_runtime(receipt):
    """Apply only the audited card edit and forbid all hard-event production."""
    run = Path(receipt["run_dir"])
    reference = four.artifact_stage_record(run / RECEIPT)
    original_card, original_manifest, original_sample = ho.herwig_card, four.build_manifest, four.run_sample
    original_popen = runtime.subprocess.Popen

    def checked_record(path, expected):
        if four.artifact_stage_record(Path(path)) != expected:
            raise ValueError(f"completed recovery input changed: {path}")

    def card(sample, lhe, seed, count, cfg):
        saved = receipt["samples"][sample.name]
        text = original_card(sample, lhe, seed, count, cfg)
        if text != saved["old_card"] or Path(lhe).resolve() != Path(saved["source"]["path"]).resolve():
            raise ValueError("recovery cannot change any card setting other than MaxErrors")
        checked_record(lhe, saved["source"])
        if change_budget(text, saved["max_errors"]) != saved["card"]:
            raise ValueError("recorded shower card change is inconsistent")
        return saved["card"]

    def production(sample, selected, cfg):
        saved = receipt["samples"][sample.name]["production"]
        path = Path(saved["path"])
        checked_record(path, {key: saved[key] for key in ("path", "size", "sha256")})
        checked_record(saved["marker"]["path"], saved["marker"])
        return path

    def manifest(cfg, samples):
        value = original_manifest(cfg, samples)
        value["shower_error_recovery"] = reference
        return value

    def run_sample(sample, samples, value, cfg):
        entry = four.manifest_sample(value, sample.name)
        entry["shower_error_recovery"] = reference
        original_sample(sample, samples, value, cfg)
        if entry.get("status") == "complete":
            entry.pop("error", None)

    def popen(argv, *args, **kwargs):
        if len(argv) > 2 and str(argv[1]) == str(Path(runtime.__file__).resolve()) and argv[2] == "--worker":
            argv = [argv[0], str(Path(__file__).resolve()), "worker", *argv[3:]]
        return original_popen(argv, *args, **kwargs)

    with patch.object(ho, "herwig_card", card), patch.object(four, "obtain_production_lhe", production), \
            patch.object(four, "build_manifest", manifest), patch.object(four, "run_sample", run_sample), \
            patch.object(runtime.subprocess, "Popen", popen):
        yield


def run_recovery(run):
    require_idle(run)
    receipt, cfg = load_receipt(run)
    with recovery_runtime(receipt):
        for profile in ("perfect", "ssc"):
            current = replace(cfg, detector_response=profile)
            if current.manifest_path.exists():
                saved = json.loads(current.manifest_path.read_text())
                check_identity(current, saved)
                if saved.get("status") == "complete":
                    from hfourlepton.reprocess import audit_completed_manifest
                    audit_completed_manifest(current.manifest_path)
                    print(f"Already complete: {profile}", flush=True)
                    continue
            # Use precisely the saved configuration, never a reconstructed set
            # of CLI defaults. four.main still checks its normal fingerprint.
            with patch.object(four, "parse_config", lambda argv: current):
                four.main([])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("audit", "prepare", "run"):
        command = sub.add_parser(action)
        command.add_argument("--run-dir", required=True, type=Path)
        if action != "run":
            command.add_argument("--max-errors", required=True, type=int)
    worker = sub.add_parser("worker", help=argparse.SUPPRESS)
    worker.add_argument("job", type=Path)
    args = parser.parse_args(argv)
    os.chdir(four.REPO_ROOT)
    if args.action == "worker":
        payload = json.loads(args.job.read_text(), object_hook=runtime._decode)
        cfg = four.Config(**payload["config"])
        receipt, _ = load_receipt(cfg.run_dir.resolve())
        with recovery_runtime(receipt):
            runtime.run_worker(args.job)
        return 0
    run = args.run_dir.resolve()
    with (run / ".shower-error-recovery.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.action == "prepare":
            prepare(run, args.max_errors)
        elif args.action == "audit":
            result = audit(run, args.max_errors)
            print(json.dumps({name: {"saved": row["failed_shower"]["saved_events"],
                "max_errors": row["max_errors"]} for name, row in result["samples"].items()}, indent=2))
        else:
            run_recovery(run)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}")
