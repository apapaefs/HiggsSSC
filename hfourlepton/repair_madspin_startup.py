#!/usr/bin/env python3
"""Migrate only the audited minimal MadSpin logging-startup failure.

No simulation is run. All hard production must be complete, and no decay or
shower may have started. Original metadata and worker logs are archived before
updating the manifest for the corrected launcher.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_pipeline as cli
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton.repair_loop_production import saved_configuration, completed_production

MIGRATION = "madspin-logging-startup-v1"
OLD_WORKER = "906d3f187d0d3293187658c6436621b40ba5f195d1c1ae1d2729041a2d761e1b"
FIXED_WORKER = "27de8ad84248942e880269a9f6f1fc8f59573aaeaa249bd7cdeac3cb41bd808e"


def fingerprint_update(cfg, samples, manifest):
    """Only this launcher transition is compatible with completed hard events."""
    current = four.campaign_tracked_inputs(cfg)
    previous = manifest["tracked_inputs"]
    if (previous["madspin_worker"]["sha256"] != OLD_WORKER
            or current["madspin_worker"]["sha256"] != FIXED_WORKER):
        raise ValueError("MadSpin launcher versions are outside this audited repair")
    old_inputs = deepcopy(current)
    old_inputs["madspin_worker"] = previous["madspin_worker"]
    versions = {name: manifest["software"][name + "_version"] for name in ("python", "mg5", "herwig", "root")}
    if (old_inputs != previous or four.configuration_fingerprint(cfg, samples,
            tracked_inputs=old_inputs, software_versions=versions) != manifest["configuration_fingerprint"]):
        raise ValueError("configuration, software or inputs changed beyond the MadSpin logging fix")
    new_fingerprint = four.configuration_fingerprint(cfg, samples, tracked_inputs=current, software_versions=versions)
    updated = deepcopy(manifest)
    updated.update(tracked_inputs=current, configuration_fingerprint=new_fingerprint)
    return updated


def check_startup_only(cfg, samples, manifest):
    """Require the exact pre-decay failure, without accepting partial output."""
    entries = {s["name"]: s for s in manifest["samples"]}
    if len(entries) != len(manifest["samples"]) or set(entries) != {s.name for s in samples}:
        raise ValueError("sample inventory differs from the original campaign")
    for name in ("samples", "analysis"):
        path = cfg.run_dir / name
        if path.exists() and any(path.iterdir()):
            raise ValueError(f"post-decay output exists: {path}; this repair does not apply")
    if (cfg.run_dir / "manifest_ssc.json").exists():
        raise ValueError("SSC postprocessing has started; this repair does not apply")
    workers = cfg.run_dir / "workers/perfect"
    if not workers.is_dir():
        raise ValueError("missing failed worker journals")
    failed = []
    artifacts = []
    for sample in samples:
        entry = entries[sample.name]
        journal = workers / (sample.name + ".json")
        decay = four.sample_decay_dir(sample, cfg)
        if entry["status"] == "planned":
            if journal.exists() or decay.exists():
                raise ValueError(f"unexpected progress in planned sample {sample.name}")
            continue
        if entry["status"] != "failed" or not sample.is_signal:
            raise ValueError(f"sample is outside the startup-failure repair: {sample.name}")
        private = json.loads(journal.read_text())
        if (private["configuration_fingerprint"] != manifest["configuration_fingerprint"]
                or private["tracked_inputs"] != manifest["tracked_inputs"]
                or next(s for s in private["samples"] if s["name"] == sample.name) != entry):
            raise ValueError(f"worker evidence differs from main manifest: {journal}")
        log = workers / (sample.name + ".log")
        text = log.read_text()
        if not all(fragment in text for fragment in ("logging.config.fileConfig", "No module named 'logging.ColorFormatter'")):
            raise ValueError(f"unrecognized worker failure: {log}")
        expected = {"stable_input_source.lhe", "stable_input.lhe", "madspin_card.dat", "fourlepton_madspin_banner_stage.json"}
        if {path.name for path in decay.iterdir()} != expected:
            raise ValueError(f"decay has progressed beyond startup or its inputs are incomplete: {decay}")
        card = four.madspin_card_text(decay / "stable_input.lhe", sample,
                                     four.seed_for(sample.name + ":madspin", cfg))
        if (decay / "madspin_card.dat").read_text() != card:
            raise ValueError(f"prepared MadSpin card changed: {decay}")
        banner = json.loads((decay / "fourlepton_madspin_banner_stage.json").read_text())
        if (banner["artifact"] != four.artifact_stage_record(decay / "stable_input.lhe")
                or banner["identity"]["source_sha256"] != four.sha256_file(decay / "stable_input_source.lhe")):
            raise ValueError(f"staged decay input changed: {decay}")
        artifacts.extend(four.artifact_stage_record(path) for path in sorted(decay.iterdir()))
        artifacts.extend(four.artifact_stage_record(path) for path in (journal, log))
        failed.append(sample.name)
    if not failed:
        raise ValueError("no logging-startup failures to repair")
    return failed, artifacts


def repair(saved, state_path, cfg, *, apply=False):
    if cfg.force or cfg.dry_run or cfg.preset != "minimal" or cfg.detector_response != "perfect":
        raise ValueError("repair requires the original minimal perfect-profile configuration")
    manifest = json.loads(cfg.manifest_path.read_text())
    backup = cfg.run_dir / "recovery" / MIGRATION
    previous = [item for item in manifest.get("software_migrations", []) if item.get("id") == MIGRATION]
    if previous:
        if (four.sha256_file(Path(__file__).with_name("madspin_worker.py")) != FIXED_WORKER
                or not (backup / "manifest_perfect.json").is_file() or not (backup / "workers").is_dir()):
            raise ValueError("existing startup repair evidence is missing or changed")
        print("MadSpin startup repair is already applied. No files changed; repeat the original run command.")
        return
    if saved["state"] != "failed":
        raise ValueError("repair requires a failed managed campaign")
    samples = four.selected_samples(cfg)
    updated = fingerprint_update(cfg, samples, manifest)
    failed, artifacts = check_startup_only(cfg, samples, manifest)
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    completed = {name: completed_production(sample, samples, cfg) for name, sample in groups.items()}
    audit = {"id": MIGRATION, "at": datetime.now(timezone.utc).isoformat(),
             "old_worker_sha256": OLD_WORKER, "new_worker_sha256": FIXED_WORKER,
             "old_configuration_fingerprint": manifest["configuration_fingerprint"],
             "new_configuration_fingerprint": updated["configuration_fingerprint"],
             "failed_workers": failed, "completed_production": completed,
             "startup_artifacts": artifacts, "backup": str(backup), "generators_launched": False}
    print(f"Verified {len(completed)} completed production groups and {len(failed)} logging-startup failures.")
    if not apply:
        print("Audit passed. No files changed; add --apply to archive worker metadata and update the manifest.")
        return audit
    updated["software_migrations"] = [*manifest.get("software_migrations", []), audit]
    updated["status"] = "ready_to_resume"
    backup.mkdir(parents=True, exist_ok=False)
    for source in (cfg.manifest_path, state_path):
        shutil.copy2(source, backup / source.name)
    (backup / "audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    workers = cfg.run_dir / "workers"
    workers.rename(backup / "workers")
    try:
        four.write_manifest(updated, cfg)
    except BaseException:
        (backup / "workers").rename(workers)
        raise
    print("Startup repair applied. All hard events and staged decay inputs are unchanged. Repeat the original run command.")
    return audit


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--work-dir", type=Path, default=cli.REPO / "inputs/pipeline")
    parser.add_argument("--apply", action="store_true")
    options = parser.parse_args(argv)
    work = options.work_dir.expanduser().resolve()
    with cli.installation_lock(work):
        saved, state_path, cfg = saved_configuration(work, options.run_tag)
        repair(saved, state_path, cfg, apply=options.apply)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"Error: {error}")
