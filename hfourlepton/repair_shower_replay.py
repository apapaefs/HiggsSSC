#!/usr/bin/env python3
"""Archive the audited ThePEG final-event replay failure; never run simulations.

Only the known shower-launcher change is compatible with these existing hard
events. Preserve all original metadata and failed showers, migrate their stage
identities explicitly, and leave the original campaign ready for a user restart.
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
from hfourlepton import higher_order as ho
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton.repair_loop_production import saved_configuration, completed_production

MIGRATION = "finite-lhe-shower-v1"
OLD_ADAPTER = "d0b58c9da6de8e7c0d4bb911e449b2b47248709f539a4d5d4267ad81a207f96e"
FIXED_ADAPTER = "8add7d6a78eb6f5b14b956089144285aee471ad7504ba17088274061871afb45"


def fingerprint_update(cfg, samples, manifest):
    current = four.campaign_tracked_inputs(cfg)
    previous = manifest["tracked_inputs"]
    if (previous["higher_order_adapter"]["sha256"] != OLD_ADAPTER
            or current["higher_order_adapter"]["sha256"] != FIXED_ADAPTER):
        raise ValueError("shower launcher versions are outside this audited repair")
    expected_old = deepcopy(current)
    expected_old["higher_order_adapter"] = previous["higher_order_adapter"]
    versions = {name: manifest["software"][name + "_version"] for name in ("python", "mg5", "herwig", "root")}
    if (expected_old != previous or four.configuration_fingerprint(cfg, samples,
            tracked_inputs=expected_old, software_versions=versions) != manifest["configuration_fingerprint"]):
        raise ValueError("configuration, software or inputs changed beyond the shower launcher fix")
    updated = deepcopy(manifest)
    updated.update(tracked_inputs=current, configuration_fingerprint=four.configuration_fingerprint(
        cfg, samples, tracked_inputs=current, software_versions=versions))
    return updated


def stage_updates(cfg, samples):
    """Verify unchanged producer/slicing identities before migrating their hash."""
    changes = {}
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    for sample in groups.values():
        group = four.production_dir(sample, cfg)
        if sample.matching == "MC@NLO":
            # The complete identity, native event count and artifact are audited
            # by completed_production before this function is called.
            marker = group / "fourlepton_production_stage.json"
            changes[marker] = json.loads(marker.read_text())
        if sample.matching == "POWHEG":
            source = cfg.signal_lhe if sample.production_mode == "ggf" else cfg.vbf_lhe
            info = ho.inspect_production(source, sample, cfg)
            marker = group / "powheg_init_preparation.json"
            destination = group / "stable_higgs_for_madspin.lhe.gz"
            expected = {"source": str(source), "source_sha256": four.sha256_file(source),
                        "native_cross_section_pb": info["cross_section_pb"],
                        "native_error_pb": info["cross_section_error_pb"], "adapter_sha256": OLD_ADAPTER}
            if not four.artifact_stage_matches(marker, expected, destination):
                raise ValueError(f"prepared POWHEG input changed: {marker}")
            changes[marker] = json.loads(marker.read_text())
    for sample in samples:
        if not sample.is_signal:
            continue
        group = four.production_dir(sample, cfg)
        destination = group / "decay_inputs" / f"{sample.channel}.lhe.gz"
        marker = destination.with_suffix(".slice.json")
        if not destination.exists() and not marker.exists():
            continue
        source = (group / "stable_higgs_for_madspin.lhe.gz" if sample.matching == "POWHEG"
                  else group / "mg5_process/Events" / cfg.run_tag / (
                      "events.lhe.gz" if sample.matching == "MC@NLO" else "unweighted_events.lhe.gz"))
        start, stop = ho.slice_bounds(sample, samples, cfg)
        expected = {"kind": "disjoint_hard_event_slice", "source": str(source.resolve()),
                    "source_sha256": four.sha256_file(source), "start": start, "stop": stop,
                    "sample": sample.name, "adapter_sha256": OLD_ADAPTER}
        if not four.artifact_stage_matches(marker, expected, destination):
            raise ValueError(f"hard-event slice changed: {marker}")
        changes[marker] = json.loads(marker.read_text())
    for marker, record in changes.items():
        if record["identity"].get("adapter_sha256") != OLD_ADAPTER:
            raise ValueError(f"unexpected adapter identity: {marker}")
        record["identity"]["adapter_sha256"] = FIXED_ADAPTER
        record.setdefault("software_migrations", []).append(MIGRATION)
    return changes


def inspect_replay(sample, cfg):
    directory = four.shower_work_dir(sample, cfg) / "herwig"
    if (directory / "matched_completion.json").exists() or (directory / "shower_execution.json").exists():
        raise ValueError(f"shower is complete or already uses the EOF fix: {directory}")
    parsed = ho.completion._parse_logs(directory, sample.name, "requested_events")
    source = Path(ho.completion._audit_source_card(directory / (sample.name + ".in")))
    weights = ho.completion.source_weight_summary(source, cfg.events_for(sample))
    count = weights["lhe_events"]
    if (count != cfg.events_for(sample) or parsed["saved_events"] != count
            or parsed["attempted_events"] != count + 1 or parsed["exception_counts"]["eventerror"] != 1):
        raise ValueError(f"not the known one-event finite-source replay: {directory}")
    card = ho.herwig_card(sample, source, four.seed_for(sample.name + ":herwig", cfg), count, cfg)
    if (directory / (sample.name + ".in")).read_text() != card:
        raise ValueError(f"shower card differs from the original campaign: {directory}")
    metadata = json.loads((directory / "root-metadata.log").read_text().splitlines()[-1])
    if ho.completion._root_summary(metadata)["entries"] != count:
        raise ValueError(f"ROOT metadata differs from finalized counts: {directory}")
    roots = sorted((directory / "events").rglob("*.root"))
    if {str(p.resolve()) for p in roots} != set(parsed["root_writers"]):
        raise ValueError(f"ROOT inventory differs from the finalized shower: {directory}")
    if {r["path"] for r in metadata} != {str(p.resolve()) for p in roots}:
        raise ValueError(f"ROOT metadata inventory differs: {directory}")
    for entry in metadata:
        stat = Path(entry["path"]).stat()
        if (stat.st_size, stat.st_mtime_ns) != (entry["size"], entry["mtime_ns"]):
            raise ValueError(f"ROOT changed since metadata inspection: {entry['path']}")
    artifacts = []
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"refusing to archive shared shower artifacts: {path}")
        if path.is_file():
            artifacts.append(four.artifact_stage_record(path))
    return {"sample": sample.name, "source_weights": weights, "counts": {
        key: parsed[key] for key in ("generated_events", "attempted_events", "saved_events", "exception_counts")},
        "artifacts": artifacts}


def check_samples(cfg, samples, manifest):
    entries = {s["name"]: s for s in manifest["samples"]}
    if len(entries) != len(manifest["samples"]) or set(entries) != {s.name for s in samples}:
        raise ValueError("sample inventory changed")
    failures, completed = [], []
    for sample in samples:
        entry = entries[sample.name]
        journal = cfg.run_dir / "workers/perfect" / (sample.name + ".json")
        if journal.exists():
            private = json.loads(journal.read_text())
            if (private["configuration_fingerprint"] != manifest["configuration_fingerprint"]
                    or private["tracked_inputs"] != manifest["tracked_inputs"]
                    or next(s for s in private["samples"] if s["name"] == sample.name) != entry):
                raise ValueError(f"worker evidence differs from manifest: {journal}")
        output = four.expected_analysis_output(sample, cfg)
        directory = four.shower_work_dir(sample, cfg) / "herwig"
        if entry["status"] == "complete":
            if not output.is_file():
                raise ValueError(f"completed analysis is missing: {output}")
            # The completed branch only checks hashes/counts and returns; no
            # simulation can start after the status and output checks above.
            four.run_sample(sample, samples, manifest, cfg)
            saved = json.loads((directory / "matched_completion.json").read_text())
            if saved != entry.get("shower_completion"):
                raise ValueError(f"completed shower differs from manifest: {sample.name}")
            ho.completion.validate_completion(saved, {"sample": sample.name,
                "nevents_requested": cfg.events_for(sample), "lhe": {"events": saved["lhe_events"]}},
                directory=directory)
            completed.append(sample.name)
        elif entry["status"] == "failed":
            if not journal.exists() or output.exists():
                raise ValueError(f"unexpected failed-sample progress: {sample.name}")
            log = journal.with_suffix(".log").read_text()
            if "inconsistent source/generated/saved shower counts" not in log:
                raise ValueError(f"unrecognized worker failure: {sample.name}")
            failures.append(inspect_replay(sample, cfg))
            decay = four.sample_decay_dir(sample, cfg)
            if sample.is_signal:
                source = four.production_dir(sample, cfg) / "decay_inputs" / f"{sample.channel}.lhe.gz"
                output = four.madspin_output_path(decay / "stable_input.lhe")
                identity = four.madspin_stage_identity(sample, source, decay / "stable_input.lhe",
                                                       decay / "madspin_card.dat", cfg)
                if not four.artifact_stage_matches(decay / "fourlepton_madspin_stage.json", identity, output):
                    raise ValueError(f"completed decay provenance changed: {sample.name}")
        elif entry["status"] != "planned" or directory.exists() or journal.exists():
            raise ValueError(f"unexpected sample state: {sample.name}")
    if not failures:
        raise ValueError("no finite-source replay failures to repair")
    return failures, completed


def repair(saved, state_path, cfg, *, apply=False):
    if (cfg.force or cfg.dry_run or cfg.preset != "minimal" or cfg.detector_response != "perfect"
            or (cfg.run_dir / "manifest_ssc.json").exists()):
        raise ValueError("repair requires the original minimal perfect-profile campaign")
    manifest = json.loads(cfg.manifest_path.read_text())
    backup = cfg.run_dir / "recovery" / MIGRATION
    if any(item.get("id") == MIGRATION for item in manifest.get("software_migrations", [])):
        if (four.sha256_file(Path(ho.__file__)) != FIXED_ADAPTER
                or not (backup / "audit.json").is_file() or not (backup / "workers").is_dir()):
            raise ValueError("existing shower repair evidence is missing or changed")
        print("Finite-shower repair is already applied. No files changed; repeat the original run command.")
        return
    if saved["state"] != "failed" or backup.exists():
        raise ValueError("repair requires a failed campaign with no previous repair archive")
    samples = four.selected_samples(cfg)
    updated = fingerprint_update(cfg, samples, manifest)
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    production = {name: completed_production(sample, samples, cfg, adapter_sha256=OLD_ADAPTER)
                  for name, sample in groups.items()}
    changes = stage_updates(cfg, samples)
    failures, completed = check_samples(cfg, samples, manifest)
    report = {"id": MIGRATION, "at": datetime.now(timezone.utc).isoformat(),
              "old_adapter_sha256": OLD_ADAPTER, "new_adapter_sha256": FIXED_ADAPTER,
              "old_configuration_fingerprint": manifest["configuration_fingerprint"],
              "new_configuration_fingerprint": updated["configuration_fingerprint"],
              "completed_production": production, "completed_samples": completed,
              "failed_showers": failures, "generators_launched": False,
              "stage_markers": [four.artifact_stage_record(p) for p in changes]}
    print(f"Verified {len(production)} hard-production groups, {len(completed)} completed analyses, "
          f"and {len(failures)} failed showers.")
    if not apply:
        print("Audit passed; no files changed. Add --apply to archive only the failed showers and update metadata.")
        return report
    updated.setdefault("software_migrations", []).append(report)
    updated["status"] = "ready_to_resume"
    failed_names = {item["sample"] for item in failures}
    for entry in updated["samples"]:
        if entry["name"] in failed_names:
            entry["status"] = "planned"
            entry.pop("error", None)
            entry.pop("last_attempt_wall_seconds", None)
    backup.mkdir(parents=True)
    metadata = [cfg.manifest_path, *changes]
    for path in metadata:
        target = backup / "metadata" / path.relative_to(cfg.run_dir)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    shutil.copy2(state_path, backup / "managed_campaign.json")
    (backup / "audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    moves = [(cfg.run_dir / "workers", backup / "workers")]
    for sample in samples:
        if sample.name in failed_names:
            moves.append((four.shower_work_dir(sample, cfg) / "herwig", backup / "showers" / sample.name))
    moved = []
    try:
        for source, destination in moves:
            destination.parent.mkdir(parents=True, exist_ok=True)
            source.rename(destination)
            moved.append((source, destination))
        for path, record in changes.items():
            temporary = path.with_suffix(path.suffix + ".tmp")
            temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            temporary.replace(path)
        four.write_manifest(updated, cfg)
    except BaseException:
        for path in metadata:
            shutil.copy2(backup / "metadata" / path.relative_to(cfg.run_dir), path)
        for source, destination in reversed(moved):
            destination.rename(source)
        raise
    print("Repair applied. Hard events, decays and completed analyses are unchanged. "
          "Repeat the original run command to shower the failed samples and finish the campaign.")
    return report


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
