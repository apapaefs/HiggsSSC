#!/usr/bin/env python3
"""Audit and repair the known first-WH MadSpin banner failure, without simulation.

This migration is deliberately restricted to a specific pair of runner hashes
and to 13 TeV campaigns with only the six POWHEG samples complete. It is not a
general override for configuration or software provenance mismatches.
"""
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import higher_order as ho

OLD_RUNNER = "23e56d1875acdac9e565abb472a5274fcdba8bb9fbbd8ede9fa8beb3c45b9747"
FIXED_RUNNER = "d8324d2b41e1d3380c52896cf931d005f941e633ba6d1d047924c413bcc3e08b"
MIGRATION = "madspin-no-b-mass-v1"
FAILED_SAMPLE = "signal_wp_h_4e"


def checked_file(path, checksum):
    path = Path(path)
    if not checksum or not path.is_file() or four.sha256_file(path) != checksum:
        raise ValueError(f"missing or changed recovery artifact: {path}")


def verify_completed(entry, sample, cfg):
    if entry["output_root"] != str(four.expected_analysis_output(sample, cfg)):
        raise ValueError("completed output path differs from the campaign configuration")
    for path_key, hash_key in (("output_root", "output_root_sha256"),
                               ("analysis_summary", "analysis_summary_sha256"),
                               ("source_lhe", "source_lhe_sha256"),
                               ("input_lhe", "input_lhe_sha256")):
        checked_file(entry[path_key], entry[hash_key])
    summary = four.read_analysis_summary(Path(entry["analysis_summary"]))
    expected = entry.get("shower_completion", {}).get("saved_events", cfg.events_for(sample))
    if int(summary.get("events_read", -1)) != expected:
        raise ValueError(f"unexpected completed event count: {sample.name}")


def repair(cfg, *, apply=False):
    if cfg.force or cfg.dry_run or cfg.order != "ho" or cfg.ebeam != 6500 or cfg.detector_response != "perfect":
        raise ValueError("repair requires the original 13 TeV HO perfect-profile configuration, without --force")
    if (cfg.run_dir / "manifest_ssc.json").exists():
        raise ValueError("this migration does not cover campaigns with an SSC manifest")
    samples = four.selected_samples(cfg)
    existing = four.load_resumable_manifest(cfg)
    if existing is None:
        raise ValueError("no existing campaign manifest to repair")
    planned = four.build_manifest(cfg, samples)
    if existing["configuration_fingerprint"] == planned["configuration_fingerprint"]:
        if any(item.get("id") == MIGRATION for item in existing.get("software_migrations", [])):
            print("MadSpin model repair is already applied; no files changed.")
            return
        raise ValueError("campaign does not need this specific software migration")
    old_inputs, new_inputs = existing["tracked_inputs"], planned["tracked_inputs"]
    if old_inputs["runner"]["sha256"] != OLD_RUNNER or new_inputs["runner"]["sha256"] != FIXED_RUNNER:
        raise ValueError("runner versions are outside this audited migration")
    previous_inputs = copy.deepcopy(new_inputs)
    previous_inputs["runner"] = old_inputs["runner"]
    if previous_inputs != old_inputs or four.configuration_fingerprint(
        cfg, samples, tracked_inputs=previous_inputs
    ) != existing["configuration_fingerprint"]:
        raise ValueError("configuration, software or inputs changed beyond the MadSpin runner fix")

    entries = {entry["name"]: entry for entry in existing["samples"]}
    if len(entries) != len(existing["samples"]) or set(entries) != {s.name for s in samples}:
        raise ValueError("sample inventory differs from the original campaign")
    completed = []
    for sample in samples:
        entry = entries[sample.name]
        if sample.is_signal and sample.matching == "POWHEG":
            if entry["status"] != "complete":
                raise ValueError("repair expects all six ggF/VBF samples complete")
            verify_completed(entry, sample, cfg)
            completed.append(sample.name)
        elif sample.name == FAILED_SAMPLE:
            if entry["status"] != "running":
                raise ValueError("first WH sample is not at the expected failed stage")
        elif entry["status"] != "planned":
            raise ValueError("other samples have progressed beyond this migration's scope")
    if len(completed) != 6:
        raise ValueError("repair expects the original full signal inventory")

    sample = next(s for s in samples if s.name == FAILED_SAMPLE)
    group = four.production_dir(sample, cfg)
    process_dir = group / "mg5_process"
    source = process_dir / "Events" / cfg.run_tag / "events.lhe.gz"
    count = four.group_event_count(sample.production_group_id, samples, cfg)
    process, launch = ho.nlo_cards(sample, cfg, process_dir, count)
    production_identity = dict(process_card=process, compile_card=ho.nlo_compile_card(launch),
        launch_card=launch, adapter_sha256=four.sha256_file(Path(ho.__file__)),
        guard_sha256=four.sha256_file(Path(ho.__file__).with_name("mg5_guard.py")))
    if not four.artifact_stage_matches(group / "fourlepton_production_stage.json", production_identity, source):
        raise ValueError("WH hard production is incomplete or its provenance differs")
    ho.validate_nlo_pole_checks(process_dir)
    if ho.inspect_production(source, sample, cfg)["events"] < count:
        raise ValueError("WH hard production has too few events")
    sliced = group / "decay_inputs" / f"{sample.channel}.lhe.gz"
    start, stop = ho.slice_bounds(sample, samples, cfg)
    slice_identity = dict(kind="disjoint_hard_event_slice", source=str(source.resolve()),
        source_sha256=four.sha256_file(source), start=start, stop=stop, sample=sample.name,
        adapter_sha256=four.sha256_file(Path(ho.__file__)))
    if not four.artifact_stage_matches(sliced.with_suffix(".slice.json"), slice_identity, sliced):
        raise ValueError("WH decay slice is missing or has changed")

    failed = four.sample_decay_dir(sample, cfg)
    raw, staged = failed / "stable_input_source.lhe", failed / "stable_input.lhe"
    if list(failed.glob("*_decayed*.lhe*")) or (failed / "fourlepton_madspin_stage.json").exists():
        raise ValueError("decayed WH output exists; this migration cannot replace it")
    with four.open_maybe_gzip(sliced) as stream:
        if stream.read() != raw.read_text():
            raise ValueError("failed MadSpin input differs from the validated WH slice")
    marker = failed / "fourlepton_madspin_banner_stage.json"
    banner_identity = json.loads(marker.read_text())["identity"]
    if (banner_identity["runner_sha256"] != OLD_RUNNER
            or banner_identity["source_sha256"] != four.sha256_file(raw)
            or not four.artifact_stage_matches(marker, banner_identity, staged)):
        raise ValueError("failed banner provenance is missing or has changed")
    raw_plan, bad_plan = four.madspin_banner_plan(raw, cfg), four.madspin_banner_plan(staged, cfg)
    if (raw_plan["model"] != "loop_sm-no_b_mass" or not raw_plan["preserved_slha"]
            or "model_requires_zero_yukawa_5" not in bad_plan["slha_validation"]["reasons"]):
        raise ValueError("failed stage is not the known bottom-Yukawa banner mismatch")

    backup = cfg.run_dir / "recovery" / MIGRATION
    audit = {"id": MIGRATION, "at": datetime.now(timezone.utc).isoformat(),
             "old_runner_sha256": OLD_RUNNER, "new_runner_sha256": FIXED_RUNNER,
             "old_configuration_fingerprint": existing["configuration_fingerprint"],
             "new_configuration_fingerprint": planned["configuration_fingerprint"],
             "completed_samples_verified": completed, "wh_production_sha256": four.sha256_file(source),
             "backup": str(backup), "generators_launched": False}
    print(json.dumps(audit, indent=2))
    if not apply:
        print("Audit passed. No files changed; use --apply to archive the failed stage and update the manifest.")
        return audit

    backup.mkdir(parents=True, exist_ok=False)
    shutil.copy2(cfg.manifest_path, backup / cfg.manifest_path.name)
    (backup / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    shutil.copy2(__file__, backup / Path(__file__).name)
    updated = four.merge_manifest_samples(existing, planned)
    updated["software_migrations"] = [*existing.get("software_migrations", []), audit]
    for entry in updated["samples"]:
        if entry["name"] in completed:
            entry["completion_software"] = copy.deepcopy(existing["software"])
        elif entry["name"] == FAILED_SAMPLE:
            entry["status"] = "planned"
    updated["status"] = "ready_to_resume"
    failed.rename(backup / "failed_madspin")
    try:
        four.write_manifest(updated, cfg)
    except BaseException:
        (backup / "failed_madspin").rename(failed)
        raise
    print("Repair applied. Six completed samples and WH hard production are unchanged. Repeat the original run command.")
    return audit


if __name__ == "__main__":
    arguments = sys.argv[1:]
    apply = "--apply" in arguments
    arguments = [arg for arg in arguments if arg != "--apply"]
    try:
        repair(four.parse_config(arguments), apply=apply)
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"Error: {error}")
