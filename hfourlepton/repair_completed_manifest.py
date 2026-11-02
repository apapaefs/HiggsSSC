#!/usr/bin/env python3
"""Re-audit the completed minimal campaign after a resume rewrote its timings.

No simulation, analysis or reporting is run. The frozen ROOT/summary inventory,
configuration, source/shower provenance and completed worker evidence must all
pass before --apply can replace only the managed manifest checksum. Original
metadata and the old checksum remain in a recovery archive.
"""
import argparse
from copy import deepcopy
from dataclasses import asdict, replace
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_pipeline as cli
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import higher_order as ho
from hfourlepton.minimal_runtime import _decode
from hfourlepton.repair_loop_production import saved_configuration, saved_options, completed_production

REPAIR = "completed-manifest-resume-v1"


def frozen_inventory(saved, manifest_path, manifest):
    """Only the manifest may differ from the independently frozen inventory."""
    recorded = saved["analysis_artifacts"]
    expected = {str(manifest_path.resolve())}
    products = {}
    for sample in manifest["samples"]:
        for path_key, hash_key in (("root_file", "output_root_sha256"),
                                   ("analysis_summary", "analysis_summary_sha256")):
            path = Path(sample[path_key]).resolve()
            name = str(path)
            if name in expected:
                raise ValueError(f"duplicate completed product: {name}")
            expected.add(name)
            digest = cli.shared.digest(path)
            if digest != recorded.get(name) or digest != sample.get(hash_key):
                raise ValueError(f"completed ROOT/summary changed: {path}")
            products[name] = digest
    if set(recorded) != expected:
        raise ValueError("frozen analysis inventory differs from the completed campaign")
    old = recorded[str(manifest_path.resolve())]
    current = cli.shared.digest(manifest_path)
    if not re.fullmatch(r"[0-9a-f]{64}", str(old)) or old == current:
        raise ValueError("no changed manifest checksum matching this repair")
    return products, old, current


def physics_fields(manifest):
    result = deepcopy(manifest)
    result.pop("status", None)
    result.pop("production_wall_seconds", None)
    for entry in result["samples"]:
        entry.pop("last_attempt_wall_seconds", None)
    return result


def worker_evidence(cfg, samples, manifest):
    """Each saved job must already have been complete before the faulty resume."""
    expected = physics_fields(manifest)
    files = []
    for sample in samples:
        directory = cfg.run_dir / "workers/perfect"
        job_path = directory / (sample.name + ".job.json")
        journal = directory / (sample.name + ".json")
        log = directory / (sample.name + ".log")
        job = json.loads(job_path.read_text(), object_hook=_decode)
        private = json.loads(journal.read_text())
        expected_cfg = replace(cfg, postprocess_only=True, manifest_override=journal)
        if job["sample"] != sample.name or job["config"] != asdict(expected_cfg):
            raise ValueError(f"saved worker configuration changed: {job_path}")
        for value in (job["manifest"], private):
            if (value.get("status") not in ("running", "complete")
                    or any(s.get("status") != "complete" for s in value["samples"])
                    or physics_fields(value) != expected):
                raise ValueError(f"worker evidence differs beyond timing/status: {sample.name}")
        if not log.read_text().rstrip().endswith(f"==> Resuming: {sample.name} is already complete"):
            raise ValueError(f"latest worker did not only resume completed output: {log}")
        files.extend((job_path, journal, log))
    return files


def audit(saved, state_path, cfg):
    manifest_path = cfg.manifest_path
    manifest = json.loads(manifest_path.read_text())
    if (cfg.force or cfg.dry_run or cfg.preset != "minimal" or cfg.detector_response != "perfect"
            or cfg.ebeam != 6500. or (cfg.run_dir / "manifest_ssc.json").exists()):
        raise ValueError("repair requires the original minimal 13 TeV perfect-profile campaign")
    expected_stage = cfg.run_tag + "-fourlepton-perfect"
    expected_error = (f"completed analysis product changed: {manifest_path}; "
                      "preserve it and use a new run tag")
    if (saved.get("state") != "failed" or saved.get("error") != expected_error
            or saved.get("stage") != expected_stage or saved.get("last_completed_stage") != expected_stage):
        raise ValueError("managed state is not the known post-resume manifest-checksum failure")
    samples = four.selected_samples(cfg)
    if (manifest.get("status") != "complete" or len(manifest["samples"]) != len(samples)
            or {e["name"] for e in manifest["samples"]} != {s.name for s in samples}
            or any(e["status"] != "complete" for e in manifest["samples"])):
        raise ValueError("all campaign samples must already be complete")
    current = four.campaign_tracked_inputs(cfg)
    versions = {name: manifest["software"][name + "_version"] for name in ("python", "mg5", "herwig", "root")}
    if (manifest["tracked_inputs"] != current or four.configuration_fingerprint(cfg, samples,
            tracked_inputs=current, software_versions=versions) != manifest["configuration_fingerprint"]):
        raise ValueError("campaign configuration, software or inputs changed")
    products, old_hash, new_hash = frozen_inventory(saved, manifest_path, manifest)
    files = worker_evidence(cfg, samples, manifest)
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    production = {name: completed_production(sample, samples, cfg) for name, sample in groups.items()}
    entries = {entry["name"]: entry for entry in manifest["samples"]}
    for sample in samples:
        entry = entries[sample.name]
        if entry["nevents_requested"] != cfg.events_for(sample):
            raise ValueError(f"event budget changed: {sample.name}")
        # This branch validates analysis, source and decay hashes and returns;
        # no producer or analyzer is invoked for an existing completed output.
        four.run_sample(sample, samples, manifest, cfg)
        directory = four.shower_work_dir(sample, cfg) / "herwig"
        completion = entry["shower_completion"]
        if (json.loads((directory / "matched_completion.json").read_text()) != completion
                or Path(completion["herwig_dir"]).resolve() != directory.resolve()):
            raise ValueError(f"completed shower marker changed: {sample.name}")
        ho.completion.validate_completion(completion, {"sample": sample.name,
            "nevents_requested": cfg.events_for(sample), "lhe_file": entry["input_lhe"],
            "lhe": {"events": completion["lhe_events"]}}, directory=directory)
    report = {"id": REPAIR, "at": datetime.now(timezone.utc).isoformat(),
              "manifest": str(manifest_path), "old_manifest_sha256": old_hash,
              "new_manifest_sha256": new_hash, "unchanged_products": products,
              "configuration_fingerprint": manifest["configuration_fingerprint"],
              "completed_production": production, "samples": [s.name for s in samples],
              "worker_evidence": [four.artifact_stage_record(p) for p in files],
              "original_manifest_bytes_available": False,
              "validation": "frozen outputs, configured inputs, shower proofs and complete worker manifests",
              "generators_launched": False, "manifest_rewritten": False,
              "managed_state": four.artifact_stage_record(state_path)}
    return report, files


def repair(saved, state_path, cfg, args, *, apply=False):
    backup = cfg.run_dir / "recovery" / REPAIR
    if any(item.get("id") == REPAIR for item in saved.get("metadata_repairs", [])):
        cli.physics_artifacts(args)
        if not (backup / "managed_campaign.json").is_file() or not (backup / "audit.json").is_file():
            raise ValueError("previous manifest repair evidence is missing")
        print("Completed-manifest repair is already applied. Repeat the original run command to rebuild the report.")
        return
    if backup.exists():
        raise ValueError(f"existing recovery archive requires inspection: {backup}")
    state_bytes = state_path.read_bytes()
    report, files = audit(saved, state_path, cfg)
    print(f"Verified {len(report['unchanged_products'])} unchanged ROOT/summary products, "
          f"{len(report['completed_production'])} hard-production groups and {len(report['samples'])} completed showers.")
    if not apply:
        print("Audit passed; no files changed. Add --apply to archive the evidence and update only the managed manifest checksum.")
        return report
    if (state_path.read_bytes() != state_bytes
            or cli.shared.digest(cfg.manifest_path) != report["new_manifest_sha256"]):
        raise ValueError("campaign metadata changed during the recovery audit")
    backup.mkdir(parents=True)
    shutil.copy2(state_path, backup / "managed_campaign.json")
    for path in [cfg.manifest_path, *files]:
        target = backup / path.relative_to(cfg.run_dir)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    (backup / "audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    updated = deepcopy(saved)
    updated["analysis_artifacts"][str(cfg.manifest_path.resolve())] = report["new_manifest_sha256"]
    updated.setdefault("metadata_repairs", []).append(report)
    updated.update(state="prepared", stage="report", error=None, updated_at=report["at"])
    try:
        cli.shared.save(state_path, updated)
        cli.physics_artifacts(args)
    except BaseException:
        shutil.copy2(backup / "managed_campaign.json", state_path)
        raise
    print("Checksum repair applied. No manifest, ROOT, LHE or report was rewritten. "
          "Repeat the original run command; it will resume directly at reporting.")
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
        _, _, args, _ = saved_options(work, options.run_tag)
        repair(saved, state_path, cfg, args, apply=options.apply)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"Error: {error}")
