#!/usr/bin/env python3
"""Audit an undersized loop-continuum sample; --run explicitly regenerates it.

Archive the original group, then retry only that group with MG5's alternative
phase-space strategy and stronger survey. Physics settings and event budgets
are unchanged. No simulation runs unless --run is supplied by the operator.
"""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_pipeline as cli
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import repair_madspin_startup as startup
from hfourlepton import checked_campaign
from hfourlepton.checked_campaign import require_events
from hfourlepton.repair_loop_production import (
    archive, completed_production, saved_configuration, saved_options,
)

RETRY = "loop-short-sample-v1"
SETTINGS = {"sde_strategy": 2, "hard_survey": 1}


def numerical_settings(text):
    values = {}
    for key in SETTINGS:
        matches = re.findall(rf"(?mi)^\s*([0-9]+)\s*=\s*{key}\b", text)
        if len(matches) != 1:
            raise ValueError(f"expected exactly one explicit {key} setting")
        values[key] = int(matches[0])
    return values


def patch_sampling(card):
    text = card.read_text()
    before = numerical_settings(text)
    if before != {"sde_strategy": 1, "hard_survey": 0}:
        raise ValueError(f"unexpected original integration settings in {card}: {before}")
    for key, value in SETTINGS.items():
        text, count = re.subn(rf"(?mi)^(\s*)[0-9]+(\s*=\s*{key}\b)",
                             lambda m: m[1] + str(value) + m[2], text)
        if count != 1:
            raise ValueError(f"ambiguous run-card setting {key}")
    card.write_text(text)


def audit(saved, cfg, sample_name):
    if (saved["state"] != "failed" or cfg.force or cfg.dry_run
            or cfg.preset != "minimal" or cfg.detector_response != "perfect"):
        raise ValueError("retry requires a failed minimal campaign with its original configuration")
    manifest = json.loads(cfg.manifest_path.read_text())
    samples = four.selected_samples(cfg)
    if manifest["tracked_inputs"]["madspin_worker"]["sha256"] == startup.OLD_WORKER:
        startup.fingerprint_update(cfg, samples, manifest)
        startup.check_startup_only(cfg, samples, manifest)
    else:
        versions = {name: manifest["software"][name + "_version"] for name in ("python", "mg5", "herwig", "root")}
        if four.configuration_fingerprint(cfg, samples, software_versions=versions) != manifest["configuration_fingerprint"]:
            raise ValueError("campaign configuration or software changed")
        if any(entry["status"] != "planned" for entry in manifest["samples"]):
            raise ValueError("postprocessing has progressed; this production-only retry does not apply")
        for name in ("workers", "decays", "samples", "analysis"):
            path = cfg.run_dir / name
            if path.exists() and any(path.iterdir()):
                raise ValueError(f"postprocessing output exists: {path}")
        if (cfg.run_dir / "manifest_ssc.json").exists():
            raise ValueError("SSC postprocessing has started")
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    if sample_name not in groups or not sample_name.startswith("bkg_gg_"):
        raise ValueError("select one loop-induced continuum production group")
    sample = groups[sample_name]
    group = four.production_dir(sample, cfg)
    process = group / "mg5_process"
    if group.is_symlink() or process.is_symlink():
        raise ValueError("refusing to archive a shared production directory")
    short = completed_production(sample, samples, cfg, require_count=False)
    if not 0 < short["events"] < short["required_events"]:
        raise ValueError("the selected sample is not an undersized event file")
    if numerical_settings((process / "Cards/run_card.dat").read_text()) != {"sde_strategy": 1, "hard_survey": 0}:
        raise ValueError("the sample already used different integration settings; investigate before retrying")
    if (group / "cards/process_card.dat").read_text() != four.process_card_text(sample, process, cfg):
        raise ValueError("the original process card changed")
    completed, queued = {}, []
    for name, other in groups.items():
        if name == sample_name:
            continue
        if other.matching != "POWHEG" and not four.production_dir(other, cfg).exists():
            queued.append(name)
        else:
            completed[name] = completed_production(other, samples, cfg)
    return dict(id=RETRY, run_tag=cfg.run_tag, sample=sample_name,
                at=datetime.now(timezone.utc).isoformat(),
                manifest=four.artifact_stage_record(cfg.manifest_path),
                short_production=short, completed_production=completed, queued_production=queued,
                previous_settings={"sde_strategy": 1, "hard_survey": 0}, retry_settings=SETTINGS,
                retry_code=four.artifact_stage_record(Path(__file__).resolve()),
                validation_code=four.artifact_stage_record(Path(checked_campaign.__file__).resolve()),
                generators_launched=False)


@contextmanager
def retry_hooks(sample, samples, cfg, group, report):
    """Override only sampling, and reject short output before writing a marker."""
    patch_card, generate = four.patch_run_card, four.run_generate_events
    expected = group / "mg5_process"
    receipt = group / "numerical_retry.json"

    def card_hook(card, **kwargs):
        if card != expected / "Cards/run_card.dat" or kwargs["sample"] != sample:
            raise ValueError("retry attempted to change another production group")
        patch_card(card, **kwargs)
        patch_sampling(card)
        receipt.write_text(json.dumps(dict(report, status="generating",
            run_card=four.artifact_stage_record(card)), indent=2, sort_keys=True) + "\n")

    def generation_hook(process, runtime, config):
        if process != expected or config != cfg:
            raise ValueError("retry attempted to launch another production group")
        generate(process, runtime, config)
        path = four.find_lhe_file(process, cfg.run_tag)
        require_events(path, four.group_event_count(sample.production_group_id, samples, cfg))
        if numerical_settings(four.read_lhe_prefix(path)) != SETTINGS:
            raise ValueError("generated LHE does not record the requested numerical retry settings")

    four.patch_run_card, four.run_generate_events = card_hook, generation_hook
    try:
        yield
    finally:
        four.patch_run_card, four.run_generate_events = patch_card, generate


def execute_retry(sample, samples, cfg, report):
    group = four.production_dir(sample, cfg)
    with retry_hooks(sample, samples, cfg, group, report):
        path = four.obtain_production_lhe(sample, samples, cfg)
    receipt = group / "numerical_retry.json"
    record = json.loads(receipt.read_text())
    record.update(status="complete", artifact=four.artifact_stage_record(path),
                  events=require_events(path, four.group_event_count(sample.production_group_id, samples, cfg)))
    receipt.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    marker = group / "fourlepton_production_stage.json"
    payload = json.loads(marker.read_text())
    payload["numerical_retry"] = four.artifact_stage_record(receipt)
    temporary = marker.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(marker)
    return completed_production(sample, samples, cfg)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--sample", required=True, choices=("bkg_gg_4e", "bkg_gg_4mu", "bkg_gg_2e2mu"))
    parser.add_argument("--work-dir", type=Path, default=cli.REPO / "inputs/pipeline")
    parser.add_argument("--run", action="store_true", help="archive the short group and generate its replacement")
    options = parser.parse_args(argv)
    work = options.work_dir.expanduser().resolve()
    with cli.installation_lock(work):
        saved, state_path, cfg = saved_configuration(work, options.run_tag)
        samples = four.selected_samples(cfg)
        sample = next(s for s in samples if s.name == options.sample)
        group = four.production_dir(sample, cfg)
        destination = cfg.run_dir / "recovery" / RETRY / options.sample
        if destination.exists():
            report = json.loads((destination / "audit.json").read_text())
            if report.get("id") != RETRY or report.get("sample") != options.sample or not (destination / "failed_production").is_dir():
                raise ValueError(f"incomplete retry archive: {destination}")
            # Never overwrite a second failed attempt or launch automatically.
            completed_production(sample, samples, cfg)
            print("The replacement sample is complete; no simulation launched. Continue with the campaign recovery.")
            return
        report = audit(saved, cfg, options.sample)
        short = report["short_production"]
        print(f"Verified {len(report['completed_production'])} complete groups; {options.sample} has "
              f"{short['events']}/{short['required_events']} events. Retry settings: {SETTINGS}.", flush=True)
        if not options.run:
            print("Audit passed. No campaign files changed or generators launched; add --run to regenerate only this group.")
            return report
        # Activate exactly the original managed runtime before archiving data.
        _saved, _state, args, environment = saved_options(work, options.run_tag)
        pipeline = cli.StudentPipeline(args, environment)
        from hfourlepton.minimal_runtime import validate_libraries
        from hfourlepton.mg5_loop_recovery import prepare_campaign
        validate_libraries(pipeline.mg5)
        prepare_campaign(pipeline.mg5, cfg.run_dir)
        archive(group, destination, report, [cfg.manifest_path, state_path])
        report = dict(report, archive=str(destination), generators_launched=True)
        # The lower runner receives the same environment as run_pipeline's child.
        previous = dict(os.environ)
        os.environ.clear()
        os.environ.update(pipeline.env)
        try:
            result = execute_retry(sample, samples, cfg, report)
        finally:
            os.environ.clear()
            os.environ.update(previous)
        (destination / "completion.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(f"Replacement verified: {result['events']} events. Original output and grids: {destination}")
        print("Continue with the MadSpin startup repair if applicable, then repeat the original run command.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"ERROR: {error}")
