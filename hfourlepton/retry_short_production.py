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
CARD_RECOVERY = "hidden-survey-card-startup-v1"
OLD_RETRY_SHA256 = "916826a949c0167a910d2f3f5119c596228fc98616b9f512b7ee62a26b1f1c01"


def numerical_settings(text, *, fresh_card=False):
    values = {}
    for key in SETTINGS:
        matches = re.findall(rf"(?mi)^[ \t]*([^\n=!#]+)=[ \t]*{key}\b", text)
        if not matches and fresh_card and key == "hard_survey":
            # MG5 3.5.15 RunCardLO declares this hidden parameter with default 0.
            values[key] = 0
            continue
        if len(matches) != 1 or not re.fullmatch(r"[0-9]+", matches[0].strip()):
            raise ValueError(f"expected exactly one explicit {key} setting")
        values[key] = int(matches[0])
    return values


def patch_sampling(card):
    text = card.read_text()
    before = numerical_settings(text, fresh_card=True)
    if before != {"sde_strategy": 1, "hard_survey": 0}:
        raise ValueError(f"unexpected original integration settings in {card}: {before}")
    for key, value in SETTINGS.items():
        text, count = re.subn(rf"(?mi)^([ \t]*)[0-9]+([ \t]*=[ \t]*{key}\b)",
                             lambda m: m[1] + str(value) + m[2], text)
        if count == 0 and key == "hard_survey":
            text = text.rstrip("\n") + f"\n {value} = {key} ! explicit numerical retry setting\n"
            count = 1
        if count != 1:
            raise ValueError(f"ambiguous run-card setting {key}")
    if numerical_settings(text) != SETTINGS:
        raise ValueError(f"failed to make numerical retry settings explicit in {card}")
    card.write_text(text)


def campaign_samples(saved, cfg):
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
    return samples


def audit(saved, cfg, sample_name):
    samples = campaign_samples(saved, cfg)
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


def same_content(path, record):
    actual = four.artifact_stage_record(path)
    return all(actual[key] == record[key] for key in ("sha256", "size"))


def audit_card_startup(saved, cfg, sample, destination, report):
    """Recognize only the old parser failure, before compilation/integration.

    This permits one additional export attempt, archiving the unfinished export
    alongside the original short sample. A numerical failure cannot use it.
    """
    samples = campaign_samples(saved, cfg)
    group = four.production_dir(sample, cfg)
    process = group / "mg5_process"
    if (sample.name != "bkg_gg_4mu" or report["retry_code"]["sha256"] != OLD_RETRY_SHA256
            or report["retry_settings"] != SETTINGS
            or report["manifest"] != four.artifact_stage_record(cfg.manifest_path)):
        raise ValueError("startup retry evidence or campaign metadata changed")
    if (destination / CARD_RECOVERY).exists():
        raise ValueError("the card-startup recovery was already attempted; preserve and inspect that attempt")
    if group.is_symlink() or process.is_symlink():
        raise ValueError("refusing to archive a shared production directory")
    if ({p.name for p in group.iterdir()} != {"cards", "mg5_process"}
            or (process / "Source/run_card.inc").exists()
            or any(p.name != ".keep" for p in (process / "Events").iterdir())
            or any(p.is_dir() for p in (process / "SubProcesses").glob("P*/G[0-9]*"))):
        raise ValueError("retry progressed beyond the run-card startup failure; preserve its outputs")
    card = process / "Cards/run_card.dat"
    text = card.read_text()
    if (numerical_settings(text, fresh_card=True) != {"sde_strategy": 1, "hard_survey": 0}
            or re.search(r"(?mi)^[^!#\n]*=\s*hard_survey\b", text)):
        raise ValueError("failed export is not the missing hidden hard_survey case")
    if (group / "cards/process_card.dat").read_text() != four.process_card_text(sample, process, cfg):
        raise ValueError("failed export process card changed")
    log = Path(saved["identity"]["environment"]["work_dir"]) / "logs" / (cfg.run_tag + "-retry-gg4mu.log")
    if "ERROR: expected exactly one explicit hard_survey setting" not in log.read_text().splitlines():
        raise ValueError(f"unrecognized startup failure in {log}")
    old = destination / "failed_production"
    source = old / "mg5_process/Events" / cfg.run_tag / "unweighted_events.lhe.gz"
    marker = old / "fourlepton_production_stage.json"
    short = report["short_production"]
    expected = four.production_stage_identity(sample, sample, samples, process, cfg)
    archived_marker = json.loads(marker.read_text())
    if (not same_content(source, short) or not same_content(marker, short["marker"])
            or archived_marker["identity"] != expected
            or archived_marker["artifact"] != {key: short[key] for key in ("path", "sha256", "size")}
            or not 0 < checked_campaign.count_events(source) == short["events"] < expected["events"]):
        raise ValueError("original short-sample archive changed")
    groups = {s.production_group_id: s for s in reversed(samples)}
    if set(report["completed_production"]) | set(report["queued_production"]) != set(groups) - {sample.production_group_id}:
        raise ValueError("production inventory differs from the archived retry audit")
    for name, previous in report["completed_production"].items():
        if completed_production(groups[name], samples, cfg) != previous:
            raise ValueError(f"previously completed production changed: {name}")
    for name in report["queued_production"]:
        if four.production_dir(groups[name], cfg).exists():
            raise ValueError(f"unexpected new production in {name}")
    return dict(report, retry_code=four.artifact_stage_record(Path(__file__).resolve()),
                card_startup_recovery=dict(id=CARD_RECOVERY,
                    original_retry_code=report["retry_code"],
                    failed_card=four.artifact_stage_record(card), log=four.artifact_stage_record(log)),
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
        card_recovery = False
        if destination.exists():
            report = json.loads((destination / "audit.json").read_text())
            if report.get("id") != RETRY or report.get("sample") != options.sample or not (destination / "failed_production").is_dir():
                raise ValueError(f"incomplete retry archive: {destination}")
            if (group / "fourlepton_production_stage.json").exists():
                completed_production(sample, samples, cfg)
                print("The replacement sample is complete; no simulation launched. Continue with the campaign recovery.")
                return
            report = audit_card_startup(saved, cfg, sample, destination, report)
            card_recovery = True
        else:
            report = audit(saved, cfg, options.sample)
        short = report["short_production"]
        summary = (f"Verified {len(report['completed_production'])} complete groups; {options.sample} has "
                   f"{short['events']}/{short['required_events']} events. Retry settings: {SETTINGS}.")
        if not options.run:
            if card_recovery:
                print("Verified the missing hidden run-card setting failure; integration has not started.")
            print(summary)
            print("Audit passed. No campaign files changed or generators launched; add --run to regenerate only this group.")
            return report
        # Activate exactly the original managed runtime before archiving data.
        _saved, _state, args, environment = saved_options(work, options.run_tag)
        pipeline = cli.StudentPipeline(args, environment)
        from hfourlepton.minimal_runtime import validate_libraries
        from hfourlepton.mg5_loop_recovery import prepare_campaign
        validate_libraries(pipeline.mg5)
        prepare_campaign(pipeline.mg5, cfg.run_dir)
        archive_destination = destination / CARD_RECOVERY if card_recovery else destination
        metadata = [cfg.manifest_path, state_path]
        if card_recovery:
            evidence = report["card_startup_recovery"]["log"]
            log = Path(evidence["path"])
            if not same_content(log, evidence):
                raise ValueError("startup log changed during the recovery audit")
            metadata.append(log)
        archive(group, archive_destination, report, metadata)
        # The launcher tees stdout to the original log; copy its evidence
        # before appending our own recovery messages.
        print(summary, flush=True)
        if card_recovery:
            print(f"Archived the pre-integration export: {archive_destination}", flush=True)
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
