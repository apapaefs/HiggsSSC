#!/usr/bin/env python3
"""Audit/archive the specific MG5 loop-combination failure; never run generators.

The original managed options, manifest, completed production markers and event
files must still agree. Only the incomplete failed group is moved to recovery/.
Repeat the original run command afterwards to regenerate that group and finish.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_pipeline as cli
from hfourlepton import higher_order as ho, workflow
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import mg5_loop_recovery as fix


def saved_options(work, tag):
    state_path = work / "campaigns/fourlepton-ho" / (tag + ".json")
    saved = json.loads(state_path.read_text())
    identity = saved["identity"]
    sampling = identity.get("sampling", {})
    if identity["workflow"] != "fourlepton-ho" or sampling.get("preset") != "minimal":
        raise ValueError("this repair covers managed minimal four-lepton campaigns only")
    options = ["run", "--workflow", "fourlepton-ho", "--work-dir", str(work), "--run-tag", tag,
               "--preset", "minimal", "--energy-tev", str(identity["energy_tev"]),
               "--jobs", str(identity["jobs"]), "--event-multiplier", sampling["event_multiplier"]]
    if sampling["continuum_m4l_gev"] is None:
        options.append("--wide-continuum")
    for name in ("signal_lhe", "vbf_lhe", "signal_normalization"):
        if identity[name]:
            options.extend(["--" + name.replace("_", "-"), identity[name]["path"]])
    args = cli.parse_args(options)
    environment = cli.shared.configuration(args)
    if cli.identity(args, environment) != identity:
        raise ValueError("saved campaign options, inputs or environment have changed")
    return saved, state_path, args, environment


def saved_configuration(work, tag):
    saved, state_path, args, _environment = saved_options(work, tag)
    pipeline = SimpleNamespace(args=args, work=work, mg5=work / "MG5_aMC_v3_5_15",
                               python=work / "venv/bin/python")
    cfg = four.parse_config([str(x) for x in workflow.campaign_arguments(pipeline, "perfect")])
    return saved, state_path, cfg


def completed_production(sample, samples, cfg, *, require_count=True, adapter_sha256=None):
    """Validate a completed group without invoking any producer or staging step."""
    needed = four.group_event_count(sample.production_group_id, samples, cfg)
    if sample.matching == "POWHEG":
        source = cfg.signal_lhe if sample.production_mode == "ggf" else cfg.vbf_lhe
        info = ho.inspect_production(source, sample, cfg)
        if info["events"] < needed:
            raise ValueError(f"too few completed events in {source}")
        return dict(four.artifact_stage_record(source), events=info["events"], required_events=needed)
    group = four.production_dir(sample, cfg)
    process_dir = group / "mg5_process"
    marker = group / "fourlepton_production_stage.json"
    if sample.matching == "MC@NLO":
        from hfourlepton.minimal_runtime import production_settings
        process, launch = ho.nlo_cards(sample, cfg, process_dir, needed)
        identity = dict(process_card=process, compile_card=ho.nlo_compile_card(launch), launch_card=launch,
                        adapter_sha256=adapter_sha256 or four.sha256_file(Path(ho.__file__)),
                        guard_sha256=four.sha256_file(Path(ho.__file__).with_name("mg5_guard.py")),
                        runtime=production_settings(sample, cfg))
        source = process_dir / "Events" / cfg.run_tag / "events.lhe.gz"
        ho.validate_nlo_pole_checks(process_dir)
    else:
        identity = four.production_stage_identity(sample, sample, samples, process_dir, cfg)
        source = process_dir / "Events" / cfg.run_tag / "unweighted_events.lhe.gz"
    if not source.is_file() or not four.artifact_stage_matches(marker, identity, source):
        raise ValueError(f"completed production provenance is missing or changed: {group}")
    from hfourlepton.checked_campaign import count_events, validate_retry_receipt
    validate_retry_receipt(marker)
    ho.inspect_production(source, sample, cfg)
    count = count_events(source)
    if require_count and count < needed:
        raise ValueError(f"{source} has {count} events, expected at least {needed}")
    return dict(four.artifact_stage_record(source), events=count, required_events=needed,
                marker=four.artifact_stage_record(marker))


def known_failure(group, tag):
    process = group / "mg5_process"
    if group.is_symlink() or process.is_symlink():
        raise ValueError("refusing to archive a shared production directory")
    if (group / "fourlepton_production_stage.json").exists() or list((process / "Events").rglob("*.lhe*")):
        raise ValueError("a final production artifact exists; this repair cannot replace it")
    source = process / "bin/internal/gen_ximprove.py"
    if fix.digest(source.read_bytes()) != fix.UPSTREAM_SHA256:
        raise ValueError("failed export does not contain the known unpatched MG5 source")
    debug = process / (tag + "_tag_1_debug.log")
    text = debug.read_text()
    if not all(word in text for word in ("gen_ximprove_share.combine_iteration", "FileNotFoundError", "/results.dat")):
        raise ValueError("failure is not the known MG5 loop-combination exception")
    discarded = [path for path in (process / "SubProcesses").glob("P*/DiscardedUnstableEvents/discarded_G*.dat")
                 if path.is_file() and path.stat().st_size]
    if not discarded:
        raise ValueError("no unstable-event recovery evidence; investigate the original error")
    return {"debug": four.artifact_stage_record(debug),
            "discarded_event_records": [four.artifact_stage_record(path) for path in sorted(discarded)],
            "mg5_source_sha256": fix.UPSTREAM_SHA256}


def audit(saved, cfg, sample_name):
    if saved["state"] != "failed":
        raise ValueError("repair requires a failed managed campaign")
    if cfg.force or cfg.dry_run or cfg.preset != "minimal" or cfg.detector_response != "perfect":
        raise ValueError("repair requires the original minimal perfect-profile configuration")
    manifest = json.loads(cfg.manifest_path.read_text())
    if any(s["status"] != "planned" for s in manifest["samples"]) or (cfg.run_dir / "workers").exists():
        raise ValueError("postprocessing has begun; this production-only repair does not apply")
    samples = four.selected_samples(cfg)
    versions = {name: manifest["software"][name + "_version"] for name in ("python", "mg5", "herwig", "root")}
    if four.configuration_fingerprint(cfg, samples, software_versions=versions) != manifest["configuration_fingerprint"]:
        raise ValueError("campaign fingerprint changed; do not bypass its provenance checks")
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    if sample_name not in groups or not sample_name.startswith("bkg_gg_"):
        raise ValueError("select the failed loop-induced continuum production group")
    sample = groups[sample_name]
    group = four.production_dir(sample, cfg)
    failure = known_failure(group, cfg.run_tag)
    card = group / "cards/process_card.dat"
    if card.read_text() != four.process_card_text(sample, group / "mg5_process", cfg):
        raise ValueError("failed process card differs from the campaign")
    completed = {}
    seen_failed = False
    for name, other in groups.items():
        if name == sample_name:
            seen_failed = True
        elif not seen_failed:
            completed[name] = completed_production(other, samples, cfg)
        elif four.production_dir(other, cfg).exists():
            raise ValueError(f"later production already exists: {name}; repair scope differs")
    return {"id": fix.PATCH_ID, "run_tag": cfg.run_tag, "failed_group": sample_name,
            "configuration_fingerprint": manifest["configuration_fingerprint"],
            "manifest": four.artifact_stage_record(cfg.manifest_path),
            "completed_production": completed, "failure": failure, "generators_launched": False}


def archive(group, destination, report, metadata):
    """A reversible rename; existing recovery data must never be overwritten."""
    destination.mkdir(parents=True, exist_ok=False)
    for path in metadata:
        shutil.copy2(path, destination / path.name)
    (destination / "audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    group.rename(destination / "failed_production")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--work-dir", type=Path, default=cli.REPO / "inputs/pipeline")
    parser.add_argument("--apply", action="store_true", help="archive only the failed group; does not generate events")
    options = parser.parse_args(argv)
    work = options.work_dir.expanduser().resolve()
    # The same lock used by run_pipeline prevents repair during a managed run.
    with cli.installation_lock(work):
        saved, state_path, cfg = saved_configuration(work, options.run_tag)
        group = cfg.run_dir / "production" / options.sample
        destination = cfg.run_dir / "recovery" / fix.PATCH_ID / options.sample
        if destination.exists():
            if not (destination / "failed_production").is_dir():
                raise ValueError(f"incomplete recovery archive: {destination}; inspect it before continuing")
            record = json.loads((destination / "audit.json").read_text())
            if record.get("id") != fix.PATCH_ID or record.get("failed_group") != options.sample:
                raise ValueError("unexpected recovery record")
            print("This failure is already archived. No files changed; repeat the original run command.")
            return 0
        report = audit(saved, cfg, options.sample)
        report.update(at=datetime.now(timezone.utc).isoformat(), backup=str(destination))
        print(f"Verified {len(report['completed_production'])} completed production groups.")
        print(f"Failed group: {options.sample}; archive: {destination}")
        if not options.apply:
            print("Audit passed. No campaign files changed; add --apply to archive the failed group.")
            return 0
        report["installation"] = fix.install(cfg.mg5_dir)
        archive(group, destination, report, [cfg.manifest_path, state_path])
        print("Failed group archived; all completed groups are unchanged. Repeat the original run command.")
        return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"Error: {error}")
