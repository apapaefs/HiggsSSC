"""Managed four-lepton HO workflow; preparation never runs generators."""
import json
from dataclasses import replace
from pathlib import Path
from hgammagamma import run_gammagamma_ho_pipeline as shared
from hgammagamma import run_ihixs_normalization as ihixs
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import higher_order
from hfourlepton import presets
from hfourlepton import run_powheg_vbf as vbf

REPO = Path(__file__).resolve().parents[1]


def profiles(args):
    return ("perfect",) if args.energy_tev == 13 else ("perfect", "ssc")


def selection(args):
    return "cms" if args.energy_tev == 13 else "ssc"


def configuration(args):
    return json.loads((REPO / f"hfourlepton/campaign_cards/{'lhc13' if args.energy_tev == 13 else 'ssc40'}.json").read_text())


def planning_config(args):
    event_options = presets.arguments(args) or ["--nevents", str(args.nevents)]
    return four.parse_config(["--order", "ho", "--signal-modes", "all", "--run-tag", args.run_tag,
        *event_options, "--nb-core", str(args.jobs), "--ebeam", str(500 * args.energy_tev),
        "--higgs-mass", "125", "--seed-base", "31122002", "--detector-response", "perfect",
        "--selection-profile", selection(args), "--dry-run"])


def prepare(args):
    """Write reviewable generator cards and sample inventory, without launching jobs."""
    destination = REPO / "hfourlepton/plans" / args.run_tag
    directory = REPO / "hfourlepton/runs" / args.run_tag
    cfg = planning_config(args)
    samples = higher_order.selected_samples(cfg)
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    plan = {"schema_version": 1, "configuration": configuration(args), "run_tag": args.run_tag,
            "events_per_sample": None if cfg.preset == "minimal" else args.nevents,
            "jobs": args.jobs, "status": "prepared_not_run", "sampling": presets.description(args),
            "higgs_events_per_production_group": None if cfg.preset == "minimal" else 3 * args.nevents,
            "production_events": {name: four.group_event_count(name, samples, cfg) for name in groups},
            "decay_sampling": "disjoint hard-event slices; independent events for 4e, 4mu, 2e2mu",
            "samples": [{"name": s.name, "production_group": s.production_group_id,
                         "process": s.process, "order": s.perturbative_order,
                         "matching": s.matching, "pdf": s.pdf_name, "nevents": cfg.events_for(s)} for s in samples]}
    shared.frozen(destination / "plan.json", plan)
    for name, sample in groups.items():
        procdir = directory / "production" / name / "mg5_process"
        if sample.production_mode == "ggf":
            continue
        if sample.production_mode == "vbf":
            ihixs.write_same_or_new(destination / name / "powheg.input",
                vbf.input_card(cfg.ebeam, cfg.higgs_mass, four.group_event_count(name, samples, cfg)))
            continue
        if sample.matching == "MC@NLO":
            process, launch = higher_order.nlo_cards(sample, cfg, procdir, four.group_event_count(name, samples, cfg))
            ihixs.write_same_or_new(destination / name / "compile.mg5", higher_order.nlo_compile_card(launch))
            ihixs.write_same_or_new(destination / name / "launch.mg5", launch)
        else:
            process = four.process_card_text(sample, procdir, cfg)
            ihixs.write_same_or_new(destination / name / "run_settings.json", json.dumps({
                "nevents": four.group_event_count(name, samples, cfg), "seed": four.seed_for(name, cfg),
                "ebeam_gev": cfg.ebeam, "lhaid": sample.lhaid, "higgs_mass_gev": cfg.higgs_mass,
                "generation_cuts": higher_order.GENERATION_CUTS,
                "minimal_runtime": presets.description(cfg),
                "applied_by": "run_four_lepton_campaign.patch_run_card"}, indent=2) + "\n")
        ihixs.write_same_or_new(destination / name / "process.mg5", process)
        if cfg.preset == "minimal":
            from hfourlepton.minimal_runtime import continuum_cut
            order = "7|6|1" if sample.name.startswith("bkg_gg_") else "6|7|1"
            ihixs.write_same_or_new(destination / name / "MadLoopParams.override", f"#MLReductionLib\n{order}\n")
            if not sample.is_signal and not cfg.wide_continuum:
                ihixs.write_same_or_new(destination / name / "dummy_cuts.f", continuum_cut(sample.matching == "MC@NLO"))
    for sample in samples:
        filename = "unweighted_events.lhe.gz" if sample.matching == "LO" else "events.lhe.gz"
        lhe = four.sample_decay_dir(sample, cfg) / "stable_input_decayed.lhe.gz" if sample.is_signal else directory / "production" / sample.production_group_id / "mg5_process/Events" / args.run_tag / filename
        ihixs.write_same_or_new(destination / sample.name / "herwig.in",
            higher_order.herwig_card(sample, lhe, four.seed_for(f"{sample.name}:herwig", cfg), cfg.events_for(sample), cfg))
        if sample.is_signal:
            source = four.sample_decay_dir(sample, cfg) / "stable_input.lhe"
            ihixs.write_same_or_new(destination / sample.name / "madspin.dat",
                four.madspin_card_text(source, sample, four.seed_for(f"{sample.name}:madspin", cfg)))
    settings, card = rate_inputs(args)
    ihixs.write_same_or_new(destination / "ihixs.json", json.dumps(settings, indent=2) + "\n")
    ihixs.write_same_or_new(destination / "powheg.input", card)
    print(f"Prepared {len(samples)} samples / {len(groups)} hard-production groups: {destination}")
    print(f"Total hard events: {sum(cfg.events_for(s) for s in samples):,}")
    if cfg.preset == "minimal":
        window = "wide continuum generation" if cfg.wide_continuum else "90–180 GeV continuum window"
        print(f"Minimal multiplier {args.event_multiplier}; {min(8, args.jobs)} postprocessing workers; {window}")
    return destination


def rate_inputs(args):
    tag = "lhc13" if args.energy_tev == 13 else "ssc40"
    settings = json.loads((shared.HO / "ihixs-ssc40.json").read_text())
    settings.update(profile=f"{tag}-ho-pure-heft-v1", sqrt_s_gev=1000. * args.energy_tev)
    card = shared.powheg.patched_input((shared.HO / "powheg-hjminnlo-ssc40-nnpdf40nnloqed.input").read_text(),
        {"ebeam1": 500 * args.energy_tev, "ebeam2": 500 * args.energy_tev,
         "ncall2": 500000, "check_bad_st2": 1})
    return settings, card


def campaign_arguments(pipeline, profile):
    """One configuration for both a managed run and its bounded metadata repair."""
    args = pipeline.args
    directory = REPO / "hfourlepton/runs" / args.run_tag
    tag = "lhc13" if args.energy_tev == 13 else "ssc40"
    return ["--order", "ho", "--signal-modes", "all",
            "--selection-profile", selection(args), "--detector-response", profile,
            "--run-tag", args.run_tag, *(presets.arguments(args) or ["--nevents", args.nevents]), "--nb-core", args.jobs,
            "--ebeam", args.energy_tev * 500, "--mg5-dir", pipeline.mg5, "--mg5-python", pipeline.python,
            "--higgs-mass", "125", "--no-herwig-env", "--no-herwig-module",
            "--signal-lhe", args.signal_lhe or directory / "powheg/powheg-hjminnlo-merged.lhe",
            "--vbf-lhe", args.vbf_lhe or directory / "powheg-vbf" / vbf.LHE_NAME,
            "--signal-normalization", args.signal_normalization or pipeline.work / f"normalization/ggf-{tag}-n3lo.json"]


def repair_madspin(pipeline):
    pipeline.call(pipeline.args.run_tag + "-repair-madspin", [pipeline.python,
        REPO / "hfourlepton/repair_madspin_banner.py", "--apply",
        *campaign_arguments(pipeline, "perfect")])


def run(pipeline):
    args = pipeline.args
    pipeline.check()
    cfg = planning_config(args)
    samples = higher_order.selected_samples(cfg)
    if cfg.preset == "minimal":
        from hfourlepton.minimal_runtime import validate_libraries
        validate_libraries(pipeline.mg5)
        prepare(args)
    directory = REPO / "hfourlepton/runs" / args.run_tag
    directory.mkdir(parents=True, exist_ok=True)
    tag = "lhc13" if args.energy_tev == 13 else "ssc40"
    settings, card_text = rate_inputs(args)
    settings_path = directory / "ihixs-settings.json"
    card = directory / "powheg-pipeline.input"
    ihixs.write_same_or_new(settings_path, json.dumps(settings, indent=2) + "\n")
    ihixs.write_same_or_new(card, card_text)
    record = args.signal_normalization or pipeline.work / f"normalization/ggf-{tag}-n3lo.json"
    if not record.exists():
        options = ["--settings", settings_path, "--powheg-input", card,
            "--work-dir", pipeline.work / f"normalization/ihixs-{tag}", "--record", record,
            "--lhapdf-dir", pipeline.prefix, "--cuba-dir", pipeline.cuba, "--boost-dir", pipeline.boost,
            "--cc", pipeline.config["cc"], "--cxx", pipeline.config["cxx"], "--jobs", args.jobs, "--resume"]
        for stage in ("build", "benchmark", "calculate"):
            pipeline.call(args.run_tag + "-ihixs-" + stage,
                [pipeline.python, REPO / "hgammagamma/run_ihixs_normalization.py", "--stage", stage,
                 *options, *(["--refine-failed"] if stage == "calculate" else [])])
    checked = shared.norm.load_record(record)
    if checked["sqrt_s_gev"] != args.energy_tev * 1000.:
        raise ValueError("ihixs record has the wrong collider energy")
    if cfg.preset == "minimal":
        # External ggF inputs must be checked against the same normalization
        # record as the downstream campaign, including the cached default.
        runtime_cfg = replace(cfg, dry_run=False, signal_normalization=record)
        for mode, source in (("ggf", args.signal_lhe), ("vbf", args.vbf_lhe)):
            if source:
                sample = next(s for s in samples if s.production_mode == mode)
                info = higher_order.inspect_production(source, sample, runtime_cfg)
                needed = four.group_event_count(sample.production_group_id, samples, cfg)
                if info["events"] < needed:
                    raise ValueError(f"{source} has {info['events']} events; minimal {mode} needs {needed} "
                                     "across its decay channels. Omit this LHE override to generate a new file.")
    if not args.signal_lhe:
        pipeline.call(args.run_tag + "-powheg", [pipeline.python, REPO / "hgammagamma/run_powheg_hjminnlo.py",
            "--nevents", four.group_event_count("signal_gg_h_stable", samples, cfg), "--jobs", args.jobs, "--ebeam", args.energy_tev * 500,
            "--powheg-dir", pipeline.powheg / "HJ/HJMiNNLO", "--run-dir", directory / "powheg",
            "--card", card, "--no-herwig-env", "--no-herwig-module", "--skip-build", "--resume",
            *(["--seed-base", cfg.seed_base] if cfg.preset == "minimal" else [])])
    if not args.vbf_lhe:
        pipeline.call(args.run_tag + '-powheg-vbf', [pipeline.python, REPO / 'hfourlepton/run_powheg_vbf.py',
            '--nevents', four.group_event_count("signal_vbf_h_stable", samples, cfg), '--jobs', args.jobs, '--ebeam', args.energy_tev * 500,
            '--powheg-dir', pipeline.powheg / 'VBF_H', '--run-dir', directory / 'powheg-vbf',
            '--seed', four.seed_for('signal_vbf_h_stable', cfg)])
    for profile in profiles(args):
        pipeline.call(args.run_tag + "-fourlepton-" + profile, [pipeline.python,
            REPO / "hfourlepton/run_four_lepton_campaign.py", *campaign_arguments(pipeline, profile)])
