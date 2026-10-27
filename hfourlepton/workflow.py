"""Managed four-lepton HO workflow; preparation never runs generators."""
import json
from pathlib import Path
from hgammagamma import run_gammagamma_ho_pipeline as shared
from hgammagamma import run_ihixs_normalization as ihixs
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import higher_order
from hfourlepton import run_powheg_vbf as vbf

REPO = Path(__file__).resolve().parents[1]


def profiles(args):
    return ("perfect",) if args.energy_tev == 13 else ("perfect", "ssc")


def selection(args):
    return "cms" if args.energy_tev == 13 else "ssc"


def configuration(args):
    return json.loads((REPO / f"hfourlepton/campaign_cards/{'lhc13' if args.energy_tev == 13 else 'ssc40'}.json").read_text())


def prepare(args):
    """Write reviewable generator cards and sample inventory, without launching jobs."""
    destination = REPO / "hfourlepton/plans" / args.run_tag
    directory = REPO / "hfourlepton/runs" / args.run_tag
    cfg = four.parse_config(["--order", "ho", "--signal-modes", "all", "--run-tag", args.run_tag,
        "--nevents", str(args.nevents), "--nb-core", str(args.jobs), "--ebeam", str(500 * args.energy_tev),
        "--higgs-mass", "125", "--seed-base", "31122002", "--detector-response", "perfect",
        "--selection-profile", selection(args), "--dry-run"])
    samples = higher_order.selected_samples(cfg)
    groups = {}
    for sample in samples:
        groups.setdefault(sample.production_group_id, sample)
    plan = {"schema_version": 1, "configuration": configuration(args), "run_tag": args.run_tag,
            "events_per_sample": args.nevents, "jobs": args.jobs, "status": "prepared_not_run",
            "higgs_events_per_production_group": 3 * args.nevents,
            "decay_sampling": "disjoint hard-event slices; independent events for 4e, 4mu, 2e2mu",
            "samples": [{"name": s.name, "production_group": s.production_group_id,
                         "process": s.process, "order": s.perturbative_order,
                         "matching": s.matching, "pdf": s.pdf_name} for s in samples]}
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
            process = four.process_card_text(sample, procdir)
            ihixs.write_same_or_new(destination / name / "run_settings.json", json.dumps({
                "nevents": four.group_event_count(name, samples, cfg), "seed": four.seed_for(name, cfg),
                "ebeam_gev": cfg.ebeam, "lhaid": sample.lhaid, "higgs_mass_gev": cfg.higgs_mass,
                "generation_cuts": higher_order.GENERATION_CUTS,
                "applied_by": "run_four_lepton_campaign.patch_run_card"}, indent=2) + "\n")
        ihixs.write_same_or_new(destination / name / "process.mg5", process)
    for sample in samples:
        filename = "unweighted_events.lhe.gz" if sample.matching == "LO" else "events.lhe.gz"
        lhe = four.sample_decay_dir(sample, cfg) / "stable_input_decayed.lhe.gz" if sample.is_signal else directory / "production" / sample.production_group_id / "mg5_process/Events" / args.run_tag / filename
        ihixs.write_same_or_new(destination / sample.name / "herwig.in",
            higher_order.herwig_card(sample, lhe, four.seed_for(f"{sample.name}:herwig", cfg), args.nevents, cfg))
        if sample.is_signal:
            source = four.sample_decay_dir(sample, cfg) / "stable_input.lhe"
            ihixs.write_same_or_new(destination / sample.name / "madspin.dat",
                four.madspin_card_text(source, sample, four.seed_for(f"{sample.name}:madspin", cfg)))
    settings, card = rate_inputs(args)
    ihixs.write_same_or_new(destination / "ihixs.json", json.dumps(settings, indent=2) + "\n")
    ihixs.write_same_or_new(destination / "powheg.input", card)
    print(f"Prepared {len(samples)} samples / {len(groups)} hard-production groups: {destination}")
    return destination


def rate_inputs(args):
    tag = "lhc13" if args.energy_tev == 13 else "ssc40"
    settings = json.loads((shared.HO / "ihixs-ssc40.json").read_text())
    settings.update(profile=f"{tag}-ho-pure-heft-v1", sqrt_s_gev=1000. * args.energy_tev)
    card = shared.powheg.patched_input((shared.HO / "powheg-hjminnlo-ssc40-nnpdf40nnloqed.input").read_text(),
        {"ebeam1": 500 * args.energy_tev, "ebeam2": 500 * args.energy_tev,
         "ncall2": 500000, "check_bad_st2": 1})
    return settings, card


def run(pipeline):
    args = pipeline.args
    pipeline.check()
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
    lhe = args.signal_lhe or directory / "powheg/powheg-hjminnlo-merged.lhe"
    if not args.signal_lhe:
        pipeline.call(args.run_tag + "-powheg", [pipeline.python, REPO / "hgammagamma/run_powheg_hjminnlo.py",
            "--nevents", 3 * args.nevents, "--jobs", args.jobs, "--ebeam", args.energy_tev * 500,
            "--powheg-dir", pipeline.powheg / "HJ/HJMiNNLO", "--run-dir", directory / "powheg",
            "--card", card, "--no-herwig-env", "--no-herwig-module", "--skip-build", "--resume"])
    vbf_lhe = args.vbf_lhe or directory / 'powheg-vbf' / vbf.LHE_NAME
    if not args.vbf_lhe:
        cfg = four.parse_config(['--order', 'ho', '--run-tag', args.run_tag, '--seed-base', '31122002', '--dry-run'])
        pipeline.call(args.run_tag + '-powheg-vbf', [pipeline.python, REPO / 'hfourlepton/run_powheg_vbf.py',
            '--nevents', 3 * args.nevents, '--jobs', args.jobs, '--ebeam', args.energy_tev * 500,
            '--powheg-dir', pipeline.powheg / 'VBF_H', '--run-dir', directory / 'powheg-vbf',
            '--seed', four.seed_for('signal_vbf_h_stable', cfg)])
    for profile in profiles(args):
        pipeline.call(args.run_tag + "-fourlepton-" + profile, [pipeline.python,
            REPO / "hfourlepton/run_four_lepton_campaign.py", "--order", "ho", "--signal-modes", "all",
            "--selection-profile", selection(args), "--detector-response", profile,
            "--run-tag", args.run_tag, "--nevents", args.nevents, "--nb-core", args.jobs,
            "--ebeam", args.energy_tev * 500, "--mg5-dir", pipeline.mg5, "--mg5-python", pipeline.python,
            "--higgs-mass", "125",
            "--no-herwig-env", "--no-herwig-module", "--signal-lhe", lhe,
            "--vbf-lhe", vbf_lhe,
            "--signal-normalization", record])
