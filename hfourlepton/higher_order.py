"""Four-lepton production and matched-shower adapters using the diphoton machinery.

No inclusive K factor is applied to non-ggF samples. NLO rates are the native
MC@NLO integrals; loop-induced components are explicitly LO. Event generation
and showering happen only when the user runs the campaign.
"""
from dataclasses import replace
import gzip
import json
import math
import re
from pathlib import Path
import subprocess
from string import Template
from types import SimpleNamespace

from hgammagamma import run_gammagamma_ho_campaign as ho
from hgammagamma import ho_signal_normalization as norm
from hgammagamma import ho_shower_completion as completion
from hfourlepton import run_four_lepton_campaign as four

NLO_PDF = "NNPDF40_nlo_as_01180_qed"
NLO_ID = 335900
GENERATION_CUTS = {"ptl_min_gev": 2., "etal_max": 5., "mmll_min_gev": 2., "drll_min": 0.}
# The VBF definition excludes s-channel VH. WH charges and ggZH are disjoint.
PRODUCTIONS = {
    "gg_h": ("ggf", "HJMiNNLO", "NNLO+PS", "POWHEG"),
    "vbf_h": ("vbf", "p p > h j j $$ w+ w- z QCD=0 [QCD]", "NLO+PS", "MC@NLO"),
    "wp_h": ("wh", "p p > h w+ [QCD]", "NLO+PS", "MC@NLO"),
    "wm_h": ("wh", "p p > h w- [QCD]", "NLO+PS", "MC@NLO"),
    "z_h": ("zh", "p p > h z [QCD]", "NLO+PS", "MC@NLO"),
    "gg_z_h": ("zh", "g g > h z [noborn=QCD]", "LO_loop_induced", "LO"),
    "tt_h": ("tth", "p p > t t~ h [QCD]", "NLO+PS", "MC@NLO"),
}


def all_samples():
    result = []
    for prefix, (mode, process, order, matching) in PRODUCTIONS.items():
        for original in four.SAMPLES[:3]:
            result.append(replace(original, name=f"signal_{prefix}_{original.channel}",
                label=f"{prefix}, H->{original.channel}", process=process,
                model="loop_sm-no_b_mass", perturbative_order=order, matching=matching,
                source_kind="powheg" if mode == "ggf" else "mg5",
                production_mode=mode, production_group_id=f"signal_{prefix}_stable",
                normalization_rule="ihixs_n3lo_times_h4l_br" if mode == "ggf" else "native_production_times_h4l_br",
                pdf_name=norm.NATIVE_PDF_SET if mode == "ggf" else NLO_PDF,
                lhaid=norm.NATIVE_PDF_ID if mode == "ggf" else NLO_ID))
    for original in four.SAMPLES[3:9]:
        nlo = original.name.startswith("bkg_qq_")
        result.append(replace(original, process=original.process + (" [QCD]" if nlo else ""),
            model="loop_sm-no_b_mass", perturbative_order="NLO+PS" if nlo else "LO_loop_induced",
            matching="MC@NLO" if nlo else "LO", pdf_name=NLO_PDF, lhaid=NLO_ID))
    return result


def selected_samples(cfg):
    if cfg.sample_set != "irreducible" or cfg.heavy_flavour_bias != "unbiased":
        raise ValueError("HO nominal samples are irreducible; run reducible controls separately with --order lo")
    modes = set(PRODUCTIONS[p][0] for p in PRODUCTIONS) if cfg.signal_modes == "all" else set(cfg.signal_modes.split(","))
    if not modes or modes - {"ggf", "vbf", "wh", "zh", "tth"}:
        raise ValueError("signal modes must be ggf,vbf,wh,zh,tth or all")
    samples = [s for s in all_samples() if not s.is_signal or s.production_mode in modes]
    if cfg.run_samples not in ("all", ""):
        aliases = {"signal": [s.name for s in samples if s.is_signal],
                   "bkg_qq": [s.name for s in samples if s.name.startswith("bkg_qq_")],
                   "bkg_gg": [s.name for s in samples if s.name.startswith("bkg_gg_")]}
        names = set()
        for token in cfg.run_samples.split(","):
            names.update(aliases.get(token, [token]))
        unknown = names - {s.name for s in samples}
        if unknown:
            raise ValueError(f"unknown/inactive HO samples: {sorted(unknown)}")
        samples = [s for s in samples if s.name in names]
    return samples


def load_signal_record(cfg):
    if cfg.signal_normalization is None:
        raise ValueError("HO ggF requires --signal-normalization; run the prepared ihixs calculation")
    record = norm.load_record(cfg.signal_normalization)
    if record["sqrt_s_gev"] != 2 * cfg.ebeam or record["higgs_mass_gev"] != cfg.higgs_mass:
        raise ValueError("ihixs energy/mass differs from the four-lepton campaign")
    return record


def powheg_higgs_mass_bounds(path, cfg):
    """Read the HJ Breit-Wigner bounds, after checking its pole inputs.

    HJ/Born_phsp.f samples mH between max(0, hmass-bwcutoff*hwidth)
    and min(sqrt(s), hmass+bwcutoff*hwidth). higgsfixedwidth changes
    the width prescription; it does not put each Higgs on its pole mass.
    """
    record = load_signal_record(cfg)
    prefix = four.read_lhe_prefix(path)
    values = {}
    for key in ("hmass", "hwidth", "bwcutoff", "minnlo"):
        matches = re.findall(rf"(?mi)^\s*{key}\s+(\S+)", prefix)
        if len(matches) != 1:
            raise ValueError(f"finite-width Higgs validation requires one explicit POWHEG {key}")
        value = ho.number(matches[0])
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"invalid POWHEG {key} for finite-width Higgs validation")
        values[key] = value
    ho.validate_signal_hard_inputs(prefix, SimpleNamespace(
        expected_higgs_mass=cfg.higgs_mass,
        expected_native_parameters=record["native_parameters"]), non_nominal_higgs_mass=True)
    offset = values["bwcutoff"] * values["hwidth"]
    if not math.isfinite(offset):
        raise ValueError("non-finite POWHEG Higgs mass window")
    return (max(0., values["hmass"] - offset),
            min(2 * cfg.ebeam, values["hmass"] + offset), values["hwidth"])


def inspect_production(path, sample, cfg):
    args = SimpleNamespace(ebeam=cfg.ebeam, expected_higgs_mass=cfg.higgs_mass,
                           nevents=cfg.events_for(sample))
    if sample.production_mode == "ggf":
        args.expected_native_parameters = load_signal_record(cfg)["native_parameters"]
    if sample.matching != "LO":
        proxy = SimpleNamespace(signal=sample.production_mode == "ggf", lhaid=sample.lhaid)
        info = ho.inspect_lhe(path, proxy, args)
    else:
        header = four.validate_lhe_collider(path, cfg)
        if header["pdf_set_ids"] != [sample.lhaid, sample.lhaid]:
            raise ValueError("loop-induced LHE has incompatible PDF IDs")
        info = {"cross_section_pb": four.parse_lhe_init_cross_section_pb(path)}
    if sample.is_signal:
        four.validate_stable_higgs_lhe(path, cfg, cfg.events_for(sample), sample=sample)
    return info


def normalization(sample, native_lhe, cfg, stable_lhe=None):
    native = four.parse_lhe_init_cross_section_pb(native_lhe)
    if not sample.is_signal:
        return dict(native_cross_section_pb=native, target_cross_section_pb=native,
                    weight_scale=1., formula="target=native; weight_scale=1", normalization_kind="native_lhe")
    if stable_lhe is None:
        raise ValueError("stable production LHE is required before four-body decay")
    info = inspect_production(stable_lhe, sample, cfg)
    production = info["cross_section_pb"]
    record = load_signal_record(cfg) if sample.production_mode == "ggf" else None
    rate = record["cross_section_pb"] if record else production
    br = four.signal_target_branching_fraction(sample, cfg)
    result = dict(native_cross_section_pb=native, stable_production_cross_section_pb=production,
        production_cross_section_pb=rate, target_branching_fraction=br,
        target_cross_section_pb=rate * br, weight_scale=rate * br / native,
        normalization_kind="ihixs_n3lo" if record else "native_lhe",
        formula="production_xsec * physical_BR / MadSpin_decayed_xsec")
    if record:
        result.update(ihixs=record, ihixs_record_sha256=record["fingerprint"])
    return result


def nlo_cards(sample, cfg, process_dir, count):
    process = "set automatic_html_opening False\nset output_dependencies internal\n" + four.process_card_text(sample, process_dir)
    settings = dict(nevents=count, iseed=four.seed_for(sample.production_group_id, cfg),
        ebeam1=cfg.ebeam, ebeam2=cfg.ebeam, lpp1=1, lpp2=1, pdlabel="lhapdf", lhaid=sample.lhaid,
        parton_shower="HERWIGPP", event_norm="average", ickkw=0,
        fixed_ren_scale=False, fixed_fac_scale=False, dynamical_scale_choice=3,
        mur_over_ref=1., muf_over_ref=1., reweight_scale=True, reweight_PDF=False,
        jetalgo=-1, jetradius=0.4, maxjetflavor=5, ptj=0., etaj=-1.,
        ptl=2., etal=5., mll=2., mll_sf=2., drll=0., drll_sf=0.)
    launch = (f"set nb_core {cfg.nb_core}\nlaunch {process_dir} aMC@NLO --parton --name={cfg.run_tag}\n" +
              "".join(f"set {k} {v}\n" for k, v in settings.items()) +
              f"set mass 25 {cfg.higgs_mass}\ndone\n")
    return process, launch


def obtain_production(sample, selected, cfg):
    if sample.production_mode == "ggf":
        if cfg.signal_lhe is None:
            raise ValueError("generate HJMiNNLO first and pass --signal-lhe")
        if cfg.dry_run:
            return cfg.signal_lhe
        info = inspect_production(cfg.signal_lhe, sample, cfg)
        path = prepare_powheg_for_madspin(cfg.signal_lhe, info, sample, cfg)
        return independent_decay_input(path, sample, selected, cfg)
    if cfg.source_backend == "external-lhe":
        path = four.resolve_external_lhe(sample, cfg).resolve()
    elif sample.matching == "LO":
        path = four.obtain_mg5_lhe(sample, selected, cfg)
    else:
        group = four.production_dir(sample, cfg)
        process_dir = group / "mg5_process"
        cards = nlo_cards(sample, cfg, process_dir, four.group_event_count(sample.production_group_id, selected, cfg))
        marker = group / "fourlepton_production_stage.json"
        identity = {"process_card": cards[0], "launch_card": cards[1], "adapter_sha256": four.sha256_file(Path(__file__))}
        path = process_dir / "Events" / cfg.run_tag / "events.lhe.gz"
        if not cfg.dry_run and path.exists():
            if not four.artifact_stage_matches(marker, identity, path):
                raise ValueError(f"existing production differs from cards: {path}; use a new tag")
        else:
            if not cfg.dry_run and process_dir.exists():
                raise ValueError(f"incomplete NLO export: {process_dir}; inspect it before using a new tag")
            four.ensure_dir(group / "cards", cfg)
            for name, content in zip(("process_card.dat", "launch.mg5"), cards):
                four.write_text(group / "cards" / name, content, cfg, name)
            four.run_runtime_command(four.mg5_command(cfg, group / "cards/process_card.dat"), cfg)
            if not cfg.dry_run:
                ho.patch_mg5_fortran(process_dir)
            four.run_runtime_command(four.mg5_command(cfg, group / "cards/launch.mg5"), cfg)
            if not cfg.dry_run:
                if not path.is_file():
                    raise ValueError(f"MG5 did not produce {path}")
                four.write_artifact_stage_marker(marker, identity, path)
    if not cfg.dry_run:
        inspect_production(path, sample, cfg)
    return independent_decay_input(path, sample, selected, cfg) if sample.is_signal else path


def slice_bounds(sample, selected, cfg):
    """Disjoint half-open source ranges, in the declared decay-channel order."""
    start = 0
    for other in selected:
        if other.production_group_id != sample.production_group_id:
            continue
        stop = start + cfg.events_for(other)
        if other.name == sample.name:
            return start, stop
        start = stop
    raise ValueError("sample missing from its production group")


def independent_decay_input(source, sample, selected, cfg):
    """Split the common integration's hard events into independent decay inputs.

    Cross sections describe the same inclusive process in each slice. No rate
    factor is applied for the slice size. Event records and optional weights
    are copied verbatim, with source ranges and checksums saved for auditing.
    """
    destination = four.production_dir(sample, cfg) / "decay_inputs" / f"{sample.channel}.lhe.gz"
    if cfg.dry_run:
        return destination
    start, stop = slice_bounds(sample, selected, cfg)
    identity = {"kind": "disjoint_hard_event_slice", "source": str(source.resolve()),
                "source_sha256": four.sha256_file(source), "start": start, "stop": stop,
                "sample": sample.name, "adapter_sha256": four.sha256_file(Path(__file__))}
    marker = destination.with_suffix(".slice.json")
    if destination.exists():
        if not four.artifact_stage_matches(marker, identity, destination):
            raise ValueError("decay input/source range changed; use a new tag")
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".tmp")
    count = 0
    in_event = False
    with four.open_maybe_gzip(source) as stream, gzip.open(temporary, "wt") as output:
        for line in stream:
            if re.search(r"<event\b", line):
                if in_event:
                    raise ValueError("nested LHE event in hard-event slicing")
                in_event = True
            if not in_event or start <= count < stop:
                output.write(line)
            if "</event>" in line:
                in_event = False
                count += 1
    if in_event or count < stop:
        temporary.unlink()
        raise ValueError(f"{source} needs at least {stop} hard events for independent {sample.channel} input; found {count}")
    temporary.replace(destination)
    four.write_artifact_stage_marker(marker, identity, destination)
    return destination


def prepare_powheg_for_madspin(source, info, sample, cfg):
    """Give MadSpin a physical init integral when POWHEG writes XSECUP=-1.

    Event records, including signs and SCALUP, are copied verbatim. The native
    signed mean is diagnostic; the physical ggF rate still comes from ihixs.
    """
    group = four.production_dir(sample, cfg)
    group.mkdir(parents=True, exist_ok=True)
    destination = group / "stable_higgs_for_madspin.lhe.gz"
    marker = group / "powheg_init_preparation.json"
    identity = {"source": str(source), "source_sha256": four.sha256_file(source),
                "native_cross_section_pb": info["cross_section_pb"],
                "native_error_pb": info["cross_section_error_pb"],
                "adapter_sha256": four.sha256_file(Path(__file__))}
    if destination.exists():
        if not four.artifact_stage_matches(marker, identity, destination):
            raise ValueError("prepared HJMiNNLO input differs; use a new tag")
        return destination
    in_init = False
    row = 0
    temporary = destination.with_suffix(".tmp")
    with four.open_maybe_gzip(source) as src, gzip.open(temporary, "wt") as out:
        for line in src:
            if line.strip() == "<init>":
                in_init = True
            elif line.strip() == "</init>":
                in_init = False
            elif in_init and line.strip() and not line.lstrip().startswith(("#", "<")):
                fields = line.split()
                if row == 0:
                    if int(fields[9]) != 1:
                        raise ValueError("HJMiNNLO preparation requires one subprocess")
                else:
                    fields[0] = f"{info['cross_section_pb']:.16e}"
                    fields[1] = f"{info['cross_section_error_pb']:.16e}"
                    line = " ".join(fields) + "\n"
                row += 1
            out.write(line)
    temporary.replace(destination)
    four.write_artifact_stage_marker(marker, identity, destination)
    return destination


def herwig_card(sample, lhe, seed, count, cfg):
    text = Template((ho.HO_DIR / "HW-LHE.in").read_text()).substitute(
        lhe_file=lhe, pdf_name=sample.pdf_name, nevents=count, seed=seed,
        sample_name=sample.name, max_errors=ho.default_herwig_error_budget(count),
        matching_settings=(ho.HO_DIR / "MCatNLO.in").read_text() if sample.matching == "MC@NLO" else "",
        decay_settings=f"set /Herwig/Particles/h0:NominalMass {cfg.higgs_mass}*GeV")
    # The h -> 4l decay is already in the LHE. Associated W/Z/top decays are inclusive.
    text = text.replace("SaveObjects No", "SaveObjects Yes")
    text = text.replace("set /Herwig/Analysis/HwSim:SavePartons Yes", "set /Herwig/Analysis/HwSim:SavePartons Yes\nset /Herwig/Analysis/HwSim:SaveOptionalWeights Yes")
    for species in ("Electron", "Muon", "Photon"):
        text = text.replace(f"PTCut{species} 10.0", f"PTCut{species} 0.4")
    return text


def write_herwig_input(sample, lhe, seed, count, output, cfg):
    text = herwig_card(sample, lhe, seed, count, cfg)
    four.write_text(output, text, cfg, "matched four-lepton Herwig card")
    return {"prompt_selection": None, "decay_bias": {"enabled": False},
            "branching_ratio_reweighter": False, "handlers": []}


def logged_command(argv, cfg, cwd, log):
    if cfg.dry_run:
        four.run_runtime_command(argv, cfg, cwd=cwd)
        return 0
    command = [str(x) for x in argv]
    if cfg.herwig_env or cfg.herwig_module:
        command = ["bash", "-lc", "; ".join([*four.runtime_setup_commands(cfg), four.format_command(argv)])]
    with log.open("w") as out:
        result = subprocess.run(command, cwd=cwd, stdout=out, stderr=subprocess.STDOUT,
                                env=four.clean_shell_env())
    return result.returncode


def run_shower(sample, lhe, cfg):
    work = four.sample_work_dir(sample, cfg)
    directory = work / "herwig"
    four.ensure_dir(directory / "events", cfg)
    card_path = directory / f"{sample.name}.in"
    card = herwig_card(sample, lhe, four.seed_for(f"{sample.name}:herwig", cfg), cfg.events_for(sample), cfg)
    marker = directory / "matched_completion.json"
    if not cfg.dry_run and card_path.exists() and card_path.read_text() != card:
        raise ValueError(f"saved shower card changed: {card_path}; use a new tag")
    metadata = write_herwig_input(sample, lhe, four.seed_for(f"{sample.name}:herwig", cfg), cfg.events_for(sample), card_path, cfg)
    if not cfg.dry_run and marker.exists():
        saved = json.loads(marker.read_text())
        completion.validate_completion(saved, {"sample": sample.name,
            "nevents_requested": cfg.events_for(sample), "lhe": {"events": saved["lhe_events"]}})
        # The completion validator verifies file provenance when a directory is supplied.
        completion.validate_completion(saved, {"sample": sample.name,
            "nevents_requested": cfg.events_for(sample), "lhe": {"events": saved["lhe_events"]}}, directory=work)
    else:
        if not cfg.dry_run and list((directory / "events").glob("*.root")):
            raise ValueError("partial shower output exists; inspect logs and use a new run tag")
        if logged_command([cfg.herwig, "read", card_path], cfg, directory, directory / "read.log"):
            raise ValueError(f"Herwig read failed: {directory / 'read.log'}")
        code = logged_command([cfg.herwig, "run", f"{sample.name}.run", f"-N{cfg.events_for(sample)}"], cfg, directory, directory / "run.log")
        if cfg.dry_run:
            return work / f"{sample.name}_hwsim_roots.input", metadata
        termination = "source_exhausted" if code else "requested_events"
        roots = sorted((directory / "events").glob("**/*.root"))
        if logged_command([cfg.mg5_python or "python3", "-c", ho.ROOT_METADATA_SCRIPT, *roots], cfg, directory, directory / "root-metadata.log"):
            raise ValueError("cannot audit shower ROOT output")
        root_info = json.loads((directory / "root-metadata.log").read_text().splitlines()[-1])
        # source_weight_summary checks complete event syntax and scans the full LHE maximum.
        with four.open_maybe_gzip(lhe) as stream:
            lhe_count = sum(line.strip() == "</event>" for line in stream)
        saved = completion.build_completion(directory, sample.name, cfg.events_for(sample), lhe_count, root_info, termination)
        saved["source_weights"] = completion.source_weight_summary(lhe, saved["attempted_events"])
        saved["fingerprint"] = completion.completion_fingerprint(saved)
        completion.validate_completion(saved, {"sample": sample.name,
            "nevents_requested": cfg.events_for(sample), "lhe": {"events": saved["lhe_events"]}})
        marker.write_text(json.dumps(saved, indent=2) + "\n")
    root_input = work / f"{sample.name}_hwsim_roots.input"
    four.write_root_input_list(directory / "events", root_input, cfg)
    metadata["shower_completion"] = saved
    return root_input, metadata
