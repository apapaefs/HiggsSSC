#!/usr/bin/env python3
"""Prepare/run POWHEG HJMiNNLO + Herwig and MG5 MC@NLO diphoton backgrounds.

The default stage only writes cards. Signal LHE production is handled by
run_powheg_hjminnlo.py; pass its merged output with --signal-lhe.
"""

from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from string import Template

try:
    from . import run_gammagamma_campaign as lo
except ImportError:
    import run_gammagamma_campaign as lo

SCRIPT_DIR = Path(__file__).resolve().parent
HO_DIR = SCRIPT_DIR / "HOAnalysis"
ANALYSIS_DIR = SCRIPT_DIR / "LOAnalysis" / "Code"
ANALYSIS_EXE = ANALYSIS_DIR / "HwSimPostAnalysis_gammagamma_SSC"


@dataclass(frozen=True)
class Sample:
    name: str
    category: str
    process: str
    response_mode: str
    matching: str

    @property
    def signal(self) -> bool:
        return self.category == "Signal"

    @property
    def pdf_name(self) -> str:
        return lo.NNLO_PDF_NAME if self.signal else lo.NLO_PDF_NAME

    @property
    def lhaid(self) -> int:
        return lo.NNLO_PDF_LHAID if self.signal else lo.NLO_PDF_LHAID


SAMPLES = (
    Sample("signal_gg_h_aa", "Signal", "HJMiNNLO; h -> gamma gamma", "genuine", "POWHEG"),
    Sample("bkg_prompt_aa", "Backgrounds", "p p > a a [QCD]", "genuine", "MC@NLO"),
    Sample("bkg_gamma_j", "Backgrounds", "p p > a j [QCD]", "gammajet", "MC@NLO"),
    Sample("bkg_dy_ee", "Backgrounds", "p p > e+ e- [QCD]", "dielectron", "MC@NLO"),
)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "build", "generate", "shower", "analyze", "all"), default="prepare")
    parser.add_argument("--run-tag", default="ho_run_01")
    parser.add_argument("--output-dir", type=Path, help="default: HOAnalysis/runs/RUN_TAG")
    parser.add_argument("--run-samples", default="all", help="all, backgrounds, or comma-separated names")
    parser.add_argument("--signal-lhe", type=Path, help="complete merged HJMiNNLO LHE (.gz accepted), with undecayed Higgs")
    parser.add_argument("--nevents", type=int, default=10000, help="events per background and maximum signal shower events")
    parser.add_argument("--ebeam", type=float, default=20000.)
    parser.add_argument("--seed-base", type=int, default=730001)
    parser.add_argument("--nb-core", type=int, default=1)
    parser.add_argument("--higgs-br", type=float, default=lo.YR4_BR_H_TO_GAMMAGAMMA)
    parser.add_argument("--mg5-dir", type=Path, default=lo.REPO_ROOT / "MG5_aMC_v3_5_15")
    parser.add_argument("--mg5-python", default=os.environ.get("MG5_PYTHON", "python3"))
    parser.add_argument("--mg5-fortran", help="override MG5 Fortran compiler (e.g. gfortran)")
    parser.add_argument("--mg5-cxx", help="override MG5 C++ compiler, compatible with LHAPDF")
    parser.add_argument("--mg5-dependencies", choices=("internal", "external"), default="internal",
                        help="internal keeps rebuildable dependency copies inside each exported process")
    parser.add_argument("--analysis-cxx", default="/usr/bin/clang++" if sys.platform == "darwin" else "c++",
                        help="C++ compiler compatible with the installed ROOT")
    parser.add_argument("--herwig", default=os.environ.get("HERWIG", "Herwig"))
    parser.add_argument("--herwig-env", type=Path, default=lo.maybe_path(os.environ.get("HERWIG_ENV")))
    parser.add_argument("--herwig-module", default=os.environ.get("HERWIG_MODULE"))
    parser.add_argument("--no-herwig-module", action="store_true")
    parser.add_argument("--gen-photon-pt-min", type=float, default=10.)
    parser.add_argument("--gen-jet-pt-min", type=float, default=10.)
    parser.add_argument("--gen-lepton-pt-min", type=float, default=10.)
    parser.add_argument("--gen-eta-max", type=float, default=6.)
    parser.add_argument("--gen-mll-min", type=float, default=30.)
    parser.add_argument("--isolation-radius", type=float, default=0.4)
    parser.add_argument("--isolation-epsilon", type=float, default=1.)
    parser.add_argument("--isolation-power", type=float, default=1.)
    parser.add_argument("--resume", action="store_true", help="reuse completed stages with identical configuration")
    parser.add_argument("--dry-run", action="store_true", help="print cards and commands without writing or running")
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.run_tag):
        parser.error("run-tag must contain only letters, digits, underscores and hyphens")
    for key in ("nevents", "nb_core", "seed_base", "ebeam", "gen_photon_pt_min", "gen_jet_pt_min",
                "gen_lepton_pt_min", "gen_eta_max", "gen_mll_min", "isolation_radius",
                "isolation_epsilon", "isolation_power"):
        if not math.isfinite(getattr(args, key)) or getattr(args, key) <= 0:
            parser.error(f"{key.replace('_', '-')} must be finite and positive")
    if not 0 < args.higgs_br <= 1:
        parser.error("higgs-br must be in (0, 1]")
    args.output_dir = (args.output_dir or HO_DIR / "runs" / args.run_tag).expanduser().resolve()
    args.mg5_dir = args.mg5_dir.expanduser().resolve()
    if args.signal_lhe:
        args.signal_lhe = args.signal_lhe.expanduser().resolve()
    if args.herwig_env:
        args.herwig_env = args.herwig_env.expanduser().resolve()
        if args.herwig_env.is_dir():
            args.herwig_env /= "bin/activate"
    if args.no_herwig_module:
        args.herwig_module = None
    elif not args.herwig_env and not args.herwig_module:
        args.herwig_module = (lo.DEFAULT_DARWIN_HERWIG_MODULE if sys.platform == "darwin"
                              else lo.DEFAULT_LINUX_HERWIG_MODULE)
    if args.herwig_module:
        args.herwig_env = None
    # MG5/ThePEG card parsers do not share shell quoting rules.
    for path in (args.output_dir, args.mg5_dir, args.signal_lhe):
        if path and (re.search(r"\s|[;#\"']", str(path))):
            parser.error(f"generator paths must not contain whitespace or card metacharacters: {path}")
    return args


def selected_samples(args):
    if args.run_samples == "all":
        return list(SAMPLES)
    if args.run_samples == "backgrounds":
        return [s for s in SAMPLES if not s.signal]
    names = set(args.run_samples.split(","))
    unknown = names - {s.name for s in SAMPLES}
    if unknown:
        raise ValueError(f"unknown samples: {', '.join(sorted(unknown))}")
    return [s for s in SAMPLES if s.name in names]


def sample_dir(args, sample):
    return args.output_dir / sample.category / "events" / sample.name


def lhe_path(args, sample):
    if sample.signal:
        return args.signal_lhe or args.output_dir / "powheg" / "powheg-hjminnlo-merged.lhe"
    return sample_dir(args, sample) / "mg5_process" / "Events" / args.run_tag / "events.lhe.gz"


def run_settings(args, sample):
    """NLO run-card keys, deliberately distinct from MadEvent LO cuts."""
    return {
        "nevents": args.nevents, "iseed": args.seed_base + SAMPLES.index(sample),
        "ebeam1": args.ebeam, "ebeam2": args.ebeam, "lpp1": 1, "lpp2": 1,
        "pdlabel": "lhapdf", "lhaid": lo.NLO_PDF_LHAID,
        "parton_shower": "HERWIGPP", "event_norm": "average", "ickkw": 0,
        "fixed_ren_scale": False, "fixed_fac_scale": False,
        "dynamical_scale_choice": 3, "mur_over_ref": 1., "muf_over_ref": 1.,
        "reweight_scale": True, "reweight_PDF": False,
        "jetalgo": -1, "jetradius": 0.4, "maxjetflavor": 5,
        # No generation cut on the extra real-emission jet in colour-singlet samples.
        "ptj": args.gen_jet_pt_min if sample.response_mode == "gammajet" else 0.,
        "etaj": args.gen_eta_max,
        "ptl": args.gen_lepton_pt_min if sample.response_mode == "dielectron" else 0.,
        "etal": args.gen_eta_max, "mll": args.gen_mll_min, "mll_sf": args.gen_mll_min,
        "drll": 0., "drll_sf": 0.,
        "gamma_is_j": False, "ptgmin": args.gen_photon_pt_min,
        "etagamma": args.gen_eta_max, "R0gamma": args.isolation_radius,
        "epsgamma": args.isolation_epsilon, "xn": args.isolation_power, "isoEM": True,
    }


def mg5_cards(args, sample):
    process_dir = sample_dir(args, sample) / "mg5_process"
    setup = "set automatic_html_opening False\n"
    for option, value in (("fortran_compiler", args.mg5_fortran), ("cpp_compiler", args.mg5_cxx)):
        if value:
            setup += f"set {option} {value} --no_save\n"
    # A massless five-flavour model and p/j definitions consistent with the NLO PDF.
    process = (setup + f"set output_dependencies {args.mg5_dependencies}\n" +
               "import model loop_sm-no_b_mass\n"
               "define p = g u c d s b u~ c~ d~ s~ b~\n"
               "define j = g u c d s b u~ c~ d~ s~ b~\n"
               f"generate {sample.process}\noutput {process_dir}\n")
    launch = (setup +
              f"set nb_core {args.nb_core}\n"
              f"launch {process_dir} aMC@NLO --parton --name={args.run_tag}\n")
    launch += "".join(f"set {key} {value}\n" for key, value in run_settings(args, sample).items())
    launch += "done\n"
    return process, launch


def herwig_card(args, sample):
    decay = ""
    if sample.signal:
        decay = ("# Unit-probability forced decay; physical BR is applied only by post-analysis.\n"
                 "set /Herwig/Particles/h0:NominalMass 125.0*GeV\n"
                 "do /Herwig/Particles/h0:SelectDecayModes h0->gamma,gamma;\n"
                 "set /Herwig/Particles/h0/h0->gamma,gamma;:BranchingRatio 1.0\n")
    matching = ((HO_DIR / "MCatNLO.in").read_text() if not sample.signal
                else "# POWHEG SCALUP veto; retain the default POWHEG recoil settings.")
    return Template((HO_DIR / "HW-LHE.in").read_text()).substitute(
        lhe_file=lhe_path(args, sample), pdf_name=sample.pdf_name,
        matching_settings=matching, decay_settings=decay, nevents=args.nevents,
        seed=args.seed_base + SAMPLES.index(sample), sample_name=sample.name,
    )


def write_text(path, value, args):
    if args.dry_run:
        print(f"+ write {path}\n{value}", end="" if value.endswith("\n") else "\n")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(value)
        temporary.replace(path)


def write_json(path, value, args):
    write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n", args)


def command(args, argv, cwd, log_name):
    print(f"+ [{cwd}] {shlex.join(str(x) for x in argv)}", flush=True)
    if args.dry_run:
        return
    setup = lo.runtime_setup_commands(args)
    # Activation scripts may export FC/CXX and thereby override MG5's
    # make_opts compiler selection. Honour explicit overrides after activation.
    for key, value in (("FC", args.mg5_fortran), ("F77", args.mg5_fortran), ("CXX", args.mg5_cxx)):
        if value:
            setup.append(f"export {key}={shlex.quote(value)}")
    shell = "\n".join([*setup, shlex.join(str(x) for x in argv)])
    log = cwd / log_name
    with log.open("w") as stream:
        result = subprocess.run(["bash", "-c", shell], cwd=cwd, env=lo.clean_shell_env(),
                                stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode or re.search(r"interrupted with error:|Error detected in", log.read_text(errors="replace")):
        raise RuntimeError(f"generator/build error (exit status {result.returncode}); see {log}")


def prepare(args, sample):
    directory = sample_dir(args, sample)
    process, launch = ("", "") if sample.signal else mg5_cards(args, sample)
    herwig = herwig_card(args, sample)
    identity = {"process_card": process, "launch_card": launch, "herwig_card": herwig,
                "higgs_br": args.higgs_br if sample.signal else 1.0,
                "response_mode": sample.response_mode,
                "herwig": args.herwig, "herwig_env": str(args.herwig_env),
                "herwig_module": args.herwig_module, "mg5_dir": str(args.mg5_dir),
                "mg5_python": args.mg5_python}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    manifest_path = directory / "campaign.json"
    manifest = {
        "schema_version": 1, "run_tag": args.run_tag, "sample": sample.name,
        "process": sample.process, "matching": sample.matching,
        "hard_accuracy": "NNLO+PS (HJMiNNLO)" if sample.signal else "NLO QCD + PS",
        "response_mode": sample.response_mode, "detector_response": "ssc",
        "pdf": sample.pdf_name, "lhaid": sample.lhaid, "ebeam_gev": args.ebeam,
        "nevents_requested": args.nevents, "seed": args.seed_base + SAMPLES.index(sample),
        "weight_scale": args.higgs_br if sample.signal else 1.0,
        "normalization_rule": "inclusive LHE cross section times physical BR once" if sample.signal
                              else "signed MC@NLO cross section; no K factor",
        "lhe_file": str(lhe_path(args, sample)), "fingerprint": fingerprint,
        "configuration": identity, "completed": [],
    }
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if old["fingerprint"] != fingerprint or old["ebeam_gev"] != args.ebeam:
            raise ValueError(f"configuration changed in {directory}; use a new --run-tag/output-dir")
        manifest = old
    elif directory.exists() and any(directory.iterdir()):
        raise ValueError(f"unmanaged non-empty sample directory: {directory}")
    if not sample.signal:
        write_text(directory / "cards/process.mg5", process, args)
        write_text(directory / "cards/launch.mg5", launch, args)
    write_text(directory / "herwig" / f"{sample.name}.in", herwig, args)
    write_json(manifest_path, manifest, args)
    return manifest


def number(value):
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError(f"non-finite LHE value: {value}")
    return result


def inspect_lhe(path, sample, args):
    """Stream full LHE: signed cross section, weights, beams/PDFs, and decay state."""
    opener = gzip.open if path.suffix == ".gz" else open
    header = []
    init = []
    in_init = in_event = closed = False
    particle_lines = 0
    higgs = stable_higgs = count = negatives = 0
    sumw = sumabs = sumw2 = 0.
    with opener(path, "rt") as stream:
        for line in stream:
            stripped = line.strip()
            if count == 0 and not in_init and not in_event:
                header.append(line)
            if stripped == "<init>":
                in_init = True
            elif stripped == "</init>":
                in_init = False
            elif in_init and stripped and not stripped.startswith(("#", "<")):
                init.append(stripped.split())
            elif re.match(r"<event(?:\s[^>]*)?>$", stripped):
                if in_event:
                    raise ValueError(f"nested/incomplete LHE event in {path}")
                in_event = True
                particle_lines = -1
                higgs = stable_higgs = 0
            elif stripped == "</event>":
                if not in_event or particle_lines != 0:
                    raise ValueError(f"incomplete LHE particle record in {path}")
                if sample.signal and (higgs != 1 or stable_higgs != 1):
                    raise ValueError("signal LHE must contain exactly one undecayed, status-1 Higgs per event")
                in_event = False
                count += 1
            elif in_event and stripped and not stripped.startswith(("#", "<")):
                fields = stripped.split()
                if particle_lines == -1:
                    particle_lines = int(fields[0])
                    weight = number(fields[2])
                    sumw += weight
                    sumabs += abs(weight)
                    sumw2 += weight * weight
                    negatives += weight < 0
                elif particle_lines > 0:
                    if len(fields) < 13:
                        raise ValueError(f"short LHE particle record in {path}")
                    higgs += int(fields[0]) == 25
                    stable_higgs += int(fields[0]) == 25 and int(fields[1]) == 1
                    particle_lines -= 1
            elif stripped == "</LesHouchesEvents>":
                closed = True
    if not closed or in_event or not count or len(init) < 2:
        raise ValueError(f"empty or incomplete LHE file: {path}")
    beam = init[0]
    if [int(x) for x in beam[:2]] != [2212, 2212]:
        raise ValueError("expected two proton beams")
    if any(not math.isclose(number(x), args.ebeam) for x in beam[2:4]):
        raise ValueError("LHE beam energy differs from --ebeam")
    if len(init) != 1 + int(beam[9]):
        raise ValueError("LHE init subprocess count is inconsistent")
    header_text = "".join(header)
    pdf_ids = [int(x) for x in beam[6:8]]
    if sample.signal and pdf_ids == [-1, -1]:
        # HJMiNNLO uses IDWTUP=-4 and leaves PDF and cross-section init fields
        # unspecified. Its complete POWHEG input is embedded in the header.
        matches = [re.search(rf"(?m)^\s*lhans{i}\s+(\d+)", header_text) for i in (1, 2)]
        if not all(matches):
            raise ValueError("unspecified LHE PDF IDs without embedded POWHEG lhans1/lhans2")
        pdf_ids = [int(match[1]) for match in matches]
    if pdf_ids != [sample.lhaid, sample.lhaid]:
        raise ValueError(f"LHE PDF IDs differ from required {sample.lhaid}")
    if not sample.signal and not re.search(r"HERWIGPP\s*['\"]?\s*=\s*parton_shower", header_text, re.I):
        raise ValueError("MC@NLO LHE lacks parton_shower=HERWIGPP provenance")
    if sample.signal and not re.search(r"\bminnlo\s+1\b", header_text, re.I):
        raise ValueError("signal LHE lacks HJMiNNLO provenance (minnlo 1)")
    if count < args.nevents:
        raise ValueError(f"only {count} LHE events, but --nevents={args.nevents}")
    xsec = sum(number(row[0]) for row in init[1:])
    error = math.sqrt(sum(number(row[1]) ** 2 for row in init[1:]))
    normalization_source = "LHE init"
    if xsec <= 0 and abs(int(beam[8])) == 4:
        # LHA IDWTUP=+/-4: average signed XWGTUP is the cross section in pb.
        xsec = sumw / count
        error = math.sqrt(max(0., sumw2 / count - xsec ** 2) / (count - 1)) if count > 1 else 0.
        normalization_source = "mean signed XWGTUP (IDWTUP=+/-4; unspecified init cross section)"
    if xsec <= 0 or sumw <= 0:
        raise ValueError("non-positive signed cross section/weight sum; inspect integration statistics")
    return {"cross_section_pb": xsec,
            "cross_section_error_pb": error, "normalization_source": normalization_source,
            "events": count, "negative_events": negatives, "negative_fraction": negatives / count,
            "sum_weight": sumw, "sum_abs_weight": sumabs, "sum_weight_squared": sumw2,
            "effective_events": sumw * sumw / sumw2 if sumw2 else 0,
            "idwtup": int(beam[8]), "source": str(path), "source_size": path.stat().st_size,
            "source_mtime_ns": path.stat().st_mtime_ns}


def save_manifest(args, sample, manifest):
    write_json(sample_dir(args, sample) / "campaign.json", manifest, args)


def run_mg5(args, directory, card_name, log_name):
    # MG5 3.5.15 can turn a bare fastjet-config from its saved configuration
    # into a nonexistent path under MG5_DIR. Resolve both tools in exactly
    # the activated runtime used for generation, and pass absolute paths.
    card = directory / "cards" / card_name
    runtime_card = card.with_name("runtime-" + card.name)
    if args.dry_run:
        print("+ resolve fastjet-config and lhapdf-config in the Herwig runtime")
    else:
        command(args, ["bash", "-c", "command -v fastjet-config && command -v lhapdf-config"],
                directory, "mg5-runtime-paths.log")
        paths = (directory / "mg5-runtime-paths.log").read_text().splitlines()
        if len(paths) != 2 or any(not Path(p).is_file() for p in paths):
            raise ValueError("could not resolve FastJet/LHAPDF executables; see mg5-runtime-paths.log")
        prefix = "".join(f"set {key} {Path(path).resolve()} --no_save\n"
                         for key, path in zip(("fastjet", "lhapdf"), paths))
        # The launch interface reloads the exported configuration and does
        # not forward every option from the parent MG5 session (FastJet in
        # particular). Repair only this campaign's process-local settings.
        configure_mg5_runtime(directory / "mg5_process", paths)
        write_text(runtime_card, prefix + card.read_text(), args)
    command(args, [args.mg5_python, args.mg5_dir / "bin/mg5_aMC", runtime_card], directory, log_name)


def configure_mg5_runtime(process_dir, paths):
    config = process_dir / "Cards/amcatnlo_configuration.txt"
    if not config.exists():
        return  # A new export inherits the runtime card's settings.
    if not config.resolve().is_relative_to(process_dir.resolve()):
        raise ValueError(f"refusing to change a shared MG5 configuration outside {process_dir}")
    text = config.read_text()
    for key, path in zip(("fastjet", "lhapdf"), paths):
        text = re.sub(rf"(?m)^\s*{key}\s*=.*\n?", "", text)
        text += f"\n{key} = {Path(path).resolve()}\n"
    config.write_text(text)


def build(args, sample, manifest):
    directory = sample_dir(args, sample)
    if not sample.signal:
        process_dir = directory / "mg5_process"
        if not process_dir.exists():
            run_mg5(args, directory, "process.mg5", "mg5-process.log")
            if not args.dry_run and not (process_dir / "bin/generate_events").exists():
                raise RuntimeError(f"MG5 process export failed; see {directory / 'mg5-process.log'}")
        elif not args.resume and "build" not in manifest["completed"]:
            raise ValueError(f"process exists: {process_dir}; use --resume to reuse it")
        if not args.dry_run:
            patch_mg5_fortran(process_dir)


def patch_mg5_fortran(process_dir):
    """Declare the MG5 3.5.15 procedure argument required by gfortran 16.

    Modify only the exported copy. This does not change the phase-space
    mapping or touch the MG5 installation used by other campaigns.
    """
    path = process_dir / "SubProcesses/genps_fks.f"
    if not path.exists():
        raise ValueError(f"incomplete MG5 process export: missing {path}")
    if not path.resolve().is_relative_to(process_dir.resolve()):
        raise ValueError(f"refusing to patch a shared MG5 source outside {process_dir}")
    original = path.read_text()
    anchor = "      external derivative,ran2,virtgranny,xinv_redvirtgranny\n"
    updated = anchor + "      external virtgranny_red\n"
    if updated in original or "derivative(virtgranny_red" not in original:
        return
    if anchor not in original:
        # A different MG5 version may already declare it on a continued line.
        return
    path.write_text(original.replace(anchor, updated, 1))


def generate(args, sample, manifest):
    directory = sample_dir(args, sample)
    if sample.signal:
        if not args.signal_lhe:
            raise ValueError("generate signal LHE with run_powheg_hjminnlo.py, then pass --signal-lhe")
    else:
        build(args, sample, manifest)
        # A complete existing LHE is recoverable even if the wrapper was interrupted.
        if not lhe_path(args, sample).exists():
            run_mg5(args, directory, "launch.mg5", "mg5-generate.log")
    if not args.dry_run:
        if not lhe_path(args, sample).exists():
            raise RuntimeError(f"expected LHE file missing; inspect logs in {directory}")
        manifest["lhe"] = inspect_lhe(lhe_path(args, sample), sample, args)


def shower(args, sample, manifest):
    directory = sample_dir(args, sample)
    herwig_dir = directory / "herwig"
    if not args.dry_run:
        metadata = inspect_lhe(lhe_path(args, sample), sample, args)
        if "lhe" in manifest and metadata != manifest["lhe"]:
            raise ValueError("LHE changed since the previous stage; use a new run tag")
        manifest["lhe"] = metadata
        # Never overwrite a recoverable partial shower.
        if list((herwig_dir / "events").glob("**/*.root")):
            raise ValueError(f"existing HwSim events in {herwig_dir}; use analyze or a new run tag")
        (herwig_dir / "events").mkdir(exist_ok=True)
        save_manifest(args, sample, manifest)
    command(args, [args.herwig, "--version"], herwig_dir, "version.log")
    command(args, [args.herwig, "read", f"{sample.name}.in"], herwig_dir, "read.log")
    command(args, [args.herwig, "run", f"{sample.name}.run", f"-N{args.nevents}"], herwig_dir, "run.log")
    if not args.dry_run:
        manifest["root_files"] = root_inventory(directory)


def root_inventory(directory):
    roots = sorted((directory / "herwig/events").glob("**/*.root"))
    if not roots or any(path.stat().st_size == 0 for path in roots):
        raise ValueError(f"missing or empty HwSim ROOT files in {directory}")
    return [{"path": str(path), "size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
            for path in roots]


def verify_completed(args, sample, stage, manifest):
    directory = sample_dir(args, sample)
    if stage == "build" and not sample.signal:
        if not (directory / "mg5_process/bin/generate_events").exists():
            raise ValueError(f"completed process export is missing in {directory}")
    if stage in ("generate", "shower", "analyze"):
        current = inspect_lhe(lhe_path(args, sample), sample, args)
        if current != manifest.get("lhe"):
            raise ValueError(f"LHE changed since completed {stage}: {directory}")
    if stage in ("shower", "analyze"):
        if root_inventory(directory) != manifest.get("root_files"):
            raise ValueError(f"ROOT products changed since completed {stage}: {directory}")
    if stage == "analyze":
        prefix = directory / f"{sample.name}_hwsim_roots-{args.run_tag}"
        for suffix in (".dat", ".top", "_var.root"):
            if not Path(str(prefix) + suffix).is_file():
                raise ValueError(f"completed analysis output is missing: {prefix}{suffix}")
        if not (directory / f"normalization-{args.run_tag}.json").is_file():
            raise ValueError(f"completed analysis normalization is missing in {directory}")


def analyze(args, sample, manifest):
    directory = sample_dir(args, sample)
    roots = sorted((directory / "herwig/events").glob("**/*.root"))
    if not roots and not args.dry_run:
        raise ValueError(f"no HwSim ROOT files in {directory}")
    if not args.dry_run:
        inventory = root_inventory(directory)
        if "root_files" in manifest and inventory != manifest["root_files"]:
            raise ValueError("HwSim ROOT inputs changed since showering")
        manifest["root_files"] = inventory
    root_input = directory / f"{sample.name}_hwsim_roots.input"
    write_text(root_input, "".join(f"{p}\n" for p in roots), args)
    weight_scale = args.higgs_br if sample.signal else 1.
    command(args, [ANALYSIS_EXE, root_input, "-t", args.run_tag, "-w", weight_scale,
                   "--response-mode", sample.response_mode,
                   "--seed", args.seed_base + SAMPLES.index(sample)], directory, "analysis.log")
    if not args.dry_run:
        dat_path = root_input.with_name(root_input.stem + f"-{args.run_tag}.dat")
        data = dict(line.split(maxsplit=1) for line in dat_path.read_text().splitlines() if line.strip())
        if int(data["events_read"]) != args.nevents:
            raise ValueError(f"HwSim event count {data['events_read']} differs from requested {args.nevents}")
        if number(data["sum_weight"]) <= 0:
            raise ValueError("non-positive signed analysis weight sum; increase the pilot statistics")
        manifest["analysis"] = data
        if "lhe" not in manifest:
            manifest["lhe"] = inspect_lhe(lhe_path(args, sample), sample, args)
        write_json(directory / f"normalization-{args.run_tag}.json", {
            "run_tag": args.run_tag, "sample": sample.name,
            "cross_section_pb": manifest["lhe"]["cross_section_pb"],
            "cross_section_error_pb": manifest["lhe"]["cross_section_error_pb"],
            "source": manifest["lhe"]["normalization_source"] + ": before detector response and BR",
            "weight_scale": weight_scale,
        }, args)


def run_campaign(args):
    samples = selected_samples(args)
    if args.stage not in ("prepare", "build") and any(s.signal for s in samples) and not args.signal_lhe:
        raise ValueError("--signal-lhe is required for the signal; use --run-samples backgrounds for MG5 only")
    if args.stage in ("analyze", "all"):
        command(args, ["make", "-C", ANALYSIS_DIR, ANALYSIS_EXE.name,
                       f"CPP={args.analysis_cxx}", f"CXX={args.analysis_cxx}",
                       "FASTJET_CPPFLAGS=", "FASTJET_LDFLAGS="],
                SCRIPT_DIR, "ho-analysis-build.log")
    stages = ("generate", "shower", "analyze") if args.stage == "all" else (args.stage,)
    for sample in samples:
        manifest = prepare(args, sample)
        for stage in stages:
            if stage == "prepare":
                continue
            if stage in manifest["completed"]:
                if args.resume:
                    if not args.dry_run:
                        verify_completed(args, sample, stage, manifest)
                    print(f"Reusing completed {stage}: {sample.name}")
                    continue
                raise ValueError(f"{sample.name}: {stage} already complete; use --resume")
            {"build": build, "generate": generate, "shower": shower, "analyze": analyze}[stage](args, sample, manifest)
            manifest["completed"].append(stage)
            save_manifest(args, sample, manifest)
    print(f"HO campaign directory: {args.output_dir}")
    return 0


def main(argv=None):
    args = parse_args(argv)
    if args.dry_run:
        return run_campaign(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / ".campaign.lock").open("a") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError(f"another campaign command is active in {args.output_dir}")
        return run_campaign(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError) as error:
        sys.exit(f"Error: {error}")
