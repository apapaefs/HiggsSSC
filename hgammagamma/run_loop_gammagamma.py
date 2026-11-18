#!/usr/bin/env python3
"""Generate the separate gg -> gamma gamma continuum at loop-induced LO+PS.

Uses the existing diphoton CLI/runtime and SSC detector response. The hard
matrix element is LO; NLO PDFs match the continuum campaign convention.
No Higgs diagrams, signal branching factors or additional K-factor enter.
"""
from __future__ import annotations

from contextlib import contextmanager
import argparse
import gzip
import json
import math
import os
from pathlib import Path
import re
import shutil
from string import Template
import sys
from types import SimpleNamespace
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hgammagamma import run_gammagamma_ho_campaign as ho
from hgammagamma import ho_shower_completion as completion
from hfourlepton import run_four_lepton_campaign as four

PROCESS = "g g > a a / h [noborn=QCD]"
MODEL = "loop_sm-no_b_mass"
SAMPLE = ho.Sample("bkg_loop_gg_aa", "Backgrounds", PROCESS, "genuine", "LO")
ORIGINAL_INSPECT = ho.inspect_lhe
ORIGINAL_PREPARE = ho.prepare
ORIGINAL_PARSE = ho.parse_args


def parse_args(argv=None):
    extra = argparse.ArgumentParser(add_help=False)
    extra.add_argument("--gen-mgg-min", type=float, default=0., help="minimum Born diphoton mass in GeV")
    extra.add_argument("--gen-mgg-max", type=float, default=-1., help="maximum Born diphoton mass; -1 disables it")
    extra.add_argument("--reuse-loop-export", type=Path, help="copy a completed private loop export into the new campaign")
    values, remaining = extra.parse_known_args(argv)
    if any(item in remaining for item in ("-h", "--help")):
        print("Additional loop-continuum options:\n" + extra.format_help())
    args = ORIGINAL_PARSE(remaining)
    if (not math.isfinite(values.gen_mgg_min) or values.gen_mgg_min < 0
            or not math.isfinite(values.gen_mgg_max)
            or (values.gen_mgg_max != -1 and values.gen_mgg_max <= values.gen_mgg_min)):
        raise ValueError("invalid Born diphoton mass interval")
    args.gen_mgg_min, args.gen_mgg_max = values.gen_mgg_min, values.gen_mgg_max
    args.reuse_loop_export = values.reuse_loop_export.resolve() if values.reuse_loop_export else None
    return args


def lhe_path(args, sample=SAMPLE):
    return ho.sample_dir(args, sample) / "mg5_process/Events" / args.run_tag / "unweighted_events.lhe.gz"


def run_settings(args):
    return dict(nevents=args.nevents, iseed=args.seed_base, ebeam1=args.ebeam, ebeam2=args.ebeam,
        pdlabel="lhapdf", pdlabel1="lhapdf", pdlabel2="lhapdf", lhaid=SAMPLE.lhaid,
        pta=args.gen_photon_pt_min, etaa=args.gen_eta_max, draa=0., draj=0., dral=0.,
        mmaa=getattr(args, "gen_mgg_min", 0.), mmaamax=getattr(args, "gen_mgg_max", -1.),
        ptj=0., etaj=-1., drjj=0., fixed_ren_scale=False, fixed_fac_scale=False,
        dynamical_scale_choice=3, scalefact=1., alpsfact=1., ickkw=0,
        nhel=1, survey_splitting=args.nb_core)


def mg5_cards(args, sample=SAMPLE):
    destination = ho.sample_dir(args, sample) / "mg5_process"
    compilers = "".join(f"set {key} {value} --no_save\n" for key, value in (
        ("fortran_compiler", args.mg5_fortran), ("cpp_compiler", args.mg5_cxx)) if value)
    process = ("set automatic_html_opening False --no_save\n"
               "set output_dependencies internal --no_save\n"
               + compilers +
               f"import model {MODEL}\ngenerate {PROCESS}\noutput {destination}\n")
    # The run card is patched explicitly before generate_events; this JSON is
    # a declarative settings artifact, not an aMC@NLO launch card.
    return process, json.dumps(run_settings(args), indent=2, sort_keys=True) + "\n"


def herwig_card(args, sample=SAMPLE):
    return Template((ho.HO_DIR / "HW-LHE.in").read_text()).substitute(
        lhe_file=lhe_path(args, sample), pdf_name=sample.pdf_name,
        matching_settings="# Loop-induced Born sample with an ordinary LO parton shower.\n",
        decay_settings="", nevents=args.nevents, seed=args.seed_base, sample_name=sample.name,
        max_errors=ho.herwig_error_budget(args, sample))


def prepare(args, sample=SAMPLE):
    manifest = ORIGINAL_PREPARE(args, sample)
    manifest.update(hard_accuracy="loop-induced LO + PS", hard_model=MODEL,
        normalization_rule="native loop-induced LO LHE cross section; no K factor or Higgs BR",
        physics_scope={"continuum_only": True, "higgs_excluded": True,
            "amplitude": "one-loop squared", "cross_section_order": "alpha_s^2 alpha_EM^2",
            "pdf_convention": "NLO PDF and alpha_s with loop-induced LO matrix element",
            "overlap": "separate loop-squared gg Born channel absent from ordinary diphoton NLO",
            "generation_isolation": "Born two-photon final state; no outgoing QCD parton",
            "born_mgg_gev": {"min": getattr(args, "gen_mgg_min", 0.),
                             "max": getattr(args, "gen_mgg_max", -1.)},
            "rare_prompt_photon_backgrounds": False},
        generator_settings=run_settings(args))
    ho.save_manifest(args, sample, manifest)
    return manifest


def patch_run_card(path, settings):
    path = Path(path)
    lines, found = [], set()
    for line in path.read_text().splitlines(keepends=True):
        match = re.match(r"^(\s*)([^=#]+?)(\s*=\s*)(\w+)\b(.*)$", line.rstrip("\r\n"))
        if match and match[4] in settings:
            key = match[4]
            if key in found:
                raise ValueError(f"duplicate loop run-card setting {key}")
            found.add(key)
            line = f"{match[1]}{settings[key]}{match[3]}{key}{match[5]}\n"
        lines.append(line)
    for key in settings:
        if key not in found:
            lines.append(f"\n{settings[key]} = {key}\n")
    path.write_text("".join(lines))


def clone_completed_export(source, destination):
    """Copy completed generated code without sharing any writable process file."""
    source, destination = Path(source).resolve(), Path(destination).resolve()
    manifest_path = source.parent / "campaign.json"
    manifest = json.loads(manifest_path.read_text())
    if (source.name != "mg5_process" or manifest.get("sample") != SAMPLE.name
            or manifest.get("process") != PROCESS or "analyze" not in manifest.get("completed", [])):
        raise ValueError("compiled reuse requires a completed continuum-loop campaign")
    completion.validate_completion(manifest["shower_completion"], manifest, directory=source.parent)
    for name in ("run_card.dat", "param_card.dat"):
        if four.file_record(source / "Cards" / name) != manifest.get("generator_artifacts", {}).get(name):
            raise ValueError("completed export cards changed after its validated generation")
    card = source / "Cards/proc_card_mg5.dat"
    commands = [line.strip() for line in card.read_text().splitlines()
                if line.strip().startswith(("generate ", "add process ", "import model "))]
    if commands != [f"import model {MODEL}", f"generate {PROCESS}"]:
        raise ValueError("compiled export process/model differs from the requested continuum")
    if destination.exists() or destination.is_relative_to(source):
        raise ValueError("compiled export destination must be a new independent directory")
    external = []
    for path in source.rglob("*"):
        if path.is_symlink() and not path.resolve().is_relative_to(source):
            if str(path.relative_to(source)) != "lib/libLHAPDF.a":
                raise ValueError(f"unexpected shared dependency in compiled export: {path}")
            external.append({"link": str(path.relative_to(source)), "target": str(path.resolve())})
    evidence = {"source": str(source), "source_manifest": completion._artifact(manifest_path),
        "process_card": completion._artifact(card), "shared_read_only_libraries": external,
        "parameter_card": completion._artifact(source / "Cards/param_card.dat"),
        "copy_policy": "independent file copies; internal absolute links rebased; no hard links"}
    shutil.copytree(source, destination, symlinks=True)
    for path in destination.rglob("*"):
        if path.is_symlink() and path.resolve().is_relative_to(source):
            target = destination / path.resolve().relative_to(source)
            path.unlink()
            path.symlink_to(os.path.relpath(target, path.parent))
    for name in ("Cards/run_card.dat", "Cards/param_card.dat", "bin/internal/banner.py",
                 "bin/internal/madevent_interface.py"):
        if os.path.samefile(source / name, destination / name):
            raise ValueError("compiled export copy shares a writable file with its source")
    return evidence


def build(args, sample, manifest):
    directory = ho.sample_dir(args, sample)
    process = directory / "mg5_process"
    if not process.exists():
        reuse = getattr(args, "reuse_loop_export", None)
        if reuse is not None and not args.dry_run:
            manifest["compiled_export_reuse"] = clone_completed_export(reuse, process)
            ho.save_manifest(args, sample, manifest)
        else:
            ho.run_mg5(args, directory, "process.mg5", "mg5-process.log")
    elif not args.resume and "build" not in manifest["completed"]:
        raise ValueError(f"process exists: {process}; use --resume")
    if args.dry_run:
        return
    if not (process / "bin/generate_events").is_file():
        raise ValueError("loop process export did not produce generate_events")
    for path in (process / "bin/internal/banner.py", process / "bin/internal/madevent_interface.py"):
        if not path.resolve().is_relative_to(process.resolve()):
            raise ValueError("loop PDF repair must modify only this exported process")
    # Reuse the checked MG5 3.5 beam-PDF repair. The order switch here selects
    # PDF ID335900, while the physical accuracy remains explicitly LO.
    cfg = SimpleNamespace(dry_run=False, order="ho")
    four.patch_mg5_hidden_pdf_defaults(process, cfg)
    four.patch_mg5_pdf_label_sync(process, cfg)
    patch_run_card(process / "Cards/run_card.dat", run_settings(args))


def inspect_lhe(path, sample, args):
    result = ORIGINAL_INSPECT(path, sample, args)
    opener = gzip.open if Path(path).suffix == ".gz" else open
    count, processes, minimum_mass, maximum_mass = 0, [], math.inf, 0.
    root = None
    with opener(path, "rb") as stream:
        for action, node in ET.iterparse(stream, events=("start", "end")):
            if root is None:
                root = node
            if action != "end":
                continue
            if node.tag == "MG5ProcCard":
                processes = [line.strip() for line in (node.text or "").splitlines()
                             if line.strip().startswith(("generate ", "add process "))]
                if processes != [f"generate {PROCESS}"]:
                    raise ValueError("loop LHE process differs or fails to exclude the Higgs")
            elif node.tag == "event":
                rows = [line.split() for line in (node.text or "").splitlines()
                        if line.strip() and not line.lstrip().startswith("#")]
                particles = rows[1:1 + int(rows[0][0])]
                incoming = sorted(int(row[0]) for row in particles if int(row[1]) == -1)
                outgoing = sorted(int(row[0]) for row in particles if int(row[1]) == 1)
                if incoming != [21, 21] or outgoing != [22, 22] or any(int(row[0]) == 25 for row in particles):
                    raise ValueError("loop LHE must contain only gg incoming and two final photons, with no Higgs")
                if ho.number(rows[0][2]) <= 0:
                    raise ValueError("loop-induced LO source requires positive hard-event weights")
                photons = [row for row in particles if int(row[0]) == 22 and int(row[1]) == 1]
                momentum = [sum(ho.number(row[index]) for row in photons) for index in (6, 7, 8, 9)]
                mass = math.sqrt(max(0., momentum[3] ** 2 - sum(value ** 2 for value in momentum[:3])))
                lower, upper = getattr(args, "gen_mgg_min", 0.), getattr(args, "gen_mgg_max", -1.)
                if mass < lower - 1e-4 or (upper >= 0 and mass > upper + 1e-4):
                    raise ValueError("loop LHE event violates the declared Born diphoton mass interval")
                minimum_mass, maximum_mass = min(minimum_mass, mass), max(maximum_mass, mass)
                count += 1
                node.clear()
                root.clear()
    if not processes or count != result["events"]:
        raise ValueError("loop LHE lacks complete process/event provenance")
    result["process_validation"] = dict(process=PROCESS, incoming=[21, 21], outgoing=[22, 22],
                                         higgs_excluded=True, events_checked=count)
    result["sha256"] = four.sha256_file(Path(path))
    result["born_mgg_validation"] = {"minimum_observed_gev": minimum_mass, "maximum_observed_gev": maximum_mass,
        "minimum_cut_gev": getattr(args, "gen_mgg_min", 0.), "maximum_cut_gev": getattr(args, "gen_mgg_max", -1.),
        "serialization_tolerance_gev": 1e-4}
    return result


def generate(args, sample, manifest):
    directory = ho.sample_dir(args, sample)
    build(args, sample, manifest)
    process = directory / "mg5_process"
    if not lhe_path(args, sample).exists():
        existing = process / "Events" / args.run_tag
        if existing.exists() and any(existing.iterdir()) and not args.dry_run:
            raise ValueError(f"unfinished native generation at {existing}; preserve and inspect before recovery")
        ho.command(args, [args.mg5_python, process / "bin/generate_events", args.run_tag,
                          "-f", f"--nb_core={args.nb_core}"], process, "loop-generate.log")
    if not args.dry_run:
        manifest["lhe"] = inspect_lhe(lhe_path(args, sample), sample, args)
        manifest["generator_artifacts"] = {name: four.file_record(process / "Cards" / name)
                                            for name in ("run_card.dat", "param_card.dat")}


def shower(args, sample, manifest):
    directory = ho.sample_dir(args, sample)
    work = directory / "herwig"
    if not args.dry_run:
        current = inspect_lhe(lhe_path(args, sample), sample, args)
        if manifest.get("lhe", current) != current:
            raise ValueError("loop source changed since generation")
        manifest["lhe"] = current
        if list((work / "events").glob("**/*.root")):
            if not args.resume:
                raise ValueError("existing loop shower products require --resume and finalized EOF evidence")
            ho.record_shower_completion(args, sample, manifest, "source_exhausted")
            return
        plan = completion.finite_shower_plan(sample.name, args.nevents, current["events"])
        ho.write_json(work / "shower_execution.json", plan, args)
        (work / "events").mkdir(exist_ok=True)
    else:
        plan = completion.finite_shower_plan(sample.name, args.nevents, args.nevents)
    ho.command(args, [args.herwig, "--version"], work, "version.log")
    ho.command(args, [args.herwig, "read", f"{sample.name}.in"], work, "read.log")
    try:
        ho.command(args, [args.herwig, "run", f"{sample.name}.run", f"-N{plan['execution_events']}"], work, "run.log")
    except RuntimeError:
        if args.dry_run or completion.EOF_MESSAGE not in (work / "run.log").read_text(errors="replace"):
            raise
        ho.record_shower_completion(args, sample, manifest, "source_exhausted")
        return
    if not args.dry_run:
        raise ValueError("finite loop shower must terminate at first source EOF")


@contextmanager
def loop_adapter():
    """Scope the extra sample to this dedicated command; nominal defaults stay fixed."""
    changes = dict(SAMPLES=(SAMPLE,), parse_args=parse_args, lhe_path=lhe_path, mg5_cards=mg5_cards,
        herwig_card=herwig_card, prepare=prepare, build=build, generate=generate,
        inspect_lhe=inspect_lhe, shower=shower)
    original = {key: getattr(ho, key) for key in changes}
    try:
        for key, value in changes.items():
            setattr(ho, key, value)
        yield
    finally:
        for key, value in original.items():
            setattr(ho, key, value)


def main(argv=None):
    with loop_adapter():
        return ho.main(argv)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as error:
        raise SystemExit(str(error))
