#!/usr/bin/env python3
"""Prepare, build, benchmark or calculate the HO signal's inclusive ggF rate.

The default stage writes a calculation plan/cards only. No rate is installed
until all calculations, numerical checks and provenance checks succeed.
The pinned external/ihixs checkout is never modified: the PDF alpha_s adapter
is applied in a separate build-source copy, alongside an upstream benchmark.
"""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import statistics
import subprocess
import sys

try:
    from . import ho_signal_normalization as norm
    from . import run_gammagamma_campaign as runtime
except ImportError:
    import ho_signal_normalization as norm
    import run_gammagamma_campaign as runtime


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_SETTINGS = SCRIPT_DIR / "HOAnalysis/ihixs-ssc40.json"
DEFAULT_WORK = SCRIPT_DIR / "HOAnalysis/normalization/ihixs-ssc40"
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?"
PRECISION_REFINEMENT_LIMIT = 3

ADAPTER_REPLACEMENTS = {
    "src/core/input_parameters.cpp": (
        "_as_over_pi = _model.alpha_strong()/consts::Pi;",
        "// HiggsSSC: direct PDF alpha_s for pure HEFT with on-shell masses.\n"
        "    _as_over_pi = _lumi->alpha_s_at_scale(_mur)/consts::Pi;"),
    "src/tools/luminosity.h": (
        "double alpha_s_at_mz(){return _pdf[0]->alphasQ(91.1876);}",
        "double alpha_s_at_mz(){return _pdf[0]->alphasQ(91.1876);}\n"
        "    double alpha_s_at_scale(double q){return _pdf[0]->alphasQ(q);}"),
    "src/tools/vegas_adaptor.cpp": (
        "_prob);",
        "_prob);\n"
        "    if (_fail != 0) {\n"
        "        std::cerr << \"HiggsSSC: Cuba integration failed to converge (fail=\"\n"
        "                  << _fail << \", neval=\" << _neval << \")\" << std::endl;\n"
        "        exit(1);\n"
        "    }"),
}
ADAPTER_COUNTS = {"src/core/input_parameters.cpp": 1, "src/tools/luminosity.h": 1,
                  "src/tools/vegas_adaptor.cpp": 2}

# Keep this separate from the physics adapter: existing integrations retain
# their original executable, while parser-repaired builds have their own hash.
PARSER_REPAIRS = {
    "src/tools/user_interface.cpp": (
        ("new char[strlen(options[i].name.c_str())]",
         "new char[strlen(options[i].name.c_str()) + 1]"),
        ("long_options[N+1].", "long_options[N]."),
    ),
}
PARSER_REPAIR_COUNTS = {"src/tools/user_interface.cpp": (1, 4)}

# Probe the linked C++ library, avoiding a second Python LHAPDF installation.
PDF_PROBE = r'''#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <cstdlib>
std::string quote(const std::string& s) {
  std::string r="\"";
  for (char c: s) {
    if(c=='\\' || c=='"') {r+='\\'; r+=c;}
    else if(c=='\n') r+="\\n";
    else if(c=='\r') r+="\\r";
    else if(c=='\t') r+="\\t";
    else r+=c;
  }
  return r+'"';
}
int main(int argc, char** argv) {
  if(argc!=4) return 2;
  try {
    LHAPDF::setVerbosity(0);
    const std::string name=argv[1];
    const int member=std::stoi(argv[2]);
    const double q=std::stod(argv[3]);
    std::unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF(name,member));
    std::cout << std::setprecision(17) << "{\"version\":" << quote(LHAPDF::version());
    std::cout << ",\"info_path\":" << quote(LHAPDF::findpdfsetinfopath(name));
    std::cout << ",\"member_path\":" << quote(LHAPDF::findpdfmempath(name,member));
    for(const std::string key: {"AlphaS_MZ","MZ","AlphaS_OrderQCD","OrderQCD","NumMembers","ErrorType","NumFlavors"})
      std::cout << "," << quote(key) << ":" << quote(pdf->info().get_entry(key));
    std::cout << ",\"alpha_s_mur\":" << pdf->alphasQ(q);
    std::cout << ",\"alpha_s_91_1876\":" << pdf->alphasQ(91.1876);
    std::cout << ",\"photon\":" << (pdf->hasFlavor(22)?"true":"false") << "}\n";
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n'; return 1;}
}
'''


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def write_same_or_new(path, text):
    path = Path(path)
    if path.exists() and path.read_text() != text:
        raise ValueError(f"existing input changed: {path}; use a new --work-dir")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "build", "benchmark", "calculate"), default="prepare")
    parser.add_argument("--settings", type=Path, default=DEFAULT_SETTINGS)
    parser.add_argument("--ihixs-source", type=Path, default=REPO_ROOT / "external/ihixs")
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--record", type=Path, default=norm.DEFAULT_RECORD)
    parser.add_argument("--powheg-input", "--powheg-card", type=Path,
                        default=SCRIPT_DIR / "HOAnalysis/powheg-hjminnlo-ssc40-nnpdf40nnloqed.input")
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--cc", default="cc")
    parser.add_argument("--cxx", default="c++")
    parser.add_argument("--jobs", type=int, default=1, help="build parallelism; integrations remain sequential")
    parser.add_argument("--lhapdf-config", default="lhapdf-config")
    parser.add_argument("--lhapdf-dir", type=Path, help="LHAPDF installation prefix")
    parser.add_argument("--cuba-dir", type=Path, help="Cuba installation prefix")
    parser.add_argument("--boost-dir", type=Path, help="Boost header directory")
    parser.add_argument("--pdf-data-dir", type=Path, help="prepend a directory containing LHAPDF set folders")
    parser.add_argument("--herwig-env", type=Path)
    parser.add_argument("--herwig-module")
    parser.add_argument("--resume", action="store_true", help="reuse verified results and archive interrupted run directories")
    parser.add_argument("--refine-failed", action="store_true",
                        help="with calculate --resume, retry inaccurate points at tighter epsrel in separate directories")
    parser.add_argument("--repair-parser", action="store_true",
                        help="with --resume, build/use a separate getopt bounds repair and preserve completed runs")
    parser.add_argument("--dry-run", action="store_true", help="print the plan without running or writing anything")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.refine_failed and (args.stage != "calculate" or not args.resume):
        parser.error("--refine-failed requires --stage calculate --resume")
    if args.repair_parser and (args.stage == "prepare" or not args.resume):
        parser.error("--repair-parser requires --resume and --stage build, benchmark or calculate")
    for key in ("settings", "ihixs_source", "work_dir", "record", "powheg_input",
                "lhapdf_dir", "cuba_dir", "boost_dir", "pdf_data_dir", "herwig_env"):
        value = getattr(args, key)
        if value is not None:
            setattr(args, key, value.expanduser().resolve())
    if args.herwig_env and args.herwig_env.is_dir():
        args.herwig_env = args.herwig_env / "bin/activate"
    if args.herwig_env and args.herwig_module:
        parser.error("choose either --herwig-env or --herwig-module")
    return args


def command(args, argv, *, cwd=None, log=None):
    argv = [str(value) for value in argv]
    setup = runtime.runtime_setup_commands(args)
    if args.pdf_data_dir:
        setup.append("export LHAPDF_DATA_PATH=" + shlex.quote(str(args.pdf_data_dir))
                     + ':"${LHAPDF_DATA_PATH:-}"')
    shell = "\n".join([*setup, shlex.join(argv)])
    print("+ " + shlex.join(argv), flush=True)
    if log:
        with Path(log).open("w") as stream:
            result = subprocess.run(["bash", "-c", shell], cwd=cwd, env=runtime.clean_shell_env(),
                                    stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"command failed ({result.returncode}); see {log}")
        return Path(log).read_text(errors="replace")
    result = subprocess.run(["bash", "-c", shell], cwd=cwd, env=runtime.clean_shell_env(),
                            capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {result.stderr.strip()}")
    return result.stdout.strip()


def load_settings(path):
    settings = json.loads(Path(path).read_text())
    fixed = {"profile": norm.PROFILE, "ihixs_commit": norm.IHIXS_COMMIT,
             "sqrt_s_gev": 40000., "higgs_mass_gev": 125., "top_mass_gev": 173.2,
             "top_scheme": "on-shell", "gf_gev_minus2": norm.GF,
             "gev_minus2_to_pb": 389379660., "ihixs_prefactor_pb": 35.0309,
             "mur_gev": 62.5, "muf_gev": 62.5, "pdf_set": norm.PDF_SET,
             "native_pdf_set": norm.NATIVE_PDF_SET, "pdf_member": 0, "alpha_s_mz": .118,
             "qcd_order": "N3LO", "hard_model": "pure_heft", "electroweak_corrections": False,
             "finite_mass_corrections": False, "resummation": False}
    for key, expected in fixed.items():
        if settings.get(key) != expected:
            raise ValueError(f"settings {key} differs from the approved HO normalization profile")
    tolerance = norm.finite_number(settings["numerical_relative_tolerance"], "numerical tolerance", positive=True)
    if tolerance > norm.NUMERICAL_TOLERANCE:
        raise ValueError("numerical tolerance cannot exceed 0.05%")
    integration = settings["integration"]
    for key in ("epsrel", "mineval", "maxeval", "nstart", "nincrease"):
        norm.finite_number(integration[key], key, positive=True)
    if integration["epsrel"] > tolerance or integration["mineval"] > integration["maxeval"]:
        raise ValueError("incompatible integration accuracy/point limits")
    for key in ("epsabs", "cuba_verbose"):
        norm.finite_number(integration[key], key, nonnegative=True)
    benchmark = settings["benchmark"]
    if any(benchmark.get(key) != value for key, value in {
        "sqrt_s_gev": 13000., "higgs_mass_gev": 125., "pdf_set": "PDF4LHC15_nnlo_100",
        "pdf_member": 0, "result_key": "eftn3lo", "expected_pb": 45.1816,
        "reference": "https://arxiv.org/html/1802.00827v1"}.items()):
        raise ValueError("benchmark settings do not match the published ihixs example")
    if not 0 < benchmark["relative_tolerance"] <= .005:
        raise ValueError("benchmark tolerance cannot exceed 0.5%")
    return settings


def native_parameters(path):
    """Freeze supplied card and source defaults, including the native HOPPET choice."""
    text = Path(path).read_text()
    values = {}
    for line in text.splitlines():
        line = re.split(r"[!#]", line, maxsplit=1)[0].strip()
        match = re.fullmatch(rf"(\w+)\s+({NUMBER})", line)
        if match:
            values[match[1].lower()] = float(match[2].replace("d", "e").replace("D", "E"))
    # The real POWHEG input keyword is topmass; missing/nonpositive values
    # request the process's fallback. Store the effective mass, not its sentinel.
    values["topmass"] = values.get("topmass", -1.) if values.get("topmass", -1.) > 0 else 173.2
    defaults = {"topmass": 173.2, "bmass": 0., "cmass": 0.,
                "renscfact": 1., "facscfact": 1., "lhapdf_in_hoppet": 0.}
    required = {"hmass": 125., "hwidth": .004152, "ebeam1": 20000., "ebeam2": 20000.,
                "ih1": 1., "ih2": 1., "lhans1": 336100., "lhans2": 336100.,
                "alphas_from_pdf": 1., "use_NNLOPS_pdfs": 1., "minlo": 1., "minnlo": 1.,
                "renscfact": 1., "facscfact": 1., "topmass": 173.2}
    for key, expected in required.items():
        actual = values.get(key.lower(), defaults.get(key))
        if actual is None or not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError(f"POWHEG {key} differs from the inclusive normalization reference")
    if values.get("lhapdf_in_hoppet", 0.) != 0.:
        raise ValueError("native reference assumes hybrid HOPPET evolution, not lhapdf_in_hoppet")
    if "tmass" in values and not math.isclose(values["tmass"], 173.2, rel_tol=1e-10):
        raise ValueError("POWHEG tmass conflicts with the pure-HEFT top-mass reference")
    for key in ("quarkmasseffects", "nnloint", "nnlo", "nnlopsreweight"):
        if values.get(key, 0.) > 0:
            raise ValueError(f"POWHEG {key} changes the pure-HEFT/native normalization reference")
    result = {**defaults, **{key: values.get(key.lower(), defaults.get(key)) for key in required}}
    result.update(tmass=173.2, quarkmasseffects=0., nnloint=0., nnlo=0., nnlopsreweight=0.,
                  gf_gev_minus2=norm.GF, hard_model="pure_heft",
                  alpha_s_evolution="HOPPET three-loop reconstructed from NNLO PDF alpha_s(MZ)",
                  top_mass_default_source="HJMiNNLO/init_couplings.f: #topmass fallback ph_topmass=173.2 GeV")
    return result


def calculation_points(settings):
    base = {"pdf_set": settings["pdf_set"], "pdf_member": 0, "qcd_order": "N3LO",
            "mur_gev": 62.5, "muf_gev": 62.5}
    points = []
    for index, (r, f) in enumerate(norm.SCALE_POINTS):
        points.append({**base, "label": "central" if index == 0 else f"scale_{index}",
                       "mur_gev": 62.5 * r, "muf_gev": 62.5 * f})
    points.extend({**base, "label": f"replica_{i:03d}", "pdf_member": i} for i in range(1, 101))
    points.extend({**base, "label": label, "pdf_set": settings["native_pdf_set"], "qcd_order": order}
                  for label, order in (("nnlo_native", "NNLO"), ("n3lo_nnlo_pdf", "N3LO")))
    return points


def production_card(settings, point):
    options = {"Etot": settings["sqrt_s_gev"], "m_higgs": settings["higgs_mass_gev"],
               "mur": point["mur_gev"], "muf": point["muf_gev"], "pdf_set": point["pdf_set"],
               "pdf_member": point["pdf_member"], "qcd_perturbative_order": point["qcd_order"],
               "qcd_order_evol": 3 if point["pdf_set"] == norm.PDF_SET else 2,
               "with_eft": True, "with_fixed_as_at_mz": 0., "top_scheme": "on-shell",
               "bottom_scheme": "on-shell", "charm_scheme": "on-shell",
               "mt_on_shell": settings["top_mass_gev"], "y_top": 1., "y_bot": 0., "y_charm": 0.,
               "with_exact_qcd_corrections": False, "with_ew_corrections": False,
               "with_mt_expansion": False, "with_delta_pdf_th": False,
               "with_scale_variation": False, "with_pdf_error": False, "with_a_s_error": False,
               "with_resummation": False, "with_scet": False, "with_indiv_mass_effects": False,
               "with_lower_ord_scale_var": False, "with_eft_channel_info": False,
               "output_filename": "ihixs.out", "verbose": "minimal", **settings["integration"]}
    return "# Pure HEFT production: read eftn3lo/eftnnlo, never Higgs XS.\n" + "".join(
        f"{key} = {str(value).lower() if isinstance(value, bool) else value}\n" for key, value in options.items())


def prefactor_multiplier(settings):
    return (settings["gf_gev_minus2"] * math.pi / (math.sqrt(2.) * 288.)
            * settings["gev_minus2_to_pb"] / settings["ihixs_prefactor_pb"])


def parse_output(text, order="N3LO"):
    key = {"N3LO": "eftn3lo", "NNLO": "eftnnlo"}[order]
    match = re.findall(rf"(?m)^\s*{key}\s*=\s*({NUMBER})\s*\[\s*({NUMBER})\s*\]\s*$", text)
    if len(match) != 1:
        raise ValueError(f"expected exactly one raw {key} result with its integration error")
    values = [float(item.replace("d", "e").replace("D", "E")) for item in match[0]]
    norm.finite_number(values[0], "raw EFT cross section", positive=True)
    norm.finite_number(values[1], "raw EFT integration error", nonnegative=True)
    def scalar(name):
        matches = re.findall(rf"(?m)^\s*{name}\s*=\s*({NUMBER})\s*$", text)
        if len(matches) != 1:
            raise ValueError(f"missing/ambiguous {name} in saved ihixs output")
        return norm.finite_number(float(matches[0].replace("d", "e").replace("D", "E")), name, positive=True)
    return {"result_key": key, "raw_cross_section_pb": values[0], "raw_cross_section_error_pb": values[1],
            "alpha_s_mur": scalar("as_at_mur"), "alpha_s_at_91_1876": scalar("as_at_mz")}


def prepare(args, settings):
    parameters = native_parameters(args.powheg_input)
    points = calculation_points(settings)
    if args.dry_run:
        print(json.dumps({"stage": args.stage, "work_dir": str(args.work_dir), "record": str(args.record),
                          "calculations": len(points), "points": points, "native_parameters": parameters}, indent=2))
        return
    args.work_dir.mkdir(parents=True, exist_ok=True)
    plan = {"settings_sha256": norm.canonical_sha256(settings), "settings": settings,
            "powheg_input_sha256": file_sha256(args.powheg_input), "native_parameters": parameters,
            "points": points}
    write_same_or_new(args.work_dir / "plan.json", json.dumps(plan, indent=2, sort_keys=True) + "\n")
    write_same_or_new(args.work_dir / "powheg.reference.input", args.powheg_input.read_text())
    for point in points:
        write_same_or_new(args.work_dir / "cards" / f"{point['label']}.card", production_card(settings, point))
    print(f"Prepared {len(points)} integrations. No cross-section record was created.")


def source_inventory(args):
    if not (args.ihixs_source / "CMakeLists.txt").is_file():
        raise ValueError("ihixs submodule is missing; initialize external/ihixs before building")
    commit = command(args, ["git", "-C", args.ihixs_source, "rev-parse", "HEAD"])
    if commit != norm.IHIXS_COMMIT:
        raise ValueError("ihixs checkout differs from its pinned commit")
    dirty = command(args, ["git", "-C", args.ihixs_source, "status", "--porcelain", "--untracked-files=no"])
    if dirty:
        raise ValueError("ihixs has tracked local changes; preserve them and use the pinned source")
    paths = command(args, ["git", "-C", args.ihixs_source, "ls-files"]).splitlines()
    hashes = {name: file_sha256(args.ihixs_source / name) for name in paths}
    return paths, hashes


def validate_cuba_prefix(prefix):
    """Explain invalid explicit prefixes before writing a build plan."""
    if prefix is None:
        return  # CMake can discover a system installation.
    headers = (prefix / "cuba.h", prefix / "include/cuba.h")
    libraries = [directory / name for directory in (prefix, prefix / "lib")
                 for name in ("libcuba.a", "libcuba.so", "libcuba.dylib")]
    if not any(path.is_file() for path in headers) or not any(path.is_file() for path in libraries):
        raise ValueError(f"--cuba-dir {prefix} is not a Cuba installation: expected cuba.h and libcuba. "
                         "Install Cuba 4.2 and use its real prefix; see HOAnalysis/README.md.")


def parser_repair_sha256():
    return norm.canonical_sha256({"replacements": PARSER_REPAIRS, "counts": PARSER_REPAIR_COUNTS})


def patched_source(name, original, variant, parser_repair=False):
    if variant == "lhapdf" and name in ADAPTER_REPLACEMENTS:
        before, after = ADAPTER_REPLACEMENTS[name]
        if original.count(before) != ADAPTER_COUNTS[name]:
            raise ValueError(f"adapter anchor changed in {name}")
        original = original.replace(before, after)
    if parser_repair and name in PARSER_REPAIRS:
        for (before, after), count in zip(PARSER_REPAIRS[name], PARSER_REPAIR_COUNTS[name]):
            if original.count(before) != count:
                raise ValueError(f"parser repair anchor changed in {name}")
            original = original.replace(before, after)
    return original


def parser_repair_args(args):
    result = copy.copy(args)
    result.work_dir = args.work_dir / "parser-repair"
    result.repair_parser = False
    return result


def build(args, settings, *, parser_repair=False):
    paths, hashes = source_inventory(args)
    source_hash = norm.canonical_sha256(hashes)
    adapter_hash = norm.canonical_sha256({"replacements": ADAPTER_REPLACEMENTS, "counts": ADAPTER_COUNTS})
    identity = {"ihixs_commit": norm.IHIXS_COMMIT, "source_sha256": source_hash,
                "adapter_sha256": adapter_hash, "cc": args.cc, "cxx": args.cxx,
                "lhapdf_dir": str(args.lhapdf_dir), "cuba_dir": str(args.cuba_dir),
                "boost_dir": str(args.boost_dir), "probe_source_sha256": hashlib.sha256(PDF_PROBE.encode()).hexdigest()}
    if parser_repair:
        identity["parser_patch_sha256"] = parser_repair_sha256()
    manifest_path = args.work_dir / "build-manifest.json"
    if manifest_path.exists():
        prior = load_build(args)
        if prior["identity"] != identity:
            raise ValueError("build inputs changed; use a new --work-dir")
        print("Reusing verified upstream/production builds.")
        return
    for variant in ("upstream", "lhapdf"):
        source = args.work_dir / f"source-{variant}"
        for name in paths:
            target = source / name
            target.parent.mkdir(parents=True, exist_ok=True)
            if (variant == "lhapdf" and name in ADAPTER_REPLACEMENTS
                    or parser_repair and name in PARSER_REPAIRS):
                original = (args.ihixs_source / name).read_text()
                write_same_or_new(target, patched_source(name, original, variant, parser_repair))
            elif target.exists():
                if file_sha256(target) != hashes[name]:
                    raise ValueError(f"build source changed: {target}; use a new --work-dir")
            else:
                shutil.copy2(args.ihixs_source / name, target)
        directory = args.work_dir / f"build-{variant}"
        directory.mkdir(parents=True, exist_ok=True)
        configure = [args.cmake, "-S", source, "-B", directory,
                     "-DCMAKE_POLICY_VERSION_MINIMUM=3.5", f"-DCMAKE_C_COMPILER={args.cc}",
                     f"-DCMAKE_CXX_COMPILER={args.cxx}"]
        for key, value in (("LHAPDF_DIR", args.lhapdf_dir), ("CUBA_DIR_USER", args.cuba_dir),
                           ("BOOST_DIR_USER", args.boost_dir)):
            if value:
                configure.append(f"-D{key}={value}")
        if not args.lhapdf_dir:
            prefix = command(args, [args.lhapdf_config, "--prefix"])
            if not prefix or "\n" in prefix:
                raise ValueError("lhapdf-config did not return one installation prefix")
            configure.append(f"-DLHAPDF_DIR={prefix}")
        command(args, configure, log=directory / "configure.log")
        command(args, [args.cmake, "--build", directory, "--target", "ihixs", "--parallel", args.jobs],
                log=directory / "build.log")
    cache = (args.work_dir / "build-lhapdf/CMakeCache.txt").read_text()
    def cache_path(key):
        match = re.search(rf"(?m)^{key}:[^=]+=(.+)$", cache)
        if not match or "NOTFOUND" in match[1]:
            raise ValueError(f"CMake did not resolve {key}")
        return Path(match[1].strip())
    includes = cache_path("LHAPDF_INC_DIR")
    library = cache_path("LHAPDF_LIB_NAMES")
    probe_source = args.work_dir / "lhapdf_probe.cpp"
    probe_source.write_text(PDF_PROBE)
    probe = args.work_dir / "lhapdf_probe"
    command(args, [args.cxx, "-std=c++11", probe_source, f"-I{includes}", library,
                   f"-Wl,-rpath,{library.parent}", "-o", probe], log=args.work_dir / "probe-build.log")
    binaries = {variant: file_sha256(args.work_dir / f"build-{variant}/ihixs") for variant in ("upstream", "lhapdf")}
    write_json(manifest_path, {"schema_version": 1, "identity": identity, "source_files": hashes,
                              "executable_sha256": binaries, "probe_sha256": file_sha256(probe),
                              "lhapdf_library_path": str(library),
                              "lhapdf_library_sha256": file_sha256(library)})


def load_build(args):
    path = args.work_dir / "build-manifest.json"
    if not path.is_file():
        raise ValueError("verified ihixs builds are missing; run --stage build first")
    data = json.loads(path.read_text())
    for variant in ("upstream", "lhapdf"):
        if file_sha256(args.work_dir / f"build-{variant}/ihixs") != data["executable_sha256"][variant]:
            raise ValueError("ihixs executable changed since its recorded build")
    if file_sha256(args.work_dir / "lhapdf_probe") != data["probe_sha256"]:
        raise ValueError("LHAPDF probe changed since its recorded build")
    if file_sha256(data["lhapdf_library_path"]) != data["lhapdf_library_sha256"]:
        raise ValueError("linked LHAPDF library changed since the recorded build")
    return data


def load_parser_repair(args, base_build):
    repaired_args = parser_repair_args(args)
    data = load_build(repaired_args)
    expected = {**base_build["identity"], "parser_patch_sha256": parser_repair_sha256()}
    if data["identity"] != expected or data["source_files"] != base_build["source_files"]:
        raise ValueError("parser repair must use the original source, physics adapter and build settings")
    for key in ("lhapdf_library_path", "lhapdf_library_sha256"):
        if data[key] != base_build[key]:
            raise ValueError("parser repair must link the original LHAPDF library")
    return repaired_args, data


def pdf_probe(args, point, expected_order=None):
    payload = json.loads(command(args, [args.work_dir / "lhapdf_probe", point["pdf_set"],
                                       point["pdf_member"], point["mur_gev"]]))
    if expected_order is not None:
        if int(payload["OrderQCD"]) != expected_order or int(payload["AlphaS_OrderQCD"]) != expected_order:
            raise ValueError("PDF/alpha_s perturbative order differs from the required PDF")
        if not math.isclose(float(payload["AlphaS_MZ"]), .118, rel_tol=1e-10):
            raise ValueError("PDF alpha_s(MZ) differs from 0.118")
        if not payload["photon"] or int(payload["NumFlavors"]) != 5:
            raise ValueError("expected five-flavour QED PDFs with a photon")
    if point["pdf_set"] == norm.PDF_SET:
        if int(payload["NumMembers"]) != 101 or payload["ErrorType"].lower() != "replicas":
            raise ValueError("expected the 100-replica NNPDF aN3LO QED MHOU set")
    payload["info_sha256"] = file_sha256(payload["info_path"])
    payload["member_sha256"] = file_sha256(payload["member_path"])
    return payload


def run_card(args, label, card_text, executable, identity):
    directory = args.work_dir / "runs" / label
    complete = directory / "complete.json"
    input_hash = hashlib.sha256(card_text.encode()).hexdigest()
    expected = {**identity, "input_sha256": input_hash, "executable_sha256": file_sha256(executable)}
    if complete.is_file():
        data = json.loads(complete.read_text())
        if data["identity"] != expected or file_sha256(directory / "input.card") != input_hash:
            raise ValueError(f"completed {label} inputs/build changed; use a new --work-dir")
        if file_sha256(directory / "ihixs.out") != data["output_sha256"]:
            raise ValueError(f"completed {label} raw output hash changed")
        if not args.resume:
            raise ValueError(f"{label} already completed; use --resume to reuse it")
        return directory, data
    if directory.exists() and any(directory.iterdir()):
        if not args.resume:
            raise ValueError(f"partial run {directory}; --resume archives it before retrying")
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        backup = args.work_dir / "resume-backups" / f"{label}-{stamp}"
        backup.parent.mkdir(parents=True, exist_ok=True)
        directory.rename(backup)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "input.card").write_text(card_text)
    command(args, [executable, "-i", "input.card"], cwd=directory, log=directory / "ihixs.log")
    if not (directory / "ihixs.out").is_file():
        raise ValueError(f"ihixs did not write its numerical result in {directory}")
    data = {"identity": expected, "output_sha256": file_sha256(directory / "ihixs.out")}
    write_json(complete, data)
    return directory, data


def benchmark(args, settings):
    build_data = load_build(args)
    card = (args.work_dir / "source-upstream/runcard/default.card").read_text()
    # Use the paper's unmodified default physics and explicitly pin saved output.
    card = re.sub(r"(?m)^output_filename\s*=.*$", "output_filename = ihixs.out", card)
    pdf = pdf_probe(args, {"pdf_set": "PDF4LHC15_nnlo_100", "pdf_member": 0, "mur_gev": 62.5})
    directory, run = run_card(args, "benchmark", card, args.work_dir / "build-upstream/ihixs",
                              {"settings_sha256": norm.canonical_sha256(settings), "variant": "upstream",
                               "pdf_info_sha256": pdf["info_sha256"], "pdf_member_sha256": pdf["member_sha256"]})
    result = parse_output((directory / "ihixs.out").read_text())
    reference = settings["benchmark"]
    delta = abs(result["raw_cross_section_pb"] / reference["expected_pb"] - 1.)
    passed = (delta <= reference["relative_tolerance"]
              and result["raw_cross_section_error_pb"] / result["raw_cross_section_pb"] <= .0005)
    payload = {**result, "passed": passed, "relative_difference": delta, **reference,
               "ihixs_commit": norm.IHIXS_COMMIT, "executable_sha256": build_data["executable_sha256"]["upstream"],
               "production_executable_sha256": build_data["executable_sha256"]["lhapdf"],
               "input_sha256": run["identity"]["input_sha256"], "output_sha256": run["output_sha256"],
               "pdf": pdf, "variant": "unmodified upstream; original internal alpha_s evolution"}
    if "parser_patch_sha256" in build_data["identity"]:
        payload["parser_patch_sha256"] = build_data["identity"]["parser_patch_sha256"]
        payload["variant"] = "upstream physics with getopt bounds repaired; original internal alpha_s evolution"
    write_json(args.work_dir / "benchmark.json", payload)
    if not passed:
        raise ValueError("published ihixs benchmark failed; inspect benchmark.json and raw output")
    print("Published raw-EFT benchmark passed. Production rates are not yet calculated.")


def load_benchmark(args, build_data):
    path = args.work_dir / "benchmark.json"
    if not path.is_file():
        raise ValueError(f"benchmark is missing in {args.work_dir}; run --stage benchmark first")
    data = json.loads(path.read_text())
    if data.get("passed") is not True or data["ihixs_commit"] != norm.IHIXS_COMMIT:
        raise ValueError("published benchmark must pass before calculating the production rate")
    for key, artifact in (("input_sha256", args.work_dir / "runs/benchmark/input.card"),
                          ("output_sha256", args.work_dir / "runs/benchmark/ihixs.out")):
        if file_sha256(artifact) != data[key]:
            raise ValueError("benchmark artifact changed")
    if (data["executable_sha256"] != build_data["executable_sha256"]["upstream"]
            or data["production_executable_sha256"] != build_data["executable_sha256"]["lhapdf"]):
        raise ValueError("build changed since benchmark validation")
    if data.get("parser_patch_sha256") != build_data["identity"].get("parser_patch_sha256"):
        raise ValueError("benchmark parser repair differs from its build")
    return data


def select_production_execution(args, label, executable, identity, repair_context):
    """Reuse a completed attempt only with its original, verified executable."""
    if repair_context is None:
        return executable, identity
    repaired_args, _ = repair_context
    complete = args.work_dir / "runs" / label / "complete.json"
    if complete.is_file():
        cached = json.loads(complete.read_text())["identity"]
        variant = cached["variant"]
        if variant == "lhapdf":
            return executable, identity
        if variant != "lhapdf-parser-repair":
            raise ValueError(f"completed {label} has an unrecognized build variant")
    # No completion marker means the process did not finish successfully.
    # run_card will archive a partial directory before using the repaired build.
    return repaired_args.work_dir / "build-lhapdf/ihixs", {**identity, "variant": "lhapdf-parser-repair"}


def run_precision_checked(args, settings, point, executable, identity, pdf, *, repair_context=None):
    """Preserve every attempt and tighten only points failing the total-rate gate."""
    multiplier = prefactor_multiplier(settings)
    enabled = getattr(args, "refine_failed", False)
    refinements = PRECISION_REFINEMENT_LIMIT if enabled else 0
    base_epsrel = settings["integration"]["epsrel"]
    for attempt in range(refinements + 1):
        effective = copy.deepcopy(settings)
        run_label = point["label"]
        if attempt:
            # EFT terms use Cuhre, with per-term accuracy relaxations up to
            # 10000. Extra Vegas statistics do not tighten these targets.
            epsrel = min(base_epsrel / 10., 1.e-6) / 10. ** (attempt - 1)
            effective["integration"]["epsrel"] = epsrel
            effective["integration"]["epsabs"] *= epsrel / base_epsrel
            run_label += f"__precision_{attempt}"
        run_identity = {**identity, "settings_sha256": norm.canonical_sha256(effective)}
        selected_executable, run_identity = select_production_execution(
            args, run_label, executable, run_identity, repair_context)
        directory, completed = run_card(
            args, run_label, production_card(effective, point), selected_executable, run_identity)
        result = parse_output((directory / "ihixs.out").read_text(), point["qcd_order"])
        if not math.isclose(result["alpha_s_mur"], pdf["alpha_s_mur"], rel_tol=1e-7, abs_tol=1e-10):
            raise ValueError(f"{point['label']}: hard alpha_s does not match LHAPDF alphasQ(muR)")
        if not math.isclose(result["alpha_s_at_91_1876"], pdf["alpha_s_91_1876"], rel_tol=1e-7, abs_tol=1e-10):
            raise ValueError(f"{point['label']}: input alpha_s(MZ) does not match the selected PDF")
        sigma = result["raw_cross_section_pb"] * multiplier
        error = result["raw_cross_section_error_pb"] * multiplier
        relative_error = error / sigma
        if relative_error <= settings["numerical_relative_tolerance"]:
            result.update(cross_section_pb=sigma, cross_section_error_pb=error,
                          integration=effective["integration"], precision_refinements=attempt,
                          run_directory=str(directory.relative_to(args.work_dir)))
            return directory, completed, result
        print(f"{point['label']}: retained {run_label} with numerical error "
              f"{100. * relative_error:.6g}% (required <= "
              f"{100. * settings['numerical_relative_tolerance']:.6g}%)", flush=True)
        if attempt < refinements:
            next_epsrel = min(base_epsrel / 10., 1.e-6) / 10. ** attempt
            print(f"{point['label']}: retrying only this point with epsrel={next_epsrel:g}", flush=True)
    hint = ("precision refinements exhausted; inspect the saved attempts"
            if enabled else "use --resume --refine-failed to retry this point with tighter epsrel")
    raise ValueError(f"{point['label']}: numerical uncertainty exceeds "
                     f"{100. * settings['numerical_relative_tolerance']:g}%; {hint}")


def check_repaired_central(args, settings, original, repair_context):
    """Require unchanged central physics before accepting mixed-build results."""
    if original["build_variant"] != "lhapdf":
        raise ValueError("parser repair requires a retained central result from the original build")
    repaired_args, repaired_build = repair_context
    effective = copy.deepcopy(settings)
    effective["integration"] = copy.deepcopy(original["integration"])
    point = {**calculation_points(settings)[0], "label": "central_check"}
    pdf = pdf_probe(repaired_args, point, 3)
    if (pdf["info_sha256"] != original["pdf_info_sha256"]
            or pdf["member_sha256"] != original["pdf_member_sha256"]):
        raise ValueError("repaired central check must use the original central PDF grids")
    identity = {"settings_sha256": norm.canonical_sha256(effective),
                "pdf_info_sha256": pdf["info_sha256"], "pdf_member_sha256": pdf["member_sha256"],
                "probe_sha256": repaired_build["probe_sha256"], "variant": "lhapdf-parser-repair"}
    directory, completed, result = run_precision_checked(
        repaired_args, effective, point, repaired_args.work_dir / "build-lhapdf/ihixs", identity, pdf)
    delta = abs(result["cross_section_pb"] / original["cross_section_pb"] - 1.)
    tolerance = settings["numerical_relative_tolerance"]
    if delta > tolerance:
        raise ValueError("parser-repaired central cross section differs from the saved central result "
                         f"by {100. * delta:.6g}% (required <= {100. * tolerance:g}%); "
                         "preserve both outputs and investigate before using the normalization")
    print(f"Parser-repaired 40 TeV central check passed (difference {100. * delta:.6g}%).", flush=True)
    return {**point, **result, "passed": True, "relative_difference": delta, "relative_tolerance": tolerance,
            "original_cross_section_pb": original["cross_section_pb"],
            "original_cross_section_error_pb": original["cross_section_error_pb"],
            "pdf_alpha_s_mur": pdf["alpha_s_mur"],
            "input_sha256": completed["identity"]["input_sha256"],
            "output_sha256": completed["output_sha256"],
            "pdf_info_sha256": pdf["info_sha256"], "pdf_member_sha256": pdf["member_sha256"],
            "executable_sha256": completed["identity"]["executable_sha256"],
            "build_variant": completed["identity"]["variant"],
            "run_directory": str(directory.relative_to(args.work_dir))}


def calculate(args, settings):
    build_data = load_build(args)
    benchmark_path = args.work_dir / "benchmark.json"
    benchmark_data = load_benchmark(args, build_data)
    repair_context = None
    if getattr(args, "repair_parser", False):
        repair_context = load_parser_repair(args, build_data)
        repaired_args, repaired_build = repair_context
        repaired_benchmark = load_benchmark(repaired_args, repaired_build)
        central_complete = args.work_dir / "runs/central/complete.json"
        if not central_complete.is_file():
            raise ValueError("parser repair requires a verified original central completion for comparison")
        central_identity = json.loads(central_complete.read_text())["identity"]
        if (central_identity["variant"] != "lhapdf"
                or central_identity["executable_sha256"] != build_data["executable_sha256"]["lhapdf"]):
            raise ValueError("parser repair central baseline must come from the original production build")
    result_runs = []
    pdfs = {}
    multiplier = prefactor_multiplier(settings)
    for point in calculation_points(settings):
        pdf = pdf_probe(args, point, 3 if point["pdf_set"] == norm.PDF_SET else 2)
        pdf_key = f"{point['pdf_set']}/{point['pdf_member']}"
        if pdf_key not in pdfs:
            pdfs[pdf_key] = {**pdf, "alpha_s_by_q_gev": {}}
        elif any(pdfs[pdf_key][key] != pdf[key] for key in ("info_sha256", "member_sha256", "version")):
            raise ValueError("PDF grid/library provenance changed during the calculation")
        pdfs[pdf_key]["alpha_s_by_q_gev"][str(point["mur_gev"])] = pdf["alpha_s_mur"]
        identity = {"settings_sha256": norm.canonical_sha256(settings),
                    "pdf_info_sha256": pdf["info_sha256"], "pdf_member_sha256": pdf["member_sha256"],
                    "probe_sha256": build_data["probe_sha256"], "variant": "lhapdf"}
        directory, completed, result = run_precision_checked(
            args, settings, point, args.work_dir / "build-lhapdf/ihixs", identity, pdf,
            repair_context=repair_context)
        result_runs.append({**point, **result,
                            "pdf_alpha_s_mur": pdf["alpha_s_mur"],
                            "input_sha256": completed["identity"]["input_sha256"],
                            "output_sha256": completed["output_sha256"],
                            "executable_sha256": completed["identity"]["executable_sha256"],
                            "build_variant": completed["identity"]["variant"],
                            "pdf_info_sha256": pdf["info_sha256"], "pdf_member_sha256": pdf["member_sha256"]})
        if point["label"] == "central" and repair_context is not None:
            repaired_central = check_repaired_central(args, settings, result_runs[0], repair_context)
    central = result_runs[0]
    scale_rates = [run["cross_section_pb"] for run in result_runs[:7]]
    replicas = [run["cross_section_pb"] for run in result_runs if run["label"].startswith("replica_")]
    artifact_paths = [args.work_dir / "plan.json", args.work_dir / "powheg.reference.input",
                      args.work_dir / "build-manifest.json", benchmark_path,
                      args.work_dir / "runs/benchmark/input.card", args.work_dir / "runs/benchmark/ihixs.out",
                      args.work_dir / "runs/benchmark/ihixs.log"]
    for run in result_runs:
        artifact_paths.extend(args.work_dir / run["run_directory"] / name
                              for name in ("input.card", "ihixs.out", "ihixs.log", "complete.json"))
    if repair_context is not None:
        artifact_paths.extend(repaired_args.work_dir / name for name in
                              ("build-manifest.json", "benchmark.json", "runs/benchmark/input.card",
                               "runs/benchmark/ihixs.out", "runs/benchmark/ihixs.log", "runs/benchmark/complete.json"))
        artifact_paths.extend(args.work_dir / repaired_central["run_directory"] / name
                              for name in ("input.card", "ihixs.out", "ihixs.log", "complete.json"))
    payload = {"schema_version": 1, "profile": norm.PROFILE, "normalization_kind": "ihixs_n3lo",
               "process": "ggF", "created_utc": datetime.now(timezone.utc).isoformat(),
               **{key: settings[key] for key in ("sqrt_s_gev", "higgs_mass_gev", "top_mass_gev", "top_scheme",
                  "gf_gev_minus2", "mur_gev", "muf_gev", "pdf_set", "pdf_member", "alpha_s_mz",
                  "qcd_order", "hard_model", "electroweak_corrections", "finite_mass_corrections", "resummation",
                  "numerical_relative_tolerance")},
               "branching_fraction_included": False, "alpha_s_source": "LHAPDF::PDF::alphasQ",
               "cross_section_pb": central["cross_section_pb"], "cross_section_error_pb": central["cross_section_error_pb"],
               "prefactor_correction": {"multiplier": multiplier, "upstream_prefactor_pb": 35.0309,
                                        "gev_minus2_to_pb": settings["gev_minus2_to_pb"],
                                        "formula": "[GF*pi/(sqrt(2)*288)*GeV^-2_to_pb]/35.0309"},
               "native_parameters": native_parameters(args.powheg_input),
               "uncertainties": {"scale_up_pb": max(scale_rates) - central["cross_section_pb"],
                                 "scale_down_pb": central["cross_section_pb"] - min(scale_rates),
                                 "pdf_stddev_pb": statistics.stdev(replicas),
                                 "pdf_replica_mean_pb": statistics.mean(replicas),
                                 "prescription": "7-point scale envelope; sample standard deviation of members 1-100; separate numerical error"},
               "runs": result_runs, "pdfs": pdfs,
               "validations": {"benchmark": benchmark_data, "alpha_s": {"passed": True,
                    "relative_tolerance": 1e-7, "checked_runs": len(result_runs),
                    "prescription": "compare saved as_at_mur against the linked LHAPDF library for every member and scale"}},
               "provenance": {"ihixs_commit": norm.IHIXS_COMMIT, "settings_sha256": norm.canonical_sha256(settings),
                    "source_sha256": build_data["identity"]["source_sha256"],
                    "adapter_sha256": build_data["identity"]["adapter_sha256"],
                    "executable_sha256": build_data["executable_sha256"]["lhapdf"],
                    "build_manifest_sha256": file_sha256(args.work_dir / "build-manifest.json"),
                    "powheg_input_sha256": file_sha256(args.powheg_input),
                    "artifact_root": os.path.relpath(args.work_dir, args.record.parent),
                    "artifacts": [{"path": str(path.relative_to(args.work_dir)), "sha256": file_sha256(path)}
                                  for path in artifact_paths],
                    "adapter": ADAPTER_REPLACEMENTS, "settings": settings}}
    if repair_context is not None:
        payload["validations"]["parser_repair_benchmark"] = repaired_benchmark
        payload["validations"]["parser_repair_central"] = repaired_central
        payload["provenance"]["parser_repair"] = {
            "ihixs_commit": norm.IHIXS_COMMIT,
            "source_sha256": repaired_build["identity"]["source_sha256"],
            "adapter_sha256": repaired_build["identity"]["adapter_sha256"],
            "parser_patch_sha256": repaired_build["identity"]["parser_patch_sha256"],
            "executable_sha256": repaired_build["executable_sha256"]["lhapdf"],
            "build_manifest_sha256": file_sha256(repaired_args.work_dir / "build-manifest.json"),
            "patches": PARSER_REPAIRS,
        }
    payload = norm.seal_record(payload)
    norm.validate_record(payload)
    if args.record.exists():
        previous = norm.load_record(args.record)
        if previous["cross_section_pb"] != payload["cross_section_pb"] or previous["provenance"] != payload["provenance"]:
            raise ValueError("an existing normalization record differs; preserve it and select a new --record path")
        print("Existing verified normalization record already matches these completed calculations.")
        return
    write_json(args.record, payload)
    print(f"Installed validated inclusive production normalization: {args.record}")


def main(argv=None):
    args = parse_args(argv)
    settings = load_settings(args.settings)
    if args.dry_run:
        prepare(args, settings)
        return 0
    if args.stage == "build":
        validate_cuba_prefix(args.cuba_dir)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    with (args.work_dir / ".normalization.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("another ihixs normalization command is using this work directory") from exc
        prepare(args, settings)
        if args.stage == "build":
            if args.repair_parser:
                base_build = load_build(args)
                for key in ("cc", "cxx", "lhapdf_dir", "cuba_dir", "boost_dir"):
                    value = getattr(args, key)
                    if (str(value) if isinstance(value, Path) or value is None else value) != base_build["identity"][key]:
                        raise ValueError(f"parser repair {key} must match the original build")
                repaired_args = parser_repair_args(args)
                repaired_args.work_dir.mkdir(parents=True, exist_ok=True)
                build(repaired_args, settings, parser_repair=True)
                load_parser_repair(args, base_build)
            else:
                build(args, settings)
        elif args.stage == "benchmark":
            if args.repair_parser:
                repaired_args, _ = load_parser_repair(args, load_build(args))
                benchmark(repaired_args, settings)
            else:
                benchmark(args, settings)
        elif args.stage == "calculate":
            calculate(args, settings)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
