#!/usr/bin/env python3
"""Run the LO SSC H -> ZZ* -> four-lepton campaign.

The runner keeps hard production, forced decays, shower/detector simulation,
and post-analysis normalization as separate, auditable stages.  Generated
products are resumable and intentionally live below ``hfourlepton/runs``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_MG5_DIR = REPO_ROOT / "MG5_aMC_v3_5_15"
DEFAULT_HERWIG_PREFIX = Path.home() / "Projects/Herwig/Herwig-REAL-stable-gcc-full"
DEFAULT_LINUX_HERWIG_MODULE = "herwig/stable"
DEFAULT_DARWIN_HERWIG_MODULE = "herwig/730"

MANIFEST_SCHEMA_VERSION = 1
HIGGS_MASS_GEV = 125.09
BR_H_TO_4L = 1.251e-4
SIGNAL_FRACTION_PROVENANCE = "MG5_aMC_3.5.15_sm_direct_1to4_integration_mH125.09"
DEFAULT_SIGNAL_CHANNEL_FRACTIONS: dict[str, float] = {
    "4e": 0.263127188,
    "4mu": 0.263127188,
    "2e2mu": 0.473745624,
}

LO_PDF_NAME = "NNPDF40_lo_as_01180"
LO_PDF_LHAID = 331900
DETECTOR_RESPONSES = ("perfect", "ssc")
SAMPLE_SETS = ("irreducible", "reducible", "all")
SOURCE_BACKENDS = ("mg5", "external-lhe")
HEAVY_FLAVOUR_BIASES = ("unbiased", "semileptonic")
PENDING_SIGNAL_WEIGHT = "TARGET_XSEC_OVER_MADSPIN_DECAYED_XSEC"
MG5_PDF_LABEL_SYNC_MARKER = "# four-lepton campaign PDF label sync"
MG5_POST_TREATCARDS_PDF_PATCH_MARKER = (
    "# four-lepton campaign post-treatcards PDF include patch"
)

GENERATION_CUTS = {
    "ptl_min_gev": 5.0,
    "etal_max": 3.0,
    "mmll_min_gev": 4.0,
    "drll_min": 0.0,
}

CUT_MASK_DEFINITION = (
    "bit0:at_least_four_reconstructed;"
    "bit1:pt_eta_acceptance;"
    "bit2:electron_crack_acceptance;"
    "bit3:isolation;"
    "bit4:disjoint_ossf_pairing;"
    "bit5:all_ossf_masses_gt_4;"
    "bit6:z1_70_100;"
    "bit7:z2_10_100;"
    "bit8:trigger_weight_applied;"
    "bit9:selected"
)

CUT_STAGES = [
    ("four_reconstructed", "At least four reconstructed leptons", 0),
    ("fiducial", "Lepton pT and eta", 1),
    ("electron_crack", "Electron crack veto", 2),
    ("isolation", "Lepton isolation", 3),
    ("ossf_pairs", "Two disjoint OSSF pairs", 4),
    ("low_mass_veto", "All OSSF masses above 4 GeV", 5),
    ("z1_window", "70 < mZ1 < 100 GeV", 6),
    ("z2_window", "10 < mZ2 < 100 GeV", 7),
    ("trigger", "Trigger efficiency", 8),
    ("selected", "Final selection", 9),
]


@dataclass(frozen=True)
class SampleSpec:
    """Declarative hard-process and normalization description."""

    name: str
    label: str
    category: str
    channel: str
    perturbative_order: str
    source_kind: str
    model: str
    process: str
    decay_strategy: str
    normalization_rule: str
    response_mode: str
    production_group_id: str
    madspin_decay: str = ""
    prompt_decay_profile: str = ""
    heavy_flavour_candidate: bool = False
    bias_of: str | None = None

    @property
    def is_signal(self) -> bool:
        return self.category == "signal"

    @property
    def is_hf_biased(self) -> bool:
        return self.bias_of is not None


def _signal(name: str, label: str, channel: str, decay: str) -> SampleSpec:
    return SampleSpec(
        name=name,
        label=label,
        category="signal",
        channel=channel,
        perturbative_order="LO",
        source_kind="mg5",
        model="loop_sm",
        process="g g > h [noborn=QCD]",
        decay_strategy="madspin_direct_1to4_spinmode_none",
        normalization_rule="stable_ggh_xsec_times_external_h4l_channel_br",
        response_mode="fourlepton",
        production_group_id="signal_gg_h_stable",
        madspin_decay=decay,
    )


def _qq(name: str, label: str, channel: str, final_state: str) -> SampleSpec:
    return SampleSpec(
        name=name,
        label=label,
        category="irreducible",
        channel=channel,
        perturbative_order="LO",
        source_kind="mg5",
        model="sm",
        process=f"p p > {final_state} / h",
        decay_strategy="matrix_element_final_state_higgs_excluded",
        normalization_rule="native_lhe_cross_section",
        response_mode="fourlepton",
        production_group_id=name,
    )


def _gg(name: str, label: str, channel: str, final_state: str) -> SampleSpec:
    return SampleSpec(
        name=name,
        label=label,
        category="irreducible",
        channel=channel,
        perturbative_order="LO_loop_induced",
        source_kind="mg5",
        model="loop_sm",
        process=f"g g > {final_state} / h [noborn=QCD]",
        decay_strategy="matrix_element_final_state_higgs_excluded",
        normalization_rule="native_lhe_cross_section",
        response_mode="fourlepton",
        production_group_id=name,
    )


def _ttz(
    z_flavour: str,
    wplus_flavour: str,
    wminus_flavour: str,
    channel: str,
) -> SampleSpec:
    stratum = f"z{z_flavour}_wp{wplus_flavour}_wm{wminus_flavour}"
    return SampleSpec(
        name=f"bkg_ttz_{stratum}",
        label=f"ttZ ({stratum})",
        category="reducible",
        channel=channel,
        perturbative_order="LO",
        source_kind="mg5",
        model="sm",
        process="p p > t t~ z",
        decay_strategy="herwig_disjoint_prompt_lepton_stratum_plus_heavy_hadron_decays",
        normalization_rule="native_lhe_cross_section_with_br_reweighted_prompt_decays",
        response_mode="fourlepton",
        production_group_id="bkg_ttz_stable",
        prompt_decay_profile=stratum,
        heavy_flavour_candidate=True,
    )


SAMPLES: tuple[SampleSpec, ...] = (
    _signal("signal_gg_h_4e", "ggH, H->4e", "4e", "h > e+ e- e+ e-"),
    _signal("signal_gg_h_4mu", "ggH, H->4mu", "4mu", "h > mu+ mu- mu+ mu-"),
    _signal(
        "signal_gg_h_2e2mu",
        "ggH, H->2e2mu",
        "2e2mu",
        "h > e+ e- mu+ mu-",
    ),
    _qq("bkg_qq_4e", "qq->4e", "4e", "e+ e- e+ e-"),
    _qq("bkg_qq_4mu", "qq->4mu", "4mu", "mu+ mu- mu+ mu-"),
    _qq("bkg_qq_2e2mu", "qq->2e2mu", "2e2mu", "e+ e- mu+ mu-"),
    _gg("bkg_gg_4e", "gg->4e (no H)", "4e", "e+ e- e+ e-"),
    _gg("bkg_gg_4mu", "gg->4mu (no H)", "4mu", "mu+ mu- mu+ mu-"),
    _gg("bkg_gg_2e2mu", "gg->2e2mu (no H)", "2e2mu", "e+ e- mu+ mu-"),
    SampleSpec(
        name="bkg_zbb_ee",
        label="Z/gamma*+bb, ee",
        category="reducible",
        channel="inclusive",
        perturbative_order="LO",
        source_kind="mg5",
        model="sm",
        process="p p > e+ e- b b~",
        decay_strategy="matrix_element_zll_plus_herwig_heavy_hadron_decays",
        normalization_rule="native_lhe_cross_section",
        response_mode="fourlepton",
        production_group_id="bkg_zbb_ee",
        heavy_flavour_candidate=True,
    ),
    SampleSpec(
        name="bkg_zbb_mumu",
        label="Z/gamma*+bb, mumu",
        category="reducible",
        channel="inclusive",
        perturbative_order="LO",
        source_kind="mg5",
        model="sm",
        process="p p > mu+ mu- b b~",
        decay_strategy="matrix_element_zll_plus_herwig_heavy_hadron_decays",
        normalization_rule="native_lhe_cross_section",
        response_mode="fourlepton",
        production_group_id="bkg_zbb_mumu",
        heavy_flavour_candidate=True,
    ),
    _ttz("ee", "e", "e", "4e"),
    _ttz("ee", "e", "mu", "inclusive"),
    _ttz("ee", "mu", "e", "inclusive"),
    _ttz("ee", "mu", "mu", "2e2mu"),
    _ttz("mumu", "e", "e", "2e2mu"),
    _ttz("mumu", "e", "mu", "inclusive"),
    _ttz("mumu", "mu", "e", "inclusive"),
    _ttz("mumu", "mu", "mu", "4mu"),
    SampleSpec(
        name="bkg_ttbar",
        label="ttbar",
        category="reducible",
        channel="inclusive",
        perturbative_order="LO",
        source_kind="mg5",
        model="sm",
        process="p p > t t~",
        decay_strategy="herwig_prompt_emu_plus_heavy_hadron_decays",
        normalization_rule="native_lhe_cross_section_with_br_reweighted_prompt_decays",
        response_mode="fourlepton",
        production_group_id="bkg_ttbar",
        prompt_decay_profile="w_emu",
        heavy_flavour_candidate=True,
    ),
)

SAMPLES_BY_NAME = {sample.name: sample for sample in SAMPLES}
RUN_SAMPLE_ALIASES: dict[str, set[str]] = {
    "signal": {sample.name for sample in SAMPLES if sample.category == "signal"},
    "bkg_qq": {sample.name for sample in SAMPLES if sample.name.startswith("bkg_qq_")},
    "bkg_gg": {sample.name for sample in SAMPLES if sample.name.startswith("bkg_gg_")},
    "bkg_zbb": {sample.name for sample in SAMPLES if sample.name.startswith("bkg_zbb_")},
    "bkg_ttz": {sample.name for sample in SAMPLES if sample.name.startswith("bkg_ttz_")},
}


@dataclass
class Config:
    nevents: int
    sample_nevents: dict[str, int]
    sample_set: str
    run_samples: str
    detector_response: str
    source_backend: str
    external_lhe: dict[str, Path]
    external_cross_sections_pb: dict[str, float]
    external_perturbative_order: str
    signal_channel_fractions: dict[str, float]
    signal_channel_fraction_source: str
    signal_channel_fraction_input: str
    ebeam: float
    luminosity_fb: float
    higgs_mass: float
    mg5_dir: Path
    mg5_python: str | None
    herwig: str
    herwig_env: Path | None
    herwig_module: str | None
    herwig_cxx: str | None
    herwig_pdf: str
    herwig_b_decays: Path | None
    heavy_flavour_bias: str
    bias_control_events: int
    run_tag: str
    dry_run: bool
    force: bool
    nb_core: int
    seed_base: int
    collier_library: Path | None
    muon_resolution_scale: float
    pileup_noise: bool

    @property
    def template_in(self) -> Path:
        return SCRIPT_DIR / "HW-template.in"

    @property
    def run_dir(self) -> Path:
        return SCRIPT_DIR / "runs" / self.run_tag

    @property
    def manifest_path(self) -> Path:
        return self.run_dir / f"manifest_{self.detector_response}.json"

    @property
    def analysis_code_dir(self) -> Path:
        return SCRIPT_DIR / "LOAnalysis" / "Code"

    @property
    def analysis_target(self) -> str:
        return "HwSimPostAnalysis_fourlepton"

    @property
    def analysis_exe(self) -> Path:
        return self.analysis_code_dir / self.analysis_target

    @property
    def analysis_source(self) -> Path:
        return self.analysis_code_dir / f"{self.analysis_target}.cc"

    @property
    def decay_reweighter_dir(self) -> Path:
        return SCRIPT_DIR / "Herwig"

    @property
    def decay_reweighter_source(self) -> Path:
        return self.decay_reweighter_dir / "LHEBranchingRatioReweighter.cc"

    @property
    def decay_reweighter_plugin(self) -> Path:
        return self.decay_reweighter_dir / "LHEBranchingRatioReweighter.so"

    def events_for(self, sample: SampleSpec) -> int:
        if sample.name in self.sample_nevents:
            return self.sample_nevents[sample.name]
        if (
            self.heavy_flavour_bias == "semileptonic"
            and sample.heavy_flavour_candidate
            and not sample.is_hf_biased
        ):
            return self.bias_control_events
        return self.nevents


def env_text(name: str, default: str) -> str:
    return os.environ.get(name, default)


def env_int(name: str, default: int) -> int:
    return int(os.environ.get(name, str(default)))


def env_float(name: str, default: float) -> float:
    return float(os.environ.get(name, str(default)))


def env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    return default if value is None else value.lower() in {"1", "true", "yes", "on"}


def clean_optional_text(value: str | None) -> str | None:
    value = value.strip() if value is not None else ""
    return value or None


def maybe_path(value: str | None) -> Path | None:
    return Path(value).expanduser() if value else None


def parse_key_values(
    values: Iterable[str],
    *,
    value_name: str,
    converter: type[int] | type[float] | type[str] = str,
) -> dict[str, int | float | str]:
    parsed: dict[str, int | float | str] = {}
    for item in values:
        if "=" not in item:
            raise ValueError(f"expected NAME={value_name}, got {item!r}")
        key, raw_value = item.split("=", 1)
        key = key.strip()
        raw_value = raw_value.strip()
        if not key or not raw_value:
            raise ValueError(f"expected NAME={value_name}, got {item!r}")
        if key in parsed:
            raise ValueError(f"duplicate assignment for {key!r}")
        parsed[key] = converter(raw_value)
    return parsed


def parse_signal_channel_fractions(text: str) -> dict[str, float]:
    parsed = parse_key_values(
        [part.strip() for part in text.split(",") if part.strip()],
        value_name="FRACTION",
        converter=float,
    )
    expected = set(DEFAULT_SIGNAL_CHANNEL_FRACTIONS)
    if set(parsed) != expected:
        raise ValueError(
            "signal fractions must define exactly " + ",".join(sorted(expected))
        )
    fractions = {key: float(value) for key, value in parsed.items()}
    if any(not math.isfinite(value) or value <= 0.0 for value in fractions.values()):
        raise ValueError("all signal channel fractions must be finite and positive")
    if not math.isclose(sum(fractions.values()), 1.0, rel_tol=0.0, abs_tol=1e-6):
        raise ValueError("signal channel fractions must sum to one within 1e-6")
    return fractions


def default_signal_fraction_text() -> str:
    return ",".join(
        f"{channel}={fraction:.9f}"
        for channel, fraction in DEFAULT_SIGNAL_CHANNEL_FRACTIONS.items()
    )


def guess_default_herwig_env() -> Path | None:
    explicit = maybe_path(os.environ.get("DEFAULT_HERWIG_ENV"))
    if explicit and explicit.exists():
        return explicit
    candidate = DEFAULT_HERWIG_PREFIX / "bin" / "activate"
    return candidate if candidate.exists() else None


def guess_herwig_b_decays() -> Path | None:
    explicit = maybe_path(os.environ.get("HERWIG_B_DECAYS"))
    if explicit:
        return explicit
    candidate = DEFAULT_HERWIG_PREFIX / "share" / "Herwig" / "defaults" / "HerwigBDecays.in"
    return candidate if candidate.exists() else None


def guess_collier_library(herwig_env: Path | None) -> Path | None:
    explicit = os.environ.get("COLLIER_LIBRARY") or os.environ.get("COLLIER_DYLIB")
    if explicit:
        return Path(explicit).expanduser()
    prefix = DEFAULT_HERWIG_PREFIX
    if herwig_env and herwig_env.name == "activate":
        prefix = herwig_env.parent.parent
    system = platform.system()
    names = ("libcollier.2.dylib", "libcollier.dylib") if system == "Darwin" else (
        "libcollier.so",
        "libcollier.so.2",
    )
    for name in names:
        for candidate in sorted((prefix / "opt").glob(f"OpenLoops-*/lib/{name}")):
            if candidate.exists():
                return candidate
    return None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nevents", type=int, default=env_int("NEVENTS", 10000))
    parser.add_argument(
        "--sample-nevents",
        action="append",
        default=[],
        metavar="NAME=N",
        help="per-sample event override; repeat as needed",
    )
    parser.add_argument("--sample-set", choices=SAMPLE_SETS, default=env_text("SAMPLE_SET", "irreducible"))
    parser.add_argument("--run-samples", default=env_text("RUN_SAMPLES", "all"))
    parser.add_argument(
        "--detector-response",
        "--response-profile",
        dest="detector_response",
        choices=DETECTOR_RESPONSES,
        default=env_text("DETECTOR_RESPONSE", "perfect"),
    )
    parser.add_argument(
        "--source-backend",
        choices=SOURCE_BACKENDS,
        default=env_text("SOURCE_BACKEND", "mg5"),
    )
    parser.add_argument(
        "--external-lhe",
        action="append",
        default=[],
        metavar="SAMPLE_OR_GROUP=PATH",
    )
    parser.add_argument(
        "--external-cross-section",
        action="append",
        default=[],
        metavar="SAMPLE_OR_GROUP=PB",
        help=(
            "explicit cross section for external LHE files whose XSECUP is "
            "missing or a placeholder; signal group keys set stable production "
            "and signal sample keys set the MadSpin-decayed native cross section"
        ),
    )
    parser.add_argument(
        "--external-order",
        choices=("LO", "NLO"),
        default=env_text("EXTERNAL_ORDER", "LO").upper(),
        help="declared hard-process order for external LHE provenance",
    )
    parser.add_argument(
        "--signal-channel-fractions",
        default=env_text("SIGNAL_CHANNEL_FRACTIONS", default_signal_fraction_text()),
        metavar="4e=X,4mu=Y,2e2mu=Z",
    )
    parser.add_argument("--ebeam", type=float, default=env_float("EBEAM", 20000.0))
    parser.add_argument("--luminosity-fb", type=float, default=env_float("LUMINOSITY_FB", 100.0))
    parser.add_argument("--higgs-mass", type=float, default=env_float("HIGGS_MASS", HIGGS_MASS_GEV))
    parser.add_argument("--mg5-dir", type=Path, default=Path(env_text("MG5_DIR", str(DEFAULT_MG5_DIR))))
    parser.add_argument("--mg5-python", default=clean_optional_text(os.environ.get("MG5_PYTHON")))
    parser.add_argument("--herwig", default=env_text("HERWIG", "Herwig"))
    parser.add_argument("--herwig-env", type=Path, default=maybe_path(os.environ.get("HERWIG_ENV")))
    parser.add_argument("--herwig-module", default=os.environ.get("HERWIG_MODULE"))
    parser.add_argument(
        "--herwig-cxx",
        default=clean_optional_text(os.environ.get("HERWIG_CXX")),
        help="C++ compiler executable for the LHE decay-weight plugin",
    )
    parser.add_argument("--no-herwig-module", action="store_true", default=env_bool("NO_HERWIG_MODULE", False))
    parser.add_argument("--herwig-pdf", default=env_text("HERWIG_PDF", LO_PDF_NAME))
    parser.add_argument("--herwig-b-decays", type=Path, default=guess_herwig_b_decays())
    parser.add_argument(
        "--heavy-flavour-bias",
        choices=HEAVY_FLAVOUR_BIASES,
        default=env_text("HEAVY_FLAVOUR_BIAS", "unbiased"),
    )
    parser.add_argument(
        "--bias-control-events",
        type=int,
        default=env_int("BIAS_CONTROL_EVENTS", 1000),
        help="unbiased nominal events retained when semileptonic strata are enabled",
    )
    parser.add_argument("--run-tag", default=env_text("RUN_TAG", "run_01"))
    parser.add_argument("--dry-run", action="store_true", default=env_bool("DRY_RUN", False))
    parser.add_argument("--force", action="store_true", default=env_bool("FORCE", False))
    parser.add_argument("--nb-core", type=int, default=env_int("NB_CORE", 1))
    parser.add_argument("--seed-base", type=int, default=env_int("SEED_BASE", 31122002))
    parser.add_argument("--collier-library", type=Path, default=None)
    parser.add_argument(
        "--muon-resolution-scale",
        type=float,
        default=env_float("MUON_RESOLUTION_SCALE", 1.0),
    )
    parser.add_argument(
        "--no-pileup-noise",
        action="store_true",
        default=env_bool("NO_PILEUP_NOISE", False),
        help="disable the GEM pileup-noise term while retaining thermal noise",
    )
    return parser


def parse_config(argv: Sequence[str]) -> Config:
    parser = build_parser()
    args = parser.parse_args(argv)
    for option, value, choices in (
        ("--sample-set", args.sample_set, SAMPLE_SETS),
        ("--detector-response", args.detector_response, DETECTOR_RESPONSES),
        ("--source-backend", args.source_backend, SOURCE_BACKENDS),
        ("--heavy-flavour-bias", args.heavy_flavour_bias, HEAVY_FLAVOUR_BIASES),
        ("--external-order", args.external_order, ("LO", "NLO")),
    ):
        if value not in choices:
            parser.error(
                f"argument {option}: invalid choice: {value!r} "
                f"(choose from {', '.join(choices)})"
            )
    try:
        sample_nevents_raw = parse_key_values(
            args.sample_nevents, value_name="N", converter=int
        )
        external_lhe_raw = parse_key_values(
            args.external_lhe, value_name="PATH", converter=str
        )
        external_cross_section_raw = parse_key_values(
            args.external_cross_section,
            value_name="PB",
            converter=float,
        )
        fractions = parse_signal_channel_fractions(args.signal_channel_fractions)
    except ValueError as error:
        parser.error(str(error))

    sample_nevents = {name: int(value) for name, value in sample_nevents_raw.items()}
    known_names = set(SAMPLES_BY_NAME) | {
        f"{sample.name}_hfbiased" for sample in SAMPLES if sample.heavy_flavour_candidate
    }
    unknown_overrides = sorted(set(sample_nevents) - known_names)
    if unknown_overrides:
        parser.error("unknown --sample-nevents sample(s): " + ",".join(unknown_overrides))
    if args.nevents <= 0 or any(value <= 0 for value in sample_nevents.values()):
        parser.error("event counts must be positive")
    if any(
        not math.isfinite(float(value)) or float(value) <= 0.0
        for value in external_cross_section_raw.values()
    ):
        parser.error("--external-cross-section values must be finite and positive")
    if args.bias_control_events <= 0:
        parser.error("--bias-control-events must be positive")
    if not math.isclose(args.higgs_mass, HIGGS_MASS_GEV, abs_tol=1e-9):
        parser.error(f"this campaign fixes --higgs-mass to {HIGGS_MASS_GEV}")
    if args.ebeam <= 0 or args.luminosity_fb <= 0:
        parser.error("beam energy and luminosity must be positive")
    if args.muon_resolution_scale <= 0:
        parser.error("--muon-resolution-scale must be positive")
    if re.fullmatch(r"[A-Za-z0-9_.-]+", args.run_tag) is None:
        parser.error("--run-tag may contain only letters, digits, '.', '_', and '-'")

    explicit_herwig_env = args.herwig_env.expanduser() if args.herwig_env else None
    if explicit_herwig_env is not None and explicit_herwig_env.is_dir():
        explicit_herwig_env = None
    herwig_module = clean_optional_text(args.herwig_module)
    if herwig_module is None and explicit_herwig_env is None and not args.no_herwig_module:
        if platform.system() == "Linux":
            herwig_module = DEFAULT_LINUX_HERWIG_MODULE
        elif platform.system() == "Darwin":
            herwig_module = DEFAULT_DARWIN_HERWIG_MODULE
    if args.no_herwig_module:
        herwig_module = None
    herwig_env = None if herwig_module else explicit_herwig_env
    if herwig_env is None and herwig_module is None:
        herwig_env = guess_default_herwig_env()

    collier = args.collier_library or guess_collier_library(herwig_env)
    fraction_cli_override = any(
        argument == "--signal-channel-fractions"
        or argument.startswith("--signal-channel-fractions=")
        for argument in argv
    )
    fraction_source = (
        "cli_override"
        if fraction_cli_override
        else "environment_override"
        if "SIGNAL_CHANNEL_FRACTIONS" in os.environ
        else "default_mg5_integration"
    )
    return Config(
        nevents=args.nevents,
        sample_nevents=sample_nevents,
        sample_set=args.sample_set,
        run_samples=args.run_samples,
        detector_response=args.detector_response,
        source_backend=args.source_backend,
        external_lhe={key: Path(str(value)).expanduser() for key, value in external_lhe_raw.items()},
        external_cross_sections_pb={
            key: float(value) for key, value in external_cross_section_raw.items()
        },
        external_perturbative_order=args.external_order,
        signal_channel_fractions=fractions,
        signal_channel_fraction_source=fraction_source,
        signal_channel_fraction_input=args.signal_channel_fractions,
        ebeam=args.ebeam,
        luminosity_fb=args.luminosity_fb,
        higgs_mass=args.higgs_mass,
        mg5_dir=args.mg5_dir.expanduser(),
        mg5_python=clean_optional_text(args.mg5_python),
        herwig=args.herwig,
        herwig_env=herwig_env,
        herwig_module=herwig_module,
        herwig_cxx=clean_optional_text(args.herwig_cxx),
        herwig_pdf=args.herwig_pdf,
        herwig_b_decays=args.herwig_b_decays.expanduser() if args.herwig_b_decays else None,
        heavy_flavour_bias=args.heavy_flavour_bias,
        bias_control_events=args.bias_control_events,
        run_tag=args.run_tag,
        dry_run=args.dry_run,
        force=args.force,
        nb_core=args.nb_core,
        seed_base=args.seed_base,
        collier_library=collier.expanduser() if collier else None,
        muon_resolution_scale=args.muon_resolution_scale,
        pileup_noise=not args.no_pileup_noise,
    )


def q(value: object) -> str:
    return shlex.quote(str(value))


def format_command(args: Sequence[object]) -> str:
    return " ".join(q(arg) for arg in args)


def log(message: str) -> None:
    print(f"\n==> {message}")


def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def clean_shell_env(env: Mapping[str, str] | None = None) -> dict[str, str]:
    cleaned = dict(os.environ if env is None else env)
    for key in list(cleaned):
        if key.startswith("BASH_FUNC_") and key.endswith("%%"):
            cleaned.pop(key, None)
    return cleaned


def ensure_dir(path: Path, cfg: Config) -> None:
    if cfg.dry_run:
        print(f"+ mkdir -p {q(path)}")
    else:
        path.mkdir(parents=True, exist_ok=True)


def run(
    args: Sequence[object],
    cfg: Config,
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
) -> None:
    prefix = f"+ cd {q(cwd)} && " if cwd else "+ "
    print(prefix + format_command(args))
    if not cfg.dry_run:
        subprocess.run(
            [str(arg) for arg in args],
            cwd=cwd,
            env=dict(env) if env is not None else None,
            check=True,
        )


def module_setup_snippet() -> str:
    return (
        "if ! type module >/dev/null 2>&1; then "
        "for module_init in "
        "/opt/homebrew/opt/modules/init/bash "
        "/etc/profile.d/modules.sh "
        "/usr/share/Modules/init/bash "
        "/usr/share/lmod/lmod/init/bash; do "
        "[ -r \"$module_init\" ] && source \"$module_init\" && break; "
        "done; fi"
    )


def sanitize_toolchain_exports_snippet() -> str:
    """Discard stale compiler exports left by a relocated runtime prefix.

    Herwig activation scripts can legitimately select a compiler toolchain,
    so valid exports are retained.  If the command token no longer resolves,
    unsetting it lets Make/MG5 fall back to the compilers currently on PATH.
    """

    return (
        "for compiler_var in CC CXX FC F77 F90; do "
        "eval \"compiler_value=\\${${compiler_var}-}\"; "
        "compiler_command=${compiler_value%% *}; "
        "if [ -n \"$compiler_command\" ] && "
        "! command -v \"$compiler_command\" >/dev/null 2>&1; then "
        "unset \"$compiler_var\"; "
        "fi; "
        "done"
    )


def runtime_setup_commands(cfg: Config) -> list[str]:
    commands = ["set -e"]
    if cfg.herwig_env:
        commands.append(f"source {q(cfg.herwig_env)}")
    if cfg.herwig_module:
        commands.extend([module_setup_snippet(), f"module load {q(cfg.herwig_module)}"])
    if cfg.herwig_env or cfg.herwig_module:
        commands.append(sanitize_toolchain_exports_snippet())
    return commands


def run_runtime_command(
    args: Sequence[object],
    cfg: Config,
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
) -> None:
    if not cfg.herwig_env and not cfg.herwig_module:
        run(args, cfg, cwd=cwd, env=env)
        return
    command = "; ".join([*runtime_setup_commands(cfg), format_command(args)])
    run(
        ["bash", "-lc", command],
        cfg,
        cwd=cwd,
        env=clean_shell_env(env),
    )


def capture_runtime_stdout(args: Sequence[object], cfg: Config) -> str | None:
    """Run a small runtime query and return stripped stdout on success."""

    if cfg.dry_run:
        return None
    if not cfg.herwig_env and not cfg.herwig_module:
        command = [str(arg) for arg in args]
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    else:
        shell_command = "; ".join(
            [*runtime_setup_commands(cfg), format_command(args)]
        )
        result = subprocess.run(
            ["bash", "-lc", shell_command],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            env=clean_shell_env(),
            check=False,
        )
    if result.returncode != 0:
        return None
    value = result.stdout.strip()
    return value or None


def _compiler_version_key(path: Path) -> tuple[int, ...]:
    match = re.search(r"-(\d+(?:\.\d+)*)$", path.name)
    if not match:
        return ()
    return tuple(int(part) for part in match.group(1).split("."))


def resolve_compiler_executable(command_text: str) -> Path | None:
    """Resolve a compiler command, including a newer versioned replacement."""

    try:
        words = shlex.split(command_text)
    except ValueError:
        return None
    if not words:
        return None
    token = words[0]
    configured = Path(token).expanduser()
    if configured.is_absolute():
        if configured.is_file() and os.access(configured, os.X_OK):
            return configured
        search_dirs = [configured.parent]
        name = configured.name
    else:
        found = shutil.which(token)
        if found:
            return Path(found)
        search_dirs = [
            Path(part)
            for part in os.environ.get("PATH", "").split(os.pathsep)
            if part
        ]
        name = token

    match = re.fullmatch(r"(.+)-\d+(?:\.\d+)*", name)
    if match is None:
        return None
    stem = match.group(1)
    replacements: list[Path] = []
    for directory in search_dirs:
        for candidate in directory.glob(f"{stem}-*"):
            if (
                _compiler_version_key(candidate)
                and candidate.is_file()
                and os.access(candidate, os.X_OK)
            ):
                replacements.append(candidate)
    return max(replacements, key=_compiler_version_key) if replacements else None


def herwig_plugin_compiler(cfg: Config) -> Path:
    requested = cfg.herwig_cxx
    if requested is None:
        requested = capture_runtime_stdout(["lhapdf-config", "--cxx"], cfg)
    if requested:
        resolved = resolve_compiler_executable(requested)
        if resolved is not None:
            return resolved
        die(
            "cannot resolve the Herwig/LHAPDF ABI compiler "
            f"{requested!r}; set --herwig-cxx to a compatible executable"
        )
    fallback = shutil.which("c++")
    if fallback is None:
        die("no C++ compiler found for the LHE branching-ratio plugin")
    return Path(fallback)


def _runtime_command_path(command: str, cfg: Config) -> Path | None:
    """Resolve a command in the same environment used for Herwig and MG5."""

    resolved = resolve_compiler_executable(command)
    if resolved is not None:
        return resolved
    try:
        words = shlex.split(command)
    except ValueError:
        return None
    if not words:
        return None
    located = capture_runtime_stdout(
        ["sh", "-c", f"command -v {q(words[0])}"],
        cfg,
    )
    return resolve_compiler_executable(located) if located else None


def _compiler_works(compiler: Path, cfg: Config) -> bool:
    return capture_runtime_stdout([compiler, "--version"], cfg) is not None


def mg5_fortran_compiler(cfg: Config, cxx_compiler: Path) -> Path:
    """Choose a working gfortran, preferring the C++ toolchain's version."""

    candidates: list[str] = []
    exported_fc = capture_runtime_stdout(
        ["sh", "-c", 'printf "%s" "${FC-}"'],
        cfg,
    )
    if exported_fc:
        candidates.append(exported_fc)

    cxx_version = _compiler_version_key(cxx_compiler)
    if cxx_version:
        version = ".".join(str(part) for part in cxx_version)
        candidates.extend(
            [
                str(cxx_compiler.parent / f"gfortran-{version}"),
                f"gfortran-{version}",
            ]
        )
    candidates.append("gfortran")

    seen: set[Path] = set()
    for candidate in candidates:
        resolved = _runtime_command_path(candidate, cfg)
        if resolved is None or resolved in seen:
            continue
        seen.add(resolved)
        if _compiler_works(resolved, cfg):
            return resolved
    die(
        "cannot resolve a working gfortran for the generated MG5 process; "
        "install gfortran in the Herwig runtime environment"
    )


def rewrite_mg5_make_opts(
    text: str,
    *,
    cxx_compiler: Path | str,
    fortran_compiler: Path | str,
    system: str,
) -> str:
    """Rewrite MG5's generated compiler defaults without touching other flags."""

    updates = {
        "DEFAULT_CPP_COMPILER": str(cxx_compiler),
        "DEFAULT_F_COMPILER": str(fortran_compiler),
    }
    if system == "Darwin":
        # Generated MG5 3.x templates default to clang/libc++ on macOS.
        # A GCC-built Herwig/LHAPDF stack must instead link its matching
        # libstdc++, and GCC does not understand ``-stdlib=libc++``.
        updates.update(
            {
                "STDLIB": "-lstdc++",
                "STDLIB_FLAG": "",
            }
        )

    found: set[str] = set()
    output: list[str] = []
    for line in text.splitlines(keepends=True):
        bare = line.rstrip("\r\n")
        newline = line[len(bare) :]
        rewritten = line
        for key, value in updates.items():
            match = re.match(rf"^(\s*){re.escape(key)}(\s*=\s*).*$", bare)
            if match is None:
                continue
            rewritten = (
                f"{match.group(1)}{key}{match.group(2)}{value}{newline}"
            )
            found.add(key)
            break
        output.append(rewritten)

    missing = sorted(set(updates) - found)
    if missing:
        raise ValueError(
            "generated MG5 Source/make_opts is missing "
            + ", ".join(missing)
        )
    return "".join(output)


def patch_mg5_toolchain(process_dir: Path, cfg: Config) -> None:
    """Patch generated MG5 make options before any Source target is built."""

    make_opts = process_dir / "Source" / "make_opts"
    if cfg.dry_run:
        cxx = cfg.herwig_cxx or "<Herwig/LHAPDF ABI C++ compiler>"
        print(
            f"+ patch MG5 toolchain in {make_opts}: "
            f"CXX={cxx} FC=<working matching gfortran>"
        )
        if platform.system() == "Darwin":
            print("  use libstdc++ and remove the generated libc++ flag")
        return
    if not make_opts.exists():
        die(f"missing generated MG5 make options: {make_opts}")

    cxx = herwig_plugin_compiler(cfg)
    fortran = mg5_fortran_compiler(cfg, cxx)
    try:
        rewritten = rewrite_mg5_make_opts(
            make_opts.read_text(),
            cxx_compiler=cxx,
            fortran_compiler=fortran,
            system=platform.system(),
        )
    except ValueError as error:
        die(str(error))
    make_opts.write_text(rewritten)
    log(f"MG5 toolchain: CXX={cxx}, FC={fortran}")


def build_decay_reweighter(cfg: Config) -> None:
    compiler: Path | str = (
        cfg.herwig_cxx or "c++" if cfg.dry_run else herwig_plugin_compiler(cfg)
    )
    run_runtime_command(
        [
            "make",
            "-C",
            cfg.decay_reweighter_dir,
            f"HERWIG_CXX={compiler}",
            cfg.decay_reweighter_plugin.name,
        ],
        cfg,
    )
    if not cfg.dry_run and not cfg.decay_reweighter_plugin.exists():
        die(
            "decay reweighter build did not create "
            f"{cfg.decay_reweighter_plugin}"
        )


def mg5_command(cfg: Config, *extra: object) -> list[object]:
    launcher = cfg.mg5_dir / "bin" / "mg5_aMC"
    return [cfg.mg5_python, launcher, *extra] if cfg.mg5_python else [launcher, *extra]


def madspin_command(cfg: Config, *extra: object) -> list[object]:
    launcher = cfg.mg5_dir / "MadSpin" / "madspin"
    return [cfg.mg5_python, launcher, *extra] if cfg.mg5_python else [launcher, *extra]


def stable_offset(name: str, modulus: int = 800_000_000) -> int:
    return int(hashlib.sha256(name.encode("utf-8")).hexdigest()[:12], 16) % modulus


def seed_for(name: str, cfg: Config) -> int:
    return 1 + (cfg.seed_base + stable_offset(name)) % 900_000_000


def requested_sample_names(cfg: Config) -> set[str] | None:
    if cfg.run_samples == "all":
        return None
    requested = {part.strip() for part in cfg.run_samples.split(",") if part.strip()}
    if not requested:
        die("--run-samples selected no names")
    unknown = requested - set(SAMPLES_BY_NAME) - set(RUN_SAMPLE_ALIASES)
    if unknown:
        die("unknown sample(s) in --run-samples: " + ",".join(sorted(unknown)))
    names: set[str] = set()
    for token in requested:
        names.update(RUN_SAMPLE_ALIASES.get(token, {token}))
    return names


def sample_in_set(sample: SampleSpec, sample_set: str) -> bool:
    if sample_set == "all":
        return True
    if sample_set == "irreducible":
        return sample.category in {"signal", "irreducible"}
    if sample_set == "reducible":
        return sample.category == "reducible"
    raise ValueError(f"unknown sample set: {sample_set}")


def biased_variant(sample: SampleSpec) -> SampleSpec:
    if not sample.heavy_flavour_candidate or sample.is_hf_biased:
        raise ValueError(f"sample is not eligible for heavy-flavour bias: {sample.name}")
    return replace(
        sample,
        name=f"{sample.name}_hfbiased",
        label=f"{sample.label} (semileptonic HF stratum)",
        normalization_rule=(
            sample.normalization_rule
            + "_times_herwig_selected_decay_branching_fraction_weights"
        ),
        bias_of=sample.name,
    )


def selected_samples(cfg: Config) -> list[SampleSpec]:
    requested = requested_sample_names(cfg)
    nominal = [
        sample
        for sample in SAMPLES
        if sample_in_set(sample, cfg.sample_set)
        and (requested is None or sample.name in requested)
    ]
    selected = list(nominal)
    if cfg.heavy_flavour_bias == "semileptonic":
        selected.extend(
            biased_variant(sample)
            for sample in nominal
            if sample.heavy_flavour_candidate
        )
    return selected


def resolve_external_lhe(sample: SampleSpec, cfg: Config) -> Path:
    for key in (sample.name, sample.bias_of, sample.production_group_id):
        if key and key in cfg.external_lhe:
            return cfg.external_lhe[key]
    die(
        "external LHE backend needs "
        f"--external-lhe {sample.name}=PATH or "
        f"--external-lhe {sample.production_group_id}=PATH"
    )


def external_cross_section_override(
    sample: SampleSpec,
    cfg: Config,
    *,
    stable_signal_production: bool = False,
) -> tuple[float | None, str | None]:
    if stable_signal_production:
        keys = (sample.production_group_id,)
    else:
        keys: tuple[str | None, ...] = (sample.name, sample.bias_of)
        if not sample.is_signal:
            keys = (*keys, sample.production_group_id)
    for key in keys:
        if key and key in cfg.external_cross_sections_pb:
            return cfg.external_cross_sections_pb[key], key
    return None, None


def group_event_count(
    production_group_id: str,
    samples: Sequence[SampleSpec],
    cfg: Config,
) -> int:
    counts = [
        cfg.events_for(sample)
        for sample in samples
        if sample.production_group_id == production_group_id
    ]
    if not counts:
        raise ValueError(f"empty production group: {production_group_id}")
    return max(counts)


def require_inputs(cfg: Config, samples: Sequence[SampleSpec]) -> None:
    if not cfg.template_in.exists():
        die(f"missing Herwig template: {cfg.template_in}")
    needs_mg5 = cfg.source_backend == "mg5" or any(sample.is_signal for sample in samples)
    if needs_mg5 and not (cfg.mg5_dir / "bin" / "mg5_aMC").exists():
        die(f"missing MG5 executable: {cfg.mg5_dir / 'bin' / 'mg5_aMC'}")
    if cfg.mg5_python and shutil.which(cfg.mg5_python) is None:
        die(f"MG5 Python executable not found: {cfg.mg5_python}")
    if cfg.herwig_env and not cfg.herwig_env.exists():
        die(f"missing Herwig activation script: {cfg.herwig_env}")
    if not cfg.herwig_env and not cfg.herwig_module and shutil.which(cfg.herwig) is None:
        die("Herwig is unavailable; set --herwig-env, --herwig-module, or --herwig")
    if cfg.collier_library and not cfg.collier_library.exists():
        die(f"missing COLLIER library: {cfg.collier_library}")
    if any(
        sample.prompt_decay_profile or sample.is_hf_biased for sample in samples
    ) and not cfg.decay_reweighter_source.exists():
        die(
            "missing LHE branching-ratio reweighter source: "
            f"{cfg.decay_reweighter_source}"
        )
    if cfg.source_backend == "external-lhe":
        for sample in samples:
            path = resolve_external_lhe(sample, cfg)
            if not cfg.dry_run and not path.exists():
                die(f"external LHE does not exist: {path}")
    if any(sample.is_hf_biased for sample in samples):
        if cfg.herwig_b_decays is None:
            die(
                "--heavy-flavour-bias semileptonic requires "
                "--herwig-b-decays PATH"
            )
        if not cfg.herwig_b_decays.exists():
            die(f"missing HerwigBDecays.in: {cfg.herwig_b_decays}")


def write_text(path: Path, text: str, cfg: Config, description: str) -> None:
    if cfg.dry_run:
        print(f"+ write {description} {path}")
        print(text, end="" if text.endswith("\n") else "\n")
    else:
        path.write_text(text)


def process_card_text(sample: SampleSpec, output_dir: Path) -> str:
    return f"""import model {sample.model}
define p = g u c d s b u~ c~ d~ s~ b~
define j = g u c d s b u~ c~ d~ s~ b~
define l+ = e+ mu+
define l- = e- mu-
define vl = ve vm
define vl~ = ve~ vm~
generate {sample.process}
output {output_dir} -f
"""


def write_mg5_process_card(
    card: Path,
    sample: SampleSpec,
    output_dir: Path,
    cfg: Config,
) -> None:
    write_text(card, process_card_text(sample, output_dir), cfg, "MG5 process card")


def madspin_card_text(
    input_lhe: Path,
    sample: SampleSpec,
    seed: int = 1,
) -> str:
    if not sample.madspin_decay:
        raise ValueError(f"sample has no MadSpin decay: {sample.name}")
    return f"""set spinmode none
import {input_lhe}
set seed {seed}
decay {sample.madspin_decay}
launch
"""


def write_madspin_card(
    card: Path,
    input_lhe: Path,
    sample: SampleSpec,
    cfg: Config,
) -> None:
    text = madspin_card_text(
        input_lhe,
        sample,
        seed_for(f"{sample.name}:madspin", cfg),
    )
    write_text(card, text, cfg, "MadSpin card")


def patch_run_card(
    card: Path,
    *,
    sample: SampleSpec,
    nevents: int,
    seed: int,
    cfg: Config,
) -> None:
    updates = {
        "nevents": str(nevents),
        "ebeam1": str(cfg.ebeam),
        "ebeam2": str(cfg.ebeam),
        "iseed": str(seed),
        "pdlabel": "lhapdf",
        "pdlabel1": "lhapdf",
        "pdlabel2": "lhapdf",
        "lhaid": str(LO_PDF_LHAID),
    }
    has_me_leptons = any(
        token in sample.process for token in ("e+", "e-", "mu+", "mu-")
    )
    if has_me_leptons:
        updates.update(
            {
                "ptl": str(GENERATION_CUTS["ptl_min_gev"]),
                "etal": str(GENERATION_CUTS["etal_max"]),
                "mmll": str(GENERATION_CUTS["mmll_min_gev"]),
                "drll": str(GENERATION_CUTS["drll_min"]),
            }
        )
    if cfg.dry_run:
        print(f"+ patch run card {card}: " + " ".join(f"{k}={v}" for k, v in updates.items()))
        return
    pattern = re.compile(r"^(\s*)([^=]+?)(\s*=\s*)([A-Za-z0-9_]+)(.*)$")
    output: list[str] = []
    found: set[str] = set()
    for line in card.read_text().splitlines(keepends=True):
        match = pattern.match(line)
        if match:
            indent, _value, separator, key, tail = match.groups()
            if key in updates:
                newline = "\n" if line.endswith("\n") else ""
                line = f"{indent}{updates[key]:>12}{separator}{key}{tail.rstrip()}{newline}"
                found.add(key)
        output.append(line)
    card.write_text("".join(output))
    required = {"nevents", "ebeam1", "ebeam2", "iseed", "pdlabel", "lhaid"}
    if has_me_leptons:
        required.update({"ptl", "etal", "mmll"})
    missing = required - found
    if missing:
        die(f"MG5 run card is missing required settings: {','.join(sorted(missing))}")


def patch_param_card_higgs_mass(card: Path, cfg: Config) -> float | None:
    if cfg.dry_run:
        print(f"+ set MASS 25 = {cfg.higgs_mass} in {card}")
        return None
    text = card.read_text()
    pattern = re.compile(
        r"(?im)^(\s*25\s+)([-+0-9.eEdD]+)(\s*(?:#.*)?)$"
    )
    matches = list(pattern.finditer(text))
    if not matches:
        die(f"could not locate Higgs MASS 25 in {card}")
    updated = pattern.sub(
        lambda match: f"{match.group(1)}{cfg.higgs_mass:.9e}{match.group(3)}",
        text,
        count=1,
    )
    card.write_text(updated)
    width_match = re.search(
        r"(?im)^\s*DECAY\s+25\s+([-+0-9.eEdD]+)",
        updated,
    )
    return (
        float(width_match.group(1).replace("D", "E").replace("d", "e"))
        if width_match
        else None
    )


def rewrite_run_card_include(text: str) -> str:
    updates = {
        "PDLABEL": "'lhapdf'",
        "PDSUBLABEL(1)": "'lhapdf'",
        "PDSUBLABEL(2)": "'lhapdf'",
        "LHAID": str(LO_PDF_LHAID),
    }
    output: list[str] = []
    for line in text.splitlines(keepends=True):
        for key, value in updates.items():
            pattern = re.compile(rf"^(\s*){re.escape(key)}(\s*=\s*).*$")
            bare = line.rstrip("\n")
            if pattern.match(bare):
                newline = "\n" if line.endswith("\n") else ""
                line = pattern.sub(
                    lambda match: f"{match.group(1)}{key}{match.group(2)}{value}",
                    bare,
                ) + newline
                break
        output.append(line)
    return "".join(output)


def generate_and_patch_run_card_include(process_dir: Path, cfg: Config) -> None:
    source_dir = process_dir / "Source"
    include = source_dir / "run_card.inc"
    if cfg.dry_run:
        print(f"+ make -C {source_dir} run_card.inc")
        print(f"+ enforce LHAPDF settings in {include}")
        return
    run_runtime_command(["make", "-C", source_dir, "run_card.inc"], cfg)
    if not include.exists():
        die(f"MG5 did not generate {include}")
    include.write_text(rewrite_run_card_include(include.read_text()))
    setrun = source_dir / "setrun.o"
    if setrun.exists():
        setrun.unlink()


def patch_mg5_hidden_pdf_defaults(process_dir: Path, cfg: Config) -> None:
    """Patch MG5 3.5 hidden beam-PDF defaults used by generated drivers."""

    banner = process_dir / "bin" / "internal" / "banner.py"
    if cfg.dry_run:
        print(f"+ patch hidden MG5 PDF defaults in {banner}")
        return
    if not banner.exists():
        die(f"missing generated MG5 banner helper: {banner}")
    replacements = {
        'self.add_param("pdlabel", "nn23lo1", hidden=True, allowed=valid_pdf)': (
            'self.add_param("pdlabel", "lhapdf", hidden=True, allowed=valid_pdf)'
        ),
        'self.add_param("pdlabel1", "nn23lo1", hidden=True, allowed=valid_pdf, fortran_name="pdsublabel(1)")': (
            'self.add_param("pdlabel1", "lhapdf", hidden=True, allowed=valid_pdf, fortran_name="pdsublabel(1)")'
        ),
        'self.add_param("pdlabel2", "nn23lo1", hidden=True, allowed=valid_pdf, fortran_name="pdsublabel(2)")': (
            'self.add_param("pdlabel2", "lhapdf", hidden=True, allowed=valid_pdf, fortran_name="pdsublabel(2)")'
        ),
        'self.add_param("lhaid", 230000, hidden=True)': (
            f'self.add_param("lhaid", {LO_PDF_LHAID}, hidden=True)'
        ),
    }
    text = banner.read_text()
    updated = text
    for old, new in replacements.items():
        if old in updated:
            updated = updated.replace(old, new, 1)
        elif new not in updated:
            die(f"could not patch MG5 hidden PDF default in {banner}: {old}")
    if updated != text:
        banner.write_text(updated)


def patch_mg5_pdf_label_sync(process_dir: Path, cfg: Config) -> None:
    """Keep generated two-beam PDF labels synchronized in MG5 3.5.x."""

    interface = process_dir / "bin" / "internal" / "madevent_interface.py"
    if cfg.dry_run:
        print(f"+ patch MG5 beam PDF synchronization in {interface}")
        return
    if not interface.exists():
        die(f"missing generated MG5 interface: {interface}")
    updated = interface.read_text()
    if MG5_PDF_LABEL_SYNC_MARKER not in updated:
        needle = "        # set  lhapdf.\n"
        if needle not in updated:
            die(f"could not locate MG5 LHAPDF block in {interface}")
        block = (
            f"        {MG5_PDF_LABEL_SYNC_MARKER}\n"
            "        if self.run_card['pdlabel1'] == \"lhapdf\" or "
            "self.run_card['pdlabel2'] == \"lhapdf\":\n"
            "            self.run_card['pdlabel'] = \"lhapdf\"\n"
        )
        updated = updated.replace(needle, needle + block, 1)
    if MG5_POST_TREATCARDS_PDF_PATCH_MARKER not in updated:
        needle = "        self.do_treatcards('')\n"
        if needle not in updated:
            die(f"could not locate MG5 treatcards call in {interface}")
        block = (
            f"        {MG5_POST_TREATCARDS_PDF_PATCH_MARKER}\n"
            "        run_card_inc = pjoin(self.me_dir, 'Source', 'run_card.inc')\n"
            "        if os.path.exists(run_card_inc):\n"
            "            with open(run_card_inc) as run_card_stream:\n"
            "                run_card_lines = run_card_stream.readlines()\n"
            "            patched_run_card_lines = []\n"
            "            for run_card_line in run_card_lines:\n"
            "                upper_line = run_card_line.lstrip().upper()\n"
            "                if upper_line.startswith('PDLABEL'):\n"
            "                    run_card_line = \"      PDLABEL = 'lhapdf'\\n\"\n"
            "                elif upper_line.startswith('PDSUBLABEL(1)'):\n"
            "                    run_card_line = \"      PDSUBLABEL(1) = 'lhapdf'\\n\"\n"
            "                elif upper_line.startswith('PDSUBLABEL(2)'):\n"
            "                    run_card_line = \"      PDSUBLABEL(2) = 'lhapdf'\\n\"\n"
            "                elif upper_line.startswith('LHAID'):\n"
            f"                    run_card_line = \"      LHAID = {LO_PDF_LHAID}\\n\"\n"
            "                patched_run_card_lines.append(run_card_line)\n"
            "            with open(run_card_inc, 'w') as run_card_stream:\n"
            "                run_card_stream.writelines(patched_run_card_lines)\n"
            "            setrun_object = pjoin(self.me_dir, 'Source', 'setrun.o')\n"
            "            if os.path.exists(setrun_object):\n"
            "                os.remove(setrun_object)\n"
        )
        updated = updated.replace(needle, needle + block, 1)
    interface.write_text(updated)


def shared_library_glob(system: str) -> str:
    return "*.dylib" if system == "Darwin" else "*.so*"


def copy_library(source: Path, destination: Path) -> None:
    temporary = destination.with_name(f"{destination.name}.tmp.{os.getpid()}")
    shutil.copy2(source, temporary, follow_symlinks=True)
    temporary.replace(destination)


def prepare_collier_runtime(runtime_dir: Path, cfg: Config) -> None:
    if cfg.collier_library is None:
        return
    ensure_dir(runtime_dir, cfg)
    pattern = shared_library_glob(platform.system())
    if cfg.dry_run:
        print(f"+ copy {cfg.collier_library.parent / pattern} to {runtime_dir}")
        return
    for library in sorted(cfg.collier_library.parent.glob(pattern)):
        if library.is_file() or library.is_symlink():
            copy_library(library, runtime_dir / library.name)
    alias = runtime_dir / (
        "libcollier.dylib" if platform.system() == "Darwin" else "libcollier.so"
    )
    if not alias.exists():
        copy_library(cfg.collier_library, alias)
    if platform.system() == "Darwin":
        patch_openloops_dylibs(runtime_dir)


def patch_openloops_dylibs(runtime_dir: Path) -> None:
    """Make copied OpenLoops libraries resolve siblings via @loader_path."""

    if shutil.which("install_name_tool") is None or shutil.which("otool") is None:
        return
    dylibs = sorted(runtime_dir.glob("*.dylib"))
    for dylib in dylibs:
        subprocess.run(
            ["install_name_tool", "-id", f"@rpath/{dylib.name}", str(dylib)],
            stderr=subprocess.DEVNULL,
            check=False,
        )
    for dylib in dylibs:
        result = subprocess.run(
            ["otool", "-L", str(dylib)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        for raw_line in result.stdout.splitlines()[1:]:
            dependency = raw_line.strip().split(maxsplit=1)[0]
            if not dependency.startswith("lib/") or not dependency.endswith(".dylib"):
                continue
            sibling = runtime_dir / Path(dependency).name
            if sibling.exists():
                subprocess.run(
                    [
                        "install_name_tool",
                        "-change",
                        dependency,
                        f"@loader_path/{sibling.name}",
                        str(dylib),
                    ],
                    stderr=subprocess.DEVNULL,
                    check=False,
                )


def prepare_loop_runtime(process_dir: Path, runtime_dir: Path | None, cfg: Config) -> None:
    if runtime_dir is None:
        return
    make_opts = process_dir / "Source" / "make_opts"
    marker = "# four-lepton campaign runtime library path"
    block = (
        f"\n{marker}\n"
        f"LDFLAGS += -Wl,-rpath,{runtime_dir}\n"
        f"RPATH_LIBS += -Wl,-rpath,{runtime_dir}\n"
    )
    if cfg.dry_run:
        print(f"+ add runtime rpath {runtime_dir} to {make_opts}")
        return
    if not make_opts.exists():
        die(f"missing MG5 make options: {make_opts}")
    text = make_opts.read_text()
    if marker not in text:
        make_opts.write_text(text.rstrip() + block)
    for subprocess_dir in sorted((process_dir / "SubProcesses").glob("P*")):
        if not subprocess_dir.is_dir():
            continue
        link = subprocess_dir / "lib"
        if link.is_symlink() and link.resolve() != runtime_dir.resolve():
            link.unlink()
        if not link.exists():
            link.symlink_to(runtime_dir, target_is_directory=True)


def restore_mg5_symmetry_factors(process_dir: Path, cfg: Config) -> None:
    subprocess_root = process_dir / "SubProcesses"
    if cfg.dry_run:
        print(f"+ restore missing symfact.dat below {subprocess_root}")
        return
    for subprocess_dir in sorted(subprocess_root.glob("P*")):
        symfact = subprocess_dir / "symfact.dat"
        original = subprocess_dir / "symfact_orig.dat"
        if not symfact.exists() and original.exists():
            shutil.copy2(original, symfact)


def generation_environment(runtime_dir: Path | None) -> dict[str, str] | None:
    if runtime_dir is None:
        return None
    env = dict(os.environ)
    keys = (
        ("DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH")
        if platform.system() == "Darwin"
        else ("LD_LIBRARY_PATH",)
    )
    for key in keys:
        current = env.get(key)
        env[key] = str(runtime_dir) if not current else f"{runtime_dir}{os.pathsep}{current}"
    return env


def find_lhe_file(process_dir: Path) -> Path:
    candidates = sorted((process_dir / "Events").glob("*/unweighted_events.lhe*"))
    if not candidates:
        die(f"could not locate MG5 LHE output below {process_dir / 'Events'}")
    undecayed = [path for path in candidates if "decayed" not in str(path)]
    return undecayed[-1] if undecayed else candidates[-1]


def open_maybe_gzip(path: Path):
    return gzip.open(path, "rt", errors="replace") if path.suffix == ".gz" else path.open(
        "rt", errors="replace"
    )


def read_lhe_prefix(path: Path, closing_tag: str = "</init>") -> str:
    lines: list[str] = []
    with open_maybe_gzip(path) as stream:
        for line in stream:
            lines.append(line)
            if closing_tag in line:
                break
            if len(lines) > 200_000:
                raise ValueError(f"could not find {closing_tag} near the start of {path}")
    return "".join(lines)


def parse_lhe_init_records(path: Path) -> list[str]:
    """Return the init header and the declared subprocess records.

    MG5 3.5 may append XML metadata such as ``<generator ...>`` inside the
    ``<init>`` element.  Those elements are not Les Houches numeric records
    and must not be interpreted as subprocess cross sections.
    """

    prefix = read_lhe_prefix(path)
    init_match = re.search(r"<init\b[^>]*>(.*?)</init>", prefix, flags=re.DOTALL)
    if not init_match:
        raise ValueError(f"LHE file has no complete <init> block: {path}")
    records = [
        line.strip()
        for line in init_match.group(1).splitlines()
        if line.strip()
        and not line.lstrip().startswith("#")
        and not line.lstrip().startswith("<")
    ]
    if not records:
        raise ValueError(f"LHE init block is empty: {path}")
    header_fields = records[0].split()
    if len(header_fields) < 10:
        raise ValueError(f"malformed LHE init header in {path}: {records[0]}")
    try:
        n_subprocesses = int(header_fields[9])
    except ValueError as exc:
        raise ValueError(
            f"invalid NPRUP in LHE init header in {path}: {records[0]}"
        ) from exc
    if n_subprocesses < 0:
        raise ValueError(f"negative NPRUP={n_subprocesses} in {path}")
    expected_records = 1 + n_subprocesses
    if len(records) < expected_records:
        raise ValueError(
            f"LHE init block in {path} declares {n_subprocesses} subprocesses "
            f"but contains only {len(records) - 1}"
        )
    return records[:expected_records]


def parse_lhe_init_cross_section_pb(path: Path) -> float:
    """Return the sum of XSECUP entries in an LHE init block."""

    records = parse_lhe_init_records(path)
    if len(records) < 2:
        raise ValueError(f"LHE init block has no subprocess records: {path}")
    total = 0.0
    for record in records[1:]:
        fields = record.split()
        if len(fields) < 4:
            raise ValueError(f"malformed LHE subprocess record in {path}: {record}")
        total += float(fields[0].replace("D", "E").replace("d", "e"))
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError(f"non-positive LHE cross section in {path}: {total}")
    return total


def parse_lhe_init_metadata(path: Path) -> dict[str, object]:
    records = parse_lhe_init_records(path)
    fields = records[0].split()
    return {
        "beam_ids": [int(fields[0]), int(fields[1])],
        "beam_energies_gev": [
            float(fields[2].replace("D", "E").replace("d", "e")),
            float(fields[3].replace("D", "E").replace("d", "e")),
        ],
        "pdf_group_ids": [int(fields[4]), int(fields[5])],
        "pdf_set_ids": [int(fields[6]), int(fields[7])],
        "idwtup": int(fields[8]),
        "n_subprocesses": int(fields[9]),
    }


def validate_lhe_collider(path: Path, cfg: Config) -> dict[str, object]:
    metadata = parse_lhe_init_metadata(path)
    if metadata["beam_ids"] != [2212, 2212]:
        raise ValueError(
            f"{path} has LHE beam IDs {metadata['beam_ids']}, expected [2212, 2212]"
        )
    energies = metadata["beam_energies_gev"]
    assert isinstance(energies, list)
    if any(
        not math.isclose(float(energy), cfg.ebeam, rel_tol=0.0, abs_tol=1.0e-6)
        for energy in energies
    ):
        raise ValueError(
            f"{path} has LHE beam energies {energies}, expected {cfg.ebeam} GeV"
        )
    return metadata


def parse_lhe_higgs_parameters(path: Path) -> tuple[float | None, float | None]:
    prefix = read_lhe_prefix(path)
    mass_match = re.search(
        r"(?im)^\s*25\s+([-+0-9.eEdD]+)\s*(?:#.*)?$",
        prefix,
    )
    width_match = re.search(
        r"(?im)^\s*DECAY\s+25\s+([-+0-9.eEdD]+)",
        prefix,
    )

    def number(match: re.Match[str] | None) -> float | None:
        return (
            float(match.group(1).replace("D", "E").replace("d", "e"))
            if match
            else None
        )

    return number(mass_match), number(width_match)


def parse_lhe_optional_weight_definitions(
    path: Path,
) -> list[dict[str, str]]:
    """Return ordered LHE ``<initrwgt>`` definitions."""

    lines: list[str] = []
    with open_maybe_gzip(path) as stream:
        for line in stream:
            lines.append(line)
            if re.search(r"<init\b", line, flags=re.IGNORECASE):
                break
            if len(lines) > 200_000:
                raise ValueError(f"could not find <init> near the start of {path}")
    prefix = "".join(lines)
    match = re.search(
        r"<initrwgt\b[^>]*>(.*?)</initrwgt>",
        prefix,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if match is None:
        return []
    definitions: list[dict[str, str]] = []
    seen_ids: set[str] = set()
    for weight_match in re.finditer(
        r"<weight\b([^>]*)>(.*?)</weight>",
        match.group(1),
        flags=re.IGNORECASE | re.DOTALL,
    ):
        attributes = " ".join(weight_match.group(1).split())
        id_match = re.search(
            r"\bid\s*=\s*(['\"])(.*?)\1",
            weight_match.group(1),
            flags=re.IGNORECASE | re.DOTALL,
        )
        if id_match is None:
            raise ValueError(f"LHE optional weight without an id in {path}")
        weight_id = id_match.group(2).strip()
        if not weight_id or weight_id in seen_ids:
            raise ValueError(
                f"empty or duplicate LHE optional-weight id {weight_id!r} "
                f"in {path}"
            )
        seen_ids.add(weight_id)
        description = " ".join(
            re.sub(r"<[^>]+>", " ", weight_match.group(2)).split()
        )
        label = (
            f"lhe_id={weight_id}:{description}"
            if description
            else f"lhe_id={weight_id}"
        )
        definitions.append(
            {
                "id": weight_id,
                "label": label,
                "description": description,
                "attributes": attributes,
            }
        )
    return definitions


def parse_lhe_weight_metadata(path: Path) -> tuple[int, dict[int, float]]:
    """Return ``IDWTUP`` and the banner ``LPRUP -> XMAXUP`` map."""

    records = parse_lhe_init_records(path)
    if len(records) < 2:
        raise ValueError(f"LHE init block has no subprocess records: {path}")
    header_fields = records[0].split()
    idwtup = int(header_fields[8])
    result: dict[int, float] = {}
    for record in records[1:]:
        fields = record.split()
        if len(fields) < 4:
            raise ValueError(f"malformed LHE subprocess record in {path}: {record}")
        maximum = float(fields[2].replace("D", "E").replace("d", "e"))
        process_id = int(fields[3])
        if not math.isfinite(maximum) or maximum == 0.0:
            raise ValueError(
                f"invalid XMAXUP={maximum} for LPRUP={process_id} in {path}"
            )
        result[process_id] = maximum
    return idwtup, result


def parse_lhe_process_max_weights(path: Path) -> dict[int, float]:
    """Return the banner ``LPRUP -> XMAXUP`` map."""

    return parse_lhe_weight_metadata(path)[1]


def predecay_lhe_weight_summary(
    path: Path, event_limit: int
) -> dict[str, int | float | str]:
    """Summarize the weights ThePEG assigns before forced decays.

    ``LesHouchesEventHandler`` with ``VarNegWeight`` feeds each event to
    HwSim with ``XWGTUP/XMAXUP`` (for the event's ``IDPRUP``), not merely
    the sign of ``XWGTUP``.  Retaining the magnitude is required for
    variable-weight NLO LHE input.
    """

    if event_limit <= 0:
        raise ValueError("event_limit must be positive")
    idwtup, banner_maximum_by_process = parse_lhe_weight_metadata(path)
    raw_events: list[tuple[int, float]] = []
    scanned_maximum_by_process: dict[int, float] = defaultdict(float)
    awaiting_header = False
    with open_maybe_gzip(path) as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if line.startswith("<event"):
                awaiting_header = True
                continue
            if not awaiting_header or not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"malformed LHE event header in {path}: {line}")
            process_id = int(fields[1])
            xwgtup = float(fields[2].replace("D", "E").replace("d", "e"))
            if process_id not in banner_maximum_by_process:
                raise ValueError(
                    f"LHE event IDPRUP={process_id} has no init record in {path}"
                )
            if not math.isfinite(xwgtup):
                raise ValueError(f"non-finite XWGTUP in {path}")
            scanned_maximum_by_process[process_id] = max(
                scanned_maximum_by_process[process_id],
                abs(xwgtup),
            )
            if len(raw_events) < event_limit:
                raw_events.append((process_id, xwgtup))
            awaiting_header = False
    if len(raw_events) != event_limit:
        raise ValueError(
            f"{path} contains only {len(raw_events)} events; {event_limit} were requested"
        )
    # LesHouchesReader scans the complete file and replaces XMAXUP whenever
    # abs(IDWTUP) != 1 (its default MaxScan is -1).  Reproduce that behavior
    # before normalizing the subset that Herwig will actually run.
    if abs(idwtup) != 1:
        maximum_by_process = scanned_maximum_by_process
        maximum_source = "full_file_max_abs_xwgtup"
    else:
        maximum_by_process = banner_maximum_by_process
        maximum_source = "lhe_init_xmaxup"
    weights: list[float] = []
    for process_id, xwgtup in raw_events:
        maximum = maximum_by_process.get(process_id, 0.0)
        if not math.isfinite(maximum) or maximum <= 0.0:
            raise ValueError(
                f"invalid effective XMAXUP={maximum} for IDPRUP={process_id} in {path}"
            )
        weights.append(xwgtup / maximum)
    npositive = sum(weight > 0.0 for weight in weights)
    nnegative = sum(weight < 0.0 for weight in weights)
    nzero = len(weights) - npositive - nnegative
    return {
        "denominator_kind": "predecay_lhe_handler_weight_sum",
        "denominator_formula": "XWGTUP / effective_XMAXUP[IDPRUP]",
        "lhe_idwtup": idwtup,
        "effective_xmaxup_source": maximum_source,
        "denominator_sumw": float(sum(weights)),
        "denominator_sumabsw": float(sum(abs(weight) for weight in weights)),
        "denominator_sumw2": float(sum(weight * weight for weight in weights)),
        "denominator_event_count": len(weights),
        "denominator_npositive": npositive,
        "denominator_nnegative": nnegative,
        "denominator_nzero": nzero,
    }


# Backward-compatible import name for callers of the first implementation.
predecay_lhe_sign_summary = predecay_lhe_weight_summary


def validate_prompt_decay_parents(
    path: Path,
    sample: SampleSpec,
    event_limit: int,
) -> dict[str, object]:
    """Require stable hard parents before applying Herwig decay selectors.

    This prevents an externally pre-decayed ttbar/ttZ LHE from receiving the
    prompt-decay branching-ratio correction a second time.
    """

    if not sample.prompt_decay_profile:
        return {"required_status1_pdgs": [], "events_checked": 0}
    if event_limit <= 0:
        raise ValueError("event_limit must be positive")
    required = (
        (6, -6)
        if sample.prompt_decay_profile == "w_emu"
        else (6, -6, 23)
    )
    if sample.prompt_decay_profile not in PROMPT_DECAY_COMMANDS:
        raise ValueError(
            f"unknown prompt decay profile: {sample.prompt_decay_profile}"
        )

    events_checked = 0
    awaiting_header = False
    particles_remaining = 0
    status_one_pdgs: list[int] = []
    with open_maybe_gzip(path) as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if line.startswith("<event"):
                awaiting_header = True
                particles_remaining = 0
                status_one_pdgs = []
                continue
            if awaiting_header:
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) < 1:
                    raise ValueError(f"malformed LHE event header in {path}")
                particles_remaining = int(fields[0])
                if particles_remaining <= 0:
                    raise ValueError(f"non-positive NUP in {path}: {line}")
                awaiting_header = False
                continue
            if particles_remaining <= 0:
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"malformed LHE particle record in {path}: {line}")
            pdg_id = int(fields[0])
            status = int(fields[1])
            if status == 1:
                status_one_pdgs.append(pdg_id)
            particles_remaining -= 1
            if particles_remaining:
                continue

            events_checked += 1
            missing = [pdg_id for pdg_id in required if pdg_id not in status_one_pdgs]
            if missing:
                raise ValueError(
                    f"{path} event {events_checked} lacks undecayed status-1 "
                    f"parent(s) {missing} required by prompt profile "
                    f"{sample.prompt_decay_profile}; do not pass pre-decayed "
                    "LHE events to Herwig forced-decay strata"
                )
            if events_checked == event_limit:
                break

    if events_checked != event_limit:
        raise ValueError(
            f"{path} contains only {events_checked} complete events; "
            f"{event_limit} were requested"
        )
    return {
        "required_status1_pdgs": list(required),
        "events_checked": events_checked,
    }


def validate_stable_higgs_lhe(
    path: Path,
    cfg: Config,
    event_limit: int,
) -> tuple[float, float | None]:
    """Require exactly one undecayed, on-mass-shell Higgs in each used event."""

    if event_limit <= 0:
        raise ValueError("event_limit must be positive")
    banner_mass, width = parse_lhe_higgs_parameters(path)
    banner_mass_is_valid = banner_mass is not None and math.isclose(
        banner_mass,
        cfg.higgs_mass,
        rel_tol=0.0,
        abs_tol=1.0e-4,
    )

    events_checked = 0
    awaiting_header = False
    particles_remaining = 0
    stable_higgs_masses: list[float] = []
    first_mass: float | None = None
    with open_maybe_gzip(path) as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if line.startswith("<event"):
                awaiting_header = True
                particles_remaining = 0
                stable_higgs_masses = []
                continue
            if awaiting_header:
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if not fields:
                    raise ValueError(f"malformed LHE event header in {path}")
                particles_remaining = int(fields[0])
                if particles_remaining <= 0:
                    raise ValueError(f"non-positive NUP in {path}: {line}")
                awaiting_header = False
                continue
            if particles_remaining <= 0:
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 11:
                raise ValueError(
                    f"malformed LHE particle record in {path}: {line}"
                )
            pdg_id = int(fields[0])
            status = int(fields[1])
            if pdg_id == 25 and status == 1:
                particle_mass = float(
                    fields[10].replace("D", "E").replace("d", "e")
                )
                if not math.isfinite(particle_mass):
                    raise ValueError(f"non-finite status-1 Higgs mass in {path}")
                stable_higgs_masses.append(particle_mass)
            particles_remaining -= 1
            if particles_remaining:
                continue

            events_checked += 1
            if len(stable_higgs_masses) != 1:
                raise ValueError(
                    f"{path} event {events_checked} contains "
                    f"{len(stable_higgs_masses)} status-1 Higgs bosons; "
                    "exactly one undecayed Higgs is required"
                )
            event_mass = stable_higgs_masses[0]
            if not math.isclose(
                event_mass,
                cfg.higgs_mass,
                rel_tol=0.0,
                abs_tol=1.0e-3,
            ):
                raise ValueError(
                    f"{path} event {events_checked} has status-1 Higgs "
                    f"mass {event_mass} GeV, expected {cfg.higgs_mass} GeV"
                )
            if first_mass is None:
                first_mass = event_mass
            if events_checked == event_limit:
                break

    if events_checked != event_limit:
        raise ValueError(
            f"{path} contains only {events_checked} complete stable-H events; "
            f"{event_limit} were requested"
        )
    assert first_mass is not None
    # The event record is authoritative for a stable-H production LHE.  A
    # missing or stale banner MASS entry is repaired in the private MadSpin
    # staging copy; it must not make an otherwise valid external event file
    # unusable.
    return banner_mass if banner_mass_is_valid else first_mass, width


def signal_target_branching_fraction(sample: SampleSpec, cfg: Config) -> float:
    if not sample.is_signal:
        raise ValueError(f"not a signal sample: {sample.name}")
    return BR_H_TO_4L * cfg.signal_channel_fractions[sample.channel]


def normalization_from_lhe(
    sample: SampleSpec,
    native_lhe: Path,
    cfg: Config,
    *,
    stable_signal_lhe: Path | None = None,
) -> dict[str, object]:
    native_override, native_override_key = (
        external_cross_section_override(sample, cfg)
        if cfg.source_backend == "external-lhe"
        else (None, None)
    )
    native_xsec = (
        native_override
        if native_override is not None
        else parse_lhe_init_cross_section_pb(native_lhe)
    )
    native_source = (
        f"cli_external_cross_section:{native_override_key}"
        if native_override_key is not None
        else "lhe_init_xsecup"
    )
    if not sample.is_signal:
        return {
            "native_cross_section_pb": native_xsec,
            "native_cross_section_source": native_source,
            "stable_production_cross_section_pb": None,
            "target_branching_fraction": None,
            "target_cross_section_pb": native_xsec,
            "weight_scale": 1.0,
            "formula": "target=native; weight_scale=1",
        }
    if stable_signal_lhe is None:
        raise ValueError("signal normalization requires its stable-H LHE")
    stable_override, stable_override_key = (
        external_cross_section_override(
            sample,
            cfg,
            stable_signal_production=True,
        )
        if cfg.source_backend == "external-lhe"
        else (None, None)
    )
    stable_xsec = (
        stable_override
        if stable_override is not None
        else parse_lhe_init_cross_section_pb(stable_signal_lhe)
    )
    stable_source = (
        f"cli_external_cross_section:{stable_override_key}"
        if stable_override_key is not None
        else "lhe_init_xsecup"
    )
    target_br = signal_target_branching_fraction(sample, cfg)
    target_xsec = stable_xsec * target_br
    weight_scale = target_xsec / native_xsec
    if not math.isfinite(weight_scale) or weight_scale <= 0.0:
        raise ValueError(f"invalid signal normalization weight: {weight_scale}")
    return {
        "native_cross_section_pb": native_xsec,
        "native_cross_section_source": native_source,
        "stable_production_cross_section_pb": stable_xsec,
        "stable_production_cross_section_source": stable_source,
        "target_branching_fraction": target_br,
        "target_cross_section_pb": target_xsec,
        "weight_scale": weight_scale,
        "formula": "stable_ggh_xsec * external_target_BR_i / madspin_decayed_xsec",
    }


HEAVY_MESON_PARENTS = ("B0", "Bbar0", "B+", "B-", "B_s0", "B_sbar0")
HEAVY_FLAVOUR_SEMILEPTONIC_BIAS_FACTOR = 4.0

PROMPT_DECAY_COMMANDS: dict[str, tuple[str, ...]] = {
    "w_emu": (
        "do /Herwig/Particles/W-:SelectDecayModes "
        "/Herwig/Particles/W-/W-->nu_mubar,mu-; "
        "/Herwig/Particles/W-/W-->nu_ebar,e-;",
        "do /Herwig/Particles/W+:SelectDecayModes "
        "/Herwig/Particles/W+/W+->nu_mu,mu+; "
        "/Herwig/Particles/W+/W+->nu_e,e+;",
    ),
}

_Z_MODE = {
    "ee": "/Herwig/Particles/Z0/Z0->e-,e+;",
    "mumu": "/Herwig/Particles/Z0/Z0->mu-,mu+;",
}
_WPLUS_MODE = {
    "e": "/Herwig/Particles/W+/W+->nu_e,e+;",
    "mu": "/Herwig/Particles/W+/W+->nu_mu,mu+;",
}
_WMINUS_MODE = {
    "e": "/Herwig/Particles/W-/W-->nu_ebar,e-;",
    "mu": "/Herwig/Particles/W-/W-->nu_mubar,mu-;",
}
for _z_flavour in ("ee", "mumu"):
    for _wp_flavour in ("e", "mu"):
        for _wm_flavour in ("e", "mu"):
            PROMPT_DECAY_COMMANDS[
                f"z{_z_flavour}_wp{_wp_flavour}_wm{_wm_flavour}"
            ] = (
                "do /Herwig/Particles/Z0:SelectDecayModes "
                + _Z_MODE[_z_flavour],
                "do /Herwig/Particles/W+:SelectDecayModes "
                + _WPLUS_MODE[_wp_flavour],
                "do /Herwig/Particles/W-:SelectDecayModes "
                + _WMINUS_MODE[_wm_flavour],
            )


def parse_semileptonic_decay_modes(
    text: str,
) -> tuple[dict[str, list[tuple[str, float]]], dict[str, float]]:
    """Extract active direct e/mu modes and total B-meson branching sums.

    The selected modes receive a finite enhancement, while every other
    active mode remains enabled with unit bias.  Returning the complete
    parent sums lets the manifest state the exact proposal normalization
    without relying on the decay table being rounded to precisely one.
    """

    selected: dict[str, list[tuple[str, float]]] = {
        parent: [] for parent in HEAVY_MESON_PARENTS
    }
    totals = {parent: 0.0 for parent in HEAVY_MESON_PARENTS}
    pattern = re.compile(
        r"^\s*decaymode\s+([^;]+);\s*([-+0-9.eEdD]+)\s+([01])\b",
        flags=re.IGNORECASE,
    )
    for line in text.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        mode = match.group(1).strip()
        branching_fraction = float(
            match.group(2).replace("D", "E").replace("d", "e")
        )
        active = int(match.group(3))
        if branching_fraction <= 0.0 or active == 0:
            continue
        if "->" not in mode:
            continue
        parent, products = mode.split("->", 1)
        parent = parent.strip()
        if parent not in selected:
            continue
        totals[parent] += branching_fraction
        # Restrict the enhancement to simple one-hadron semileptonic modes.
        # ThePEG canonicalizes the repository tag to hadron, neutrino,
        # charged lepton and retains the trailing semicolon in the object
        # name.
        product_list = [product.strip() for product in products.split(",")]
        if len(product_list) != 3 or any("=" in product for product in product_list):
            continue
        leptons = [
            product
            for product in product_list
            if product in {"e+", "e-", "mu+", "mu-"}
        ]
        neutrinos = [
            product
            for product in product_list
            if product in {"nu_e", "nu_ebar", "nu_mu", "nu_mubar"}
        ]
        hadrons = [
            product
            for product in product_list
            if product not in set(leptons + neutrinos)
        ]
        if len(leptons) == len(neutrinos) == len(hadrons) == 1:
            selected[parent].append(
                (
                    f"{parent}->{hadrons[0]},{neutrinos[0]},{leptons[0]}",
                    branching_fraction,
                )
            )
    return selected, totals


def semileptonic_decay_commands(
    path: Path,
) -> tuple[list[str], dict[str, int], dict[str, dict[str, float | int]]]:
    modes, totals = parse_semileptonic_decay_modes(path.read_text())
    missing = [parent for parent, parent_modes in modes.items() if not parent_modes]
    if missing:
        raise ValueError(
            "no direct e/mu decay modes found for " + ",".join(sorted(missing))
        )
    commands: list[str] = []
    counts: dict[str, int] = {}
    proposal: dict[str, dict[str, float | int]] = {}
    for parent in HEAVY_MESON_PARENTS:
        parent_modes = modes[parent]
        counts[parent] = len(parent_modes)
        original_total = totals[parent]
        if not math.isfinite(original_total) or original_total <= 0.0:
            raise ValueError(f"invalid active branching-ratio sum for {parent}")
        selected_sum = sum(branching for _, branching in parent_modes)
        proposal_normalization = (
            original_total
            + (HEAVY_FLAVOUR_SEMILEPTONIC_BIAS_FACTOR - 1.0)
            * selected_sum
        )
        for mode, branching_fraction in parent_modes:
            biased_branching = (
                branching_fraction * HEAVY_FLAVOUR_SEMILEPTONIC_BIAS_FACTOR
            )
            if not (0.0 < biased_branching <= 1.0):
                raise ValueError(
                    f"biased branching ratio for {mode} is outside (0,1]: "
                    f"{biased_branching}"
                )
            commands.append(
                f"set /Herwig/Particles/{parent}/{mode};:BranchingRatio "
                f"{biased_branching:.17g}"
            )
        proposal[parent] = {
            "selected_mode_count": len(parent_modes),
            "original_active_branching_sum": original_total,
            "selected_branching_sum": selected_sum,
            "bias_factor": HEAVY_FLAVOUR_SEMILEPTONIC_BIAS_FACTOR,
            "proposal_normalization": proposal_normalization,
        }
    return commands, counts, proposal


def decay_reweighter_handlers(sample: SampleSpec) -> list[dict[str, str]]:
    handlers: list[dict[str, str]] = []
    if sample.prompt_decay_profile:
        handlers.append(
            {
                "class": "HiggsSSC::LHEBranchingRatioReweighter",
                "stage": "PostHadronizationHandlers",
                "scope": "prompt_abs_pdg_23_24",
            }
        )
    if sample.is_hf_biased:
        handlers.append(
            {
                "class": "HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter",
                "stage": "PostDecayHandlers",
                "scope": "terminal_ground_state_B_abs_pdg_511_521_531",
            }
        )
    return handlers


def decay_selection_block(
    sample: SampleSpec,
    cfg: Config,
    *,
    reweighter_library: Path | None = None,
) -> tuple[str, dict[str, object]]:
    prompt_commands: list[str] = []
    if sample.prompt_decay_profile:
        try:
            prompt_commands.extend(PROMPT_DECAY_COMMANDS[sample.prompt_decay_profile])
        except KeyError as error:
            raise ValueError(
                f"unknown prompt decay profile: {sample.prompt_decay_profile}"
            ) from error

    heavy_commands: list[str] = []
    mode_counts: dict[str, int] = {}
    heavy_proposal: dict[str, dict[str, float | int]] = {}
    if sample.is_hf_biased:
        if cfg.herwig_b_decays is None:
            raise ValueError("semileptonic heavy-flavour bias needs HerwigBDecays.in")
        heavy_commands, mode_counts, heavy_proposal = (
            semileptonic_decay_commands(cfg.herwig_b_decays)
        )

    all_commands = [*prompt_commands, *heavy_commands]
    if not all_commands:
        return (
            "# No forced Herwig decay selection for this sample.",
            {
                "prompt_selection": sample.prompt_decay_profile or None,
                "branching_ratio_reweighter": False,
                "branching_ratio_reweighter_implementations": [],
                "handlers": [],
                "decay_bias": {
                    "enabled": False,
                    "scheme": "unbiased",
                    "closure_validated": None,
                },
            },
        )

    # Herwig::BranchingRatioReweighter only updates StandardEventHandler
    # events.  LHE showering uses LesHouchesEventHandler, so the campaign
    # loads LHE-compatible handlers.  Prompt W/Z factors are evaluated after
    # the cascade, before weak hadron decays can create virtual W/Z history
    # lines.  Selected ground-state B factors are evaluated separately after
    # recursive decays, when B mesons from excited-hadron decays are present.
    handler_lines: list[str] = []
    handlers = decay_reweighter_handlers(sample)
    implementations = [handler["class"] for handler in handlers]
    if prompt_commands:
        handler_lines.extend(
            [
                "create HiggsSSC::LHEBranchingRatioReweighter "
                "/Herwig/Generators/PromptBRReweighter",
                "insert /Herwig/Generators/theGenerator:EventHandler:"
                "PostHadronizationHandlers 0 "
                "/Herwig/Generators/PromptBRReweighter",
            ]
        )
    if heavy_commands:
        handler_lines.extend(
            [
                "create HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter "
                "/Herwig/Generators/HeavyFlavorBRReweighter",
                "insert /Herwig/Generators/theGenerator:EventHandler:"
                "PostDecayHandlers 0 "
                "/Herwig/Generators/HeavyFlavorBRReweighter",
            ]
        )
    lines = [
        *(
            [
                "cd /Herwig/Particles",
                f"read {cfg.herwig_b_decays}",
                "cd /Herwig/Generators",
            ]
            if heavy_commands
            else []
        ),
        f"globallibrary {reweighter_library or cfg.decay_reweighter_plugin}",
        *handler_lines,
        *all_commands,
    ]
    return (
        "\n".join(lines),
        {
            "prompt_selection": sample.prompt_decay_profile or None,
            "branching_ratio_reweighter": True,
            "branching_ratio_reweighter_implementations": implementations,
            "handlers": handlers,
            "decay_bias": {
                "enabled": sample.is_hf_biased,
                "scheme": (
                    "full_support_direct_b_semileptonic_importance_sampling"
                    if sample.is_hf_biased
                    else "unbiased"
                ),
                "bias_of": sample.bias_of,
                "selected_parents": mode_counts,
                "parent_proposals": heavy_proposal or None,
                "proposal": (
                    "q_m = b_m*p_m/sum_n(b_n*p_n)"
                    if sample.is_hf_biased
                    else None
                ),
                "bias_function": (
                    f"b_m={HEAVY_FLAVOUR_SEMILEPTONIC_BIAS_FACTOR:g} for "
                    "selected direct one-hadron e/mu modes; b_m=1 otherwise"
                    if sample.is_hf_biased
                    else None
                ),
                "weight_formula": (
                    "product_over_decayed_B_parents("
                    "sum_n(b_n*p_n)/(sum_n(p_n)*b_selected))"
                    if sample.is_hf_biased
                    else None
                ),
                "full_tail_support": True if sample.is_hf_biased else None,
                "closure_validated": False if sample.is_hf_biased else None,
            },
        },
    )


def write_herwig_input(
    sample: SampleSpec,
    lhe_file: Path,
    seed: int,
    nevents: int,
    output: Path,
    cfg: Config,
    *,
    reweighter_library: Path | None = None,
) -> dict[str, object]:
    block, metadata = decay_selection_block(
        sample,
        cfg,
        reweighter_library=reweighter_library,
    )
    if "SelectDecayModes" in block and "BranchingRatioReweighter" not in block:
        raise AssertionError("forced decay selection lacks BranchingRatioReweighter")
    if cfg.dry_run:
        print(f"+ write Herwig input {output}")
        print(f"  LHE={lhe_file} events={nevents} seed={seed} PDF={cfg.herwig_pdf}")
        if block:
            print(block)
        return metadata

    text = cfg.template_in.read_text()
    text = text.replace("LHEFILE.lhe.gz", str(lhe_file))
    text = text.replace("# FOURLEPTON_DECAY_SELECTION", block)
    text = text.replace("FOURLEPTON_RUN", sample.name)
    text = re.sub(
        r"set theGenerator:NumberOfEvents\s+\S+",
        f"set theGenerator:NumberOfEvents {nevents}",
        text,
    )
    text = re.sub(
        r"set theGenerator:RandomNumberGenerator:Seed\s+\S+",
        f"set theGenerator:RandomNumberGenerator:Seed {seed}",
        text,
    )
    text = re.sub(
        r"set /Herwig/Partons/thePDFset:PDFName\s+\S+",
        f"set /Herwig/Partons/thePDFset:PDFName {cfg.herwig_pdf}",
        text,
    )
    output.write_text(text)
    return metadata


def write_root_input_list(
    events_dir: Path,
    output: Path,
    cfg: Config,
) -> list[Path]:
    if cfg.dry_run:
        print(f"+ list {events_dir}/**/*.root in {output}")
        return [events_dir / "HwSim.root"]
    roots = sorted(events_dir.glob("**/*.root"))
    if not roots:
        die(f"no HwSim ROOT files found below {events_dir}")
    output.write_text("".join(f"{path}\n" for path in roots))
    return roots


def expected_analysis_output(sample: SampleSpec, cfg: Config) -> Path:
    return (
        cfg.run_dir
        / "analysis"
        / cfg.detector_response
        / sample.name
        / f"{sample.name}_{cfg.run_tag}_{cfg.detector_response}.root"
    )


def analysis_command(
    sample: SampleSpec,
    root_input: Path,
    seed: int,
    weight_scale: float | str,
    cfg: Config,
    optional_weight_names: Path | None = None,
) -> list[object]:
    output_dir = expected_analysis_output(sample, cfg).parent
    command: list[object] = [
        cfg.analysis_exe,
        "--input-list",
        root_input,
        "--response-profile",
        cfg.detector_response,
        "--seed",
        seed,
        "--weight-scale",
        weight_scale,
        "--muon-resolution-scale",
        cfg.muon_resolution_scale,
        "--tag",
        cfg.run_tag,
        "--output-dir",
        output_dir,
        "--sample",
        sample.name,
        "--category",
        sample.category,
        "--channel",
        sample.channel,
    ]
    if not cfg.pileup_noise:
        command.append("--no-pileup-noise")
    if optional_weight_names is not None:
        command.extend(
            ["--optional-weight-names-file", optional_weight_names]
        )
    return command


_SHA256_CACHE: dict[tuple[str, int, int], str] = {}


def sha256_file(path: Path) -> str:
    stat = path.stat()
    key = (str(path.resolve()), stat.st_size, stat.st_mtime_ns)
    cached = _SHA256_CACHE.get(key)
    if cached is not None:
        return cached
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    value = digest.hexdigest()
    _SHA256_CACHE[key] = value
    return value


def artifact_stage_record(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def artifact_stage_matches(
    marker: Path,
    identity: dict[str, object],
    artifact: Path,
) -> bool:
    try:
        payload = json.loads(marker.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return (
        isinstance(payload, dict)
        and payload.get("identity") == identity
        and payload.get("artifact") == artifact_stage_record(artifact)
    )


def write_artifact_stage_marker(
    marker: Path,
    identity: dict[str, object],
    artifact: Path,
) -> None:
    temporary = marker.with_name(f"{marker.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(
            {
                "identity": identity,
                "artifact": artifact_stage_record(artifact),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    temporary.replace(marker)


def mg5_version(mg5_dir: Path) -> str | None:
    version_file = mg5_dir / "VERSION"
    if not version_file.exists():
        return None
    match = re.search(r"(?im)^\s*version\s*=\s*(\S+)", version_file.read_text())
    return match.group(1) if match else None


def herwig_version(cfg: Config) -> str | None:
    candidates: list[Path] = []
    configured = Path(cfg.herwig).expanduser()
    if configured.is_absolute():
        candidates.append(configured)
    if cfg.herwig_env and cfg.herwig_env.name == "activate":
        candidates.append(cfg.herwig_env.parent / "Herwig")
    found = shutil.which(cfg.herwig)
    if found:
        candidates.append(Path(found))
    for candidate in candidates:
        if not candidate.exists():
            continue
        result = subprocess.run(
            [str(candidate), "--version"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        match = re.search(r"(?m)^Herwig\s+(\S+)", result.stdout)
        if match:
            return match.group(1)
    if cfg.herwig_env or cfg.herwig_module:
        command = "; ".join(
            [*runtime_setup_commands(cfg), f"{q(cfg.herwig)} --version"]
        )
        result = subprocess.run(
            ["bash", "-lc", command],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            env=clean_shell_env(),
            check=False,
        )
        match = re.search(r"(?m)^Herwig\s+(\S+)", result.stdout)
        if match:
            return match.group(1)
    return None


def command_version(command: str, *arguments: str) -> str | None:
    executable = shutil.which(command)
    if executable is None:
        return None
    result = subprocess.run(
        [executable, *arguments],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    first_line = result.stdout.strip().splitlines()
    return first_line[0] if first_line else None


def runtime_command_version(
    cfg: Config,
    command: str,
    *arguments: str,
) -> str | None:
    direct = command_version(command, *arguments)
    if direct is not None:
        return direct
    if not cfg.herwig_env and not cfg.herwig_module:
        return None
    shell_command = "; ".join(
        [
            *runtime_setup_commands(cfg),
            format_command([command, *arguments]),
        ]
    )
    result = subprocess.run(
        ["bash", "-lc", shell_command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        env=clean_shell_env(),
        check=False,
    )
    lines = result.stdout.strip().splitlines()
    return lines[0] if lines else None


def file_record(path: Path, *, dry_run: bool = False) -> dict[str, object]:
    exists = path.exists()
    return {
        "path": str(path),
        "exists": exists,
        "sha256": None if dry_run or not exists else sha256_file(path),
        "checksum_status": (
            "pending_dry_run"
            if dry_run
            else "available"
            if exists
            else "file_not_available"
        ),
    }


def sample_work_dir(sample: SampleSpec, cfg: Config) -> Path:
    return cfg.run_dir / "samples" / cfg.detector_response / sample.name


def sample_decay_dir(sample: SampleSpec, cfg: Config) -> Path:
    """Profile-independent forced-decay products shared by perfect and SSC."""

    return cfg.run_dir / "decays" / sample.name / "madspin"


def production_dir(sample: SampleSpec, cfg: Config) -> Path:
    return cfg.run_dir / "production" / sample.production_group_id


def sample_manifest_entry(sample: SampleSpec, cfg: Config) -> dict[str, object]:
    work_dir = sample_work_dir(sample, cfg)
    process_dir = production_dir(sample, cfg) / "mg5_process"
    process_card = production_dir(sample, cfg) / "cards" / "process_card.dat"
    madspin_card = sample_decay_dir(sample, cfg) / "madspin_card.dat"
    if sample.is_signal:
        target_br = signal_target_branching_fraction(sample, cfg)
        weight: float | str | None = PENDING_SIGNAL_WEIGHT
        formula = "stable_ggh_xsec * external_target_BR_i / madspin_decayed_xsec"
    else:
        target_br = None
        weight = 1.0
        formula = "target=native; weight_scale=1"
    me_lepton_cuts = any(
        token in sample.process for token in ("e+", "e-", "mu+", "mu-")
    )
    return {
        "name": sample.name,
        "label": sample.label,
        "category": sample.category,
        "channel": sample.channel,
        "perturbative_order": (
            cfg.external_perturbative_order
            if cfg.source_backend == "external-lhe"
            else sample.perturbative_order
        ),
        "source_kind": cfg.source_backend,
        "declared_source_kind": sample.source_kind,
        "model": sample.model,
        "process": sample.process,
        "decay_strategy": sample.decay_strategy,
        "normalization_rule": sample.normalization_rule,
        "response_mode": sample.response_mode,
        "cut_mask_definition": CUT_MASK_DEFINITION,
        "production_group_id": sample.production_group_id,
        "nevents_requested": cfg.events_for(sample),
        "seed": seed_for(f"{sample.name}:analysis", cfg),
        "seeds": {
            "production": seed_for(sample.production_group_id, cfg),
            "madspin": (
                seed_for(f"{sample.name}:madspin", cfg)
                if sample.is_signal
                else None
            ),
            "herwig": seed_for(f"{sample.name}:herwig", cfg),
            "analysis": seed_for(f"{sample.name}:analysis", cfg),
        },
        "status": "planned",
        "matrix_element_generation_cuts": (
            dict(GENERATION_CUTS)
            if me_lepton_cuts
            else {
                "ptl_min_gev": None,
                "etal_max": None,
                "mmll_min_gev": None,
                "drll_min": None,
                "reason": "no charged leptons in the hard-process final state",
            }
        ),
        "cards": {
            "process": {"path": str(process_card), "sha256": None},
            "madspin": (
                {"path": str(madspin_card), "sha256": None}
                if sample.madspin_decay
                else None
            ),
            "run": {
                "path": str(process_dir / "Cards" / "run_card.dat"),
                "sha256": None,
            },
            "parameters": {
                "path": str(process_dir / "Cards" / "param_card.dat"),
                "sha256": None,
            },
            "herwig": {
                "path": str(
                    work_dir / "herwig" / f"{sample.name}.in"
                ),
                "sha256": None,
            },
        },
        "normalization": {
            "native_cross_section_pb": None,
            "stable_production_cross_section_pb": None,
            "external_total_h4l_branching_fraction": (
                BR_H_TO_4L if sample.is_signal else None
            ),
            "channel_fraction": (
                cfg.signal_channel_fractions[sample.channel]
                if sample.is_signal
                else None
            ),
            "target_branching_fraction": target_br,
            "target_cross_section_pb": None,
            "weight_scale": weight,
            "generated_sumw": None,
            "generated_sumw2": None,
            "denominator_kind": "predecay_lhe_handler_weight_sum",
            "denominator_formula": "XWGTUP / effective_XMAXUP[IDPRUP]",
            "denominator_sumw": None,
            "denominator_sumabsw": None,
            "denominator_sumw2": None,
            "denominator_event_count": None,
            "denominator_npositive": None,
            "denominator_nnegative": None,
            "denominator_nzero": None,
            "formula": formula,
            "applied_exactly_once_by": "HwSimPostAnalysis_fourlepton",
        },
        "prompt_decay_selection": sample.prompt_decay_profile or None,
        "decay_bias": {
            "enabled": sample.is_hf_biased,
            "scheme": (
                "full_support_direct_b_semileptonic_importance_sampling"
                if sample.is_hf_biased
                else "unbiased"
            ),
            "bias_of": sample.bias_of,
            "proposal": (
                "q_m = b_m*p_m/sum_n(b_n*p_n)"
                if sample.is_hf_biased
                else None
            ),
            "full_tail_support": True if sample.is_hf_biased else None,
            "closure_validated": False if sample.is_hf_biased else None,
        },
        "branching_ratio_reweighter": bool(
            sample.prompt_decay_profile or sample.is_hf_biased
        ),
        "decay": {
            "prompt_selection": sample.prompt_decay_profile or None,
            "branching_ratio_reweighter_enabled": bool(
                sample.prompt_decay_profile or sample.is_hf_biased
            ),
            "heavy_flavour_bias_enabled": sample.is_hf_biased,
            "handlers": decay_reweighter_handlers(sample),
        },
        "external_lhe": (
            str(resolve_external_lhe(sample, cfg))
            if cfg.source_backend == "external-lhe"
            else None
        ),
        "input_lhe": None,
        "input_lhe_init": None,
        "input_lhe_sha256": None,
        "input_lhe_checksum_status": "pending_generation",
        "source_lhe": None,
        "source_lhe_init": None,
        "source_lhe_sha256": None,
        "source_lhe_checksum_status": "pending_generation",
        "higgs_mass_gev": cfg.higgs_mass if sample.is_signal else None,
        "generator_higgs_width_gev": None,
        "madspin_banner_augmentation": None,
        "optional_weight_definitions_status": "pending_generation",
        "optional_weight_definitions": None,
        "optional_weight_names_file": None,
        "output_root": str(expected_analysis_output(sample, cfg)),
        "output_root_sha256": None,
        "output_root_checksum_status": "pending_analysis",
        "root_file": str(expected_analysis_output(sample, cfg)),
        "analysis_summary": str(
            expected_analysis_output(sample, cfg).with_suffix(".summary.json")
        ),
        "analysis_summary_sha256": None,
        "weight_sums": None,
        "analyzer_cli": None,
    }


def campaign_tracked_inputs(cfg: Config) -> dict[str, dict[str, object]]:
    paths = {
        "runner": Path(__file__).resolve(),
        "herwig_template": cfg.template_in,
        "analyzer_source": cfg.analysis_source,
        "analyzer_makefile": cfg.analysis_code_dir / "Makefile",
        "decay_reweighter_source": cfg.decay_reweighter_source,
        "decay_reweighter_makefile": cfg.decay_reweighter_source.parent
        / "Makefile",
    }
    if cfg.herwig_b_decays is not None:
        paths["herwig_b_decays"] = cfg.herwig_b_decays
    return {
        name: file_record(path, dry_run=False)
        for name, path in sorted(paths.items())
    }


def campaign_software_versions(cfg: Config) -> dict[str, str | None]:
    return {
        "python": platform.python_version(),
        "mg5": mg5_version(cfg.mg5_dir),
        "herwig": herwig_version(cfg),
        "root": runtime_command_version(cfg, "root-config", "--version"),
    }


def configuration_fingerprint(
    cfg: Config,
    samples: Sequence[SampleSpec],
    *,
    tracked_inputs: dict[str, dict[str, object]] | None = None,
    software_versions: dict[str, str | None] | None = None,
) -> str:
    if tracked_inputs is None:
        tracked_inputs = campaign_tracked_inputs(cfg)
    if software_versions is None:
        software_versions = campaign_software_versions(cfg)
    payload = {
        "sample_set": cfg.sample_set,
        "samples": [
            {
                "name": sample.name,
                "model": sample.model,
                "process": sample.process,
                "madspin_decay": sample.madspin_decay,
                "prompt_decay_profile": sample.prompt_decay_profile,
                "nevents": cfg.events_for(sample),
            }
            for sample in samples
        ],
        "detector_response": cfg.detector_response,
        "muon_resolution_scale": cfg.muon_resolution_scale,
        "pileup_noise": cfg.pileup_noise,
        "source_backend": cfg.source_backend,
        "external_lhe": {
            key: {
                "path": str(path.resolve()),
                "sha256": (
                    sha256_file(path)
                    if path.exists() and not cfg.dry_run
                    else None
                ),
            }
            for key, path in sorted(cfg.external_lhe.items())
        },
        "external_cross_sections_pb": dict(
            sorted(cfg.external_cross_sections_pb.items())
        ),
        "external_perturbative_order": cfg.external_perturbative_order,
        "ebeam": cfg.ebeam,
        "higgs_mass": cfg.higgs_mass,
        "signal_channel_fractions": cfg.signal_channel_fractions,
        "signal_channel_fraction_source": cfg.signal_channel_fraction_source,
        "signal_channel_fraction_input": cfg.signal_channel_fraction_input,
        "pdf_lhaid": LO_PDF_LHAID,
        "seed_base": cfg.seed_base,
        "heavy_flavour_bias": cfg.heavy_flavour_bias,
        "analysis_target": cfg.analysis_target,
        "tracked_inputs": tracked_inputs,
        "software_versions": software_versions,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def build_manifest(
    cfg: Config,
    samples: Sequence[SampleSpec],
) -> dict[str, object]:
    tracked_inputs = campaign_tracked_inputs(cfg)
    software_versions = campaign_software_versions(cfg)
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "analysis": "H_to_ZZstar_to_4l",
        "status": "planned",
        "run_tag": cfg.run_tag,
        "configuration_fingerprint": configuration_fingerprint(
            cfg,
            samples,
            tracked_inputs=tracked_inputs,
            software_versions=software_versions,
        ),
        "collider": {
            "beam": "pp",
            "beam_energy_gev": cfg.ebeam,
            "sqrt_s_tev": 2.0 * cfg.ebeam / 1000.0,
        },
        "luminosity_fb": cfg.luminosity_fb,
        "detector_profile": cfg.detector_response,
        "response_mode": cfg.detector_response,
        "sample_set": cfg.sample_set,
        "source_backend": cfg.source_backend,
        "external_cross_sections_pb": dict(
            sorted(cfg.external_cross_sections_pb.items())
        ),
        "pdf": {
            "name": (
                None if cfg.source_backend == "external-lhe" else LO_PDF_NAME
            ),
            "lhaid": (
                None if cfg.source_backend == "external-lhe" else LO_PDF_LHAID
            ),
            "perturbative_order": (
                cfg.external_perturbative_order
                if cfg.source_backend == "external-lhe"
                else "LO"
            ),
            "hard_process_provenance": (
                "LHE_init_header"
                if cfg.source_backend == "external-lhe"
                else "campaign_MG5_run_card"
            ),
            "herwig_shower_pdf_name": cfg.herwig_pdf,
        },
        "higgs": {
            "mass_gev": cfg.higgs_mass,
            "external_total_h4l_branching_fraction": BR_H_TO_4L,
            "channel_fractions": dict(cfg.signal_channel_fractions),
            "channel_fraction_provenance": (
                SIGNAL_FRACTION_PROVENANCE
                if cfg.signal_channel_fraction_source
                == "default_mg5_integration"
                else cfg.signal_channel_fraction_source
            ),
            "channel_fraction_input": cfg.signal_channel_fraction_input,
            "channel_fraction_note": (
                "Reference LO direct-1to4 integration; approximately 1% "
                "integration precision, not a PROPHECY4F-exact partition."
                if cfg.signal_channel_fraction_source
                == "default_mg5_integration"
                else "User-supplied normalized channel partition; no default "
                "MG5 integration provenance is claimed."
            ),
            "branching_fraction_reference": (
                "https://cds.cern.ch/record/2273853/files/"
                "ATLAS-CONF-2017-046.pdf"
            ),
        },
        "generation_cuts": dict(GENERATION_CUTS),
        "cut_stages": [
            {
                "name": name,
                "label": label,
                "bit": bit,
                "required_mask": (1 << (bit + 1)) - 1,
            }
            for name, label, bit in CUT_STAGES
        ],
        "detector": {
            "profile": cfg.detector_response,
            "muon_resolution_scale": cfg.muon_resolution_scale,
            "pileup_noise_configured": cfg.pileup_noise,
            "pileup_noise_enabled": (
                cfg.pileup_noise and cfg.detector_response == "ssc"
            ),
            "active_efficiencies": (
                {"electron": 1.0, "muon": 1.0, "trigger": {"all": 1.0}}
                if cfg.detector_response == "perfect"
                else {
                    "electron": 0.90,
                    "muon": 0.85 * 0.95,
                    "trigger": {"4e": 0.98, "4mu": 0.99, "2e2mu": 0.99},
                }
            ),
            "electron": {
                "reconstruction_efficiency": 0.90,
                "resolution_stochastic_barrel": 0.060,
                "resolution_stochastic_endcap": 0.085,
                "resolution_constant": 0.004,
                "barrel_boundary_abs_eta": 1.01,
                "thermal_noise_et_gev": {
                    "barrel": 0.100,
                    "outside_barrel": 0.175,
                },
                "pileup_noise_et_gev": {
                    "abs_eta_below_1p4": 0.120,
                    "abs_eta_above_1p4_formula": "0.120 + 0.366*(abs_eta-1.4)",
                },
                "noise_combination": (
                    "thermal_and_pileup_in_quadrature_then_divide_by_pt"
                ),
            },
            "muon": {
                "reconstruction_efficiency": 0.85 * 0.95,
                "smearing_variable": "q_over_pt",
                "table_identity": "GEM_TDR_Fig4-20_approx_digitization",
                "eta_nodes_abs": [0.1, 0.6, 1.0, 1.3, 1.6, 1.85, 2.0, 2.2, 2.5],
                "relative_resolution_rows": {
                    "pt10": [0.028, 0.0265, 0.0255, 0.0215, 0.019, 0.0195, 0.020, 0.0205, 0.022],
                    "pt25": [0.0167, 0.0155, 0.015, 0.0175, 0.0182, 0.0195, 0.014, 0.0135, 0.0125],
                    "pt50": [0.0148, 0.0138, 0.0138, 0.0185, 0.0198, 0.020, 0.0145, 0.0145, 0.0135],
                    "pt100": [0.0135, 0.0125, 0.0125, 0.0195, 0.023, 0.025, 0.0205, 0.0235, 0.0245],
                },
                "high_pt_anchor": {
                    "pt_gev": 500.0,
                    "relative_resolution_eta0": 0.05,
                    "relative_resolution_eta2p5": 0.12,
                },
                "interpolation": "linear_abs_eta_and_logarithmic_pt",
                "clamp_pt_gev": [10.0, 500.0],
                "variation": {
                    "parameter": "multiplicative_muon_resolution_scale",
                    "nominal": 1.0,
                    "down": 0.8,
                    "up": 1.2,
                },
            },
            "trigger_efficiency": {"4e": 0.98, "other": 0.99},
            "declared_unmodeled_effects": [
                "charge_misidentification",
                "explicit_pileup_event_overlay",
                "conversions",
                "impact_parameter_response",
            ],
            "references": {
                "four_lepton_study": (
                    "https://lss.fnal.gov/archive/other/calt-68-1856.pdf"
                ),
                "gem_tdr": (
                    "https://lss.fnal.gov/archive/other/ssc/"
                    "ssc-gem-tn-93-262.pdf"
                ),
            },
        },
        "heavy_flavour_bias_policy": {
            "requested": cfg.heavy_flavour_bias,
            "nominal_samples_are_unbiased": True,
            "accelerated_scheme": (
                "full_support_direct_b_semileptonic_importance_sampling"
                if cfg.heavy_flavour_bias == "semileptonic"
                else None
            ),
            "proposal": (
                "q_m = b_m*p_m/sum_n(b_n*p_n)"
                if cfg.heavy_flavour_bias == "semileptonic"
                else None
            ),
            "weighted_closure_against_unbiased_control_required": True,
            "closure_can_be_claimed_by_runner": False,
        },
        "software": {
            "runner": str(Path(__file__).resolve()),
            "runner_sha256": sha256_file(Path(__file__).resolve()),
            "python_version": software_versions["python"],
            "mg5_dir": str(cfg.mg5_dir),
            "mg5_version": software_versions["mg5"],
            "herwig_command": cfg.herwig,
            "herwig_version": software_versions["herwig"],
            "herwig_module": cfg.herwig_module,
            "herwig_pdf": cfg.herwig_pdf,
            "root_version": software_versions["root"],
            "analyzer_target": cfg.analysis_target,
            "analyzer_schema_version": 1,
        },
        "tracked_inputs": tracked_inputs,
        "samples": [sample_manifest_entry(sample, cfg) for sample in samples],
    }


def manifest_sample(
    manifest: dict[str, object],
    sample_name: str,
) -> dict[str, object]:
    samples = manifest["samples"]
    assert isinstance(samples, list)
    for entry in samples:
        assert isinstance(entry, dict)
        if entry.get("name") == sample_name:
            return entry
    raise KeyError(sample_name)


def write_manifest(manifest: dict[str, object], cfg: Config) -> None:
    text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if cfg.dry_run:
        print(f"+ would write manifest {cfg.manifest_path}")
        return
    cfg.manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = cfg.manifest_path.with_name(
        f"{cfg.manifest_path.name}.tmp.{os.getpid()}"
    )
    temporary.write_text(text)
    temporary.replace(cfg.manifest_path)


def load_resumable_manifest(cfg: Config) -> dict[str, object] | None:
    if cfg.dry_run or cfg.force or not cfg.manifest_path.exists():
        return None
    try:
        manifest = json.loads(cfg.manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        die(f"cannot read existing manifest {cfg.manifest_path}: {error}")
    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        die(
            f"manifest schema mismatch in {cfg.manifest_path}; "
            "use a new --run-tag or explicit --force"
        )
    if manifest.get("detector_profile") != cfg.detector_response:
        die(f"detector profile mismatch in {cfg.manifest_path}")
    return manifest


def run_generate_events(
    process_dir: Path,
    runtime_dir: Path | None,
    cfg: Config,
) -> None:
    command = [
        process_dir / "bin" / "generate_events",
        cfg.run_tag,
        "-f",
        f"--nb_core={cfg.nb_core}",
    ]
    run_runtime_command(
        command,
        cfg,
        cwd=process_dir,
        env=generation_environment(runtime_dir),
    )


def production_stage_identity(
    sample: SampleSpec,
    representative: SampleSpec,
    selected: Sequence[SampleSpec],
    process_dir: Path,
    cfg: Config,
) -> dict[str, object]:
    process_card = process_card_text(representative, process_dir)
    return {
        "schema_version": 1,
        "production_group_id": sample.production_group_id,
        "model": representative.model,
        "process": representative.process,
        "process_card_sha256": hashlib.sha256(
            process_card.encode("utf-8")
        ).hexdigest(),
        "events": group_event_count(
            sample.production_group_id,
            selected,
            cfg,
        ),
        "seed": seed_for(sample.production_group_id, cfg),
        "beam_energy_gev": cfg.ebeam,
        "higgs_mass_gev": cfg.higgs_mass,
        "pdf_lhaid": LO_PDF_LHAID,
        "generator_cuts": dict(GENERATION_CUTS),
        "mg5_version": mg5_version(cfg.mg5_dir),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "collier_library_sha256": (
            sha256_file(cfg.collier_library)
            if cfg.collier_library is not None
            and cfg.collier_library.exists()
            else None
        ),
    }


def obtain_mg5_lhe(
    sample: SampleSpec,
    selected: Sequence[SampleSpec],
    cfg: Config,
) -> Path:
    group_dir = production_dir(sample, cfg)
    process_dir = group_dir / "mg5_process"
    marker = group_dir / "fourlepton_production_stage.json"
    representative = next(
        candidate
        for candidate in selected
        if candidate.production_group_id == sample.production_group_id
    )
    identity = production_stage_identity(
        sample,
        representative,
        selected,
        process_dir,
        cfg,
    )
    if not cfg.dry_run and not cfg.force and process_dir.exists():
        try:
            existing_lhe = find_lhe_file(process_dir)
        except SystemExit:
            existing_lhe = None
        if existing_lhe is not None:
            if artifact_stage_matches(marker, identity, existing_lhe):
                return existing_lhe
            die(
                f"existing production artifacts for "
                f"{sample.production_group_id} do not match the current "
                "process/card/software provenance; rerun with --force"
            )

    card_dir = group_dir / "cards"
    card = card_dir / "process_card.dat"
    ensure_dir(card_dir, cfg)
    write_mg5_process_card(card, representative, process_dir, cfg)
    run_runtime_command(mg5_command(cfg, card), cfg)

    patch_mg5_toolchain(process_dir, cfg)
    patch_mg5_hidden_pdf_defaults(process_dir, cfg)
    patch_mg5_pdf_label_sync(process_dir, cfg)
    patch_run_card(
        process_dir / "Cards" / "run_card.dat",
        sample=representative,
        nevents=group_event_count(sample.production_group_id, selected, cfg),
        seed=seed_for(sample.production_group_id, cfg),
        cfg=cfg,
    )
    generator_width = patch_param_card_higgs_mass(
        process_dir / "Cards" / "param_card.dat",
        cfg,
    )
    if generator_width is not None:
        print(f"  generator Higgs width={generator_width:.9g} GeV")
    generate_and_patch_run_card_include(process_dir, cfg)

    runtime_dir: Path | None = None
    if cfg.collier_library and representative.model.startswith("loop_"):
        runtime_dir = group_dir / "runtime_lib"
        prepare_collier_runtime(runtime_dir, cfg)
        prepare_loop_runtime(process_dir, runtime_dir, cfg)
    if representative.model.startswith("loop_"):
        restore_mg5_symmetry_factors(process_dir, cfg)
    run_generate_events(process_dir, runtime_dir, cfg)
    if cfg.dry_run:
        return process_dir / "Events" / cfg.run_tag / "unweighted_events.lhe.gz"
    generated_lhe = find_lhe_file(process_dir)
    write_artifact_stage_marker(marker, identity, generated_lhe)
    return generated_lhe


def obtain_production_lhe(
    sample: SampleSpec,
    selected: Sequence[SampleSpec],
    cfg: Config,
) -> Path:
    if cfg.source_backend == "external-lhe":
        return resolve_external_lhe(sample, cfg).resolve()
    return obtain_mg5_lhe(sample, selected, cfg)


def madspin_output_path(input_lhe: Path) -> Path:
    name = input_lhe.name.replace(".lhe", "_decayed.lhe")
    if not name.endswith(".gz"):
        name += ".gz"
    return input_lhe.with_name(name)


def stage_madspin_input(source: Path, destination: Path, cfg: Config) -> None:
    """Stage an uncompressed LHE because MG5 3.5 cannot import gzip banners."""

    def logical_sha256(path: Path) -> str:
        digest = hashlib.sha256()
        opener = gzip.open if path.name.endswith(".gz") else Path.open
        with opener(path, "rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    if cfg.dry_run:
        print(f"+ stage uncompressed MadSpin input {source} -> {destination}")
        return
    if destination.exists() and not cfg.force:
        if logical_sha256(source) == logical_sha256(destination):
            return
        raise ValueError(
            f"staged MadSpin input {destination} differs from {source}; "
            "rerun with --force"
        )
    temporary = destination.with_name(f"{destination.name}.tmp.{os.getpid()}")
    if source.name.endswith(".gz"):
        with gzip.open(source, "rb") as input_stream, temporary.open(
            "wb"
        ) as output_stream:
            shutil.copyfileobj(input_stream, output_stream)
    else:
        shutil.copy2(source, temporary)
    temporary.replace(destination)


_SLHA_BLOCK_RE = re.compile(
    r"<slha\b[^>]*>(?P<body>.*?)</slha\s*>",
    flags=re.IGNORECASE | re.DOTALL,
)
_SLHA_NUMBER_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?"
_MADSPIN_REQUIRED_BLOCK_ENTRIES: dict[str, tuple[int, ...]] = {
    "sminputs": (1, 2, 3),
    "mass": (5, 6, 15, 23, 24, 25),
    "yukawa": (5, 6, 15),
    "decay": (6, 23, 24, 25),
}


def _slha_block_values(slha: str, block_name: str) -> dict[int, list[float]]:
    if block_name.lower() == "decay":
        matches = re.findall(
            rf"(?im)^\s*DECAY\s+(\d+)\s+({_SLHA_NUMBER_RE})(?:\s|#|$)",
            slha,
        )
    else:
        block = re.search(
            rf"(?ims)^\s*BLOCK\s+{re.escape(block_name)}\b"
            r"(?P<body>.*?)(?=^\s*(?:BLOCK|DECAY)\b|\Z)",
            slha,
        )
        matches = (
            re.findall(
                rf"(?im)^\s*(\d+)\s+({_SLHA_NUMBER_RE})(?:\s|#|$)",
                block.group("body"),
            )
            if block is not None
            else []
        )
    values: dict[int, list[float]] = {}
    for raw_code, raw_value in matches:
        values.setdefault(int(raw_code), []).append(
            float(raw_value.replace("D", "E").replace("d", "e"))
        )
    return values


def canonical_madspin_slha(cfg: Config) -> str:
    card = cfg.mg5_dir / "models" / "sm" / "restrict_default.dat"
    if not card.exists():
        raise ValueError(f"missing canonical SM parameter card: {card}")
    text = card.read_text()
    pattern = re.compile(
        r"(?im)^(\s*25\s+)([-+0-9.eEdD]+)(\s*(?:#.*)?)$"
    )
    replacement = rf"\g<1>{cfg.higgs_mass:.9e}\g<3>"
    updated, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise ValueError(f"cannot patch MASS 25 in canonical SM card: {card}")

    # ``restrict_default.dat`` is a UFO restriction card rather than a
    # standalone-complete param_card.  In particular the W mass is a
    # dependent parameter and is omitted, while MadSpin's generated process
    # checks MASS(24) when it reads DECAY(24).  Reconstruct the same on-shell
    # value emitted by MG5's generic SM param-card writer.
    sminputs = _slha_block_values(updated, "sminputs")
    masses = _slha_block_values(updated, "mass")
    try:
        alpha_ew = 1.0 / sminputs[1][0]
        fermi_constant = sminputs[2][0]
        z_mass = masses[23][0]
        radicand = (
            z_mass**4 / 4.0
            - alpha_ew * math.pi * z_mass**2
            / (fermi_constant * math.sqrt(2.0))
        )
        w_mass = math.sqrt(z_mass**2 / 2.0 + math.sqrt(radicand))
    except (KeyError, IndexError, ValueError, ZeroDivisionError) as error:
        raise ValueError(
            f"cannot derive dependent MASS 24 from canonical SM card: {card}"
        ) from error
    if not math.isfinite(w_mass) or w_mass <= 0.0:
        raise ValueError(
            f"canonical SM card gives an invalid dependent MASS 24: {card}"
        )
    mass_block = re.search(
        r"(?ims)^\s*BLOCK\s+MASS\b(?P<body>.*?)(?=^\s*(?:BLOCK|DECAY)\b|\Z)",
        updated,
    )
    if mass_block is None:
        raise ValueError(f"canonical SM card has no MASS block: {card}")
    mass_body = mass_block.group("body")
    if 24 not in _slha_block_values(updated, "mass"):
        replacement_body = (
            mass_body.rstrip()
            + f"\n   24 {w_mass:.9e} # MW (derived by MG5 SM model)\n\n"
        )
        updated = (
            updated[: mass_block.start("body")]
            + replacement_body
            + updated[mass_block.end("body") :]
        )

    width_match = re.search(
        r"(?im)^\s*DECAY\s+25\s+([-+0-9.eEdD]+)",
        updated,
    )
    if width_match is None or float(
        width_match.group(1).replace("D", "E").replace("d", "e")
    ) <= 0.0:
        raise ValueError(
            f"canonical SM card has no positive Higgs width: {card}"
        )
    completed = updated.rstrip() + "\n"
    validation = validate_madspin_slha(completed, cfg)
    if not validation["valid"]:
        raise ValueError(
            "canonical SM card is not standalone-complete for MadSpin: "
            + ", ".join(str(reason) for reason in validation["reasons"])
        )
    return completed


def validate_madspin_slha(
    slha: str,
    cfg: Config,
) -> dict[str, object]:
    """Validate all SM inputs required by standalone MadSpin."""

    block_values = {
        name: _slha_block_values(slha, name)
        for name in _MADSPIN_REQUIRED_BLOCK_ENTRIES
    }
    mass_matches = block_values["mass"].get(25, [])
    width_matches = block_values["decay"].get(25, [])

    mass: float | None = None
    width: float | None = None
    reasons: list[str] = []
    if len(mass_matches) != 1:
        reasons.append("missing_or_ambiguous_mass_25")
    else:
        mass = mass_matches[0]
        if not math.isfinite(mass) or not math.isclose(
            mass,
            cfg.higgs_mass,
            rel_tol=0.0,
            abs_tol=1.0e-4,
        ):
            reasons.append("invalid_mass_25")
    if len(width_matches) != 1:
        reasons.append("missing_or_ambiguous_decay_25")
    else:
        width = width_matches[0]
        if not math.isfinite(width) or width <= 0.0:
            reasons.append("nonpositive_decay_25")

    for block_name, required_codes in _MADSPIN_REQUIRED_BLOCK_ENTRIES.items():
        values = block_values[block_name]
        for code in required_codes:
            entries = values.get(code, [])
            reason = f"missing_or_ambiguous_{block_name}_{code}"
            if len(entries) != 1:
                if not (
                    block_name == "mass"
                    and code == 25
                    and "missing_or_ambiguous_mass_25" in reasons
                ) and not (
                    block_name == "decay"
                    and code == 25
                    and "missing_or_ambiguous_decay_25" in reasons
                ):
                    reasons.append(reason)
                continue
            value = entries[0]
            if not math.isfinite(value):
                reasons.append(f"nonfinite_{block_name}_{code}")
            elif block_name in {"sminputs", "mass", "decay"} and value <= 0.0:
                if not (
                    block_name == "decay"
                    and code == 25
                    and "nonpositive_decay_25" in reasons
                ):
                    reasons.append(f"nonpositive_{block_name}_{code}")

    return {
        "valid": not reasons,
        "mass_gev": mass,
        "width_gev": width,
        "reasons": reasons,
    }


def madspin_banner_plan(path: Path, cfg: Config) -> dict[str, object]:
    lines: list[str] = []
    with path.open("rt", errors="replace") as stream:
        for line in stream:
            lines.append(line)
            if re.search(r"<init\b", line, flags=re.IGNORECASE):
                break
            if len(lines) > 200_000:
                raise ValueError(f"could not find <init> near the start of {path}")
    prefix = "".join(lines)
    if not re.search(r"<init\b", prefix, flags=re.IGNORECASE):
        raise ValueError(f"LHE file has no <init> block: {path}")
    has_header_open = bool(
        re.search(r"<header\b", prefix, flags=re.IGNORECASE)
    )
    has_header_close = bool(
        re.search(r"</header>", prefix, flags=re.IGNORECASE)
    )
    if has_header_open != has_header_close:
        raise ValueError(f"LHE file has a malformed <header> block: {path}")
    run_card_open_count = len(
        re.findall(r"<MGRunCard\b", prefix, flags=re.IGNORECASE)
    )
    run_card_close_count = len(
        re.findall(r"</MGRunCard\s*>", prefix, flags=re.IGNORECASE)
    )
    if (
        run_card_open_count != run_card_close_count
        or run_card_open_count > 1
    ):
        raise ValueError(
            f"LHE file has a malformed or ambiguous <MGRunCard> block: {path}"
        )
    slha_blocks = list(_SLHA_BLOCK_RE.finditer(prefix))
    slha_open_count = len(
        re.findall(r"<slha\b", prefix, flags=re.IGNORECASE)
    )
    slha_close_count = len(
        re.findall(r"</slha\s*>", prefix, flags=re.IGNORECASE)
    )
    if slha_open_count != slha_close_count or slha_open_count > 1:
        raise ValueError(f"LHE file has a malformed or ambiguous <slha> block: {path}")
    slha_validation: dict[str, object]
    if slha_blocks:
        slha_validation = validate_madspin_slha(
            slha_blocks[0].group("body"),
            cfg,
        )
    else:
        slha_validation = {
            "valid": False,
            "mass_gev": None,
            "width_gev": None,
            "reasons": ["missing_slha"],
        }
    has_slha = bool(slha_blocks)
    valid_slha = bool(slha_validation["valid"])
    return {
        "added_mg5_process_card": not bool(
            re.search(r"<MG5ProcCard\b", prefix, flags=re.IGNORECASE)
        ),
        "added_mg_run_card": run_card_open_count == 0,
        "added_slha": not has_slha,
        "replaced_slha": has_slha and not valid_slha,
        "preserved_slha": has_slha and valid_slha,
        "slha_validation": slha_validation,
        "created_header": not has_header_open,
        "strategy": "standalone_complete_sm_madspin_metadata_v3",
    }


def prepare_madspin_banner_input(
    source: Path,
    destination: Path,
    cfg: Config,
    *,
    nevents: int | None = None,
) -> dict[str, object]:
    """Add only the metadata that standalone MadSpin needs.

    Generic POWHEG-style LHE files need not contain MG5's process card, run
    card, or an SLHA block.  The injected SM metadata is used solely to
    generate and serialize the direct Higgs decay; production momenta and
    weights remain untouched.
    """

    plan = madspin_banner_plan(source, cfg)
    needs_canonical_slha = bool(
        plan["added_slha"] or plan["replaced_slha"]
    )
    slha = canonical_madspin_slha(cfg) if needs_canonical_slha else ""
    proc_card = """<MG5ProcCard>
<![CDATA[
# Injected by HiggsSSC for standalone MadSpin decay of a stable Higgs.
import model sm
define p = g u c d s b u~ c~ d~ s~ b~
generate p p > h
]]>
</MG5ProcCard>
"""
    requested_events = cfg.nevents if nevents is None else nevents
    run_card = f"""<MGRunCard>
<![CDATA[
# Minimal source metadata required by standalone MadSpin bridge mode.
 {requested_events} = nevents
 1 = lpp1
 1 = lpp2
 {cfg.ebeam:.9e} = ebeam1
 {cfg.ebeam:.9e} = ebeam2
 15.0 = bwcutoff
 0 = ickkw
 3.0 = lhe_version
 sum = event_norm
]]>
</MGRunCard>
"""
    additions = ""
    if plan["added_mg5_process_card"]:
        additions += proc_card
    if plan["added_mg_run_card"]:
        additions += run_card
    if plan["added_slha"]:
        additions += f"<slha>\n{slha}</slha>\n"

    canonical_card = cfg.mg5_dir / "models" / "sm" / "restrict_default.dat"
    identity = {
        "schema_version": 1,
        "source_sha256": sha256_file(source),
        "augmentation": plan,
        "higgs_mass_gev": cfg.higgs_mass,
        "injected_run_card_nevents": (
            requested_events if plan["added_mg_run_card"] else None
        ),
        "injected_run_card_beam_energy_gev": (
            cfg.ebeam if plan["added_mg_run_card"] else None
        ),
        "canonical_sm_card_sha256": (
            sha256_file(canonical_card)
            if needs_canonical_slha and canonical_card.exists()
            else None
        ),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
    }
    marker = destination.with_name("fourlepton_madspin_banner_stage.json")
    if destination.exists() and not cfg.force:
        if artifact_stage_matches(marker, identity, destination):
            return plan
        raise ValueError(
            f"prepared MadSpin input {destination} does not match the current "
            "source/banner augmentation; rerun with --force"
        )
    if cfg.force:
        if destination.exists():
            destination.unlink()
        if marker.exists():
            marker.unlink()

    temporary = destination.with_name(f"{destination.name}.tmp.{os.getpid()}")
    with source.open("rt", errors="replace") as input_stream:
        preamble_lines: list[str] = []
        for line in input_stream:
            preamble_lines.append(line)
            if re.search(r"<init\b", line, flags=re.IGNORECASE):
                break
            if len(preamble_lines) > 200_000:
                raise ValueError(f"could not find <init> near the start of {source}")
        preamble = "".join(preamble_lines)
        if not re.search(r"<init\b", preamble, flags=re.IGNORECASE):
            raise ValueError(f"LHE file has no <init> block: {source}")

        if plan["replaced_slha"]:
            replacement = f"<slha>\n{slha}</slha>"
            preamble, replacement_count = _SLHA_BLOCK_RE.subn(
                replacement,
                preamble,
                count=1,
            )
            if replacement_count != 1:
                raise ValueError(f"could not replace invalid SLHA block in {source}")

        if plan["created_header"]:
            preamble, insertion_count = re.subn(
                r"(?i)(<init\b)",
                f"<header>\n{additions}</header>\n\\1",
                preamble,
                count=1,
            )
        else:
            preamble, insertion_count = re.subn(
                r"(?i)(</header>)",
                additions + r"\1",
                preamble,
                count=1,
            )
        if insertion_count != 1:
            raise ValueError(f"could not augment MadSpin banner in {source}")

        with temporary.open("wt") as output_stream:
            output_stream.write(preamble)
            shutil.copyfileobj(input_stream, output_stream)
    prepared_plan = madspin_banner_plan(temporary, cfg)
    if (
        prepared_plan["added_mg5_process_card"]
        or prepared_plan["added_mg_run_card"]
        or not bool(prepared_plan["slha_validation"]["valid"])
    ):
        temporary.unlink(missing_ok=True)
        raise ValueError(
            f"prepared MadSpin banner still lacks valid model metadata: {source}"
        )
    temporary.replace(destination)
    write_artifact_stage_marker(marker, identity, destination)
    return plan


def madspin_stage_identity(
    sample: SampleSpec,
    stable_lhe: Path,
    staged_lhe: Path,
    card: Path,
    cfg: Config,
) -> dict[str, object]:
    frontend = cfg.mg5_dir / "MadSpin" / "madspin"
    return {
        "schema_version": 1,
        "sample": sample.name,
        "decay": sample.madspin_decay,
        "spinmode": "none",
        "seed": seed_for(f"{sample.name}:madspin", cfg),
        "stable_lhe_sha256": sha256_file(stable_lhe),
        "staged_lhe_sha256": sha256_file(staged_lhe),
        "card_sha256": sha256_file(card),
        "madspin_frontend_sha256": (
            sha256_file(frontend) if frontend.exists() else None
        ),
        "mg5_version": mg5_version(cfg.mg5_dir),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
    }


def run_signal_madspin(
    sample: SampleSpec,
    stable_lhe: Path,
    cfg: Config,
) -> tuple[Path, Path, dict[str, object]]:
    if not sample.is_signal or not sample.madspin_decay:
        raise ValueError(f"not a MadSpin signal sample: {sample.name}")
    decay_dir = sample_decay_dir(sample, cfg)
    ensure_dir(decay_dir, cfg)
    raw_input = decay_dir / "stable_input_source.lhe"
    local_input = decay_dir / "stable_input.lhe"
    output = madspin_output_path(local_input)
    card = decay_dir / "madspin_card.dat"
    marker = decay_dir / "fourlepton_madspin_stage.json"
    try:
        stage_madspin_input(stable_lhe, raw_input, cfg)
    except ValueError as error:
        die(str(error))
    if cfg.dry_run:
        augmentation: dict[str, object] = {
            "status": "pending_dry_run",
            "strategy": "standalone_complete_sm_madspin_metadata_v3",
        }
    else:
        try:
            augmentation = prepare_madspin_banner_input(
                raw_input,
                local_input,
                cfg,
                nevents=cfg.events_for(sample),
            )
        except ValueError as error:
            die(str(error))
    write_madspin_card(card, local_input, sample, cfg)
    identity: dict[str, object] | None = None
    if not cfg.dry_run:
        identity = madspin_stage_identity(
            sample,
            stable_lhe,
            local_input,
            card,
            cfg,
        )
        if not cfg.force and output.exists():
            if artifact_stage_matches(marker, identity, output):
                return output, card, augmentation
            die(
                f"existing MadSpin output for {sample.name} does not match "
                "the current stable LHE/card/software provenance; "
                "rerun with --force"
            )
        if cfg.force:
            if output.exists():
                output.unlink()
            if marker.exists():
                marker.unlink()
    run_runtime_command(madspin_command(cfg, card), cfg, cwd=decay_dir)
    if cfg.dry_run:
        return output, card, augmentation
    produced: Path
    if output.exists():
        produced = output
    else:
        candidates = sorted(decay_dir.glob("*_decayed*.lhe*"))
        if not candidates:
            die(
                f"MadSpin produced no decayed LHE for {sample.name} "
                f"in {decay_dir}"
            )
        produced = candidates[-1]
    assert identity is not None
    write_artifact_stage_marker(marker, identity, produced)
    return produced, card, augmentation


def run_herwig(
    sample: SampleSpec,
    lhe_file: Path,
    cfg: Config,
) -> tuple[Path, dict[str, object]]:
    work_dir = sample_work_dir(sample, cfg)
    herwig_dir = work_dir / "herwig"
    events_dir = herwig_dir / "events"
    root_input = work_dir / f"{sample.name}_hwsim_roots.input"
    stage_marker = herwig_dir / "fourlepton_herwig_stage.json"
    ensure_dir(events_dir, cfg)

    reweighter_library: Path | None = None
    if sample.prompt_decay_profile or sample.is_hf_biased:
        source_plugin = cfg.decay_reweighter_plugin
        reweighter_library = herwig_dir / source_plugin.name
        if cfg.dry_run:
            print(
                f"+ stage {source_plugin} as {reweighter_library} "
                "for Herwig run-file reload"
            )
        else:
            if not source_plugin.exists():
                die(f"missing forced-decay reweighter plugin: {source_plugin}")
            shutil.copy2(source_plugin, reweighter_library)

    herwig_input = herwig_dir / f"{sample.name}.in"
    metadata = write_herwig_input(
        sample,
        lhe_file,
        seed_for(f"{sample.name}:herwig", cfg),
        cfg.events_for(sample),
        herwig_input,
        cfg,
        reweighter_library=reweighter_library,
    )
    if not cfg.dry_run:
        identity = {
            "schema_version": 1,
            "sample": sample.name,
            "detector_response": cfg.detector_response,
            "events": cfg.events_for(sample),
            "seed": seed_for(f"{sample.name}:herwig", cfg),
            "lhe_sha256": sha256_file(lhe_file),
            "herwig_card_sha256": sha256_file(herwig_input),
            "reweighter_plugin_sha256": (
                sha256_file(reweighter_library)
                if reweighter_library is not None
                else None
            ),
        }
        existing_roots = sorted(events_dir.glob("**/*.root"))
        if cfg.force:
            for path in existing_roots:
                path.unlink()
            if stage_marker.exists():
                stage_marker.unlink()
            existing_roots = []
        elif existing_roots:
            if herwig_stage_matches(stage_marker, identity, existing_roots):
                root_input.write_text(
                    "".join(f"{path}\n" for path in existing_roots)
                )
                return root_input, metadata
            die(
                f"partial Herwig artifacts for {sample.name} do not match "
                "the current LHE/card/plugin provenance; rerun with --force"
            )

    run_runtime_command(
        [cfg.herwig, "read", herwig_input],
        cfg,
        cwd=herwig_dir,
    )
    run_runtime_command(
        [cfg.herwig, "run", f"{sample.name}.run", f"-N{cfg.events_for(sample)}"],
        cfg,
        cwd=herwig_dir,
    )
    roots = write_root_input_list(events_dir, root_input, cfg)
    if not cfg.dry_run:
        root_records = [
            {
                "path": str(path.resolve()),
                "size": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
            }
            for path in roots
        ]
        temporary = stage_marker.with_name(
            f"{stage_marker.name}.tmp.{os.getpid()}"
        )
        temporary.write_text(
            json.dumps(
                {"identity": identity, "roots": root_records},
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        temporary.replace(stage_marker)
    return root_input, metadata


def herwig_stage_matches(
    marker: Path,
    identity: dict[str, object],
    roots: Sequence[Path],
) -> bool:
    try:
        payload = json.loads(marker.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    if not isinstance(payload, dict) or payload.get("identity") != identity:
        return False
    recorded = payload.get("roots")
    if not isinstance(recorded, list) or len(recorded) != len(roots):
        return False
    actual = [
        {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
        }
        for path in roots
    ]
    return recorded == actual


def card_metadata(path: Path, cfg: Config) -> dict[str, str | None]:
    return {
        "path": str(path),
        "sha256": None if cfg.dry_run or not path.exists() else sha256_file(path),
    }


def read_analysis_summary(path: Path) -> dict[str, object]:
    try:
        summary = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read analyzer summary {path}: {error}") from error
    for key in ("generated_raw", "generated_scaled", "selected_scaled"):
        value = summary.get(key)
        if not isinstance(value, dict):
            raise ValueError(f"analyzer summary {path} is missing object {key}")
    return summary


def run_sample(
    sample: SampleSpec,
    selected: Sequence[SampleSpec],
    manifest: dict[str, object],
    cfg: Config,
) -> None:
    entry = manifest_sample(manifest, sample.name)
    output_root = expected_analysis_output(sample, cfg)
    summary_path = output_root.with_suffix(".summary.json")
    if (
        not cfg.force
        and entry.get("status") == "complete"
        and output_root.exists()
    ):
        if not summary_path.exists():
            die(
                f"completed sample {sample.name} lacks {summary_path}; "
                "rerun with --force"
            )
        try:
            resumed_summary = read_analysis_summary(summary_path)
        except ValueError as error:
            die(str(error))
        if int(resumed_summary.get("events_read", -1)) != cfg.events_for(sample):
            die(
                f"completed sample {sample.name} has an unexpected event count; "
                "rerun with --force"
            )
        recorded_checksum = entry.get("output_root_sha256")
        if recorded_checksum and sha256_file(output_root) != recorded_checksum:
            die(
                f"completed sample {sample.name} ROOT checksum changed; "
                "rerun with --force"
            )
        recorded_summary_checksum = entry.get("analysis_summary_sha256")
        if (
            recorded_summary_checksum
            and sha256_file(summary_path) != recorded_summary_checksum
        ):
            die(
                f"completed sample {sample.name} analyzer-summary checksum "
                "changed; rerun with --force"
            )
        for label, path_key, checksum_key in (
            ("source LHE", "source_lhe", "source_lhe_sha256"),
            ("input LHE", "input_lhe", "input_lhe_sha256"),
        ):
            recorded_path = entry.get(path_key)
            recorded_lhe_checksum = entry.get(checksum_key)
            if recorded_path and recorded_lhe_checksum:
                path = Path(str(recorded_path))
                if not path.exists():
                    die(
                        f"completed sample {sample.name} {label} is missing; "
                        "rerun with --force"
                    )
                if sha256_file(path) != recorded_lhe_checksum:
                    die(
                        f"completed sample {sample.name} {label} checksum "
                        "changed; rerun with --force"
                    )
        log(f"Resuming: {sample.name} is already complete")
        return

    log(f"Running {sample.name}")
    entry["status"] = "running"
    write_manifest(manifest, cfg)

    stable_or_final_lhe = obtain_production_lhe(sample, selected, cfg)
    source_lhe_metadata: dict[str, object] | None = None
    if not cfg.dry_run:
        try:
            source_lhe_metadata = validate_lhe_collider(
                stable_or_final_lhe,
                cfg,
            )
        except ValueError as error:
            die(str(error))
    madspin_card: Path | None = None
    madspin_banner_augmentation: dict[str, object] | None = None
    if sample.is_signal:
        if not cfg.dry_run:
            try:
                _mass, width = validate_stable_higgs_lhe(
                    stable_or_final_lhe,
                    cfg,
                    cfg.events_for(sample),
                )
            except ValueError as error:
                die(str(error))
            entry["generator_higgs_width_gev"] = width
        (
            analysis_lhe,
            madspin_card,
            madspin_banner_augmentation,
        ) = run_signal_madspin(
            sample,
            stable_or_final_lhe,
            cfg,
        )
    else:
        analysis_lhe = stable_or_final_lhe
    input_lhe_metadata: dict[str, object] | None = None
    if not cfg.dry_run:
        try:
            input_lhe_metadata = validate_lhe_collider(analysis_lhe, cfg)
            if sample.prompt_decay_profile:
                entry["prompt_parent_validation"] = validate_prompt_decay_parents(
                    analysis_lhe,
                    sample,
                    cfg.events_for(sample),
                )
        except ValueError as error:
            die(str(error))

    if cfg.dry_run:
        normalization: dict[str, object] = dict(entry["normalization"])  # type: ignore[arg-type]
        weight_scale: float | str = (
            PENDING_SIGNAL_WEIGHT if sample.is_signal else 1.0
        )
    else:
        try:
            normalization = normalization_from_lhe(
                sample,
                analysis_lhe,
                cfg,
                stable_signal_lhe=stable_or_final_lhe if sample.is_signal else None,
            )
        except ValueError as error:
            die(str(error))
        weight_scale_value = normalization["weight_scale"]
        assert isinstance(weight_scale_value, float)
        weight_scale = weight_scale_value
        entry["normalization"] = {
            **entry["normalization"],  # type: ignore[arg-type]
            **normalization,
        }

    if not cfg.dry_run:
        try:
            sign_summary = predecay_lhe_weight_summary(
                analysis_lhe,
                cfg.events_for(sample),
            )
        except ValueError as error:
            die(str(error))
        entry["normalization"] = {
            **entry["normalization"],  # type: ignore[arg-type]
            **sign_summary,
        }

    optional_weight_names: Path | None = None
    optional_weight_definitions: list[dict[str, str]] | None = None
    if cfg.dry_run:
        entry["optional_weight_definitions_status"] = "pending_dry_run"
    else:
        try:
            optional_weight_definitions = (
                parse_lhe_optional_weight_definitions(analysis_lhe)
            )
        except ValueError as error:
            die(str(error))
        if optional_weight_definitions:
            optional_weight_names = (
                sample_work_dir(sample, cfg)
                / "lhe_optional_weight_names.txt"
            )
            optional_weight_names.parent.mkdir(parents=True, exist_ok=True)
            optional_weight_names.write_text(
                "".join(
                    f"{definition['label']}\n"
                    for definition in optional_weight_definitions
                )
            )
        entry["optional_weight_definitions_status"] = (
            "available" if optional_weight_definitions else "not_present"
        )
    entry["optional_weight_definitions"] = optional_weight_definitions
    entry["optional_weight_names_file"] = (
        card_metadata(optional_weight_names, cfg)
        if optional_weight_names is not None
        else None
    )

    root_input, decay_metadata = run_herwig(sample, analysis_lhe, cfg)
    output_dir = output_root.parent
    ensure_dir(output_dir, cfg)
    command = analysis_command(
        sample,
        root_input,
        seed_for(f"{sample.name}:analysis", cfg),
        weight_scale,
        cfg,
        optional_weight_names,
    )
    reuse_analysis = False
    if (
        not cfg.dry_run
        and not cfg.force
        and output_root.exists()
        and summary_path.exists()
    ):
        try:
            partial_summary = read_analysis_summary(summary_path)
            reuse_analysis = (
                int(partial_summary.get("events_read", -1))
                == cfg.events_for(sample)
            )
        except ValueError:
            reuse_analysis = False
    if not reuse_analysis:
        run_runtime_command(command, cfg)
    if not cfg.dry_run and not output_root.exists():
        die(f"analyzer did not create expected output: {output_root}")
    analysis_summary: dict[str, object] | None = None
    if not cfg.dry_run:
        if not summary_path.exists():
            die(f"analyzer did not create expected summary: {summary_path}")
        try:
            analysis_summary = read_analysis_summary(summary_path)
        except ValueError as error:
            die(str(error))
        if int(analysis_summary.get("events_read", -1)) != cfg.events_for(sample):
            die(
                f"analyzer read {analysis_summary.get('events_read')} events for "
                f"{sample.name}, expected {cfg.events_for(sample)}"
            )

    process_card = production_dir(sample, cfg) / "cards" / "process_card.dat"
    generated_process_dir = production_dir(sample, cfg) / "mg5_process"
    herwig_card = (
        sample_work_dir(sample, cfg) / "herwig" / f"{sample.name}.in"
    )
    entry["cards"] = {
        "process": card_metadata(process_card, cfg),
        "madspin": card_metadata(madspin_card, cfg) if madspin_card else None,
        "run": card_metadata(
            generated_process_dir / "Cards" / "run_card.dat",
            cfg,
        ),
        "parameters": card_metadata(
            generated_process_dir / "Cards" / "param_card.dat",
            cfg,
        ),
        "herwig": card_metadata(herwig_card, cfg),
    }
    entry["input_lhe"] = str(analysis_lhe)
    entry["input_lhe_init"] = input_lhe_metadata
    entry["input_lhe_sha256"] = (
        None if cfg.dry_run else sha256_file(analysis_lhe)
    )
    entry["input_lhe_checksum_status"] = (
        "pending_dry_run" if cfg.dry_run else "available"
    )
    entry["source_lhe"] = str(stable_or_final_lhe)
    entry["source_lhe_init"] = source_lhe_metadata
    entry["source_lhe_sha256"] = (
        None if cfg.dry_run else sha256_file(stable_or_final_lhe)
    )
    entry["source_lhe_checksum_status"] = (
        "pending_dry_run" if cfg.dry_run else "available"
    )
    entry["madspin_banner_augmentation"] = madspin_banner_augmentation
    entry["output_root_sha256"] = (
        None if cfg.dry_run else sha256_file(output_root)
    )
    entry["output_root_checksum_status"] = (
        "pending_dry_run" if cfg.dry_run else "available"
    )
    entry["analysis_summary_sha256"] = (
        None if cfg.dry_run else sha256_file(summary_path)
    )
    if analysis_summary is not None:
        entry["weight_sums"] = {
            "generated_raw": analysis_summary["generated_raw"],
            "generated_scaled": analysis_summary["generated_scaled"],
            "selected_scaled": analysis_summary["selected_scaled"],
        }
        normalization_entry = entry["normalization"]
        assert isinstance(normalization_entry, dict)
        generated_raw = analysis_summary["generated_raw"]
        assert isinstance(generated_raw, dict)
        normalization_entry.update(
            {
                # These are post-Herwig raw event-weight sums.  They must stay
                # distinct from the pre-decay LHE denominator whenever a
                # forced-decay importance weight is present.
                "generated_sumw": generated_raw["sumw"],
                "generated_sumabsw": generated_raw["sumabsw"],
                "generated_sumw2": generated_raw["sumw2"],
                "generated_npositive": generated_raw["positive_entries"],
                "generated_nnegative": generated_raw["negative_entries"],
                "generated_nzero": generated_raw["zero_entries"],
                "generated_sumw_positive": generated_raw["sumw_positive"],
                "generated_sumw_negative": generated_raw["sumw_negative"],
                "generated_negative_fraction": generated_raw[
                    "negative_fraction"
                ],
            }
        )
    entry["prompt_decay_selection"] = decay_metadata["prompt_selection"]
    entry["decay_bias"] = decay_metadata["decay_bias"]
    entry["branching_ratio_reweighter"] = decay_metadata[
        "branching_ratio_reweighter"
    ]
    entry["decay"] = {
        "prompt_selection": decay_metadata["prompt_selection"],
        "branching_ratio_reweighter_enabled": decay_metadata[
            "branching_ratio_reweighter"
        ],
        "heavy_flavour_bias_enabled": bool(
            isinstance(decay_metadata["decay_bias"], dict)
            and decay_metadata["decay_bias"].get("enabled")
        ),
        "handlers": decay_metadata["handlers"],
    }
    entry["analyzer_cli"] = [str(value) for value in command]
    entry["status"] = "planned" if cfg.dry_run else "complete"
    write_manifest(manifest, cfg)


def merge_manifest_samples(
    existing: dict[str, object] | None,
    planned: dict[str, object],
) -> dict[str, object]:
    if existing is None:
        return planned
    old_entries = {
        str(entry["name"]): entry
        for entry in existing.get("samples", [])
        if isinstance(entry, dict) and "name" in entry
    }
    planned_entries = planned["samples"]
    assert isinstance(planned_entries, list)
    merged_entries: list[dict[str, object]] = []
    for planned_entry in planned_entries:
        assert isinstance(planned_entry, dict)
        name = str(planned_entry["name"])
        merged_entries.append(old_entries.get(name, planned_entry))
    planned["samples"] = merged_entries
    planned["status"] = existing.get("status", "running")
    return planned


def main(argv: Sequence[str]) -> int:
    cfg = parse_config(argv)
    samples = selected_samples(cfg)
    if not samples:
        die(
            f"no samples selected by sample-set={cfg.sample_set} "
            f"run-samples={cfg.run_samples}"
        )
    require_inputs(cfg, samples)

    planned = build_manifest(cfg, samples)
    existing = load_resumable_manifest(cfg)
    if (
        existing is not None
        and existing.get("configuration_fingerprint")
        != planned["configuration_fingerprint"]
    ):
        die(
            "existing manifest configuration differs; use a new --run-tag "
            "or explicitly rerun with --force"
        )
    manifest = merge_manifest_samples(existing, planned)
    manifest["status"] = "planned" if cfg.dry_run else "running"
    write_manifest(manifest, cfg)

    log(
        f"Building four-lepton analyzer "
        f"(profile={cfg.detector_response}, target={cfg.analysis_target})"
    )
    run_runtime_command(
        ["make", "-C", cfg.analysis_code_dir, cfg.analysis_target],
        cfg,
    )
    software = manifest["software"]
    assert isinstance(software, dict)
    software["analyzer_executable"] = file_record(
        cfg.analysis_exe,
        dry_run=cfg.dry_run,
    )
    if any(
        sample.prompt_decay_profile or sample.is_hf_biased for sample in samples
    ):
        log("Building LHE forced-decay branching-ratio reweighter")
        build_decay_reweighter(cfg)
        software["decay_reweighter_plugin"] = file_record(
            cfg.decay_reweighter_plugin,
            dry_run=cfg.dry_run,
        )
    write_manifest(manifest, cfg)

    for sample in samples:
        run_sample(sample, samples, manifest, cfg)

    manifest["status"] = "planned" if cfg.dry_run else "complete"
    write_manifest(manifest, cfg)
    if cfg.dry_run:
        print(json.dumps(manifest, indent=2, sort_keys=True))
    else:
        print(f"\nManifest: {cfg.manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
