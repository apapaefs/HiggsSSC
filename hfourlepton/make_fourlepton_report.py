#!/usr/bin/env python3
"""Build the cut-based report for the SSC H -> ZZ* -> four-lepton analysis.

The module deliberately keeps ROOT at the I/O boundary.  All bookkeeping,
normalisation, cut-flow, and response-closure functions operate on ordinary
Python mappings and are therefore testable on machines without ROOT.

Normalisation convention
------------------------
The campaign's ``normalization.denominator_sumw`` is the preferred
pre-decay hard-event denominator.  It is mandatory when Herwig's branching
ratio reweighter is active; otherwise dividing by the reweighted
``SourceEvents.generator_weight`` would cancel the forced-decay branching
factor.  Unforced samples may fall back to the source-tree generator sum.
A selected ``FourLepton.event_weight`` contains the generator, sample,
detector-response, and trigger weights.  Consequently

    selected cross section =
        native hard-process cross section
        * sum(selected event_weight)
        / pre-decay denominator sumw.

No entry-count fallback is made when the signed denominator is zero.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import html
import json
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_LUMINOSITY_FB = 100.0
DEFAULT_SIGNAL_REGION_GEV = (120.0, 130.0)
DEFAULT_MASS_RANGE_GEV = (70.0, 200.0)
DEFAULT_MASS_BIN_WIDTH_GEV = 5.0
REPORT_SCHEMA_VERSION = 1
SUPPORTED_MANIFEST_MAJOR = 1
CHANNEL_CODES = {0: "unknown", 1: "4e", 2: "4mu", 3: "2e2mu"}


@dataclass(frozen=True)
class CutStage:
    """One cumulative cut-flow stage."""

    name: str
    label: str
    bit: int
    required_mask: int


def _cumulative_stage(name: str, label: str, bit: int) -> CutStage:
    return CutStage(name=name, label=label, bit=bit, required_mask=(1 << (bit + 1)) - 1)


# This is the public analyzer/report contract.  A manifest may replace the
# labels or masks through ``analysis.cut_stages``.
DEFAULT_CUT_STAGES: tuple[CutStage, ...] = (
    _cumulative_stage("four_reconstructed", "At least four reconstructed leptons", 0),
    _cumulative_stage("fiducial", "Lepton pT and eta", 1),
    _cumulative_stage("electron_crack", "Electron crack veto", 2),
    _cumulative_stage("isolation", "Lepton isolation", 3),
    _cumulative_stage("ossf_pairs", "Two disjoint OSSF pairs", 4),
    _cumulative_stage("low_mass_veto", "All OSSF masses above 4 GeV", 5),
    _cumulative_stage("z1_window", "70 < mZ1 < 100 GeV", 6),
    _cumulative_stage("z2_window", "10 < mZ2 < 100 GeV", 7),
    _cumulative_stage("trigger", "Trigger efficiency", 8),
    _cumulative_stage("selected", "Final selection", 9),
)


@dataclass
class Audit:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def extend(self, other: "Audit", prefix: str = "") -> None:
        self.errors.extend(f"{prefix}{message}" for message in other.errors)
        self.warnings.extend(f"{prefix}{message}" for message in other.warnings)

    def as_dict(self) -> dict[str, Any]:
        return {
            "ok": not self.errors,
            "errors": list(self.errors),
            "warnings": list(self.warnings),
        }


@dataclass
class WeightStats:
    """Signed-weight summary.

    ``entries`` counts independent source-event contributions.  When several
    response hypotheses originate from one source event, callers should first
    group the weights by source and then add that one combined contribution.
    """

    entries: int = 0
    positive_entries: int = 0
    negative_entries: int = 0
    zero_entries: int = 0
    sumw: float = 0.0
    sumabsw: float = 0.0
    sumw2: float = 0.0
    sumw_positive: float = 0.0
    sumw_negative: float = 0.0

    def add(self, weight: float) -> None:
        value = finite_float(weight, "weight")
        self.entries += 1
        self.sumw += value
        self.sumabsw += abs(value)
        self.sumw2 += value * value
        if value > 0.0:
            self.positive_entries += 1
            self.sumw_positive += value
        elif value < 0.0:
            self.negative_entries += 1
            self.sumw_negative += value
        else:
            self.zero_entries += 1

    def merge(self, other: "WeightStats") -> None:
        self.entries += other.entries
        self.positive_entries += other.positive_entries
        self.negative_entries += other.negative_entries
        self.zero_entries += other.zero_entries
        self.sumw += other.sumw
        self.sumabsw += other.sumabsw
        self.sumw2 += other.sumw2
        self.sumw_positive += other.sumw_positive
        self.sumw_negative += other.sumw_negative

    def scaled(self, factor: float) -> "WeightStats":
        scale = finite_float(factor, "scale")
        positive = self.sumw_positive * scale
        negative = self.sumw_negative * scale
        if scale < 0.0:
            positive, negative = negative, positive
        if scale > 0.0:
            positive_entries = self.positive_entries
            negative_entries = self.negative_entries
            zero_entries = self.zero_entries
        elif scale < 0.0:
            positive_entries = self.negative_entries
            negative_entries = self.positive_entries
            zero_entries = self.zero_entries
        else:
            positive_entries = 0
            negative_entries = 0
            zero_entries = self.entries
        return WeightStats(
            entries=self.entries,
            positive_entries=positive_entries,
            negative_entries=negative_entries,
            zero_entries=zero_entries,
            sumw=self.sumw * scale,
            sumabsw=self.sumabsw * abs(scale),
            sumw2=self.sumw2 * scale * scale,
            sumw_positive=positive,
            sumw_negative=negative,
        )

    def as_dict(self) -> dict[str, Any]:
        negative_abs = abs(self.sumw_negative)
        signed_components = self.sumw_positive + negative_abs
        return {
            "entries": self.entries,
            "positive_entries": self.positive_entries,
            "negative_entries": self.negative_entries,
            "zero_entries": self.zero_entries,
            "sumw": self.sumw,
            "sumabsw": self.sumabsw,
            "sumw2": self.sumw2,
            "stat_uncertainty": math.sqrt(max(0.0, self.sumw2)),
            "sumw_positive": self.sumw_positive,
            "sumw_negative": self.sumw_negative,
            "negative_weight_fraction": (
                self.negative_entries / self.entries
                if self.entries != 0
                else None
            ),
            "negative_absolute_weight_fraction": (
                negative_abs / signed_components if signed_components != 0.0 else None
            ),
        }


@dataclass
class RootRecords:
    source_events: list[dict[str, Any]]
    four_lepton: list[dict[str, Any]]
    audit: Audit = field(default_factory=Audit)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class ManifestSample:
    name: str
    label: str
    category: str
    channel: str
    detector_profile: str | None
    output_root: Path | None
    native_cross_section_pb: float | None
    target_cross_section_pb: float | None
    weight_scale: float | None
    generated_sumw: float | None
    generated_sumw2: float | None
    denominator_sumw: float | None
    denominator_sumw2: float | None
    perturbative_order: str | None
    source_kind: str | None
    production_group_id: str | None
    detector_settings: dict[str, Any]
    decay_bias_enabled: bool
    decay_bias_scheme: str | None
    decay_bias_closure_validated: bool | None
    decay_bias_of: str | None
    branching_ratio_reweighter: bool
    prompt_decay_profile: str | None
    manifest_path: Path
    raw: dict[str, Any]
    audit: Audit = field(default_factory=Audit)


@dataclass
class ManifestBundle:
    path: Path
    schema_version: str | int | float | None
    run_tag: str | None
    sqrt_s_tev: float | None
    detector_profile: str | None
    luminosity_fb: float | None
    cut_stages: tuple[CutStage, ...]
    samples: list[ManifestSample]
    audit: Audit = field(default_factory=Audit)


def finite_float(value: Any, name: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return parsed


def optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return finite_float(value, "value")
    except (TypeError, ValueError):
        return None


def optional_bool(value: Any) -> bool | None:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "yes", "on", "enabled", "validated"}:
        return True
    if text in {"false", "no", "off", "disabled", "unvalidated"}:
        return False
    return None


def first_present(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping and mapping[key] is not None:
            return mapping[key]
    return None


def nested_mapping(mapping: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = mapping.get(key, {})
    return value if isinstance(value, Mapping) else {}


def normalize_channel(value: Any, fallback: str = "unknown") -> str:
    if value is None:
        return fallback
    if isinstance(value, bool):
        return fallback
    if isinstance(value, (int, float)):
        code = int(value)
        if code == 0:
            return fallback
        return CHANNEL_CODES.get(code, f"code_{code}")
    text = str(value).strip().lower().replace("μ", "mu")
    aliases = {
        "eeee": "4e",
        "4electron": "4e",
        "4electrons": "4e",
        "mumumumu": "4mu",
        "4m": "4mu",
        "eemumu": "2e2mu",
        "2e2m": "2e2mu",
        "inclusive": "inclusive",
        "all": "inclusive",
    }
    return aliases.get(text, text or fallback)


def parse_manifest_major(version: Any) -> int | None:
    if version is None:
        return None
    try:
        return int(str(version).split(".", maxsplit=1)[0])
    except ValueError:
        return None


def parse_cut_stages(manifest: Mapping[str, Any], audit: Audit) -> tuple[CutStage, ...]:
    analysis = nested_mapping(manifest, "analysis")
    raw_stages = first_present(analysis, "cut_stages", "cut_definitions")
    if raw_stages is None:
        raw_stages = first_present(manifest, "cut_stages", "cut_definitions")
    if not isinstance(raw_stages, Sequence) or isinstance(raw_stages, (str, bytes)):
        return DEFAULT_CUT_STAGES

    stages: list[CutStage] = []
    for index, raw in enumerate(raw_stages):
        if not isinstance(raw, Mapping):
            audit.warnings.append(f"cut stage {index} is not an object; using default stages")
            return DEFAULT_CUT_STAGES
        name = str(raw.get("name", f"cut_{index}"))
        label = str(raw.get("label", name.replace("_", " ")))
        bit_value = raw.get("bit", index)
        try:
            bit = int(bit_value)
            required_mask = int(raw.get("required_mask", (1 << (bit + 1)) - 1))
        except (TypeError, ValueError):
            audit.warnings.append(f"invalid cut stage {name!r}; using default stages")
            return DEFAULT_CUT_STAGES
        if bit < 0:
            audit.warnings.append(f"negative cut bit for {name!r}; using default stages")
            return DEFAULT_CUT_STAGES
        stages.append(CutStage(name=name, label=label, bit=bit, required_mask=required_mask))
    return tuple(stages) if stages else DEFAULT_CUT_STAGES


def _sample_list(raw: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    samples = raw.get("samples")
    if isinstance(samples, Mapping):
        converted: list[Mapping[str, Any]] = []
        for name, value in samples.items():
            if isinstance(value, Mapping):
                item = dict(value)
                item.setdefault("name", name)
                converted.append(item)
        return converted
    if isinstance(samples, Sequence) and not isinstance(samples, (str, bytes)):
        return [value for value in samples if isinstance(value, Mapping)]
    sample = raw.get("sample")
    if isinstance(sample, Mapping):
        return [sample]
    if any(key in raw for key in ("name", "sample_name", "output_root", "root_file")):
        return [raw]
    return []


def _output_root(sample: Mapping[str, Any], manifest_dir: Path) -> Path | None:
    outputs = nested_mapping(sample, "outputs")
    paths = nested_mapping(sample, "paths")
    value = first_present(
        sample,
        "output_root",
        "root_file",
        "analysis_root",
        "four_lepton_root",
    )
    if value is None:
        value = first_present(outputs, "root", "root_file", "analysis_root", "four_lepton")
    if value is None:
        value = first_present(paths, "output_root", "root_file", "analysis_root")
    if value is None:
        return None
    path = Path(os.path.expandvars(os.path.expanduser(str(value))))
    return path if path.is_absolute() else (manifest_dir / path).resolve()


def load_manifest(path: Path) -> ManifestBundle:
    path = path.resolve()
    raw_value = json.loads(path.read_text())
    if not isinstance(raw_value, Mapping):
        raise ValueError(f"{path}: manifest root must be a JSON object")
    raw = dict(raw_value)
    audit = Audit()
    version = first_present(raw, "schema_version", "manifest_version", "version")
    major = parse_manifest_major(version)
    if major is None:
        audit.errors.append("missing or invalid schema_version")
    elif major != SUPPORTED_MANIFEST_MAJOR:
        audit.errors.append(
            f"unsupported manifest schema major {major}; expected {SUPPORTED_MANIFEST_MAJOR}"
        )

    collider = nested_mapping(raw, "collider")
    top_profile = first_present(raw, "detector_profile", "detector_response", "response_profile")
    top_detector = dict(nested_mapping(raw, "detector"))
    top_decay_bias = nested_mapping(raw, "decay_bias")
    luminosity = optional_float(first_present(raw, "luminosity_fb", "integrated_luminosity_fb"))
    samples: list[ManifestSample] = []
    raw_samples = _sample_list(raw)
    if not raw_samples:
        audit.errors.append("manifest contains no samples")

    for index, raw_sample_value in enumerate(raw_samples):
        raw_sample = dict(raw_sample_value)
        sample_audit = Audit()
        name_value = first_present(raw_sample, "name", "sample_name", "id")
        name = str(name_value) if name_value is not None else f"unnamed_{index}"
        if name_value is None:
            sample_audit.errors.append("missing sample name")

        category_value = first_present(raw_sample, "category", "sample_category", "kind")
        category = str(category_value).strip().lower() if category_value is not None else "unknown"
        category_aliases = {
            "background": "irreducible",
            "bkg": "irreducible",
            "irreducible_background": "irreducible",
            "reducible_background": "reducible",
        }
        category = category_aliases.get(category, category)
        if category_value is None:
            sample_audit.errors.append("missing sample category")

        channel_value = first_present(raw_sample, "channel", "final_state")
        channel = normalize_channel(channel_value, "unknown")
        if channel_value is None:
            sample_audit.errors.append("missing sample channel")

        normalization = nested_mapping(raw_sample, "normalization")
        native = optional_float(
            first_present(
                normalization,
                "native_cross_section_pb",
                "cross_section_pb",
                "generator_cross_section_pb",
            )
        )
        if native is None:
            native = optional_float(
                first_present(
                    raw_sample,
                    "native_cross_section_pb",
                    "cross_section_pb",
                    "generator_cross_section_pb",
                )
            )
        target = optional_float(
            first_present(normalization, "target_cross_section_pb", "effective_cross_section_pb")
        )
        if target is None:
            target = optional_float(
                first_present(raw_sample, "target_cross_section_pb", "effective_cross_section_pb")
            )
        weight_scale = optional_float(
            first_present(normalization, "weight_scale", "sample_weight_scale", "branching_weight")
        )
        if weight_scale is None:
            weight_scale = optional_float(
                first_present(raw_sample, "weight_scale", "sample_weight_scale", "branching_weight")
            )
        if native is None and target is not None and weight_scale not in (None, 0.0):
            native = target / weight_scale
            sample_audit.warnings.append(
                "native cross section derived from target_cross_section_pb / weight_scale"
            )
        if target is None and native is not None and weight_scale is not None:
            target = native * weight_scale
        if native is None:
            sample_audit.errors.append(
                "missing native cross section; selected yields cannot be normalized"
            )

        generated_sumw = optional_float(
            first_present(normalization, "generated_sumw", "sumw", "sum_generator_weight")
        )
        generated_sumw2 = optional_float(
            first_present(normalization, "generated_sumw2", "sumw2", "sum_generator_weight2")
        )
        denominator_sumw = optional_float(
            first_present(
                normalization,
                "denominator_sumw",
                "pre_decay_sumw",
                "predecay_sumw",
                "lhe_sum_sign",
                "sum_sign",
            )
        )
        denominator_sumw2 = optional_float(
            first_present(
                normalization,
                "denominator_sumw2",
                "pre_decay_sumw2",
                "predecay_sumw2",
                "lhe_sum_sign2",
                "sum_sign2",
            )
        )
        profile = first_present(
            raw_sample, "detector_profile", "detector_response", "response_profile"
        )
        profile_text = str(profile if profile is not None else top_profile).strip() or None
        if profile_text is None:
            sample_audit.errors.append("missing detector response profile")
        root_path = _output_root(raw_sample, path.parent)
        if root_path is None:
            sample_audit.errors.append("missing ROOT output path")
        sample_decay_bias = nested_mapping(raw_sample, "decay_bias")
        effective_decay_bias = sample_decay_bias if sample_decay_bias else top_decay_bias
        detector_settings = dict(top_detector)
        detector_settings.update(nested_mapping(raw_sample, "detector"))
        bias_enabled = optional_bool(effective_decay_bias.get("enabled"))
        bias_scheme = first_present(effective_decay_bias, "scheme", "method")
        bias_of = first_present(
            effective_decay_bias,
            "bias_of",
            "control_sample",
            "unbiased_sample",
        )
        bias_validated = optional_bool(
            first_present(
                effective_decay_bias,
                "closure_validated",
                "validated",
                "weighted_closure_validated",
            )
        )
        decay_metadata = nested_mapping(raw_sample, "decay_metadata")
        decay_config = nested_mapping(raw_sample, "decay")
        reweighter = optional_bool(
            first_present(
                raw_sample,
                "branching_ratio_reweighter",
                "br_reweighter_enabled",
            )
        )
        if reweighter is None:
            reweighter = optional_bool(
                first_present(
                    decay_metadata,
                    "branching_ratio_reweighter",
                    "br_reweighter_enabled",
                )
            )
        if reweighter is None:
            reweighter = optional_bool(
                first_present(
                    decay_config,
                    "branching_ratio_reweighter",
                    "br_reweighter_enabled",
                )
            )
        prompt_profile = first_present(
            raw_sample, "prompt_decay_profile", "prompt_decay_selection"
        )
        if prompt_profile is None:
            prompt_profile = first_present(
                decay_metadata, "prompt_decay_profile", "prompt_selection"
            )
        if prompt_profile is None:
            prompt_profile = first_present(
                decay_config, "prompt_decay_profile", "prompt_selection"
            )
        # Heavy-flavour bias and any forced prompt profile are implemented by
        # the same Herwig BR reweighter in the campaign.
        br_reweighter = bool(reweighter) or bool(bias_enabled) or bool(prompt_profile)

        samples.append(
            ManifestSample(
                name=name,
                label=str(raw_sample.get("label", name.replace("_", " "))),
                category=category,
                channel=channel,
                detector_profile=profile_text,
                output_root=root_path,
                native_cross_section_pb=native,
                target_cross_section_pb=target,
                weight_scale=weight_scale,
                generated_sumw=generated_sumw,
                generated_sumw2=generated_sumw2,
                denominator_sumw=denominator_sumw,
                denominator_sumw2=denominator_sumw2,
                perturbative_order=(
                    str(raw_sample["perturbative_order"])
                    if raw_sample.get("perturbative_order") is not None
                    else None
                ),
                source_kind=(
                    str(raw_sample["source_kind"])
                    if raw_sample.get("source_kind") is not None
                    else None
                ),
                production_group_id=(
                    str(raw_sample["production_group_id"])
                    if raw_sample.get("production_group_id") is not None
                    else None
                ),
                detector_settings=detector_settings,
                decay_bias_enabled=bool(bias_enabled),
                decay_bias_scheme=str(bias_scheme) if bias_scheme is not None else None,
                decay_bias_closure_validated=bias_validated,
                decay_bias_of=str(bias_of) if bias_of is not None else None,
                branching_ratio_reweighter=br_reweighter,
                prompt_decay_profile=(
                    str(prompt_profile) if prompt_profile is not None else None
                ),
                manifest_path=path,
                raw=raw_sample,
                audit=sample_audit,
            )
        )

    return ManifestBundle(
        path=path,
        schema_version=version,
        run_tag=str(raw["run_tag"]) if raw.get("run_tag") is not None else None,
        sqrt_s_tev=optional_float(
            first_present(collider, "sqrt_s_tev", "energy_tev", "centre_of_mass_energy_tev")
        ),
        detector_profile=str(top_profile) if top_profile is not None else None,
        luminosity_fb=luminosity,
        cut_stages=parse_cut_stages(raw, audit),
        samples=samples,
        audit=audit,
    )


SOURCE_BRANCH_ALIASES: dict[str, tuple[str, ...]] = {
    "source_index": ("source_index", "source_event", "event_index"),
    "generator_weight": ("generator_weight", "nominal_weight", "evweight"),
    "sample_weight": ("sample_weight", "sample_weight_scale"),
    "base_event_weight": ("base_event_weight",),
    "response_failure_probability": ("response_failure_probability", "failure_probability"),
    "response_probability_in_tree": ("response_probability_in_tree", "pass_probability"),
    "selected_probability": ("selected_probability",),
    "selected_event_weight": ("selected_event_weight",),
    "response_closure_delta": ("response_closure_delta",),
    "raw_lepton_count": ("raw_lepton_count",),
    "accepted_lepton_count": ("accepted_lepton_count",),
    "subset_hypothesis_count": ("subset_hypothesis_count",),
    "cut_probabilities": ("cut_probabilities",),
    "optional_weight_names": ("optional_weight_names", "optional_weights_names"),
    "optional_weights": ("optional_weights", "variation_weights"),
    "optional_selected_event_weights": (
        "optional_selected_event_weights",
        "selected_optional_weights",
    ),
}

FOUR_LEPTON_BRANCH_ALIASES: dict[str, tuple[str, ...]] = {
    "source_index": ("source_index", "source_event", "event_index"),
    "hypothesis_probability": ("hypothesis_probability", "hypothesis_weight"),
    "generator_weight": ("generator_weight", "nominal_weight"),
    "sample_weight": ("sample_weight", "sample_weight_scale"),
    "response_weight": ("response_weight",),
    "trigger_weight": ("trigger_weight",),
    "pretrigger_event_weight": ("pretrigger_event_weight",),
    "event_weight": ("event_weight", "weight"),
    "cut_mask": ("cut_mask", "selection_mask"),
    "channel": ("channel", "channel_code"),
    "m4l": ("m4l", "mass4l", "four_lepton_mass"),
    "mZ1": ("mZ1", "mz1", "z1_mass"),
    "mZ2": ("mZ2", "mz2", "z2_mass"),
    "pt4l": ("pt4l", "four_lepton_pt"),
    "y4l": ("y4l", "four_lepton_rapidity"),
    "min_delta_r": ("min_delta_r",),
    "n_jets": ("n_jets",),
    "n_bjets": ("n_bjets",),
    "met_pt": ("met_pt", "met"),
}

METADATA_BRANCH_ALIASES: dict[str, tuple[str, ...]] = {
    "schema_version": ("schema_version",),
    "response_profile": ("response_profile", "detector_profile", "detector_response"),
    "sample": ("sample", "sample_name"),
    "category": ("category",),
    "requested_channel": ("requested_channel", "channel"),
    "cut_mask_definition": ("cut_mask_definition",),
    "cut_bit_zero_semantics": ("cut_bit_zero_semantics",),
    "channel_definition": ("channel_definition",),
    "detector_parameter_source": ("detector_parameter_source",),
    "optional_weight_names_file": ("optional_weight_names_file",),
    "optional_weight_names_source": ("optional_weight_names_source",),
    "weight_scale": ("weight_scale", "sample_weight_scale"),
    "muon_resolution_scale": ("muon_resolution_scale",),
    "electron_efficiency": ("electron_efficiency",),
    "muon_efficiency": ("muon_efficiency",),
    "trigger_4e_efficiency": ("trigger_4e_efficiency",),
    "trigger_other_efficiency": ("trigger_other_efficiency",),
    "em_pileup_noise_enabled": ("em_pileup_noise_enabled",),
    "seed": ("seed",),
}

REQUIRED_REPORT_BRANCHES: dict[str, frozenset[str]] = {
    "SourceEvents": frozenset({"source_index", "generator_weight"}),
    "FourLepton": frozenset(
        {
            "source_index",
            "hypothesis_probability",
            "event_weight",
            "cut_mask",
            "channel",
            "m4l",
        }
    ),
}


def _python_value(value: Any) -> Any:
    if hasattr(value, "to_list"):
        value = value.to_list()
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return [_python_value(item) for item in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return value


def _chosen_branches(
    available: Iterable[str], aliases: Mapping[str, Sequence[str]]
) -> dict[str, str]:
    names = {str(name).split(";", maxsplit=1)[0] for name in available}
    chosen: dict[str, str] = {}
    for canonical, candidates in aliases.items():
        for candidate in candidates:
            if candidate in names:
                chosen[canonical] = candidate
                break
    return chosen


def _read_with_uproot(path: Path) -> RootRecords:
    import uproot  # type: ignore[import-not-found]

    audit = Audit()
    records: dict[str, list[dict[str, Any]]] = {"SourceEvents": [], "FourLepton": []}
    aliases_by_tree = {
        "SourceEvents": SOURCE_BRANCH_ALIASES,
        "FourLepton": FOUR_LEPTON_BRANCH_ALIASES,
    }
    metadata: dict[str, Any] = {}
    with uproot.open(path) as root_file:
        for tree_name, aliases in aliases_by_tree.items():
            if tree_name not in root_file:
                audit.errors.append(f"missing ROOT tree {tree_name}")
                continue
            tree = root_file[tree_name]
            chosen = _chosen_branches(tree.keys(), aliases)
            if not chosen:
                audit.errors.append(f"ROOT tree {tree_name} contains no recognized branches")
                continue
            missing = REQUIRED_REPORT_BRANCHES[tree_name] - chosen.keys()
            if missing:
                audit.errors.append(
                    f"ROOT tree {tree_name} lacks required report branches: "
                    + ", ".join(sorted(missing))
                )
            # Awkward handles both scalar leaves and the optional-weight
            # vectors without forcing jagged arrays into an object ndarray.
            arrays = tree.arrays(list(chosen.values()), library="ak")
            count = int(tree.num_entries)
            for index in range(count):
                record = {
                    canonical: _python_value(arrays[actual][index])
                    for canonical, actual in chosen.items()
                }
                records[tree_name].append(record)
        metadata_tree_present = "AnalysisMetadata" in root_file
        if metadata_tree_present:
            tree = root_file["AnalysisMetadata"]
            chosen = _chosen_branches(tree.keys(), METADATA_BRANCH_ALIASES)
            if chosen and int(tree.num_entries) > 0:
                arrays = tree.arrays(list(chosen.values()), library="np")
                metadata.update(
                    {
                        canonical: _python_value(arrays[actual][0])
                        for canonical, actual in chosen.items()
                    }
                )
        else:
            audit.warnings.append("missing optional ROOT tree AnalysisMetadata")
    metadata["analysis_metadata_present"] = metadata_tree_present
    metadata["analysis_metadata"] = {
        key: value
        for key, value in metadata.items()
        if key not in {"analysis_metadata_present", "analysis_metadata"}
    }
    return RootRecords(records["SourceEvents"], records["FourLepton"], audit, metadata)


def _read_with_pyroot(path: Path) -> RootRecords:
    import ROOT  # type: ignore[import-not-found]

    audit = Audit()
    output: dict[str, list[dict[str, Any]]] = {"SourceEvents": [], "FourLepton": []}
    aliases_by_tree = {
        "SourceEvents": SOURCE_BRANCH_ALIASES,
        "FourLepton": FOUR_LEPTON_BRANCH_ALIASES,
    }
    metadata: dict[str, Any] = {}
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"could not open ROOT file {path}")
    try:
        for tree_name, aliases in aliases_by_tree.items():
            tree = root_file.Get(tree_name)
            if not tree:
                audit.errors.append(f"missing ROOT tree {tree_name}")
                continue
            available = [branch.GetName() for branch in tree.GetListOfBranches()]
            chosen = _chosen_branches(available, aliases)
            if not chosen:
                audit.errors.append(f"ROOT tree {tree_name} contains no recognized branches")
                continue
            missing = REQUIRED_REPORT_BRANCHES[tree_name] - chosen.keys()
            if missing:
                audit.errors.append(
                    f"ROOT tree {tree_name} lacks required report branches: "
                    + ", ".join(sorted(missing))
                )
            for entry in tree:
                record: dict[str, Any] = {}
                for canonical, actual in chosen.items():
                    value = getattr(entry, actual)
                    if canonical == "cut_probabilities":
                        value = [float(value[i]) for i in range(10)]
                    elif canonical in {
                        "optional_weights",
                        "optional_selected_event_weights",
                    }:
                        value = list(value)
                    elif canonical == "optional_weight_names":
                        # Convert while the TFile is open. cppyy std::string
                        # proxies become invalid once their owning tree closes.
                        value = [str(item) for item in value]
                    record[canonical] = _python_value(value)
                output[tree_name].append(record)
        tree = root_file.Get("AnalysisMetadata")
        metadata_tree_present = bool(tree)
        if tree:
            available = [branch.GetName() for branch in tree.GetListOfBranches()]
            chosen = _chosen_branches(available, METADATA_BRANCH_ALIASES)
            if tree.GetEntries() > 0:
                tree.GetEntry(0)
                metadata.update(
                    {
                        canonical: (
                            str(getattr(tree, actual))
                            if canonical
                            in {
                                "response_profile",
                                "sample",
                                "category",
                                "requested_channel",
                                "cut_mask_definition",
                                "cut_bit_zero_semantics",
                                "channel_definition",
                                "detector_parameter_source",
                                "optional_weight_names_file",
                                "optional_weight_names_source",
                            }
                            else _python_value(getattr(tree, actual))
                        )
                        for canonical, actual in chosen.items()
                    }
                )
        else:
            audit.warnings.append("missing optional ROOT tree AnalysisMetadata")
    finally:
        root_file.Close()
    metadata["analysis_metadata_present"] = metadata_tree_present
    metadata["analysis_metadata"] = {
        key: value
        for key, value in metadata.items()
        if key not in {"analysis_metadata_present", "analysis_metadata"}
    }
    return RootRecords(output["SourceEvents"], output["FourLepton"], audit, metadata)


def _add_summary_sidecar(path: Path, records: RootRecords) -> RootRecords:
    summary_path = path.with_suffix(".summary.json")
    if not summary_path.exists():
        records.audit.warnings.append(f"analyzer summary sidecar not found: {summary_path}")
        return records
    try:
        summary = json.loads(summary_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        records.audit.warnings.append(f"could not read analyzer summary {summary_path}: {error}")
        return records
    if isinstance(summary, Mapping):
        records.metadata["summary"] = dict(summary)
        records.metadata["summary_path"] = str(summary_path)
        for key in (
            "schema_version",
            "response_profile",
            "sample",
            "category",
            "requested_channel",
            "weight_scale",
            "muon_resolution_scale",
            "seed",
        ):
            if records.metadata.get(key) in (None, "", b"") and summary.get(key) is not None:
                records.metadata[key] = summary[key]
    else:
        records.audit.warnings.append(f"analyzer summary is not a JSON object: {summary_path}")
    return records


def load_root_records(path: Path) -> RootRecords:
    """Load only the scalar report branches from a four-lepton ROOT output."""

    if not path.exists():
        return RootRecords([], [], Audit(errors=[f"ROOT output does not exist: {path}"]))
    try:
        return _add_summary_sidecar(path, _read_with_uproot(path))
    except ImportError:
        pass
    except Exception as error:
        uproot_error = error
    else:  # pragma: no cover - the return above is explicit
        uproot_error = None
    try:
        return _add_summary_sidecar(path, _read_with_pyroot(path))
    except ImportError:
        message = "neither uproot nor PyROOT is available"
        if "uproot_error" in locals():
            message += f"; uproot reader also failed: {uproot_error}"
        return RootRecords([], [], Audit(errors=[message]))
    except Exception as error:
        details = f"could not read {path}: {error}"
        if "uproot_error" in locals():
            details += f" (uproot: {uproot_error})"
        return RootRecords([], [], Audit(errors=[details]))


def source_index(record: Mapping[str, Any], fallback: int) -> Any:
    value = first_present(record, "source_index", "source_event", "event_index")
    return value if value is not None else ("row", fallback)


def record_event_weight(record: Mapping[str, Any], pretrigger: bool = False) -> float:
    if pretrigger:
        value = first_present(record, "pretrigger_event_weight")
        if value is not None:
            return finite_float(value, "pretrigger_event_weight")
    value = first_present(record, "event_weight", "weight")
    if value is not None:
        return finite_float(value, "event_weight")
    generator = optional_float(first_present(record, "generator_weight", "nominal_weight"))
    sample = optional_float(first_present(record, "sample_weight", "sample_weight_scale"))
    response = optional_float(record.get("response_weight"))
    hypothesis = optional_float(record.get("hypothesis_probability"))
    trigger = optional_float(record.get("trigger_weight"))
    result = generator if generator is not None else 1.0
    result *= sample if sample is not None else 1.0
    if response is not None:
        result *= response
    elif hypothesis is not None:
        result *= hypothesis
    if not pretrigger and trigger is not None:
        result *= trigger
    return result


def grouped_weight_stats(
    records: Iterable[Mapping[str, Any]],
    predicate: Callable[[Mapping[str, Any]], bool] | None = None,
    weight: Callable[[Mapping[str, Any]], float] = record_event_weight,
) -> WeightStats:
    """Sum response hypotheses by source before calculating ``sumw2``."""

    grouped: dict[Any, float] = defaultdict(float)
    for fallback, record in enumerate(records):
        if predicate is not None and not predicate(record):
            continue
        grouped[source_index(record, fallback)] += finite_float(weight(record), "record weight")
    stats = WeightStats()
    for value in grouped.values():
        stats.add(value)
    return stats


def source_weight_stats(source_events: Iterable[Mapping[str, Any]]) -> WeightStats:
    stats = WeightStats()
    seen: set[Any] = set()
    for fallback, record in enumerate(source_events):
        index = source_index(record, fallback)
        if index in seen:
            raise ValueError(f"duplicate SourceEvents source index {index!r}")
        seen.add(index)
        value = first_present(record, "generator_weight", "nominal_weight", "evweight")
        if value is None:
            raise ValueError(f"SourceEvents row {index!r} has no generator weight")
        stats.add(finite_float(value, "generator_weight"))
    return stats


def passes_stage(record: Mapping[str, Any], stage: CutStage) -> bool:
    try:
        cut_mask = int(record.get("cut_mask", 0))
    except (TypeError, ValueError):
        return False
    return (cut_mask & stage.required_mask) == stage.required_mask


def cutflow_stats(
    records: Sequence[Mapping[str, Any]], stages: Sequence[CutStage] = DEFAULT_CUT_STAGES
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stage in stages:
        # Trigger is bit 8 in the standard contract.  Earlier stages must not
        # be suppressed by an event-level trigger probability.
        use_pretrigger = stage.bit < 8
        stats = grouped_weight_stats(
            records,
            predicate=lambda record, current=stage: passes_stage(record, current),
            weight=lambda record, before_trigger=use_pretrigger: record_event_weight(
                record, pretrigger=before_trigger
            ),
        )
        row = {"name": stage.name, "label": stage.label, "bit": stage.bit}
        row.update(stats.as_dict())
        rows.append(row)
    return rows


def response_closure_audit(
    source_events: Sequence[Mapping[str, Any]],
    four_lepton: Sequence[Mapping[str, Any]],
    tolerance: float = 1.0e-9,
) -> dict[str, Any]:
    probability_by_source: dict[Any, float] = defaultdict(float)
    for fallback, record in enumerate(four_lepton):
        probability = optional_float(record.get("hypothesis_probability"))
        if probability is not None:
            probability_by_source[source_index(record, fallback)] += probability

    deltas: list[float] = []
    missing = 0
    direct_disagreements = 0
    for fallback, source in enumerate(source_events):
        index = source_index(source, fallback)
        failure = optional_float(source.get("response_failure_probability"))
        in_tree = optional_float(source.get("response_probability_in_tree"))
        row_probability = probability_by_source.get(index)
        direct_delta = optional_float(source.get("response_closure_delta"))

        if failure is not None and row_probability is not None:
            delta = failure + row_probability - 1.0
        elif failure is not None and in_tree is not None:
            delta = failure + in_tree - 1.0
        elif direct_delta is not None:
            delta = direct_delta
        else:
            missing += 1
            continue
        deltas.append(delta)
        if direct_delta is not None and abs(delta - direct_delta) > tolerance:
            direct_disagreements += 1

    failed = sum(abs(delta) > tolerance for delta in deltas)
    return {
        "tolerance": tolerance,
        "source_events": len(source_events),
        "evaluated": len(deltas),
        "missing": missing,
        "failed": failed,
        "ok": failed == 0 and missing == 0,
        "max_abs_delta": max((abs(delta) for delta in deltas), default=None),
        "mean_abs_delta": (
            sum(abs(delta) for delta in deltas) / len(deltas) if deltas else None
        ),
        "direct_delta_disagreements": direct_disagreements,
    }


def optional_weight_audit(
    source_events: Sequence[Mapping[str, Any]],
    expected_names: Sequence[str] | None = None,
    synthesized_name_events: int | float | None = None,
) -> dict[str, Any]:
    """Audit and sum optional generator variations without normalizing them.

    Optional variation denominators are not part of the current campaign
    contract.  The returned sums are therefore deliberately labelled raw and
    must not be interpreted as cross sections or event yields.
    """

    reference_names: list[str] | None = None
    generated_by_name: dict[str, WeightStats] = {}
    selected_by_name: dict[str, WeightStats] = {}
    errors: list[str] = []
    populated_events = 0
    branch_events = 0
    expected = [str(name) for name in (expected_names or [])]

    for fallback, source in enumerate(source_events):
        has_any_branch = any(
            key in source
            for key in (
                "optional_weight_names",
                "optional_weights",
                "optional_selected_event_weights",
            )
        )
        if not has_any_branch:
            continue
        branch_events += 1
        names_value = source.get("optional_weight_names", [])
        weights_value = source.get("optional_weights", [])
        selected_value = source.get("optional_selected_event_weights", [])
        names = [str(value) for value in names_value] if isinstance(names_value, list) else []
        weights = list(weights_value) if isinstance(weights_value, list) else []
        selected = list(selected_value) if isinstance(selected_value, list) else []
        index = source_index(source, fallback)

        if not names and not weights and not selected:
            continue
        populated_events += 1
        if len(set(names)) != len(names):
            errors.append(f"source {index!r}: duplicate optional-weight names")
        if len(names) != len(weights):
            errors.append(
                f"source {index!r}: {len(names)} optional names but "
                f"{len(weights)} generated weights"
            )
            continue
        if len(names) != len(selected):
            errors.append(
                f"source {index!r}: {len(names)} optional names but "
                f"{len(selected)} selected weights"
            )
            continue
        if reference_names is None:
            reference_names = names
            generated_by_name = {name: WeightStats() for name in names}
            selected_by_name = {name: WeightStats() for name in names}
        elif names != reference_names:
            errors.append(
                f"source {index!r}: optional-weight name/order differs from first event"
            )
            continue
        for name, generated_weight, selected_weight in zip(names, weights, selected):
            try:
                generated_by_name[name].add(
                    finite_float(generated_weight, f"optional generated weight {name}")
                )
                selected_by_name[name].add(
                    finite_float(selected_weight, f"optional selected weight {name}")
                )
            except (TypeError, ValueError) as error:
                errors.append(f"source {index!r}: {error}")

    anonymous_names = [
        name
        for name in (reference_names or [])
        if re.fullmatch(r"optional_[0-9]+", name)
    ]
    if expected:
        if reference_names is None:
            errors.append(
                f"manifest defines {len(expected)} optional variations but "
                "SourceEvents contains no populated variation vectors"
            )
        elif anonymous_names:
            errors.append(
                "manifest-defined optional variations use synthesized positional "
                "names in SourceEvents: " + ", ".join(anonymous_names)
            )
        elif reference_names != expected:
            errors.append(
                "SourceEvents optional-weight names/order disagree with manifest "
                "definitions"
            )
        synthesis_count = optional_float(synthesized_name_events)
        if synthesis_count is not None and synthesis_count > 0.0:
            errors.append(
                "analyzer synthesized anonymous optional-weight names for "
                f"{int(synthesis_count)} source events despite manifest definitions"
            )

    variations = [
        {
            "name": name,
            "normalization_status": "raw_unnormalized_no_variation_denominator",
            "generated_raw": generated_by_name[name].as_dict(),
            "selected_raw": selected_by_name[name].as_dict(),
        }
        for name in (reference_names or [])
    ]
    if errors:
        status = "inconsistent"
    elif variations:
        status = "ok_raw_unnormalized"
    else:
        status = "not_available"
    return {
        "status": status,
        "normalized": False,
        "normalization_note": (
            "Raw signed variation sums only; no variation denominator is available."
        ),
        "source_events": len(source_events),
        "events_with_optional_branches": branch_events,
        "events_with_variations": populated_events,
        "variation_count": len(variations),
        "variation_names": list(reference_names or []),
        "manifest_variation_names": expected,
        "anonymous_variation_names": anonymous_names,
        "synthesized_name_events": optional_float(synthesized_name_events),
        "variations": variations,
        "errors": errors,
    }


def signed_weighted_moments(
    records: Iterable[Mapping[str, Any]], value_key: str, weight_key: str = "event_weight"
) -> dict[str, Any]:
    entries = 0
    sumw = 0.0
    sumwx = 0.0
    sumwx2 = 0.0
    for record in records:
        value = optional_float(record.get(value_key))
        if value is None:
            continue
        weight_value = optional_float(record.get(weight_key))
        weight = weight_value if weight_value is not None else record_event_weight(record)
        entries += 1
        sumw += weight
        sumwx += weight * value
        sumwx2 += weight * value * value
    if sumw == 0.0:
        return {
            "entries": entries,
            "sumw": sumw,
            "mean": None,
            "sigma": None,
            "status": "undefined_zero_signed_weight",
        }
    mean = sumwx / sumw
    variance = sumwx2 / sumw - mean * mean
    tolerance = 1.0e-12 * max(1.0, abs(sumwx2 / sumw), mean * mean)
    if variance < -tolerance:
        return {
            "entries": entries,
            "sumw": sumw,
            "mean": mean,
            "sigma": None,
            "status": "undefined_negative_signed_variance",
        }
    return {
        "entries": entries,
        "sumw": sumw,
        "mean": mean,
        "sigma": math.sqrt(max(0.0, variance)),
        "status": "ok",
    }


def _final_stage(stages: Sequence[CutStage]) -> CutStage:
    for stage in reversed(stages):
        if stage.name == "selected":
            return stage
    return stages[-1]


def _normalization(
    sample: ManifestSample, generated: WeightStats, selected: WeightStats, audit: Audit
) -> dict[str, Any]:
    manifest_normalization = nested_mapping(sample.raw, "normalization")
    denominator_kind = first_present(
        manifest_normalization, "denominator_kind", "normalization_denominator_kind"
    )
    requires_predecay_denominator = sample.branching_ratio_reweighter
    if requires_predecay_denominator and sample.denominator_sumw is None:
        audit.errors.append(
            "Herwig forced-decay/BranchingRatioReweighter sample lacks "
            "normalization.denominator_sumw from pre-decay hard-event weights; "
            "using SourceEvents would cancel the decay branching weight"
        )
        denominator = None
        denominator_source = "missing_required_predecay_sumw"
    elif sample.denominator_sumw is not None:
        denominator = sample.denominator_sumw
        denominator_source = "manifest_predecay_sumw"
        if denominator_kind not in (
            None,
            "predecay_lhe_handler_weight_sum",
            "predecay_lhe_sign_sum",
        ):
            audit.warnings.append(
                f"unrecognized normalization denominator kind {denominator_kind!r}"
            )
    else:
        denominator = generated.sumw
        denominator_source = "SourceEvents.generator_weight"
    if sample.generated_sumw is not None:
        difference = generated.sumw - sample.generated_sumw
        scale = max(abs(generated.sumw), abs(sample.generated_sumw), 1.0)
        if abs(difference) > 1.0e-8 * scale:
            audit.errors.append(
                "manifest generated_sumw disagrees with SourceEvents "
                f"({sample.generated_sumw:.12g} vs {generated.sumw:.12g})"
            )
    if sample.generated_sumw2 is not None:
        difference2 = generated.sumw2 - sample.generated_sumw2
        scale2 = max(abs(generated.sumw2), abs(sample.generated_sumw2), 1.0)
        if abs(difference2) > 1.0e-8 * scale2:
            audit.errors.append(
                "manifest generated_sumw2 disagrees with SourceEvents "
                f"({sample.generated_sumw2:.12g} vs {generated.sumw2:.12g})"
            )

    if (
        sample.native_cross_section_pb is not None
        and sample.target_cross_section_pb is not None
        and sample.weight_scale is not None
    ):
        expected_target = sample.native_cross_section_pb * sample.weight_scale
        tolerance = 1.0e-8 * max(
            abs(expected_target), abs(sample.target_cross_section_pb), 1.0e-30
        )
        if abs(expected_target - sample.target_cross_section_pb) > tolerance:
            audit.errors.append(
                "target cross section is inconsistent with "
                "native_cross_section_pb * weight_scale"
            )

    if sample.native_cross_section_pb is None:
        return {
            "status": "unavailable_missing_native_cross_section",
            "denominator_sumw": denominator,
            "denominator_source": denominator_source,
            "denominator_kind": denominator_kind,
            "requires_predecay_denominator": requires_predecay_denominator,
            "scale_pb_per_weight": None,
            "selected_cross_section_pb": None,
            "selected_cross_section_sumw2_pb2": None,
        }
    if denominator is None:
        return {
            "status": "unavailable_missing_predecay_denominator",
            "denominator_sumw": None,
            "denominator_source": denominator_source,
            "denominator_kind": denominator_kind,
            "requires_predecay_denominator": requires_predecay_denominator,
            "scale_pb_per_weight": None,
            "selected_cross_section_pb": None,
            "selected_cross_section_sumw2_pb2": None,
        }
    if denominator == 0.0:
        audit.errors.append(
            "signed generated sumw is zero; normalization is undefined (no entry fallback used)"
        )
        return {
            "status": "unavailable_zero_signed_denominator",
            "denominator_sumw": denominator,
            "denominator_source": denominator_source,
            "denominator_kind": denominator_kind,
            "requires_predecay_denominator": requires_predecay_denominator,
            "scale_pb_per_weight": None,
            "selected_cross_section_pb": None,
            "selected_cross_section_sumw2_pb2": None,
        }
    scale = sample.native_cross_section_pb / denominator
    return {
        "status": "ok",
        "denominator_sumw": denominator,
        "denominator_sumw2": sample.denominator_sumw2,
        "denominator_source": denominator_source,
        "denominator_kind": denominator_kind,
        "denominator_sumabsw": optional_float(
            manifest_normalization.get("denominator_sumabsw")
        ),
        "denominator_event_count": optional_float(
            manifest_normalization.get("denominator_event_count")
        ),
        "denominator_npositive": optional_float(
            manifest_normalization.get("denominator_npositive")
        ),
        "denominator_nnegative": optional_float(
            manifest_normalization.get("denominator_nnegative")
        ),
        "denominator_nzero": optional_float(
            manifest_normalization.get("denominator_nzero")
        ),
        "requires_predecay_denominator": requires_predecay_denominator,
        "scale_pb_per_weight": scale,
        "selected_cross_section_pb": selected.sumw * scale,
        "selected_cross_section_sumw2_pb2": selected.sumw2 * scale * scale,
    }


def _yield_dict(stats: WeightStats, normalization: Mapping[str, Any], luminosity_fb: float) -> dict[str, Any]:
    scale_pb = normalization.get("scale_pb_per_weight")
    if scale_pb is None:
        return {
            "status": normalization.get("status", "unavailable"),
            "cross_section_pb": None,
            "cross_section_sumw2_pb2": None,
            "events": None,
            "events_sumw2": None,
            "stat_uncertainty_events": None,
            "raw": stats.as_dict(),
        }
    xsec = stats.sumw * scale_pb
    xsec_sumw2 = stats.sumw2 * scale_pb * scale_pb
    event_scale = luminosity_fb * 1000.0
    return {
        "status": "ok",
        "cross_section_pb": xsec,
        "cross_section_sumw2_pb2": xsec_sumw2,
        "events": xsec * event_scale,
        "events_sumw2": xsec_sumw2 * event_scale * event_scale,
        "stat_uncertainty_events": math.sqrt(max(0.0, xsec_sumw2)) * event_scale,
        "raw": stats.as_dict(),
    }


def normalized_cutflow(
    records: Sequence[Mapping[str, Any]],
    stages: Sequence[CutStage],
    normalization: Mapping[str, Any],
    luminosity_fb: float,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for raw in cutflow_stats(records, stages):
        stats = WeightStats(
            entries=int(raw["entries"]),
            positive_entries=int(raw["positive_entries"]),
            negative_entries=int(raw["negative_entries"]),
            zero_entries=int(raw["zero_entries"]),
            sumw=float(raw["sumw"]),
            sumabsw=float(raw["sumabsw"]),
            sumw2=float(raw["sumw2"]),
            sumw_positive=float(raw["sumw_positive"]),
            sumw_negative=float(raw["sumw_negative"]),
        )
        normalized = _yield_dict(stats, normalization, luminosity_fb)
        output.append(
            {
                **raw,
                "cross_section_pb": normalized["cross_section_pb"],
                "cross_section_sumw2_pb2": normalized["cross_section_sumw2_pb2"],
                "yield_events": normalized["events"],
                "yield_sumw2": normalized["events_sumw2"],
                "yield_stat_uncertainty": normalized["stat_uncertainty_events"],
                "normalization_status": normalized["status"],
            }
        )
    return output


def mass_spectrum(
    records: Sequence[Mapping[str, Any]],
    normalization: Mapping[str, Any],
    luminosity_fb: float,
    mass_range_gev: tuple[float, float] = DEFAULT_MASS_RANGE_GEV,
    bin_width_gev: float = DEFAULT_MASS_BIN_WIDTH_GEV,
) -> dict[str, Any]:
    """Return selected, signed m4l bins including underflow and overflow."""

    low, high = mass_range_gev
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError("mass spectrum range must have finite, increasing boundaries")
    if not math.isfinite(bin_width_gev) or bin_width_gev <= 0.0:
        raise ValueError("mass spectrum bin width must be positive")
    edges = [low]
    while edges[-1] < high:
        edges.append(min(high, edges[-1] + bin_width_gev))

    definitions: list[tuple[str, float | None, float | None]] = [
        ("underflow", None, low)
    ]
    definitions.extend(
        ("regular", bin_low, bin_high)
        for bin_low, bin_high in zip(edges[:-1], edges[1:])
    )
    definitions.append(("overflow", high, None))
    bins: list[dict[str, Any]] = []
    for bin_type, bin_low, bin_high in definitions:
        def in_bin(
            record: Mapping[str, Any],
            kind: str = bin_type,
            lower: float | None = bin_low,
            upper: float | None = bin_high,
        ) -> bool:
            mass = optional_float(record.get("m4l"))
            if mass is None:
                return False
            if kind == "underflow":
                return mass < float(upper)
            if kind == "overflow":
                return mass >= float(lower)
            return float(lower) <= mass < float(upper)

        stats = grouped_weight_stats(records, predicate=in_bin)
        normalized = _yield_dict(stats, normalization, luminosity_fb)
        bins.append(
            {
                "type": bin_type,
                "low_gev": bin_low,
                "high_gev": bin_high,
                "sumw": stats.sumw,
                "sumabsw": stats.sumabsw,
                "sumw2": stats.sumw2,
                "source_events": stats.entries,
                "cross_section_pb": normalized["cross_section_pb"],
                "cross_section_sumw2_pb2": normalized["cross_section_sumw2_pb2"],
                "events": normalized["events"],
                "events_sumw2": normalized["events_sumw2"],
                "stat_uncertainty_events": normalized["stat_uncertainty_events"],
                "normalization_status": normalized["status"],
            }
        )
    return {
        "mass_range_gev": [low, high],
        "bin_width_gev": bin_width_gev,
        "missing_mass_rows": sum(
            1 for record in records if optional_float(record.get("m4l")) is None
        ),
        "bins": bins,
    }


def _runtime_values_equal(field: str, expected: Any, actual: Any) -> bool:
    if field in {
        "weight_scale",
        "muon_resolution_scale",
        "electron_efficiency",
        "muon_efficiency",
        "trigger_4e_efficiency",
        "trigger_other_efficiency",
    }:
        expected_number = optional_float(expected)
        actual_number = optional_float(actual)
        return (
            expected_number is not None
            and actual_number is not None
            and math.isclose(
                expected_number,
                actual_number,
                rel_tol=1.0e-10,
                abs_tol=1.0e-15,
            )
        )
    if field == "em_pileup_noise_enabled":
        return optional_bool(expected) is not None and optional_bool(expected) == optional_bool(
            actual
        )
    if field == "requested_channel":
        return normalize_channel(expected) == normalize_channel(actual)
    if field == "response_profile":
        return str(expected).strip().lower() == str(actual).strip().lower()
    if field in {"schema_version", "seed"}:
        try:
            return int(expected) == int(actual)
        except (TypeError, ValueError):
            return False
    return str(expected) == str(actual)


def _root_analysis_metadata(records: RootRecords) -> tuple[dict[str, Any], bool]:
    nested = records.metadata.get("analysis_metadata")
    marker = records.metadata.get("analysis_metadata_present")
    if isinstance(nested, Mapping):
        root_metadata = dict(nested)
    else:
        root_metadata = {
            key: value
            for key, value in records.metadata.items()
            if key in METADATA_BRANCH_ALIASES
        }
    present = bool(root_metadata) if marker is None else bool(marker)
    return root_metadata, present


def _expected_runtime_metadata(sample: ManifestSample) -> dict[str, Any]:
    detector = sample.detector_settings
    active_efficiencies = nested_mapping(detector, "active_efficiencies")
    active_trigger = nested_mapping(active_efficiencies, "trigger")
    profile = str(sample.detector_profile or "").strip().lower()
    trigger_4e = first_present(active_trigger, "4e", "all")
    trigger_other = first_present(active_trigger, "4mu", "2e2mu", "other", "all")
    expected: dict[str, Any] = {
        "schema_version": 1,
        "response_profile": sample.detector_profile,
        "sample": sample.name,
        "category": sample.category,
        "requested_channel": sample.channel,
        "cut_mask_definition": sample.raw.get("cut_mask_definition"),
        "seed": first_present(
            nested_mapping(sample.raw, "seeds"),
            "analysis",
        )
        or sample.raw.get("seed"),
        "weight_scale": sample.weight_scale,
        "muon_resolution_scale": detector.get("muon_resolution_scale"),
        "detector_parameter_source": (
            "perfect_identity_response"
            if profile == "perfect"
            else "GEM_EM_response_and_GEM_TDR_muon_digitization"
            if profile == "ssc"
            else None
        ),
        "electron_efficiency": active_efficiencies.get("electron"),
        "muon_efficiency": active_efficiencies.get("muon"),
        "trigger_4e_efficiency": trigger_4e,
        "trigger_other_efficiency": trigger_other,
        "em_pileup_noise_enabled": detector.get("pileup_noise_enabled"),
    }
    return {key: value for key, value in expected.items() if value is not None}


def detector_provenance_audit(
    sample: ManifestSample,
    records: RootRecords,
) -> tuple[dict[str, Any], Audit]:
    """Cross-check configured and executed detector/provenance metadata."""

    local_audit = Audit()
    root_metadata, root_metadata_present = _root_analysis_metadata(records)
    summary_value = records.metadata.get("summary")
    summary = dict(summary_value) if isinstance(summary_value, Mapping) else {}
    summary_detector = nested_mapping(summary, "detector")
    summary_runtime = dict(summary)
    for target, source in (
        ("detector_parameter_source", "parameter_source"),
        ("electron_efficiency", "electron_efficiency"),
        ("muon_efficiency", "muon_efficiency"),
        ("trigger_4e_efficiency", "trigger_4e_efficiency"),
        ("trigger_other_efficiency", "trigger_other_efficiency"),
        ("em_pileup_noise_enabled", "em_pileup_noise_enabled"),
    ):
        if summary_detector.get(source) is not None:
            summary_runtime[target] = summary_detector[source]

    strict_contract = bool(sample.detector_settings) or any(
        sample.raw.get(key) is not None
        for key in ("analysis_summary", "analysis_summary_sha256")
    )
    configured_profile = first_present(
        sample.detector_settings,
        "profile",
        "response_profile",
        "detector_response",
    )
    if (
        configured_profile is not None
        and sample.detector_profile is not None
        and not _runtime_values_equal(
            "response_profile",
            sample.detector_profile,
            configured_profile,
        )
    ):
        local_audit.errors.append(
            "manifest detector.profile disagrees with the manifest sample/top-level "
            f"profile ({configured_profile!r} vs {sample.detector_profile!r})"
        )
    if strict_contract and not root_metadata_present:
        local_audit.errors.append(
            "strict detector provenance requires ROOT AnalysisMetadata"
        )
    if strict_contract and not summary:
        local_audit.errors.append(
            "strict detector provenance requires the analyzer summary sidecar"
        )

    checks: list[dict[str, Any]] = []
    expected_fields = _expected_runtime_metadata(sample)
    for field, expected in expected_fields.items():
        root_actual = root_metadata.get(field)
        summary_actual = summary_runtime.get(field)
        root_status = (
            "missing"
            if root_actual is None
            else "match"
            if _runtime_values_equal(field, expected, root_actual)
            else "mismatch"
        )
        summary_status = (
            "missing"
            if summary_actual is None
            else "match"
            if _runtime_values_equal(field, expected, summary_actual)
            else "mismatch"
        )
        checks.append(
            {
                "field": field,
                "expected": expected,
                "root_analysis_metadata": root_actual,
                "analyzer_summary": summary_actual,
                "root_status": root_status,
                "summary_status": summary_status,
            }
        )
        for source_name, status, actual in (
            ("ROOT AnalysisMetadata", root_status, root_actual),
            ("analyzer summary", summary_status, summary_actual),
        ):
            if status == "mismatch":
                local_audit.errors.append(
                    f"manifest field {field!r} disagrees with {source_name} "
                    f"({expected!r} vs {actual!r})"
                )
            elif status == "missing" and strict_contract:
                local_audit.errors.append(
                    f"{source_name} lacks required manifest field {field!r}"
                )

    events_expected = sample.raw.get("nevents_requested")
    events_summary = summary.get("events_read")
    if events_expected is not None and events_summary is not None:
        if int(events_expected) != int(events_summary):
            local_audit.errors.append(
                "manifest nevents_requested disagrees with analyzer summary "
                f"({events_expected!r} vs {events_summary!r})"
            )
    elif strict_contract and events_expected is not None:
        local_audit.errors.append("analyzer summary lacks required field 'events_read'")
    if events_summary is not None and int(events_summary) != len(records.source_events):
        local_audit.errors.append(
            "analyzer summary events_read disagrees with SourceEvents entries "
            f"({events_summary!r} vs {len(records.source_events)})"
        )

    rows_summary = summary.get("four_lepton_rows")
    if rows_summary is not None and int(rows_summary) != len(records.four_lepton):
        local_audit.errors.append(
            "analyzer summary four_lepton_rows disagrees with FourLepton entries "
            f"({rows_summary!r} vs {len(records.four_lepton)})"
        )

    summary_output = summary.get("output_root")
    if summary_output is not None and sample.output_root is not None:
        if Path(str(summary_output)).resolve() != sample.output_root.resolve():
            local_audit.errors.append(
                "analyzer summary output_root disagrees with manifest ROOT output "
                f"({summary_output!r} vs {str(sample.output_root)!r})"
            )
    elif strict_contract and sample.output_root is not None and summary:
        local_audit.errors.append("analyzer summary lacks required field 'output_root'")

    return (
        {
            "strict_contract": strict_contract,
            "manifest_detector_profile": configured_profile,
            "root_analysis_metadata_present": root_metadata_present,
            "analyzer_summary_present": bool(summary),
            "checks": checks,
            "events": {
                "manifest_requested": events_expected,
                "summary_read": events_summary,
                "source_tree_entries": len(records.source_events),
                "summary_four_lepton_rows": rows_summary,
                "four_lepton_tree_entries": len(records.four_lepton),
            },
            "summary_path": records.metadata.get("summary_path"),
            "ok": not local_audit.errors,
            "errors": list(local_audit.errors),
            "warnings": list(local_audit.warnings),
        },
        local_audit,
    )


def _manifest_optional_weight_names(sample: ManifestSample) -> list[str]:
    definitions = sample.raw.get("optional_weight_definitions")
    if not isinstance(definitions, Sequence) or isinstance(definitions, (str, bytes)):
        return []
    names: list[str] = []
    for definition in definitions:
        if isinstance(definition, Mapping):
            value = first_present(definition, "label", "name", "id")
        else:
            value = definition
        if value is not None:
            names.append(str(value))
    return names


def analyze_sample(
    sample: ManifestSample,
    records: RootRecords,
    stages: Sequence[CutStage],
    luminosity_fb: float,
    signal_region_gev: tuple[float, float],
    closure_tolerance: float,
    allow_unvalidated_bias: bool = False,
    mass_range_gev: tuple[float, float] = DEFAULT_MASS_RANGE_GEV,
    mass_bin_width_gev: float = DEFAULT_MASS_BIN_WIDTH_GEV,
) -> dict[str, Any]:
    audit = Audit()
    audit.extend(sample.audit)
    audit.extend(records.audit)
    detector_audit, detector_audit_messages = detector_provenance_audit(
        sample, records
    )
    audit.extend(detector_audit_messages)
    try:
        generated = source_weight_stats(records.source_events)
    except ValueError as error:
        audit.errors.append(str(error))
        generated = WeightStats()

    final_stage = _final_stage(stages)
    selected_rows = [record for record in records.four_lepton if passes_stage(record, final_stage)]
    selected = grouped_weight_stats(selected_rows)
    low, high = signal_region_gev
    signal_region_rows = [
        record
        for record in selected_rows
        if (optional_float(record.get("m4l")) is not None)
        and low < float(record["m4l"]) < high
    ]
    signal_region = grouped_weight_stats(signal_region_rows)
    normalization = _normalization(sample, generated, selected, audit)
    if any(optional_float(record.get("m4l")) is None for record in selected_rows):
        audit.errors.append("one or more selected FourLepton rows lack a finite m4l value")
    analyzer_summary = records.metadata.get("summary")
    analyzer_diagnostics = (
        analyzer_summary.get("diagnostics", {})
        if isinstance(analyzer_summary, Mapping)
        else {}
    )
    source_overflow_events = sum(
        1
        for source in records.source_events
        if (optional_float(source.get("accepted_lepton_count")) or 0.0) > 10.0
    )
    summary_overflow_events = optional_float(
        analyzer_diagnostics.get("accepted_lepton_overflow_events")
    )
    if summary_overflow_events is not None and summary_overflow_events < 0.0:
        audit.errors.append(
            "analyzer summary reports a negative accepted-lepton overflow count"
        )
    accepted_lepton_overflow = source_overflow_events > 0 or (
        summary_overflow_events is not None and summary_overflow_events > 0.0
    )
    if accepted_lepton_overflow:
        audit.errors.append(
            "accepted-lepton enumeration overflow occurred "
            f"(SourceEvents={source_overflow_events}, "
            f"analyzer summary={summary_overflow_events}); sample is excluded "
            "from physics totals"
        )
    if (
        summary_overflow_events is not None
        and int(summary_overflow_events) != source_overflow_events
    ):
        audit.errors.append(
            "accepted-lepton overflow count disagrees between SourceEvents and "
            "the analyzer summary "
            f"({source_overflow_events} vs {summary_overflow_events:g})"
        )

    unvalidated_bias = (
        sample.decay_bias_enabled and sample.decay_bias_closure_validated is not True
    )
    exclusion_reasons: list[str] = []
    if accepted_lepton_overflow:
        exclusion_reasons.append("accepted_lepton_overflow")
    if unvalidated_bias and not allow_unvalidated_bias:
        exclusion_reasons.append("unvalidated_decay_bias")
        audit.errors.append(
            "decay-bias closure is not validated; sample is excluded from combined "
            "physics totals (use --allow-unvalidated-bias only for an explicit study)"
        )
    elif unvalidated_bias:
        audit.warnings.append(
            "unvalidated decay-bias sample explicitly included in combined physics totals"
        )
    included_in_physics_totals = not exclusion_reasons
    exclusion_reason = exclusion_reasons[0] if exclusion_reasons else None

    channels = sorted(
        {
            normalize_channel(record.get("channel"), sample.channel)
            for record in records.four_lepton
        }
        | ({sample.channel} if sample.channel != "unknown" else set()),
        key=lambda value: ("4e", "4mu", "2e2mu", "inclusive", "unknown").index(value)
        if value in {"4e", "4mu", "2e2mu", "inclusive", "unknown"}
        else 99,
    )
    channel_results: dict[str, Any] = {}
    for channel in channels:
        channel_records = [
            record
            for record in records.four_lepton
            if normalize_channel(record.get("channel"), sample.channel) == channel
        ]
        channel_selected = [
            record for record in channel_records if passes_stage(record, final_stage)
        ]
        channel_sr = [
            record
            for record in channel_selected
            if optional_float(record.get("m4l")) is not None
            and low < float(record["m4l"]) < high
        ]
        channel_results[channel] = {
            "cutflow": normalized_cutflow(
                channel_records, stages, normalization, luminosity_fb
            ),
            "selected": _yield_dict(
                grouped_weight_stats(channel_selected), normalization, luminosity_fb
            ),
            "signal_region": _yield_dict(
                grouped_weight_stats(channel_sr), normalization, luminosity_fb
            ),
            "mass_spectrum": mass_spectrum(
                channel_selected,
                normalization,
                luminosity_fb,
                mass_range_gev,
                mass_bin_width_gev,
            ),
        }

    mass_resolution = None
    if sample.category == "signal":
        mass_resolution = signed_weighted_moments(selected_rows, "m4l")
    closure = response_closure_audit(
        records.source_events, records.four_lepton, closure_tolerance
    )
    if not closure["ok"]:
        audit.errors.append(
            "detector response closure failed or is incomplete "
            f"({closure['failed']} failed, {closure['missing']} missing)"
        )
    expected_optional_weight_names = _manifest_optional_weight_names(sample)
    optional_weights = optional_weight_audit(
        records.source_events,
        expected_names=expected_optional_weight_names,
        synthesized_name_events=analyzer_diagnostics.get(
            "optional_weight_name_syntheses"
        ),
    )
    optional_name_mismatches = optional_float(
        analyzer_diagnostics.get("optional_weight_name_mismatches")
    )
    if optional_name_mismatches is not None and optional_name_mismatches > 0.0:
        optional_weights["errors"].append(
            "analyzer reports optional-weight name/value mismatches for "
            f"{int(optional_name_mismatches)} source events"
        )
        optional_weights["status"] = "inconsistent"
    root_runtime_metadata, _ = _root_analysis_metadata(records)
    optional_name_source = str(
        root_runtime_metadata.get("optional_weight_names_source") or ""
    )
    if (
        expected_optional_weight_names
        and "synthesized" in optional_name_source.lower()
    ):
        optional_weights["errors"].append(
            "ROOT AnalysisMetadata records synthesized positional optional-weight names "
            "despite manifest definitions"
        )
        optional_weights["status"] = "inconsistent"
    for message in optional_weights["errors"]:
        audit.errors.append(f"optional weights: {message}")

    return {
        "name": sample.name,
        "label": sample.label,
        "category": sample.category,
        "channel": sample.channel,
        "detector_profile": sample.detector_profile,
        "perturbative_order": sample.perturbative_order,
        "source_kind": sample.source_kind,
        "production_group_id": sample.production_group_id,
        "comparison_identity": {
            "input_lhe_sha256": sample.raw.get("input_lhe_sha256"),
            "source_lhe_sha256": sample.raw.get("source_lhe_sha256"),
            "production_group_id": sample.production_group_id,
            "nevents_requested": sample.raw.get("nevents_requested"),
            "seeds": dict(nested_mapping(sample.raw, "seeds")),
            "source_events_read": generated.entries,
        },
        "decay_bias": {
            "enabled": sample.decay_bias_enabled,
            "scheme": sample.decay_bias_scheme,
            "closure_validated": sample.decay_bias_closure_validated,
            "bias_of": sample.decay_bias_of,
            "full_tail_support": optional_bool(
                nested_mapping(sample.raw, "decay_bias").get("full_tail_support")
            ),
            "branching_ratio_reweighter": sample.branching_ratio_reweighter,
            "prompt_decay_profile": sample.prompt_decay_profile,
        },
        "included_in_physics_totals": included_in_physics_totals,
        "exclusion_reason": exclusion_reason,
        "exclusion_reasons": exclusion_reasons,
        "manifest": str(sample.manifest_path),
        "root_file": str(sample.output_root) if sample.output_root is not None else None,
        "native_cross_section_pb": sample.native_cross_section_pb,
        "target_cross_section_pb": sample.target_cross_section_pb,
        "weight_scale": sample.weight_scale,
        "generated": generated.as_dict(),
        "selected": _yield_dict(selected, normalization, luminosity_fb),
        "signal_region": _yield_dict(signal_region, normalization, luminosity_fb),
        "mass_spectrum": mass_spectrum(
            selected_rows,
            normalization,
            luminosity_fb,
            mass_range_gev,
            mass_bin_width_gev,
        ),
        "normalization": normalization,
        "cutflow": normalized_cutflow(
            records.four_lepton, stages, normalization, luminosity_fb
        ),
        "channels": channel_results,
        "response_closure": closure,
        "optional_weights": optional_weights,
        "detector_audit": {
            **detector_audit,
            "manifest_profile": sample.detector_profile,
            "muon_pt_below_table": optional_float(
                analyzer_diagnostics.get("muon_pt_below_table")
            ),
            "muon_pt_above_table": optional_float(
                analyzer_diagnostics.get("muon_pt_above_table")
            ),
            "accepted_lepton_overflow_events": optional_float(
                analyzer_diagnostics.get("accepted_lepton_overflow_events")
            ),
            "source_tree_accepted_lepton_overflow_events": source_overflow_events,
            "accepted_lepton_overflow": accepted_lepton_overflow,
            "optional_weight_names_source": optional_name_source or None,
        },
        "mass_resolution": mass_resolution,
        "audit": audit.as_dict(),
    }


def _sum_optional(values: Iterable[float | None]) -> float | None:
    materialized = list(values)
    return None if any(value is None for value in materialized) else sum(
        float(value) for value in materialized
    )


def aggregate_mass_spectra(samples: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Combine independently normalized sample spectra by profile/channel/category."""

    groups: dict[tuple[str, str, str], dict[str, Any]] = {}
    for sample in samples:
        if not sample["included_in_physics_totals"]:
            continue
        profile = str(sample.get("detector_profile") or "missing")
        category = str(sample["category"])
        for channel, channel_data in sample["channels"].items():
            spectrum = channel_data["mass_spectrum"]
            if not any(bin_value["source_events"] for bin_value in spectrum["bins"]):
                continue
            key = (profile, channel, category)
            group = groups.get(key)
            if group is None:
                group = {
                    "detector_profile": profile,
                    "channel": channel,
                    "category": category,
                    "samples": [],
                    "mass_range_gev": list(spectrum["mass_range_gev"]),
                    "bin_width_gev": spectrum["bin_width_gev"],
                    "bins": [
                        {
                            "type": bin_value["type"],
                            "low_gev": bin_value["low_gev"],
                            "high_gev": bin_value["high_gev"],
                            "cross_section_pb": 0.0,
                            "cross_section_sumw2_pb2": 0.0,
                            "events": 0.0,
                            "events_sumw2": 0.0,
                            "normalization_available": True,
                        }
                        for bin_value in spectrum["bins"]
                    ],
                }
                groups[key] = group
            elif (
                group["mass_range_gev"] != spectrum["mass_range_gev"]
                or group["bin_width_gev"] != spectrum["bin_width_gev"]
            ):
                raise ValueError(f"incompatible mass binning in aggregate group {key}")
            group["samples"].append(sample["name"])
            for target, source in zip(group["bins"], spectrum["bins"]):
                if source["cross_section_pb"] is None:
                    target["normalization_available"] = False
                    target["cross_section_pb"] = None
                    target["cross_section_sumw2_pb2"] = None
                    target["events"] = None
                    target["events_sumw2"] = None
                    continue
                if not target["normalization_available"]:
                    continue
                target["cross_section_pb"] += source["cross_section_pb"]
                target["cross_section_sumw2_pb2"] += source[
                    "cross_section_sumw2_pb2"
                ]
                target["events"] += source["events"]
                target["events_sumw2"] += source["events_sumw2"]
    return [
        groups[key]
        for key in sorted(groups, key=lambda item: (item[0], item[1], item[2]))
    ]


def _comparison_numbers(ssc: Any, perfect: Any) -> dict[str, Any]:
    difference = (
        ssc - perfect if ssc is not None and perfect is not None else None
    )
    ratio = (
        ssc / perfect
        if ssc is not None and perfect is not None and perfect != 0.0
        else None
    )
    return {
        "perfect": perfect,
        "ssc": ssc,
        "difference_ssc_minus_perfect": difference,
        "ratio_ssc_over_perfect": ratio,
    }


def _comparison_identity(
    perfect: Mapping[str, Any],
    ssc: Mapping[str, Any],
) -> dict[str, Any]:
    perfect_identity = nested_mapping(perfect, "comparison_identity")
    ssc_identity = nested_mapping(ssc, "comparison_identity")
    checks: list[tuple[str, Any, Any]] = [
        (
            "input_lhe_sha256",
            perfect_identity.get("input_lhe_sha256"),
            ssc_identity.get("input_lhe_sha256"),
        ),
        (
            "source_lhe_sha256",
            perfect_identity.get("source_lhe_sha256"),
            ssc_identity.get("source_lhe_sha256"),
        ),
        (
            "production_group_id",
            perfect_identity.get("production_group_id"),
            ssc_identity.get("production_group_id"),
        ),
        (
            "nevents_requested",
            perfect_identity.get("nevents_requested"),
            ssc_identity.get("nevents_requested"),
        ),
        (
            "source_events_read",
            perfect_identity.get("source_events_read"),
            ssc_identity.get("source_events_read"),
        ),
    ]
    perfect_seeds = nested_mapping(perfect_identity, "seeds")
    ssc_seeds = nested_mapping(ssc_identity, "seeds")
    seed_names = ["production", "herwig", "analysis"]
    if perfect.get("category") == "signal" or ssc.get("category") == "signal":
        seed_names.insert(1, "madspin")
    for seed_name in seed_names:
        checks.append(
            (
                f"seed.{seed_name}",
                perfect_seeds.get(seed_name),
                ssc_seeds.get(seed_name),
            )
        )
    missing = [
        name for name, left, right in checks if left is None or right is None
    ]
    mismatched = [
        name
        for name, left, right in checks
        if left is not None and right is not None and left != right
    ]
    return {
        "status": (
            "mismatch"
            if mismatched
            else "unverified_missing_provenance"
            if missing
            else "verified_identical_source"
        ),
        "detector_only_interpretation_valid": not missing and not mismatched,
        "missing_fields": missing,
        "mismatched_fields": mismatched,
    }


def detector_profile_comparisons(
    samples: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    matched: dict[
        tuple[str, str, str], dict[str, list[Mapping[str, Any]]]
    ] = defaultdict(lambda: defaultdict(list))
    for sample in samples:
        profile = str(sample.get("detector_profile") or "").lower()
        if profile not in {"perfect", "ssc"}:
            continue
        key = (str(sample["name"]), str(sample["channel"]), str(sample["category"]))
        matched[key][profile].append(sample)

    comparisons: list[dict[str, Any]] = []
    for key, profiles in sorted(matched.items()):
        if "perfect" not in profiles or "ssc" not in profiles:
            continue
        if len(profiles["perfect"]) != 1 or len(profiles["ssc"]) != 1:
            continue
        perfect = profiles["perfect"][0]
        ssc = profiles["ssc"][0]
        perfect_resolution = perfect.get("mass_resolution") or {}
        ssc_resolution = ssc.get("mass_resolution") or {}
        comparisons.append(
            {
                "name": key[0],
                "channel": key[1],
                "category": key[2],
                "included_in_physics_totals": (
                    perfect["included_in_physics_totals"]
                    and ssc["included_in_physics_totals"]
                ),
                "source_identity": _comparison_identity(perfect, ssc),
                "selected_events": _comparison_numbers(
                    ssc["selected"]["events"], perfect["selected"]["events"]
                ),
                "signal_region_events": _comparison_numbers(
                    ssc["signal_region"]["events"],
                    perfect["signal_region"]["events"],
                ),
                "mass_mean_gev": _comparison_numbers(
                    ssc_resolution.get("mean"), perfect_resolution.get("mean")
                ),
                "mass_sigma_gev": _comparison_numbers(
                    ssc_resolution.get("sigma"), perfect_resolution.get("sigma")
                ),
            }
        )
    return comparisons


def _heavy_flavour_closure_metric(
    biased: Mapping[str, Any],
    control: Mapping[str, Any],
) -> dict[str, Any]:
    biased_value = optional_float(biased.get("cross_section_pb"))
    control_value = optional_float(control.get("cross_section_pb"))
    biased_variance = optional_float(biased.get("cross_section_sumw2_pb2"))
    control_variance = optional_float(control.get("cross_section_sumw2_pb2"))
    if biased_value is None or control_value is None:
        return {
            "status": "unavailable_normalization",
            "biased_cross_section_pb": biased_value,
            "control_cross_section_pb": control_value,
            "difference_pb": None,
            "ratio_biased_over_control": None,
            "combined_stat_uncertainty_pb": None,
            "difference_over_combined_stat_uncertainty": None,
            "uncertainty_assumption": "independent_samples_zero_covariance",
        }
    difference = biased_value - control_value
    combined_variance = (
        biased_variance + control_variance
        if biased_variance is not None and control_variance is not None
        else None
    )
    combined_uncertainty = (
        math.sqrt(max(0.0, combined_variance))
        if combined_variance is not None
        else None
    )
    return {
        "status": "diagnostic_only",
        "biased_cross_section_pb": biased_value,
        "control_cross_section_pb": control_value,
        "difference_pb": difference,
        "ratio_biased_over_control": (
            biased_value / control_value if control_value != 0.0 else None
        ),
        "combined_stat_uncertainty_pb": combined_uncertainty,
        "difference_over_combined_stat_uncertainty": (
            difference / combined_uncertainty
            if combined_uncertainty not in (None, 0.0)
            else None
        ),
        "uncertainty_assumption": "independent_samples_zero_covariance",
    }


def heavy_flavour_closure_diagnostics(
    samples: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Compare biased strata with controls without promoting them to validated."""

    by_name_profile: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for sample in samples:
        by_name_profile[
            (
                str(sample["name"]),
                str(sample.get("detector_profile") or "missing").lower(),
            )
        ].append(sample)

    diagnostics: list[dict[str, Any]] = []
    for biased in samples:
        bias = nested_mapping(biased, "decay_bias")
        if optional_bool(bias.get("enabled")) is not True:
            continue
        control_name_value = bias.get("bias_of")
        control_name = str(control_name_value) if control_name_value is not None else None
        profile = str(biased.get("detector_profile") or "missing").lower()
        controls = (
            by_name_profile.get((control_name, profile), [])
            if control_name is not None
            else []
        )
        base: dict[str, Any] = {
            "biased_sample": biased["name"],
            "control_sample": control_name,
            "detector_profile": biased.get("detector_profile"),
            "scheme": bias.get("scheme"),
            "full_tail_support": bias.get("full_tail_support"),
            "manifest_closure_validated": bias.get("closure_validated"),
            "closure_validation_claimed_by_report": False,
            "interpretation": (
                "Diagnostic weighted agreement only. This report never changes "
                "closure_validated or admits the biased sample to nominal totals."
            ),
        }
        if control_name is None:
            diagnostics.append({**base, "status": "missing_control_mapping"})
            continue
        if len(controls) != 1:
            diagnostics.append(
                {
                    **base,
                    "status": (
                        "missing_control_sample"
                        if not controls
                        else "ambiguous_control_sample"
                    ),
                    "control_matches": len(controls),
                }
            )
            continue
        control = controls[0]
        cutflow_control = {
            str(row["name"]): row for row in control.get("cutflow", [])
        }
        cutflow_diagnostics = []
        for biased_row in biased.get("cutflow", []):
            control_row = cutflow_control.get(str(biased_row["name"]))
            if control_row is None:
                continue
            cutflow_diagnostics.append(
                {
                    "name": biased_row["name"],
                    "label": biased_row.get("label"),
                    **_heavy_flavour_closure_metric(
                        biased_row,
                        control_row,
                    ),
                }
            )
        selected_metric = _heavy_flavour_closure_metric(
            nested_mapping(biased, "selected"),
            nested_mapping(control, "selected"),
        )
        signal_region_metric = _heavy_flavour_closure_metric(
            nested_mapping(biased, "signal_region"),
            nested_mapping(control, "signal_region"),
        )
        available = (
            selected_metric["status"] == "diagnostic_only"
            and signal_region_metric["status"] == "diagnostic_only"
        )
        diagnostics.append(
            {
                **base,
                "status": (
                    "diagnostic_available_not_validation"
                    if available
                    else "diagnostic_incomplete_normalization"
                ),
                "biased_source_events": nested_mapping(biased, "generated").get(
                    "entries"
                ),
                "control_source_events": nested_mapping(control, "generated").get(
                    "entries"
                ),
                "selected": selected_metric,
                "signal_region": signal_region_metric,
                "cutflow": cutflow_diagnostics,
            }
        )
    return diagnostics


def _physics_summary(samples: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    included = [sample for sample in samples if sample["included_in_physics_totals"]]
    signal = [sample for sample in included if sample["category"] == "signal"]
    backgrounds = [sample for sample in included if sample["category"] != "signal"]
    signal_yield = _sum_optional(sample["signal_region"]["events"] for sample in signal)
    background_yield = _sum_optional(
        sample["signal_region"]["events"] for sample in backgrounds
    )
    signal_variance = _sum_optional(
        sample["signal_region"]["events_sumw2"] for sample in signal
    )
    background_variance = _sum_optional(
        sample["signal_region"]["events_sumw2"] for sample in backgrounds
    )

    composition: list[dict[str, Any]] = []
    for sample in backgrounds:
        events = sample["signal_region"]["events"]
        fraction = (
            events / background_yield
            if events is not None and background_yield not in (None, 0.0)
            else None
        )
        composition.append(
            {
                "name": sample["name"],
                "label": sample["label"],
                "category": sample["category"],
                "detector_profile": sample.get("detector_profile"),
                "events": events,
                "stat_uncertainty_events": sample["signal_region"][
                    "stat_uncertainty_events"
                ],
                "fraction_of_total_background": fraction,
            }
        )

    return {
        "signal_region": {
            "signal_events": signal_yield,
            "signal_stat_uncertainty_events": (
                math.sqrt(max(0.0, signal_variance))
                if signal_variance is not None
                else None
            ),
            "background_events": background_yield,
            "background_stat_uncertainty_events": (
                math.sqrt(max(0.0, background_variance))
                if background_variance is not None
                else None
            ),
            "signal_to_background": (
                signal_yield / background_yield
                if signal_yield is not None and background_yield not in (None, 0.0)
                else None
            ),
        },
        "background_composition": composition,
    }


def summarize_samples(samples: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    profiles: dict[str, dict[str, Any]] = {}
    for sample in samples:
        profile = sample.get("detector_profile") or "missing"
        profile_info = profiles.setdefault(
            profile,
            {
                "samples": [],
                "response_sources_evaluated": 0,
                "response_sources_failed": 0,
                "response_sources_missing": 0,
                "max_abs_closure_delta": None,
            },
        )
        profile_info["samples"].append(sample["name"])
        closure = sample["response_closure"]
        profile_info["response_sources_evaluated"] += closure["evaluated"]
        profile_info["response_sources_failed"] += closure["failed"]
        profile_info["response_sources_missing"] += closure["missing"]
        delta = closure["max_abs_delta"]
        current = profile_info["max_abs_closure_delta"]
        if delta is not None and (current is None or delta > current):
            profile_info["max_abs_closure_delta"] = delta

    physics_by_profile = {
        profile: _physics_summary(
            [
                sample
                for sample in samples
                if (sample.get("detector_profile") or "missing") == profile
            ]
        )
        for profile in sorted(profiles)
    }
    if len(physics_by_profile) == 1:
        sole_physics = next(iter(physics_by_profile.values()))
        combined_signal_region = sole_physics["signal_region"]
    else:
        # Perfect and SSC are alternative detector hypotheses, not disjoint
        # event samples.  A cross-profile sum would double-count the physics.
        combined_signal_region = {
            "status": "multiple_detector_profiles_use_physics_by_profile",
            "signal_events": None,
            "signal_stat_uncertainty_events": None,
            "background_events": None,
            "background_stat_uncertainty_events": None,
            "signal_to_background": None,
        }
    composition_by_profile = {
        profile: values["background_composition"]
        for profile, values in physics_by_profile.items()
    }

    return {
        "signal_region": combined_signal_region,
        "physics_by_profile": physics_by_profile,
        "background_composition_by_profile": composition_by_profile,
        "background_composition": [
            item
            for profile in sorted(composition_by_profile)
            for item in composition_by_profile[profile]
        ],
        "detector_profiles": profiles,
        "mass_resolution": [
            {
                "name": sample["name"],
                "label": sample["label"],
                "channel": sample["channel"],
                "detector_profile": sample["detector_profile"],
                **sample["mass_resolution"],
            }
            for sample in samples
            if sample["included_in_physics_totals"]
            and sample["category"] == "signal"
            if sample.get("mass_resolution") is not None
        ],
        "excluded_samples": [
            {
                "name": sample["name"],
                "label": sample["label"],
                "reason": sample["exclusion_reason"],
            }
            for sample in samples
            if not sample["included_in_physics_totals"]
        ],
        "mass_spectra": aggregate_mass_spectra(samples),
        "perfect_ssc_comparisons": detector_profile_comparisons(samples),
        "heavy_flavour_closure_diagnostics": heavy_flavour_closure_diagnostics(
            samples
        ),
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _format_number(value: Any, digits: int = 5) -> str:
    if value is None:
        return "&mdash;"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return html.escape(str(value))
    if number == 0.0:
        return "0"
    if abs(number) >= 1.0e4 or abs(number) < 1.0e-3:
        return f"{number:.{digits - 1}e}"
    return f"{number:.{digits}g}"


def _audit_html(audit: Mapping[str, Any]) -> str:
    errors = audit.get("errors", [])
    warnings = audit.get("warnings", [])
    if not errors and not warnings:
        return '<p class="ok">All report provenance checks passed.</p>'
    parts: list[str] = []
    if errors:
        parts.append("<h3>Errors</h3><ul class=\"errors\">")
        parts.extend(f"<li>{html.escape(str(message))}</li>" for message in errors)
        parts.append("</ul>")
    if warnings:
        parts.append("<h3>Warnings</h3><ul class=\"warnings\">")
        parts.extend(f"<li>{html.escape(str(message))}</li>" for message in warnings)
        parts.append("</ul>")
    return "".join(parts)


def _spectrum_svg(spectrum: Mapping[str, Any]) -> str:
    regular = [bin_value for bin_value in spectrum["bins"] if bin_value["type"] == "regular"]
    values = [bin_value["events"] for bin_value in regular]
    if not regular or any(value is None for value in values):
        return '<p class="note">Normalized spectrum unavailable.</p>'
    numeric = [float(value) for value in values]
    width, height = 820.0, 260.0
    left, right, top, bottom = 62.0, 18.0, 18.0, 42.0
    plot_width = width - left - right
    plot_height = height - top - bottom
    minimum = min(0.0, min(numeric, default=0.0))
    maximum = max(0.0, max(numeric, default=0.0))
    if maximum == minimum:
        maximum = minimum + 1.0

    def y_position(value: float) -> float:
        return top + (maximum - value) / (maximum - minimum) * plot_height

    baseline = y_position(0.0)
    bin_width = plot_width / len(regular)
    bars: list[str] = []
    for index, value in enumerate(numeric):
        y_value = y_position(value)
        y = min(y_value, baseline)
        bar_height = max(0.8, abs(y_value - baseline))
        color = "#2474b5" if value >= 0.0 else "#b33b33"
        bars.append(
            f'<rect x="{left + index * bin_width + 0.7:.2f}" y="{y:.2f}" '
            f'width="{max(0.5, bin_width - 1.4):.2f}" height="{bar_height:.2f}" '
            f'fill="{color}"><title>{value:.8g} events</title></rect>'
        )
    mass_low, mass_high = spectrum["mass_range_gev"]
    return (
        f'<svg class="spectrum" viewBox="0 0 {width:g} {height:g}" '
        'role="img" aria-label="Signed normalized four-lepton mass spectrum">'
        f'<line x1="{left:g}" y1="{baseline:.2f}" x2="{width-right:g}" '
        f'y2="{baseline:.2f}" stroke="#263746" stroke-width="1"/>'
        f'<line x1="{left:g}" y1="{top:g}" x2="{left:g}" '
        f'y2="{height-bottom:g}" stroke="#263746" stroke-width="1"/>'
        + "".join(bars)
        + f'<text x="{left:g}" y="{height-13:g}" font-size="12">{mass_low:g}</text>'
        f'<text x="{width-right:g}" y="{height-13:g}" text-anchor="end" '
        f'font-size="12">{mass_high:g} GeV</text>'
        f'<text x="14" y="{top+5:g}" font-size="11">{maximum:.4g}</text>'
        f'<text x="14" y="{height-bottom:g}" font-size="11">{minimum:.4g}</text>'
        '<text x="18" y="135" transform="rotate(-90 18 135)" '
        'font-size="12">events / bin</text></svg>'
    )


def render_html(report: Mapping[str, Any]) -> str:
    metadata = report["metadata"]
    summary = report["summary"]
    low, high = metadata["signal_region_gev"]
    sample_rows: list[str] = []
    for sample in report["samples"]:
        selected = sample["selected"]
        signal_region = sample["signal_region"]
        sample_rows.append(
            "<tr>"
            f"<td>{html.escape(sample['label'])}<br><code>{html.escape(sample['name'])}</code></td>"
            f"<td>{html.escape(sample['category'])}</td>"
            f"<td>{html.escape(sample['channel'])}</td>"
            f"<td>{html.escape(str(sample['detector_profile'] or 'missing'))}</td>"
            f"<td>{_format_number(sample['generated']['sumw'])}</td>"
            f"<td>{_format_number(selected['cross_section_pb'])}</td>"
            f"<td>{_format_number(signal_region['events'])} &plusmn; "
            f"{_format_number(signal_region['stat_uncertainty_events'])}</td>"
            f"<td>{'included' if sample['included_in_physics_totals'] else 'EXCLUDED'}"
            f"<br>{'OK' if sample['audit']['ok'] else 'CHECK'}</td>"
            "</tr>"
        )

    cutflow_sections: list[str] = []
    for sample in report["samples"]:
        def cutflow_table(stages: Sequence[Mapping[str, Any]]) -> str:
            rows = "".join(
                "<tr>"
                f"<td>{html.escape(str(stage['label']))}</td>"
                f"<td>{stage['entries']}</td>"
                f"<td>{_format_number(stage['sumw'])}</td>"
                f"<td>{_format_number(stage['sumabsw'])}</td>"
                f"<td>{_format_number(stage['sumw2'])}</td>"
                f"<td>{_format_number(stage['cross_section_pb'])}</td>"
                f"<td>{_format_number(stage['yield_events'])} &plusmn; "
                f"{_format_number(stage['yield_stat_uncertainty'])}</td>"
                "</tr>"
                for stage in stages
            )
            return (
                '<table><thead><tr><th>Stage</th><th>Source events</th>'
                "<th>Signed sumw</th><th>sum|w|</th><th>sumw2</th>"
                "<th>Cross section [pb]</th><th>Yield</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>"
            )

        channel_tables = "".join(
            f"<h3>Channel {html.escape(channel)}</h3>{cutflow_table(values['cutflow'])}"
            for channel, values in sample["channels"].items()
        )
        cutflow_sections.append(
            f"<details><summary>{html.escape(sample['label'])}</summary>"
            f"<h3>All channels</h3>{cutflow_table(sample['cutflow'])}"
            f"{channel_tables}</details>"
        )

    composition_rows = "".join(
        "<tr>"
        f"<td>{html.escape(str(item['detector_profile'] or 'missing'))}</td>"
        f"<td>{html.escape(item['label'])}</td>"
        f"<td>{html.escape(item['category'])}</td>"
        f"<td>{_format_number(item['events'])} &plusmn; "
        f"{_format_number(item['stat_uncertainty_events'])}</td>"
        f"<td>{_format_number(item['fraction_of_total_background'])}</td>"
        "</tr>"
        for item in summary["background_composition"]
    )
    profile_yield_rows = "".join(
        "<tr>"
        f"<td>{html.escape(profile)}</td>"
        f"<td>{_format_number(values['signal_region']['signal_events'])} &plusmn; "
        f"{_format_number(values['signal_region']['signal_stat_uncertainty_events'])}</td>"
        f"<td>{_format_number(values['signal_region']['background_events'])} &plusmn; "
        f"{_format_number(values['signal_region']['background_stat_uncertainty_events'])}</td>"
        f"<td>{_format_number(values['signal_region']['signal_to_background'])}</td>"
        "</tr>"
        for profile, values in sorted(summary["physics_by_profile"].items())
    )
    resolution_rows = "".join(
        "<tr>"
        f"<td>{html.escape(item['label'])}</td>"
        f"<td>{html.escape(str(item['detector_profile']))}</td>"
        f"<td>{html.escape(item['channel'])}</td>"
        f"<td>{_format_number(item['mean'])}</td>"
        f"<td>{_format_number(item['sigma'])}</td>"
        f"<td>{html.escape(item['status'])}</td>"
        "</tr>"
        for item in summary["mass_resolution"]
    )
    detector_rows = "".join(
        "<tr>"
        f"<td>{html.escape(profile)}</td>"
        f"<td>{html.escape(', '.join(values['samples']))}</td>"
        f"<td>{values['response_sources_evaluated']}</td>"
        f"<td>{values['response_sources_failed']}</td>"
        f"<td>{values['response_sources_missing']}</td>"
        f"<td>{_format_number(values['max_abs_closure_delta'])}</td>"
        "</tr>"
        for profile, values in sorted(summary["detector_profiles"].items())
    )
    optional_weight_rows = "".join(
        "<tr>"
        f"<td>{html.escape(sample['label'])}</td>"
        f"<td>{html.escape(sample['optional_weights']['status'])}</td>"
        f"<td>{sample['optional_weights']['events_with_variations']}</td>"
        f"<td>{sample['optional_weights']['variation_count']}</td>"
        f"<td>{html.escape(', '.join(sample['optional_weights']['variation_names']) or 'none')}</td>"
        "<td>raw, unnormalized</td>"
        "</tr>"
        for sample in report["samples"]
    )
    spectrum_sections: list[str] = []
    for spectrum in summary["mass_spectra"]:
        underflow = spectrum["bins"][0]
        overflow = spectrum["bins"][-1]
        heading = (
            f"{spectrum['detector_profile']} / {spectrum['channel']} / "
            f"{spectrum['category']}"
        )
        spectrum_sections.append(
            f"<details><summary>{html.escape(heading)}</summary>"
            f"{_spectrum_svg(spectrum)}"
            f"<p class=\"note\">Samples: {html.escape(', '.join(spectrum['samples']))}. "
            f"Underflow: {_format_number(underflow['events'])} events; "
            f"overflow: {_format_number(overflow['events'])} events. "
            f"Each bin is {spectrum['bin_width_gev']:g} GeV and retains signed "
            "sumw and sumw2 in report.json.</p></details>"
        )
    comparison_rows = "".join(
        "<tr>"
        f"<td>{html.escape(item['name'])}</td>"
        f"<td>{html.escape(item['channel'])}</td>"
        f"<td>{_format_number(item['selected_events']['perfect'])}</td>"
        f"<td>{_format_number(item['selected_events']['ssc'])}</td>"
        f"<td>{_format_number(item['selected_events']['ratio_ssc_over_perfect'])}</td>"
        f"<td>{_format_number(item['mass_sigma_gev']['perfect'])}</td>"
        f"<td>{_format_number(item['mass_sigma_gev']['ssc'])}</td>"
        f"<td>{_format_number(item['mass_sigma_gev']['difference_ssc_minus_perfect'])}</td>"
        f"<td>{html.escape(item['source_identity']['status'])}</td>"
        "</tr>"
        for item in summary["perfect_ssc_comparisons"]
    )
    heavy_flavour_closure_rows = "".join(
        "<tr>"
        f"<td>{html.escape(str(item['biased_sample']))}</td>"
        f"<td>{html.escape(str(item.get('control_sample') or 'missing'))}</td>"
        f"<td>{html.escape(str(item.get('detector_profile') or 'missing'))}</td>"
        f"<td>{html.escape(str(item['status']))}</td>"
        f"<td>{_format_number(nested_mapping(item, 'selected').get('ratio_biased_over_control'))}</td>"
        f"<td>{_format_number(nested_mapping(item, 'selected').get('difference_over_combined_stat_uncertainty'))}</td>"
        f"<td>{_format_number(nested_mapping(item, 'signal_region').get('ratio_biased_over_control'))}</td>"
        f"<td>{_format_number(nested_mapping(item, 'signal_region').get('difference_over_combined_stat_uncertainty'))}</td>"
        "</tr>"
        for item in summary["heavy_flavour_closure_diagnostics"]
    )
    global_audit = report["audit"]
    sr = summary["signal_region"]
    title = html.escape(metadata["title"])
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
:root {{ color-scheme: light; --ink:#18212b; --muted:#5e6b78; --line:#d8e0e8;
  --blue:#1557a0; --pale:#f4f7fb; --red:#a12622; --amber:#815500; --green:#18733f; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; color:var(--ink); background:#fff; font:15px/1.48 system-ui,sans-serif; }}
main {{ max-width:1180px; margin:0 auto; padding:32px 24px 56px; }}
h1 {{ margin:0 0 6px; font-size:2rem; }}
h2 {{ margin-top:34px; border-bottom:2px solid var(--line); padding-bottom:7px; }}
h3 {{ margin-bottom:5px; }}
.subtitle {{ color:var(--muted); margin:0 0 24px; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:12px; }}
.card {{ background:var(--pale); border:1px solid var(--line); border-radius:8px; padding:14px; }}
.card b {{ display:block; color:var(--muted); font-size:.82rem; text-transform:uppercase; }}
.card span {{ display:block; margin-top:5px; font-size:1.35rem; font-variant-numeric:tabular-nums; }}
table {{ width:100%; border-collapse:collapse; margin:12px 0 20px; }}
th,td {{ border-bottom:1px solid var(--line); padding:8px 9px; text-align:right; vertical-align:top; }}
th:first-child,td:first-child,td:nth-child(2),td:nth-child(3),td:nth-child(4) {{ text-align:left; }}
th {{ background:var(--pale); color:#344554; }}
code {{ color:var(--muted); font-size:.82em; }}
details {{ border:1px solid var(--line); border-radius:6px; margin:8px 0; padding:8px 12px; }}
summary {{ cursor:pointer; font-weight:650; }}
.spectrum {{ display:block; width:100%; max-width:820px; height:auto; margin:10px auto; background:#fff; }}
.ok {{ color:var(--green); font-weight:650; }}
.errors {{ color:var(--red); }}
.warnings {{ color:var(--amber); }}
.note {{ color:var(--muted); font-size:.9rem; }}
@media(max-width:760px) {{ main {{ padding:20px 12px; }} table {{ font-size:.82rem; }} th,td {{ padding:6px 4px; }} }}
</style>
</head>
<body><main>
<h1>{title}</h1>
<p class="subtitle">LO simulation at {_format_number(metadata.get('sqrt_s_tev'))} TeV;
signed event weights; luminosity {_format_number(metadata['luminosity_fb'])} fb<sup>&minus;1</sup>.
Signal region: {low:g} &lt; m<sub>4l</sub> &lt; {high:g} GeV.</p>

<section class="cards">
<div class="card"><b>Single-profile signal yield</b><span>{_format_number(sr['signal_events'])} &plusmn; {_format_number(sr['signal_stat_uncertainty_events'])}</span></div>
<div class="card"><b>Single-profile background yield</b><span>{_format_number(sr['background_events'])} &plusmn; {_format_number(sr['background_stat_uncertainty_events'])}</span></div>
<div class="card"><b>Single-profile signal / background</b><span>{_format_number(sr['signal_to_background'])}</span></div>
<div class="card"><b>Detector profiles</b><span>{html.escape(', '.join(summary['detector_profiles']) or 'none')}</span></div>
</section>

<h2>Signal-region yields by detector profile</h2>
<table><thead><tr><th>Response</th><th>Signal events</th><th>Background events</th>
<th>Signal / background</th></tr></thead><tbody>{profile_yield_rows}</tbody></table>
<p class="note">Perfect and SSC are alternative response hypotheses. Their yields
are never added together.</p>

<h2>Samples and normalized yields</h2>
<table><thead><tr><th>Sample</th><th>Category</th><th>Channel</th><th>Response</th>
<th>Generated sumw</th><th>Selected cross section [pb]</th><th>Signal-region events</th><th>Audit</th>
</tr></thead><tbody>{''.join(sample_rows)}</tbody></table>
<p class="note">Uncertainties are propagated from source-grouped sumw2. Negative and exactly
zero signed totals are retained; an em dash denotes unavailable normalization.
Samples with unvalidated decay-bias closure are visibly excluded from combined totals.</p>

<h2>Background composition in the signal region</h2>
<table><thead><tr><th>Response</th><th>Sample</th><th>Category</th><th>Events</th><th>Fraction</th></tr></thead>
<tbody>{composition_rows or '<tr><td colspan="5">No background samples</td></tr>'}</tbody></table>

<h2>Mass resolution</h2>
<table><thead><tr><th>Signal sample</th><th>Response</th><th>Channel</th>
<th>Signed mean [GeV]</th><th>Signed RMS [GeV]</th><th>Status</th></tr></thead>
<tbody>{resolution_rows or '<tr><td colspan="6">No signal mass entries</td></tr>'}</tbody></table>

<h2>Full selected m<sub>4l</sub> spectra</h2>
{''.join(spectrum_sections) or '<p>No normalized spectra are available.</p>'}

<h2>Perfect versus SSC detector response</h2>
<table><thead><tr><th>Matched sample</th><th>Channel</th>
<th>Perfect selected yield</th><th>SSC selected yield</th><th>SSC / perfect</th>
<th>Perfect mass RMS [GeV]</th><th>SSC mass RMS [GeV]</th><th>RMS difference [GeV]</th>
<th>Source identity</th></tr></thead><tbody>{comparison_rows or '<tr><td colspan="9">No matched perfect/SSC sample pairs</td></tr>'}</tbody></table>

<h2>Detector-response closure</h2>
<table><thead><tr><th>Profile</th><th>Samples</th><th>Evaluated</th><th>Failed</th>
<th>Missing</th><th>Maximum |closure &minus; 1|</th></tr></thead><tbody>{detector_rows}</tbody></table>

<h2>Optional-weight audit</h2>
<table><thead><tr><th>Sample</th><th>Status</th><th>Events with variations</th>
<th>Variation count</th><th>Names</th><th>Interpretation</th></tr></thead>
<tbody>{optional_weight_rows}</tbody></table>
<p class="note">Variation sums in report.json retain signed generated and selected
weights, including sumw2. They are intentionally not normalized because the campaign
does not yet provide a denominator for each variation.</p>

<h2>Heavy-flavour bias/control diagnostics</h2>
<table><thead><tr><th>Biased stratum</th><th>Unbiased control</th><th>Response</th>
<th>Status</th><th>Selected ratio</th><th>Selected difference / stat. uncertainty</th>
<th>Signal-region ratio</th><th>Signal-region difference / stat. uncertainty</th>
</tr></thead>
<tbody>{heavy_flavour_closure_rows or '<tr><td colspan="8">No heavy-flavour biased strata</td></tr>'}</tbody></table>
<p class="note">These are transparent weighted comparisons only. They do not establish
full proposal support, do not set closure_validated, and never override nominal
sample exclusion. The displayed difference/statistical-uncertainty values assume
zero covariance between the biased and control samples.</p>

<h2>Cutflows</h2>
{''.join(cutflow_sections)}

<h2>Provenance audit</h2>
{_audit_html(global_audit)}
<p class="note">Machine-readable values and the complete per-channel cutflows are in
<a href="report.json">report.json</a>. Generated {_datetime.datetime.now(_datetime.timezone.utc).strftime('%Y-%m-%d')}.</p>
</main></body></html>
"""


def expand_manifest_paths(paths: Sequence[Path]) -> list[Path]:
    expanded: list[Path] = []
    for path_value in paths:
        path = path_value.resolve()
        if path.is_dir():
            candidates = sorted(path.glob("**/manifest.json"))
            candidates.extend(sorted(path.glob("**/*manifest*.json")))
            seen: set[Path] = set()
            for candidate in candidates:
                resolved = candidate.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    expanded.append(resolved)
        else:
            expanded.append(path)
    return expanded


def generate_report(
    manifest_paths: Sequence[Path],
    output_dir: Path,
    *,
    luminosity_fb: float | None = None,
    signal_region_gev: tuple[float, float] = DEFAULT_SIGNAL_REGION_GEV,
    mass_range_gev: tuple[float, float] = DEFAULT_MASS_RANGE_GEV,
    mass_bin_width_gev: float = DEFAULT_MASS_BIN_WIDTH_GEV,
    title: str = "LO H -> ZZ* -> 4l analysis",
    closure_tolerance: float = 1.0e-9,
    record_loader: Callable[[Path], RootRecords] = load_root_records,
    allow_unvalidated_bias: bool = False,
) -> dict[str, Any]:
    manifests = [load_manifest(path) for path in expand_manifest_paths(manifest_paths)]
    if not manifests:
        raise ValueError("no campaign manifests found")
    if luminosity_fb is None:
        declared = [bundle.luminosity_fb for bundle in manifests if bundle.luminosity_fb is not None]
        luminosity_fb = declared[0] if declared else DEFAULT_LUMINOSITY_FB
    if luminosity_fb <= 0.0:
        raise ValueError("luminosity must be positive")
    low, high = signal_region_gev
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError("signal region must have finite, increasing boundaries")
    mass_low, mass_high = mass_range_gev
    if not (math.isfinite(mass_low) and math.isfinite(mass_high) and mass_low < mass_high):
        raise ValueError("mass spectrum range must have finite, increasing boundaries")
    if not math.isfinite(mass_bin_width_gev) or mass_bin_width_gev <= 0.0:
        raise ValueError("mass spectrum bin width must be positive")

    global_audit = Audit()
    sample_reports: list[dict[str, Any]] = []
    sqrt_values = {bundle.sqrt_s_tev for bundle in manifests if bundle.sqrt_s_tev is not None}
    if len(sqrt_values) > 1:
        global_audit.errors.append("manifests mix different collider energies")
    stages = manifests[0].cut_stages
    for bundle in manifests:
        global_audit.extend(bundle.audit, prefix=f"{bundle.path.name}: ")
        if bundle.cut_stages != stages:
            global_audit.errors.append(f"{bundle.path}: cut-stage definition differs")
        for sample in bundle.samples:
            if sample.output_root is None:
                records = RootRecords([], [], Audit(errors=["ROOT path unavailable"]))
            else:
                records = record_loader(sample.output_root)
            result = analyze_sample(
                sample,
                records,
                stages,
                luminosity_fb,
                signal_region_gev,
                closure_tolerance,
                allow_unvalidated_bias,
                mass_range_gev,
                mass_bin_width_gev,
            )
            sample_reports.append(result)
            for message in result["audit"]["errors"]:
                global_audit.errors.append(f"{sample.name}: {message}")
            for message in result["audit"]["warnings"]:
                global_audit.warnings.append(f"{sample.name}: {message}")

    profiles = {sample["detector_profile"] for sample in sample_reports}
    if None in profiles:
        global_audit.errors.append("one or more samples lack a detector response profile")
    if len(profiles) > 1:
        global_audit.warnings.append(
            "report compares multiple detector profiles: "
            + ", ".join(sorted(str(profile) for profile in profiles))
        )
    sample_profile_counts: dict[tuple[str, str, str, str], int] = defaultdict(int)
    for sample in sample_reports:
        identity_key = (
            str(sample["name"]),
            str(sample["channel"]),
            str(sample["category"]),
            str(sample.get("detector_profile") or "missing").lower(),
        )
        sample_profile_counts[identity_key] += 1
    for identity_key, count in sorted(sample_profile_counts.items()):
        if count > 1:
            global_audit.errors.append(
                "duplicate sample/profile identity "
                f"name={identity_key[0]}, channel={identity_key[1]}, "
                f"category={identity_key[2]}, profile={identity_key[3]} "
                f"({count} entries); profile comparison is ambiguous"
            )
    summary = summarize_samples(sample_reports)
    for comparison in summary["perfect_ssc_comparisons"]:
        identity = comparison["source_identity"]
        if identity["status"] == "mismatch":
            global_audit.errors.append(
                f"{comparison['name']}: perfect/SSC source identity differs in "
                + ", ".join(identity["mismatched_fields"])
            )
        elif identity["status"] == "unverified_missing_provenance":
            global_audit.warnings.append(
                f"{comparison['name']}: perfect/SSC source identity could not "
                "be fully verified; missing "
                + ", ".join(identity["missing_fields"])
            )

    report = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "metadata": {
            "title": title,
            "generated_at_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
            "manifest_paths": [str(bundle.path) for bundle in manifests],
            "run_tags": sorted(
                {bundle.run_tag for bundle in manifests if bundle.run_tag is not None}
            ),
            "sqrt_s_tev": next(iter(sqrt_values)) if len(sqrt_values) == 1 else None,
            "luminosity_fb": luminosity_fb,
            "signal_region_gev": [low, high],
            "mass_spectrum_range_gev": list(mass_range_gev),
            "mass_spectrum_bin_width_gev": mass_bin_width_gev,
            "normalization_convention": (
                "native_cross_section_pb * sum(event_weight) / "
                "normalization.denominator_sumw; unforced-sample fallback: "
                "sum(SourceEvents.generator_weight)"
            ),
            "allow_unvalidated_bias": allow_unvalidated_bias,
            "cut_stages": [
                {
                    "name": stage.name,
                    "label": stage.label,
                    "bit": stage.bit,
                    "required_mask": stage.required_mask,
                }
                for stage in stages
            ],
        },
        "samples": sample_reports,
        "summary": summary,
        "audit": global_audit.as_dict(),
    }
    report = _json_safe(report)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    (output_dir / "index.html").write_text(render_html(report))
    return report


def positive_float(value: str) -> float:
    parsed = finite_float(value, "value")
    if parsed <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "manifests",
        nargs="+",
        type=Path,
        help="Campaign manifest JSON files, or directories containing them.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR / "report",
        help="Directory receiving report.json and index.html.",
    )
    parser.add_argument(
        "--luminosity-fb",
        type=positive_float,
        default=None,
        help=f"Integrated luminosity (default: manifest or {DEFAULT_LUMINOSITY_FB:g} fb^-1).",
    )
    parser.add_argument(
        "--signal-region",
        type=float,
        nargs=2,
        metavar=("LOW", "HIGH"),
        default=DEFAULT_SIGNAL_REGION_GEV,
        help="Open m4l signal-region bounds in GeV (default: 120 130).",
    )
    parser.add_argument(
        "--mass-range",
        type=float,
        nargs=2,
        metavar=("LOW", "HIGH"),
        default=DEFAULT_MASS_RANGE_GEV,
        help="Regular m4l spectrum bounds in GeV; under/overflow are retained.",
    )
    parser.add_argument(
        "--mass-bin-width",
        type=positive_float,
        default=DEFAULT_MASS_BIN_WIDTH_GEV,
        help="m4l spectrum bin width in GeV (default: 5).",
    )
    parser.add_argument("--title", default="LO H -> ZZ* -> 4l analysis")
    parser.add_argument("--closure-tolerance", type=positive_float, default=1.0e-9)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return a nonzero status after writing the report if an audit error is present.",
    )
    parser.add_argument(
        "--allow-unvalidated-bias",
        action="store_true",
        help=(
            "Explicitly include decay-biased samples whose weighted closure has not "
            "been validated. By default they are reported but excluded from totals."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report = generate_report(
        args.manifests,
        args.output_dir,
        luminosity_fb=args.luminosity_fb,
        signal_region_gev=(float(args.signal_region[0]), float(args.signal_region[1])),
        mass_range_gev=(float(args.mass_range[0]), float(args.mass_range[1])),
        mass_bin_width_gev=args.mass_bin_width,
        title=args.title,
        closure_tolerance=args.closure_tolerance,
        allow_unvalidated_bias=args.allow_unvalidated_bias,
    )
    print(f"Report JSON: {args.output_dir / 'report.json'}")
    print(f"Report HTML: {args.output_dir / 'index.html'}")
    errors = report["audit"]["errors"]
    if errors:
        print(f"Audit: {len(errors)} error(s)", file=sys.stderr)
        for message in errors:
            print(f"  - {message}", file=sys.stderr)
    return 2 if args.strict and errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
