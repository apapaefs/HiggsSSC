#!/usr/bin/env python3
"""Evaluate independent LO heavy-flavour controls with a frozen SSC classifier.

All new events are evaluation-only. Independent shards of the same process
estimate one cross section and are averaged by generated hard-event exposure.
Different processes are added. No training, split projection or threshold scan
is performed, and the nominal analysis artifacts are never changed.
"""
from __future__ import annotations

import argparse
import copy
import csv
import html
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hfourlepton import analyze_fourlepton as analysis
from hfourlepton import extend_xgboost_test as frozen
from hfourlepton import report_reader, reprocess as revision
from hfourlepton.replace_qq_backgrounds import PHYSICS_FIELDS
from ho_xgboost_core import signed_yield_summary

PROCESSES = frozenset({"bkg_zbb_ee", "bkg_zbb_mumu", "bkg_ttbar"})
CHANNELS = ("4e", "4mu", "2e2mu")
DETECTOR_FIELDS = ("collider", "selection_profile", "detector_profile", "detector", "cut_stages")
DECAY_BACKEND_FIELDS = ("heavy_flavour_decay_backend", "heavy_flavour_decay_table_sha256")
SHARD_FIELDS = (*PHYSICS_FIELDS, "decay_strategy", "decay", "branching_ratio_reweighter", *DECAY_BACKEND_FIELDS)
LEGACY_GENERATION_POLICY = {"name": "legacy_generator_defaults", "cut_overrides": {}}
BIAS_SUFFIX = "_hfbiased"


def physical_process(name):
    return name[:-len(BIAS_SUFFIX)] if name.endswith(BIAS_SUFFIX) else name


def is_biased(sample):
    return bool(sample.get("decay_bias", {}).get("enabled"))


def generation_policy(manifest):
    return manifest.get("reducible_control_partition", {}).get("generation_policy", LEGACY_GENERATION_POLICY)


def shard_proofs(index_paths, manifest_paths):
    """Read auditable half-open slices of common parent LHE populations."""
    requested = {str(path.resolve()) for path in manifest_paths}
    result, parents = {}, {}
    for index_path in index_paths:
        index = json.loads(index_path.read_text())
        for name, parent in index["samples"].items():
            parent_hash = parent["source_sha256"]
            if parent_hash not in parents:
                revision.verify(parent["source_lhe"], parent_hash)
                parents[parent_hash] = physical_process(name)
            elif parents[parent_hash] != physical_process(name):
                raise ValueError("one parent hard-event population is assigned to distinct processes")
            ranges = []
            for shard in parent["shards"]:
                start, stop, events = shard["start"], shard["stop"], shard["events"]
                if (not all(isinstance(value, int) for value in (start, stop, events))
                        or not 0 <= start < stop <= parent["source_events"] or stop-start != events):
                    raise ValueError("invalid parent hard-event slice")
                if any(start < old_stop and old_start < stop for old_start, old_stop in ranges):
                    raise ValueError("parent hard-event slices overlap")
                ranges.append((start, stop))
                manifest_path = str(Path(shard["manifest"]).resolve())
                if manifest_path not in requested:
                    continue
                key = (manifest_path, name)
                if key in result:
                    raise ValueError("duplicate shard provenance")
                revision.verify(shard["lhe"], shard["sha256"])
                result[key] = dict(parent_sha256=parent_hash, parent_source_events=parent["source_events"],
                    parent_cross_section_pb=parent["cross_section_pb"], start=start, stop=stop,
                    events=events, lhe_sha256=shard["sha256"], run_tag=shard["run_tag"],
                    generation_policy=parent.get("generation_policy", LEGACY_GENERATION_POLICY),
                    production_seed=parent.get("production_seed"),
                    actual_generation_cuts=parent.get("actual_generation_cuts"),
                    index=revision.artifact(index_path))
    return result


def validate_manifests(base, manifests, proofs=None, *, study_decay_bias=False):
    """Allow deliberately different LO hard physics, but identical response/cuts."""
    if not manifests:
        raise ValueError("at least one reducible manifest is required")
    base_samples = base["samples"]
    base_names = {sample["name"] for sample in base_samples}
    versions = {sample.get("selection_version") for sample in base_samples}
    if versions != {"ssc_v1"}:
        raise ValueError("frozen samples must consistently use ssc_v1 selection")
    if base.get("collider", {}).get("sqrt_s_tev") != 40.:
        raise ValueError("this first-pass evaluation requires 40 TeV")
    occupied_hashes = {sample.get(key) for sample in base_samples for key in
                       ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256")}
    occupied_hashes.discard(None)
    frozen_hashes = set(occupied_hashes)
    reference, seeds, parent_ranges = {}, {}, {}
    proofs = [{} for _ in manifests] if proofs is None else proofs
    if len(proofs) != len(manifests):
        raise ValueError("shard proof inventory differs from manifests")
    for manifest, manifest_proofs in zip(manifests, proofs):
        if manifest.get("status") != "complete" or manifest.get("order") != "lo":
            raise ValueError("reducible manifests must be complete LO campaigns")
        for field in DETECTOR_FIELDS:
            if field not in base or manifest.get(field) != base[field]:
                raise ValueError(f"reducible response/selection differs: {field}")
        names = [sample["name"] for sample in manifest["samples"]]
        allowed = PROCESSES | ({name + BIAS_SUFFIX for name in PROCESSES} if study_decay_bias else set())
        if not names or len(set(names)) != len(names) or set(names) - allowed:
            raise ValueError("declare only unique Zbb ee/mumu and ttbar controls in each manifest")
        if set(names) & base_names:
            raise ValueError("new background process already occurs in the frozen model inventory")
        for sample in manifest["samples"]:
            name = sample["name"]
            proof = manifest_proofs.get(name)
            if manifest.get("reducible_control_partition") and proof is None:
                raise ValueError("partitioned controls require matching --shard-index provenance")
            if (sample.get("status") != "complete" or sample.get("category") != "reducible"
                    or sample.get("perturbative_order") != "LO"
                    or sample.get("selection_version") != "ssc_v1"):
                raise ValueError(f"incomplete or incompatible reducible sample: {name}")
            bias = sample.get("decay_bias", {})
            biased = is_biased(sample)
            if biased:
                if not study_decay_bias:
                    raise ValueError("decay-biased controls require --study-decay-bias diagnostic mode")
                if (name != bias.get("bias_of", "") + BIAS_SUFFIX or bias.get("bias_of") not in PROCESSES
                        or bias.get("full_tail_support") is not True
                        or bias.get("scheme") != "full_support_direct_b_semileptonic_importance_sampling"
                        or sample.get("branching_ratio_reweighter") is not True):
                    raise ValueError("weighted controls require the declared full-support proposal and p/q reweighter")
                if sample.get("normalization", {}).get("denominator_kind") not in (
                        "predecay_lhe_handler_weight_sum", "predecay_lhe_sign_sum"):
                    raise ValueError("weighted controls require a pre-decay hard-event normalization denominator")
                if (sample.get("heavy_flavour_decay_backend") != "herwig_native_table"
                        or not sample.get("heavy_flavour_decay_table_sha256")
                        or not sample.get("heavy_flavour_decay_table")):
                    raise ValueError("weighted controls require explicit native Herwig decay-table provenance")
            elif (name.endswith(BIAS_SUFFIX) or sample.get("decay", {}).get("heavy_flavour_bias_enabled")
                  or manifest.get("decay_bias", {}).get("enabled")):
                raise ValueError("inconsistent weighted-stratum flags")
            hashes = [sample.get(key) for key in
                      ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256")]
            if any(not isinstance(value, str) or not value for value in hashes):
                raise ValueError(f"missing source provenance: {name}")
            if set(hashes) & occupied_hashes:
                raise ValueError(f"hard-event or ROOT overlap with another sample: {name}")
            occupied_hashes.update(hashes)
            seed = (proof.get("production_seed") if proof and proof.get("production_seed") is not None
                    else sample.get("seeds", {}).get("production"))
            parent = proof["parent_sha256"] if proof else None
            if proof:
                if (parent in frozen_hashes or proof["lhe_sha256"] not in hashes
                        or proof["run_tag"] != manifest.get("run_tag")
                        or proof["events"] != sample.get("normalization", {}).get("denominator_event_count")):
                    raise ValueError("shard proof differs from manifest source provenance")
                if proof.get("generation_policy", LEGACY_GENERATION_POLICY) != generation_policy(manifest):
                    raise ValueError("shard index and manifest generation policies differ")
                native = sample.get("normalization", {}).get("native_cross_section_pb")
                if native is None or not math.isclose(native, proof["parent_cross_section_pb"], rel_tol=1.e-10):
                    raise ValueError("partitioned shard changes the parent cross section")
                ranges = parent_ranges.setdefault(parent, [])
                start, stop = proof["start"], proof["stop"]
                if any(start < end and begin < stop for begin, end in ranges):
                    raise ValueError("supplied shards overlap in parent hard events")
                ranges.append((start, stop))
            if (not isinstance(seed, int) or isinstance(seed, bool) or seed <= 0
                    or ((name, seed) in seeds and (parent is None or seeds[name, seed] != parent))):
                raise ValueError(f"missing or reused production seed: {name}")
            seeds[name, seed] = parent
            if name in reference:
                old, old_manifest, old_proof = reference[name]
                for field in SHARD_FIELDS:
                    if old.get(field) != sample.get(field):
                        raise ValueError(f"independent shards change {name} {field}")
                for field in ("generation_cuts", "decay_bias", "heavy_flavour_bias"):
                    if old_manifest.get(field) != manifest.get(field):
                        raise ValueError(f"independent shards change {name} {field}")
                if generation_policy(old_manifest) != generation_policy(manifest):
                    raise ValueError(f"independent shards change {name} generation_policy")
                if proof and old_proof:
                    if proof.get("actual_generation_cuts") != old_proof.get("actual_generation_cuts"):
                        raise ValueError(f"independent shards change {name} actual generation cuts")
                    if proof["parent_sha256"] != old_proof["parent_sha256"] and (
                            proof.get("production_seed") is None or old_proof.get("production_seed") is None):
                        raise ValueError("independent parent runs require original production seeds in the shard indices")
            else:
                reference[name] = sample, manifest, proof
    for name, (weighted, weighted_manifest, weighted_proof) in reference.items():
        if not is_biased(weighted):
            continue
        control_name = physical_process(name)
        if control_name not in reference:
            continue
        control, control_manifest, control_proof = reference[control_name]
        fields = (*tuple(field for field in PHYSICS_FIELDS if field != "normalization_rule"), *DECAY_BACKEND_FIELDS)
        for field in fields:
            if weighted.get(field) != control.get(field):
                raise ValueError(f"weighted/control physics differs: {control_name} {field}")
        if generation_policy(weighted_manifest) != generation_policy(control_manifest):
            raise ValueError("weighted/control generation domains differ")
        if weighted_proof and control_proof:
            if weighted_proof.get("actual_generation_cuts") != control_proof.get("actual_generation_cuts"):
                raise ValueError("weighted/control actual generation cuts differ")
        # Shared parent input is supported only through the disjoint interval
        # checks above. Otherwise require distinct original production seeds.
        weighted_seed = (weighted_proof or {}).get("production_seed", weighted.get("seeds", {}).get("production"))
        control_seed = (control_proof or {}).get("production_seed", control.get("seeds", {}).get("production"))
        same_parent = (weighted_proof and control_proof
                       and weighted_proof["parent_sha256"] == control_proof["parent_sha256"])
        if weighted_seed == control_seed and not same_parent:
            raise ValueError("weighted/control independence requires different hard-event seeds or disjoint parent slices")


def rows_from_records(sample, records, stages, scale, luminosity, namespace):
    """Keep source identities and reconstructed channel; normalization uses all sources."""
    sources = [int(row["source_index"]) for row in records.source_events]
    if len(sources) != len(set(sources)) or any(source < 0 for source in sources):
        raise ValueError("duplicate or invalid SourceEvents identifiers")
    population = set(sources)
    rows = []
    group = sample.raw.get("analysis_source_group", sample.production_group_id) or sample.name
    for row in records.four_lepton:
        source = int(row["source_index"])
        if source not in population:
            raise ValueError("hypothesis has no SourceEvents entry")
        if not analysis.report.passes_stage(row, stages[-1]):
            continue
        pts = sorted((float(value) for value in row.get("lepton_pt", [])), reverse=True)
        if len(pts) != 4:
            raise ValueError("selected hypothesis requires four ordered lepton momenta")
        values = {name: float(row[name]) for name in ("m4l", "mZ1", "mZ2", "pt4l", "min_delta_r")}
        values.update(abs_y4l=abs(float(row["y4l"])), **{f"pt{i+1}": pt for i, pt in enumerate(pts)})
        weight = float(row["event_weight"]) * scale * luminosity * 1000.
        if not all(math.isfinite(value) for value in (*values.values(), weight)):
            raise ValueError("nonfinite reducible kinematics or signed weight")
        channel = analysis.report.normalize_channel(row.get("channel"), sample.channel)
        if channel not in CHANNELS:
            raise ValueError("selected hypothesis has an unknown reconstructed channel")
        rows.append(dict(sample=sample.name, shard=namespace, group=f"{namespace}/{group}",
                         source=source, channel=channel, weight=weight, **values))
    return rows


def load_shard(path, luminosity, proofs=None, *, study_decay_bias=False):
    bundle = analysis.report.load_manifest(path)
    errors = list(bundle.audit.errors)
    namespace = analysis.digest(path)
    rows, metadata = [], []
    for sample in bundle.samples:
        print(f"Auditing {path.name}: {sample.name}", flush=True)
        for key, artifact in (("output_root_sha256", sample.output_root),
                              ("analysis_summary_sha256", Path(sample.raw["analysis_summary"]))):
            if not artifact or not artifact.is_file() or analysis.digest(artifact) != sample.raw.get(key):
                raise ValueError(f"missing or changed analysis product: {sample.name} {key}")
        if sample.raw.get("heavy_flavour_decay_backend") == "herwig_native_table":
            revision.verify(sample.raw["heavy_flavour_decay_table"], sample.raw["heavy_flavour_decay_table_sha256"])
        records = report_reader.load_report_records(sample.output_root,
            selected_mask=analysis.report._final_stage(bundle.cut_stages).required_mask, progress=True)
        result = analysis.report.analyze_sample(sample, records, bundle.cut_stages, luminosity,
                                               (120., 130.), 1.e-9,
                                               allow_unvalidated_bias=study_decay_bias and sample.decay_bias_enabled)
        errors.extend(f"{sample.name}: {error}" for error in result["audit"]["errors"])
        normalization = result["normalization"]
        scale = normalization.get("scale_pb_per_weight")
        exposure = normalization.get("denominator_event_count")
        if (scale is None or exposure is None or not math.isfinite(scale + exposure)
                or scale <= 0 or exposure <= 0 or int(exposure) != exposure):
            raise ValueError(f"invalid signed normalization or generated exposure: {sample.name}")
        rows.extend(rows_from_records(sample, records, bundle.cut_stages, scale, luminosity, namespace))
        cutflow = copy.deepcopy(result["cutflow"])
        for row, stage in zip(cutflow, bundle.cut_stages):
            row["independent_source_events"] = len({int(record["source_index"]) for record in records.four_lepton
                if analysis.report.passes_stage(record, stage)})
        generated = result["generated"]
        source_scale = scale * luminosity * 1000. * (sample.weight_scale if sample.weight_scale is not None else 1.)
        cutflow.insert(0, dict(name="source_normalization", label="Inclusive source weights before response", bit=-1,
            yield_events=generated["sumw"] * source_scale,
            yield_sumw2=generated["sumw2"] * source_scale**2,
            independent_source_events=generated["entries"]))
        metadata.append(dict(sample=sample.name, shard=namespace, manifest=str(path),
            source_group=sample.raw.get("analysis_source_group", sample.production_group_id) or sample.name,
            generated_exposure=exposure, observed_source_events=len(records.source_events),
            normalization=normalization, native_cross_section_pb=sample.native_cross_section_pb,
            input_lhe_sha256=sample.raw["input_lhe_sha256"],
            source_lhe_sha256=sample.raw["source_lhe_sha256"],
            root_sha256=sample.raw["output_root_sha256"],
            analysis_summary_sha256=sample.raw["analysis_summary_sha256"],
            parent_slice=(proofs or {}).get(sample.name),
            decay_bias=result["decay_bias"],
            sampling_stratum="weighted_diagnostic" if sample.decay_bias_enabled else "unbiased",
            decay_backend={key: sample.raw.get(key) for key in (*DECAY_BACKEND_FIELDS, "heavy_flavour_decay_table")},
            source_weight_diagnostic=result["generated"],
            cutflow=cutflow,
            audit_warnings=result["audit"].get("warnings", [])))
    if errors:
        raise ValueError("reducible audit failed:\n" + "\n".join(errors))
    return rows, metadata


def combine_shards(rows, metadata):
    """Exposure-average repeated processes; namespace source IDs across shards."""
    totals, parent_exposures, parent_denominators, seen = {}, {}, {}, set()
    for item in metadata:
        key = item["sample"], item["shard"]
        if key in seen:
            raise ValueError("duplicate shard metadata")
        seen.add(key)
        exposure = item["generated_exposure"]
        if not math.isfinite(exposure) or exposure <= 0:
            raise ValueError("generated exposure must be positive")
        totals[item["sample"]] = totals.get(item["sample"], 0.) + exposure
        proof = item.get("parent_slice")
        if proof:
            parent = (item["sample"], proof["parent_sha256"])
            denominator = item["normalization"]["denominator_sumw"]
            if not math.isfinite(denominator) or denominator <= 0:
                raise ValueError("partitioned controls require a positive signed normalization denominator")
            parent_exposures[parent] = parent_exposures.get(parent, 0.) + exposure
            parent_denominators[parent] = parent_denominators.get(parent, 0.) + denominator
    coefficients = {}
    for item in metadata:
        proof = item.get("parent_slice")
        coefficient = item["generated_exposure"] / totals[item["sample"]]
        if proof:
            parent = (item["sample"], proof["parent_sha256"])
            coefficient = (parent_exposures[parent] / totals[item["sample"]]
                           * item["normalization"]["denominator_sumw"] / parent_denominators[parent])
        coefficients[item["sample"], item["shard"]] = coefficient
    combined = []
    for row in rows:
        key = row["sample"], row["shard"]
        if key not in coefficients:
            raise ValueError("row has no declared shard exposure")
        combined.append({**row, "weight": row["weight"] * coefficients[key]})
    dtypes = {"sample": str, "shard": str, "group": str, "source": np.int64, "channel": str,
              "weight": float, **{name: float for name in analysis.VARIABLES}}
    data = {key: np.asarray([row[key] for row in combined], dtype=dtype) for key, dtype in dtypes.items()}
    receipt = [{**item, "combination_coefficient": coefficients[item["sample"], item["shard"]],
                "total_process_generated_exposure": totals[item["sample"]]} for item in metadata]
    return data, receipt


def pooled_cutflow(metadata):
    result = {}
    for item in metadata:
        coefficient = item["combination_coefficient"]
        process = result.setdefault(item["sample"], {})
        for row in item.get("cutflow", []):
            pooled = process.setdefault(row["name"], dict(name=row["name"], label=row["label"],
                bit=row["bit"], yield_events=0., yield_sumw2=0., independent_source_events=0))
            if row["label"] != pooled["label"] or row["bit"] != pooled["bit"]:
                raise ValueError("shard cutflow definitions differ")
            pooled["yield_events"] += coefficient * row["yield_events"]
            pooled["yield_sumw2"] += coefficient**2 * row["yield_sumw2"]
            pooled["independent_source_events"] += row["independent_source_events"]
    for process in result.values():
        for row in process.values():
            row["mc_error"] = math.sqrt(row["yield_sumw2"])
            row["effective_count"] = row["yield_events"]**2 / row["yield_sumw2"] if row["yield_sumw2"] else 0.
    return {name: list(process.values()) for name, process in result.items()}


def selected_yields(data, mask, processes):
    def summary(selection):
        return signed_yield_summary(data["weight"], data["group"], data["source"], selected=selection)
    return dict(total=summary(mask),
        by_process={name: summary(mask & (data["sample"] == name)) for name in sorted(processes)},
        by_channel={channel: summary(mask & (data["channel"] == channel)) for channel in CHANNELS},
        by_process_channel={name: {channel: summary(mask & (data["sample"] == name)
                                                      & (data["channel"] == channel))
                                  for channel in CHANNELS} for name in sorted(processes)})


def closure_metric(weighted, control, luminosity, minimum_effective=25.):
    """Diagnostic comparison of independent strata, never a closure approval."""
    scale = luminosity * 1000.
    def as_cross_section(item):
        return dict(cross_section_pb=item["sum_weight"] / scale,
                    cross_section_sumw2_pb2=item["variance"] / scale**2)
    result = analysis.report._heavy_flavour_closure_metric(as_cross_section(weighted), as_cross_section(control))
    support = all(item["effective_count"] >= minimum_effective and item["selected_source_events"] > 0
                  for item in (weighted, control))
    result.update(weighted_effective_count=weighted["effective_count"],
                  unbiased_effective_count=control["effective_count"],
                  weighted_source_events=weighted["selected_source_events"],
                  unbiased_source_events=control["selected_source_events"],
                  statistically_supported=support, minimum_effective_events=minimum_effective,
                  closure_validated=False,
                  interpretation="Independent-stratum weighted comparison only; sparse or empty selections do not validate the proposal.")
    result["difference_95pct_gaussian_interval_pb"] = (
        [result["difference_pb"] + direction * 1.959963984540054 * result["combined_stat_uncertainty_pb"]
         for direction in (-1., 1.)] if support else None)
    result["interval_definition"] = "Approximate marginal Gaussian interval for the difference of independent MC estimates; omitted when either stratum lacks effective-event support."
    result["moment_based_difference_diagnostic"] = result["difference_over_combined_stat_uncertainty"]
    if not support:
        result["status"] = "insufficient_mc_statistics_for_closure"
        result["difference_over_combined_stat_uncertainty"] = None
    return result


def closure_diagnostics(result, cutflows, luminosity, minimum_effective=25.):
    comparisons = {}
    weighted = result.get("weighted_diagnostics")
    if weighted is None:
        return comparisons
    for name in weighted["cuts"]["by_process"]:
        control = physical_process(name)
        if control not in result["cuts"]["by_process"]:
            comparisons[name] = dict(status="missing_unbiased_control", closure_validated=False)
            continue
        output = {selection: closure_metric(weighted[selection]["by_process"][name],
                    result[selection]["by_process"][control], luminosity, minimum_effective)
                  for selection in ("preselection", "cuts", "frozen_classifier")}
        output["by_channel"] = {selection: {channel: closure_metric(
            weighted[selection]["by_process_channel"][name][channel],
            result[selection]["by_process_channel"][control][channel], luminosity, minimum_effective)
            for channel in CHANNELS} for selection in ("cuts", "frozen_classifier")}
        control_flow = {row["name"]: row for row in cutflows.get(control, [])}
        def cut_yield(row):
            return dict(sum_weight=row["yield_events"], variance=row["yield_sumw2"],
                        effective_count=row["effective_count"], selected_source_events=row["independent_source_events"])
        output["cutflow"] = [{"name": row["name"], **closure_metric(cut_yield(row), cut_yield(control_flow[row["name"]]),
                              luminosity, minimum_effective)}
                             for row in cutflows.get(name, []) if row["name"] in control_flow]
        output.update(control_sample=control, closure_validated=False,
                      covariance_assumption="disjoint_hard_event_populations_zero_covariance")
        comparisons[name] = output
    return comparisons


def sideband_diagnostics(metadata, luminosity, minimum_effective=25., reader=None):
    """Normalize pretrigger sideband moments separately for each sampling arm."""
    if reader is None:
        from hfourlepton.reducible_isolation_sidebands import shard_moments
        reader = shard_moments
    results = {}
    for item in metadata:
        name, count, cross_section, moments = reader(item["manifest"])
        if (name != item["sample"] or count != item["generated_exposure"]
                or not math.isclose(cross_section, item["native_cross_section_pb"], rel_tol=1.e-10)):
            raise ValueError("sideband reader differs from audited sample normalization")
        scale = item["normalization"]["scale_pb_per_weight"] * luminosity * 1000.
        scale *= item["combination_coefficient"]
        dest = results.setdefault(name, {})
        for region, (sumw, sumw2, count) in moments.items():
            if not all(math.isfinite(value) for value in (sumw, sumw2, count)) or sumw2 < 0 or count < 0:
                raise ValueError("invalid source-group sideband moments")
            pooled = dest.setdefault(region, dict(sum_weight=0., variance=0., selected_source_events=0))
            pooled["sum_weight"] += scale * sumw
            pooled["variance"] += scale**2 * sumw2
            pooled["selected_source_events"] += int(count)
    for regions in results.values():
        for value in regions.values():
            value["mc_error"] = math.sqrt(value["variance"])
            value["effective_count"] = value["sum_weight"]**2 / value["variance"] if value["variance"] else 0.
            value["rate_fb"] = value["sum_weight"] / luminosity
            value["mc_error_fb"] = value["mc_error"] / luminosity
    comparisons = {}
    for name, regions in results.items():
        if not name.endswith(BIAS_SUFFIX):
            continue
        control = physical_process(name)
        comparisons[name] = {region: closure_metric(weighted, results[control][region], luminosity, minimum_effective)
                             for region, weighted in regions.items() if control in results and region in results[control]}
    return dict(definition="Pretrigger four leading crack-accepted leptons in failed-isolation hypotheses; massless nearest-Z pairing. No nominal isolation or classifier selection is replaced.",
                strata=results, weighted_control_closure=comparisons, closure_validated=False,
                covariance_assumption="disjoint_hard_event_populations_zero_covariance")


def evaluate_data(data, metadata, card, features, threshold, model, *, study_decay_bias=False):
    if features != card.get("features") or not features or set(features) - set(analysis.FEATURES):
        raise ValueError("frozen feature order differs from the saved analysis card")
    if not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("invalid frozen classifier threshold")
    eligible = analysis.cut_mask(data, card.get("preselection", {}))
    baseline = eligible & analysis.cut_mask(data, card.get("cuts", {}))
    scores = np.full(len(data["weight"]), np.nan)
    if np.any(eligible):
        prediction = np.asarray(model.predict_proba(np.column_stack([data[name][eligible] for name in features])))
        if prediction.shape != (int(eligible.sum()), 2) or not np.all(np.isfinite(prediction)):
            raise ValueError("invalid frozen-model probabilities")
        if np.any((prediction < 0) | (prediction > 1)):
            raise ValueError("frozen-model probabilities are outside [0,1]")
        scores[eligible] = prediction[:, 1]
    processes = {item["sample"] for item in metadata}
    weighted_processes = {name for name in processes if name.endswith(BIAS_SUFFIX)}
    if weighted_processes and not study_decay_bias:
        raise ValueError("weighted samples require explicit diagnostic-only evaluation")
    def selected(process_set):
        stratum = np.isin(data["sample"], sorted(process_set))
        return dict(preselection=selected_yields(data, eligible & stratum, process_set),
                    cuts=selected_yields(data, baseline & stratum, process_set),
                    frozen_classifier=selected_yields(data, baseline & (scores >= threshold) & stratum, process_set))
    result = selected(processes - weighted_processes)
    if study_decay_bias:
        result["weighted_diagnostics"] = selected(weighted_processes)
        result["weighted_diagnostics"].update(included_in_nominal_totals=False, closure_validated=False)
    return result, scores, eligible, baseline


def run(manifest_paths, frozen_directory, output, index_paths=(), *, study_decay_bias=False, isolation_sidebands=False):
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory is nonempty; choose a new output directory")
    paths = [path.resolve() for path in manifest_paths]
    if len(paths) != len(set(paths)):
        raise ValueError("the same manifest was supplied more than once")
    summary, base = frozen.frozen_analysis(frozen_directory)
    raw = [json.loads(path.read_text()) for path in paths]
    indexed = shard_proofs(index_paths, paths)
    proofs = [{name: proof for (manifest, name), proof in indexed.items() if manifest == str(path)} for path in paths]
    validate_manifests(base, raw, proofs, study_decay_bias=study_decay_bias)
    luminosity = float(summary["luminosity_fb"])
    if not math.isfinite(luminosity) or luminosity <= 0:
        raise ValueError("invalid frozen luminosity")
    rows, metadata = [], []
    for path, proof in zip(paths, proofs):
        shard_rows, shard_metadata = load_shard(path, luminosity, proof, study_decay_bias=study_decay_bias)
        rows.extend(shard_rows)
        metadata.extend(shard_metadata)
    data, metadata = combine_shards(rows, metadata)
    del rows
    from xgboost import XGBClassifier
    model = XGBClassifier()
    model.load_model(frozen_directory / "model.json")
    card = copy.deepcopy(summary["card"])
    threshold = summary["threshold"]["threshold"]
    result, scores, eligible, baseline = evaluate_data(data, metadata, card, summary["features"], threshold, model,
                                                     study_decay_bias=study_decay_bias)
    cutflows = pooled_cutflow(metadata)
    if study_decay_bias:
        result["weighted_control_closure"] = closure_diagnostics(result, cutflows, luminosity,
            float(card.get("min_background_effective_events", 25.)))
    if isolation_sidebands:
        result["isolation_sidebands"] = sideband_diagnostics(metadata, luminosity,
            float(card.get("min_background_effective_events", 25.)))
    result.update(schema_version=1, mode="independent_reducible_evaluation", luminosity_fb=luminosity,
        threshold=threshold, features=summary["features"], frozen_card=card,
        scope="First-pass heavy-flavour estimate from unbiased LO Zbb and ttbar controls; not a complete experimental reducible background.",
        limitations=["Generic jet-to-lepton misidentification, material conversions and impact-parameter response are not modeled.",
                     "Zero selected events do not establish a zero background or a statistical upper bound."],
        evaluation_policy="All new source events are evaluation-only; no fit, threshold scan, partition projection, or nominal-artifact changes.",
        combination_policy="Disjoint slices of one parent are pooled by their signed normalization denominators. Independent same-process parents are averaged by all generated hard-event counts; distinct processes are summed. Source-group weights are summed before squaring for MC variance.",
        source_manifests=[revision.artifact(path) for path in paths], shards=metadata,
        cutflow_by_process=cutflows,
        missing_planned_processes=sorted(PROCESSES - {item["sample"] for item in metadata}),
        frozen_artifacts={name: revision.artifact(frozen_directory / name) for name in
                          ("completion.json", "summary.json", "model.json")},
        software={name: revision.artifact(Path(module.__file__).resolve()) for name, module in
                  (("analysis", analysis), ("report", analysis.report), ("reader", report_reader),
                   ("frozen", frozen))})
    result["software"]["evaluator"] = revision.artifact(Path(__file__).resolve())
    if isolation_sidebands:
        from hfourlepton import reducible_isolation_sidebands as sidebands_module
        result["software"]["isolation_sidebands"] = revision.artifact(Path(sidebands_module.__file__).resolve())
    if study_decay_bias:
        result["scope"] = "First-pass heavy-flavour controls and separate importance-weighted closure diagnostics; top-level yields contain unbiased samples only."
        result["combination_policy"] += " Unbiased and importance-weighted strata are never averaged or summed together."
    output.mkdir(parents=True, exist_ok=True)
    revision.save(output / "summary.json", result)
    with (output / "yields.csv").open("w", newline="") as stream:
        fields = ["stratum", "selection", "process", "channel", *result["cuts"]["total"]]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        sections_to_write = [("unbiased", result)]
        if study_decay_bias:
            sections_to_write.append(("weighted_diagnostic", result["weighted_diagnostics"]))
        for stratum, sections in sections_to_write:
            for selection in ("preselection", "cuts", "frozen_classifier"):
                section = sections[selection]
                writer.writerow(dict(stratum=stratum, selection=selection, process="all", channel="all", **section["total"]))
                for name, value in section["by_process"].items():
                    writer.writerow(dict(stratum=stratum, selection=selection, process=name, channel="all", **value))
                    for channel, channel_value in section["by_process_channel"][name].items():
                        writer.writerow(dict(stratum=stratum, selection=selection, process=name, channel=channel, **channel_value))
    with (output / "scores.csv").open("w", newline="") as stream:
        fields = ("sample", "shard", "group", "source", "channel", "m4l", "weight")
        writer = csv.writer(stream)
        writer.writerow([*fields, "score", "baseline"])
        for index in np.flatnonzero(eligible):
            writer.writerow([*(data[field][index] for field in fields), scores[index], bool(baseline[index])])
    with (output / "cutflow.csv").open("w", newline="") as stream:
        fields = ("process", "name", "label", "bit", "yield_events", "yield_sumw2",
                  "independent_source_events", "mc_error", "effective_count")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for name, cutflow in result["cutflow_by_process"].items():
            for row in cutflow:
                writer.writerow(dict(process=name, **row))
    brief = {key: result[key] for key in ("scope", "limitations", "threshold", "preselection", "cuts", "frozen_classifier")}
    if study_decay_bias:
        brief.update(weighted_diagnostics=result["weighted_diagnostics"], weighted_control_closure=result["weighted_control_closure"])
    if isolation_sidebands:
        brief["isolation_sidebands"] = result["isolation_sidebands"]
    (output / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>Reducible backgrounds</title>'
        '<h1>Heavy-flavour background estimate</h1><p><a href="summary.json">Provenance and yields</a> · '
        '<a href="yields.csv">CSV</a> · <a href="cutflow.csv">Cutflow</a> · <a href="scores.csv">Frozen scores</a></p><pre>'
        + html.escape(json.dumps(brief, indent=2)) + '</pre>')
    revision.save(output / "completion.json", dict(result="complete", files={name: analysis.digest(output / name)
        for name in ("summary.json", "yields.csv", "cutflow.csv", "scores.csv", "index.html")}))
    print(output / "index.html")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, nargs="+", action="extend", required=True,
                        help="completed independent LO control manifests; repeat or supply several paths")
    parser.add_argument("--frozen-analysis", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=Path, nargs="+", action="extend", default=[],
                        help="orchestrator indices proving disjoint slices of common parent LHEs")
    parser.add_argument("--study-decay-bias", action="store_true",
                        help="evaluate separate full-support importance-weighted strata and unbiased closure controls; never promote weighted estimates")
    parser.add_argument("--isolation-sidebands", action="store_true",
                        help="also compare stored pretrigger loose-isolation moments within each separate sampling stratum")
    args = parser.parse_args(argv)
    run(args.manifest, args.frozen_analysis.resolve(), args.output_dir.resolve(), args.shard_index,
        study_decay_bias=args.study_decay_bias, isolation_sidebands=args.isolation_sidebands)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(str(error))
