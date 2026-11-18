#!/usr/bin/env python3
"""Add a validated loop-induced continuum sample to a frozen SSC cut summary."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyze_lo_varfiles as analysis
from read_root_varfiles import read_ho_ROOT_varfile
from hgammagamma import ho_shower_completion as completion
from hgammagamma.run_loop_gammagamma import SAMPLE, PROCESS


def frozen_nominal(path):
    record = json.loads(Path(path).read_text())
    metadata, samples = record["metadata"], record["samples"]
    if metadata.get("mode") != "cuts" or metadata.get("detector_response") != "ssc":
        raise ValueError("the frozen input must be an SSC cut summary")
    if not metadata.get("cuts"):
        raise ValueError("the frozen input must specify its cuts")
    luminosity = float(metadata["luminosity_fb"])
    if not math.isfinite(luminosity) or luminosity <= 0:
        raise ValueError("invalid frozen luminosity")
    names = [item["sample"] for item in samples]
    if len(set(names)) != len(names) or SAMPLE.name in names:
        raise ValueError("the loop contribution is already included or sample names repeat")
    if "bkg_prompt_aa" not in names or not any(item["category"] == "Signal" for item in samples):
        raise ValueError("the frozen summary must include nominal diphoton continuum and signal")
    for item in samples:
        rate = float(item["selected_cross_section_pb"])
        if not math.isfinite(rate) or rate < 0 or item.get("detector_response") != "ssc":
            raise ValueError("invalid frozen sample rate or detector response")
        if not math.isclose(float(item["expected_events"]), rate * luminosity * 1000, rel_tol=1e-10):
            raise ValueError("frozen event yield differs from its cross section and luminosity")
    return record


def source_mc_error(selected, weights, source_events, *, attempted, cross_section, denominator):
    """Finite-source MC error, grouping correlated response rows by hard event.

    Failed shower events contribute zero. The generator integration error is
    reported separately; this quantity covers the finite acceptance sample.
    """
    grouped = defaultdict(float)
    for accept, weight, source in zip(selected, weights, source_events):
        if accept:
            grouped[source] += weight
    if attempted < 2 or len(grouped) > attempted:
        raise ValueError("invalid attempted source count for the MC error")
    total = math.fsum(grouped.values())
    second = math.fsum(weight * weight for weight in grouped.values())
    variance = attempted / (attempted - 1) * max(0., second - total * total / attempted)
    return abs(cross_section / denominator) * math.sqrt(variance)


def positive_source_upper_bound(selected, weights, source_events, *, attempted, cross_section):
    """Conservative one-sided 95% bound for a response weight in [0, 1].

    Its mean cannot exceed the probability of a positive selected response.
    Bound that probability using the exact binomial construction, counting
    source events rather than correlated detector-response rows.
    """
    from scipy.stats import beta
    total, accepted = defaultdict(list), defaultdict(list)
    for accept, weight, source in zip(selected, weights, source_events):
        if not math.isfinite(weight) or weight < 0:
            raise ValueError("positive-source bound requires finite nonnegative response weights")
        total[source].append(weight)
        if accept:
            accepted[source].append(weight)
    if len(total) > attempted or attempted < 1 or any(
            not math.isclose(math.fsum(values), 1., rel_tol=1e-9, abs_tol=1e-12) for values in total.values()):
        raise ValueError("positive-source bound requires a unit total response weight for each saved hard event")
    count = sum(math.fsum(values) > 0 for values in accepted.values())
    probability = 1. if count == attempted else float(beta.ppf(.95, count + 1, attempted - count))
    return {"confidence_level": .95, "selected_positive_sources": count,
        "attempted_sources": attempted, "maximum_response_weight_per_source": 1.,
        "positive_response_probability_upper": probability,
        "selected_cross_section_upper_pb": cross_section * probability,
        "method": "one-sided Clopper-Pearson bound on positive-source probability; selected weight <= 1",
        "scope": "finite independent unweighted hard events, with native generator rate fixed; no theory or detector uncertainty"}


def evaluate(analysis_root, run_tag, nominal_path, frozen_classifier=None):
    nominal = frozen_nominal(nominal_path)
    directory = Path(analysis_root).resolve() / "Backgrounds/events" / SAMPLE.name
    manifest_path = directory / "campaign.json"
    manifest = json.loads(manifest_path.read_text())
    if (manifest.get("sample") != SAMPLE.name or manifest.get("run_tag") != run_tag
            or manifest.get("hard_accuracy") != "loop-induced LO + PS"
            or manifest.get("process") != PROCESS or manifest.get("ebeam_gev") != 20000
            or manifest.get("weight_scale") != 1 or "analyze" not in manifest.get("completed", [])):
        raise ValueError("a completed 40 TeV continuum-only loop sample is required")
    proof = manifest.get("lhe", {}).get("process_validation", {})
    if proof.get("process") != PROCESS or proof.get("higgs_excluded") is not True:
        raise ValueError("the source lacks audited continuum process provenance")
    sample, = analysis.discover_samples(Path(analysis_root).resolve(), run_tag, {SAMPLE.name})
    if not sample.requires_full_sample or sample.detector_response != "ssc" or sample.response_mode != "genuine":
        raise ValueError("loop evaluation requires the validated full SSC genuine-photon response")
    source = manifest["shower_completion"].get("source_weights", {})
    attempts = manifest["shower_completion"]["attempted_events"]
    if any(not math.isclose(float(source.get(key, -1)), attempts, rel_tol=1e-10)
           for key in ("sum_weight", "sum_abs_weight", "sum_weight_squared")):
        raise ValueError("the finite-source MC error requires positive unweighted hard events")
    rows, weights, source_events, tree_entries = read_ho_ROOT_varfile(sample.var_file)
    if len(set(source_events)) != sample.events_read or len(weights) != len(source_events):
        raise ValueError("loop source-event population differs from the completed analysis")
    cuts = [analysis.Cut.from_mapping(item) for item in nominal["metadata"]["cuts"]]
    result = analysis._summarize_loaded_sample(sample, rows, weights, cuts, nominal["metadata"]["luminosity_fb"])
    result["selected_cross_section_mc_error_pb"] = source_mc_error(
        analysis.apply_cuts(rows, cuts), weights, source_events,
        attempted=attempts,
        cross_section=sample.cross_section_pb, denominator=sample.normalization_sum_weight)
    result["generator_cross_section_error_pb"] = sample.cross_section_error_pb
    result["finite_sample_upper_bound"] = positive_source_upper_bound(
        analysis.apply_cuts(rows, cuts), weights, source_events, attempted=attempts,
        cross_section=sample.cross_section_pb)
    combined = [*nominal["samples"], result]
    output = {
        "metadata": {"mode": "frozen_cuts_with_loop_continuum", "detector_response": "ssc",
            "luminosity_fb": nominal["metadata"]["luminosity_fb"], "cuts": nominal["metadata"]["cuts"],
            "loop_run_tag": run_tag, "hard_accuracy": manifest["hard_accuracy"],
            "generator_phase_space": manifest.get("generator_settings", {}),
            "born_mgg_gev": manifest.get("physics_scope", {}).get("born_mgg_gev", {"min": 0., "max": -1.}),
            "inclusive_loop_component": (manifest.get("generator_settings", {}).get("mmaa", 0.) == 0.
                and manifest.get("generator_settings", {}).get("mmaamax", -1.) == -1.),
            "normalization": "native loop-induced LO cross section times selected weight / all attempted source weight",
            "nominal_summary": completion._artifact(nominal_path),
            "loop_manifest": completion._artifact(manifest_path),
            "nominal_totals": analysis.compute_analysis_totals(nominal["samples"]),
            "totals": analysis.compute_analysis_totals(combined),
            "mc_error_scope": "finite independent hard events; generator integration error is separate",
            "classifier_evaluated": False},
        "samples": combined, "loop_contribution": result,
    }
    if frozen_classifier is not None:
        from hgammagamma.frozen_gammagamma_classifier import load_frozen_classifier, score_rows, FEATURE_NAMES
        model, metadata, ml_nominal, replay = load_frozen_classifier(frozen_classifier)
        if metadata["luminosity_fb"] != nominal["metadata"]["luminosity_fb"]:
            raise ValueError("frozen cut and classifier luminosities differ")
        if any(item["sample"] == SAMPLE.name for item in ml_nominal["samples"]):
            raise ValueError("the frozen classifier summary already includes the loop contribution")
        eligible, scores = score_rows(model, rows)
        selected = eligible & (scores >= metadata["best_threshold"])
        total = math.fsum(weight for weight, accept in zip(weights, selected) if accept)
        efficiency = total / sample.normalization_sum_weight
        rate = sample.cross_section_pb * efficiency
        contribution = {**result, "selected_entries": int(selected.sum()), "sum_selected_weight": total,
            "analysis_efficiency": efficiency, "selected_cross_section_pb": rate,
            "expected_events": rate * metadata["luminosity_fb"] * 1000,
            "mc_events_after_analysis": int(selected.sum()),
            "evaluation_partition": "new independent full sample", "partition_fraction": 1.,
            "selected_cross_section_mc_error_pb": source_mc_error(selected, weights, source_events,
                attempted=attempts, cross_section=sample.cross_section_pb, denominator=sample.normalization_sum_weight)}
        contribution["finite_sample_upper_bound"] = positive_source_upper_bound(
            selected, weights, source_events, attempted=attempts, cross_section=sample.cross_section_pb)
        output["metadata"]["classifier_evaluated"] = True
        output["classifier"] = {"feature_names": FEATURE_NAMES, "threshold": metadata["best_threshold"],
            "candidate_selection": "n_selected_photons >= 2", "additional_mass_window": False,
            "normalization": "full exposure; no training/test split or partition multiplier for new loop events",
            "replay_diagnostics": replay, "loop_contribution": contribution,
            "nominal_totals": analysis.compute_analysis_totals(ml_nominal["samples"]),
            "totals": analysis.compute_analysis_totals([*ml_nominal["samples"], contribution]),
            "nominal_evaluation": "saved test partition projected to full exposure"}
        output["_score_rows"] = [{"sourceevent": source_events[i], "tree_entry": tree_entries[i],
            "score": float(scores[i]), "selected": bool(selected[i]), "raw_weight": weights[i],
            "physical_weight": weights[i] * sample.cross_section_pb / sample.normalization_sum_weight
                * metadata["luminosity_fb"] * 1000,
            **{name: rows[i][name] for name in FEATURE_NAMES}} for i in range(len(rows)) if eligible[i]]
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", required=True, type=Path)
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--nominal-summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--frozen-classifier", type=Path, help="validated frozen model and saved-score reference bundle")
    args = parser.parse_args(argv)
    if args.output.exists():
        raise ValueError("output exists; choose a new evaluation artifact")
    score_path = args.output.with_name(args.output.stem + "_scores.csv")
    if args.frozen_classifier and score_path.exists():
        raise ValueError("score output exists; choose a new evaluation artifact")
    result = evaluate(args.analysis_root, args.run_tag, args.nominal_summary, args.frozen_classifier)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    score_rows = result.pop("_score_rows", None)
    if score_rows is not None:
        from hgammagamma.frozen_gammagamma_classifier import FEATURE_NAMES
        with score_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["sourceevent", "tree_entry", "score", "selected",
                "raw_weight", "physical_weight", *FEATURE_NAMES])
            writer.writeheader()
            writer.writerows(score_rows)
        result["classifier"]["scores"] = completion._artifact(score_path)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps(result["loop_contribution"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as error:
        raise SystemExit(str(error))
