#!/usr/bin/env python3
"""Combine independent loop diphoton samples in disjoint Born-mass regions.

Inputs must be completed campaign reports and separately validated Born maps.
The native generator rates stay fixed for the finite-source confidence bound.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

if str(Path(__file__).resolve().parents[1]) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hgammagamma import ho_shower_completion as completion


def source_statistics(selected, weights, sources, *, attempted, cross_section,
                      cross_section_error=0., confidence=.975):
    """Count independent sources, including zero-response unsaved attempts."""
    from scipy.stats import beta
    selected, weights, sources = np.asarray(selected, bool), np.asarray(weights), np.asarray(sources)
    if (selected.shape != weights.shape or sources.shape != weights.shape
            or sources.dtype.kind not in "iu" or np.any(sources < 0)
            or np.any(sources >= attempted) or attempted < 2
            or not np.all(np.isfinite(weights)) or np.any(weights < 0)
            or not math.isfinite(cross_section) or cross_section <= 0
            or not math.isfinite(cross_section_error) or cross_section_error < 0
            or not 0 < confidence < 1):
        raise ValueError("invalid independent positive-source sample")
    all_weights = np.bincount(sources, weights=weights, minlength=attempted)
    if not np.allclose(all_weights[np.unique(sources)], 1., rtol=1.e-9, atol=1.e-12):
        raise ValueError("saved source hypotheses do not sum to unit response")
    grouped = np.bincount(sources[selected], weights=weights[selected], minlength=attempted)
    if np.any(grouped > 1 + 1.e-10):
        raise ValueError("selected source response exceeds one")
    total, second = math.fsum(grouped), math.fsum(grouped * grouped)
    positive = int(np.count_nonzero(grouped > 0))
    variance = attempted / (attempted - 1) * max(0., second - total * total / attempted)
    probability = 1. if positive == attempted else float(beta.ppf(confidence, positive + 1, attempted - positive))
    return {
        "attempted_sources": attempted, "saved_sources": len(np.unique(sources)),
        "selected_positive_sources": positive, "selected_response_rows": int(np.count_nonzero(selected)),
        "sum_selected_response": total, "sum_squared_source_response": second,
        "effective_selected_sources": total * total / second if second else 0.,
        "normalization_denominator": attempted,
        "native_cross_section_pb": cross_section,
        "native_cross_section_error_pb": cross_section_error,
        "selected_cross_section_pb": cross_section * total / attempted,
        "selected_cross_section_mc_error_pb": cross_section / attempted * math.sqrt(variance),
        "selected_generator_integration_error_pb": cross_section_error * total / attempted,
        "positive_source_upper_bound": {
            "confidence_level": confidence, "selected_positive_sources": positive,
            "attempted_sources": attempted, "positive_response_probability_upper": probability,
            "selected_cross_section_upper_pb": cross_section * probability,
            "method": "exact one-sided Clopper-Pearson bound on P(selected response > 0), with response in [0,1]",
        },
    }


def combine_statistics(broad, focused):
    components = {"broad_outside": broad, "focused_inside": focused}
    if any(row["positive_source_upper_bound"]["confidence_level"] != .975 for row in components.values()):
        raise ValueError("combined 95% bound requires the declared 97.5% component bounds")
    return {
        "components": components,
        "selected_cross_section_pb": math.fsum(row["selected_cross_section_pb"] for row in components.values()),
        "selected_cross_section_mc_error_pb": math.sqrt(math.fsum(row["selected_cross_section_mc_error_pb"]**2 for row in components.values())),
        "selected_cross_section_upper_95_pb": math.fsum(row["positive_source_upper_bound"]["selected_cross_section_upper_pb"] for row in components.values()),
        "finite_source_covariance_pb2": 0.,
        "finite_source_covariance_reason": "independent generation seeds and source files in disjoint Born domains",
        "generator_integration_error_linear_envelope_pb": math.fsum(row["selected_generator_integration_error_pb"] for row in components.values()),
        "generator_integration_error_treatment": "component errors shown separately; linear envelope makes no independence assumption and is not a confidence bound",
        "upper_bound_coverage": "at least 95% coverage for this selection by the union bound on two one-sided 97.5% component bounds; not joint coverage across cuts and classifier",
    }


def verify_artifact(record):
    completion._validate_artifact(record, verify_file=True)
    return Path(record["path"]).resolve()


def normalized_herwig_card(text):
    # Only location, statistics, RNG and diagnostic controls are nonphysical.
    ignored = (":FileName ", ":NumberOfEvents ", ":RandomNumberGenerator:Seed ",
               ":MaxErrors ", ":PrintEvent ", ":DumpPeriod ")
    return [line for raw in text.splitlines()
            if (line := " ".join(raw.split("#", 1)[0].split()))
            and not any(token in line for token in ignored)]


def herwig_rng_seed(text):
    seeds = re.findall(r"^\s*set\s+\S+:RandomNumberGenerator:Seed\s+(\d+)\s*(?:#.*)?$", text, re.MULTILINE)
    if len(seeds) != 1:
        raise ValueError("the Herwig card must declare exactly one unambiguous RNG seed")
    return int(seeds[0])


def parse_run_card_text(text):
    result = {}
    for raw in text.splitlines():
        line = raw.split("!", 1)[0].split("#", 1)[0].strip()
        if "=" in line:
            value, name = line.rsplit("=", 1)
            if re.fullmatch(r"\w+", name.strip()):
                result[name.strip().lower()] = " ".join(value.lower().split())
    return result


def parse_run_card(path):
    return parse_run_card_text(path.read_text())


def embedded_run_card(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    header = []
    with opener(path, "rt") as stream:
        for line in stream:
            header.append(line)
            if "</MGRunCard>" in line:
                break
            if "<event>" in line or len(header) > 100000:
                raise ValueError("source has no bounded embedded MG run card")
    matches = re.findall(r"<MGRunCard>(.*?)</MGRunCard>", "".join(header), re.DOTALL)
    if len(matches) != 1:
        raise ValueError("source must contain exactly one embedded MG run card")
    return parse_run_card_text(matches[0]), hashlib.sha256(matches[0].encode()).hexdigest()


def validate_domains(broad, focused, low=100., high=150.):
    """Allow a mass restriction and independent sampling; reject physics drift."""
    for key in ("sample", "hard_accuracy", "hard_model", "process", "pdf", "lhaid", "ebeam_gev",
                "matching", "detector_response", "response_mode", "weight_scale"):
        if broad.get(key) != focused.get(key):
            raise ValueError(f"the two campaigns have different {key}")
    if broad["seed"] == focused["seed"] or broad["lhe"]["sha256"] == focused["lhe"]["sha256"]:
        raise ValueError("independent generation seeds and LHE files are required")
    allow = {"nevents", "iseed", "survey_splitting", "mmaa", "mmaamax"}
    settings = []
    for row in (broad, focused):
        settings.append({k: v for k, v in row["generator_settings"].items() if k not in allow})
    if settings[0] != settings[1]:
        raise ValueError("generation physics settings differ beyond the Born mass window")
    if (broad["generator_settings"].get("mmaa", 0.) != 0.
            or broad["generator_settings"].get("mmaamax", -1.) != -1.
            or focused["generator_settings"].get("mmaa") != low
            or focused["generator_settings"].get("mmaamax") != high):
        raise ValueError("native broad/focused Born mass domains differ from the partition")
    cards, embedded_proofs = [], []
    for row in (broad, focused):
        artifacts = row["generator_artifacts"]
        for item in artifacts.values():
            # Generator-artifact metadata uses a smaller schema than completion.
            actual = completion._artifact(Path(item["path"]))
            if actual["sha256"] != item["sha256"]:
                raise ValueError("generator card changed since campaign completion")
        cards.append(parse_run_card(Path(artifacts["run_card.dat"]["path"])))
        actual = cards[-1]
        def numeric(name):
            return float(actual[name].replace("d", "e"))
        for name, expected in (("mmaa", row["generator_settings"].get("mmaa", 0.)),
                               ("mmaamax", row["generator_settings"].get("mmaamax", -1.)),
                               ("nevents", row["nevents_requested"])):
            if numeric(name) != expected:
                raise ValueError(f"effective native {name} differs from the declared generation policy")
        # MG5 clears the mutable run-card seed to0 after consuming it. The
        # completed LHE embeds the actual initial run card with the input seed.
        embedded, digest = embedded_run_card(Path(row["lhe_file"]))
        if (float(embedded["iseed"]) != row["seed"] or numeric("iseed") not in (0, row["seed"])
                or {k: v for k, v in actual.items() if k != "iseed"}
                != {k: v for k, v in embedded.items() if k != "iseed"}):
            raise ValueError("embedded LHE run card does not prove the declared independent generation seed and physics")
        embedded_proofs.append({"seed": int(embedded["iseed"]), "card_sha256": digest,
                                "source_lhe_sha256": row["lhe"]["sha256"]})
    if {k: v for k, v in cards[0].items() if k not in allow} != {k: v for k, v in cards[1].items() if k not in allow}:
        raise ValueError("effective native run cards differ beyond permitted sampling settings")
    if broad["generator_artifacts"]["param_card.dat"]["sha256"] != focused["generator_artifacts"]["param_card.dat"]["sha256"]:
        raise ValueError("generator model parameter cards differ")
    shower_cards = []
    for row in (broad, focused):
        path = verify_artifact(row["shower_completion"]["artifacts"]["herwig_card"])
        actual = path.read_text()
        if actual != row["configuration"]["herwig_card"]:
            raise ValueError("actual Herwig card differs from the recorded configuration")
        if herwig_rng_seed(actual) != row["seed"]:
            raise ValueError("actual Herwig RNG seed differs from the independent declared seed")
        shower_cards.append(actual)
    if normalized_herwig_card(shower_cards[0]) != normalized_herwig_card(shower_cards[1]):
        raise ValueError("Herwig physics cards differ")
    response_sources, response_binaries = [], []
    for row in (broad, focused):
        directory = Path(row["analysis"]["input"]).parent.parents[5] / "LOAnalysis/Code"
        source = directory / "HwSimPostAnalysis_gammagamma_SSC.cc"
        includes = re.findall(r'^\s*#include\s+"([^"]+)"', source.read_text(), re.MULTILINE)
        response_sources.append({name: completion._artifact(directory/name)
                                 for name in [source.name, *includes]})
        response_binaries.append(completion._artifact(directory / "HwSimPostAnalysis_gammagamma_SSC"))
    if ({name: item["sha256"] for name, item in response_sources[0].items()}
            != {name: item["sha256"] for name, item in response_sources[1].items()}):
        raise ValueError("the archived detector-response sources differ")
    return {"same_hard_process_model_pdf_and_shower": True, "effective_native_run_cards_checked": True,
            "independent_seeds": [broad["seed"], focused["seed"]],
            "independent_actual_herwig_seeds": [herwig_rng_seed(text) for text in shower_cards],
            "embedded_lhe_run_cards": embedded_proofs,
            "matching_detector_response_sources": response_sources,
            "compiled_detector_response_artifacts": response_binaries,
            "independent_lhe_sha256": [broad["lhe"]["sha256"], focused["lhe"]["sha256"]],
            "allowed_run_card_differences": sorted(allow)}


def load_inputs(report_path, mapping_path):
    report, mapping = json.loads(report_path.read_text()), json.loads(mapping_path.read_text())
    if mapping.get("status") != "validated" or mapping.get("synthetic_skip_recovery_passed") is not True:
        raise ValueError("an independently validated Born map is required")
    for record in [*mapping["inputs"].values(), *mapping["raw_roots"], mapping["mapping_artifact"]]:
        verify_artifact(record)
    manifest_path = verify_artifact(report["metadata"]["loop_manifest"])
    if report["metadata"]["loop_manifest"]["sha256"] != mapping["inputs"]["manifest"]["sha256"]:
        raise ValueError("the final report and Born map refer to different campaigns")
    manifest = json.loads(manifest_path.read_text())
    proof = completion.validate_completion(manifest["shower_completion"], manifest, directory=manifest_path.parent)
    if (proof["attempted_events"] != mapping["attempted_events"]
            or proof["saved_events"] != mapping["saved_events"]
            or proof["source_weights"]["sum_weight"] != mapping["attempted_events"]):
        raise ValueError("Born-map populations differ from the positive unit source proof")
    partition = mapping["partition"]
    if (partition["low_gev"], partition["high_gev"], partition["high_edge_exclusive"]) != (100., 150., True):
        raise ValueError("unsupported Born-map partition")
    with np.load(mapping["mapping_artifact"]["path"], allow_pickle=False) as file:
        arrays = {key: file[key] for key in file.files}
    return report, mapping, manifest, arrays


def evaluate_component(inputs, model, model_metadata, model_sha, *, inside):
    import analyze_lo_varfiles as analysis
    from read_root_varfiles import read_ho_ROOT_varfile
    from hgammagamma.frozen_gammagamma_classifier import score_rows
    report, mapping, manifest, arrays = inputs
    directory = Path(mapping["inputs"]["manifest"]["path"]).parent
    sample, = analysis.discover_samples(directory.parents[2], manifest["run_tag"], {manifest["sample"]})
    if (not sample.requires_full_sample or sample.detector_response != "ssc" or sample.response_mode != "genuine"
            or Path(sample.var_file).resolve() != Path(mapping["inputs"]["response"]["path"]).resolve()
            or sample.normalization_sum_weight != mapping["attempted_events"]):
        raise ValueError("response or source normalization differs from the audited Born map")
    rows, weights, sources, entries = read_ho_ROOT_varfile(sample.var_file)
    if not np.array_equal(sources, arrays["response_sourceevent"]) or entries != list(range(len(rows))):
        raise ValueError("response rows do not match the audited Born-source association")
    cuts = [analysis.Cut.from_mapping(row) for row in report["metadata"]["cuts"]]
    eligible, scores = score_rows(model, rows)
    threshold = model_metadata["best_threshold"]
    classifier = report["classifier"]
    if (classifier["threshold"] != threshold or classifier["additional_mass_window"]
            or classifier["replay_diagnostics"]["files"]["signal_background_xgboost.json"] != model_sha):
        raise ValueError("the report uses a different frozen classifier or selection")
    selections = {"cuts": np.asarray(analysis.apply_cuts(rows, cuts), bool),
                  "frozen_classifier": eligible & (scores >= threshold)}
    original = {"cuts": report["loop_contribution"], "frozen_classifier": classifier["loop_contribution"]}
    region_mask = arrays["response_inside_window"] if inside else ~arrays["response_inside_window"]
    output = {}
    for name, selected in selections.items():
        common = dict(attempted=mapping["attempted_events"], cross_section=sample.cross_section_pb,
                      cross_section_error=sample.cross_section_error_pb or 0.)
        full = source_statistics(selected, weights, sources, **common)
        for key in ("selected_cross_section_pb", "selected_cross_section_mc_error_pb"):
            if not math.isclose(full[key], original[name][key], rel_tol=1.e-9, abs_tol=1.e-12):
                raise ValueError(f"the independently replayed {name} {key} does not reproduce its final report")
        output[name] = source_statistics(selected & region_mask, weights, sources, **common)
        output[name].update(run_tag=manifest["run_tag"], born_domain="100 <= m_gg < 150 GeV" if inside else "m_gg < 100 or m_gg >= 150 GeV",
                            discarded_sources=mapping["discarded_events"],
                            discarded_sources_in_domain=mapping["partition"]["discarded_inside"] if inside else mapping["discarded_events"]-mapping["partition"]["discarded_inside"],
                            original_report_reproduced=True, original_unpartitioned_selected_cross_section_pb=full["selected_cross_section_pb"])
    return output


def combine(broad_report, focused_report, broad_mapping, focused_mapping, frozen_classifier):
    from hgammagamma.frozen_gammagamma_classifier import load_frozen_classifier
    broad, focused = load_inputs(broad_report, broad_mapping), load_inputs(focused_report, focused_mapping)
    if broad[0]["metadata"]["cuts"] != focused[0]["metadata"]["cuts"]:
        raise ValueError("the frozen cut selections differ")
    if focused[1]["partition"]["attempted_inside"] != focused[1]["attempted_events"]:
        raise ValueError("a focused source event falls outside the exact common Born partition")
    domain_proof = validate_domains(broad[2], focused[2])
    model, metadata, _, replay = load_frozen_classifier(frozen_classifier)
    model_sha = replay["files"]["signal_background_xgboost.json"]
    outer = evaluate_component(broad, model, metadata, model_sha, inside=False)
    inner = evaluate_component(focused, model, metadata, model_sha, inside=True)
    inputs = {key: completion._artifact(path) for key, path in {
        "broad_report": broad_report, "focused_report": focused_report,
        "broad_mapping": broad_mapping, "focused_mapping": focused_mapping, "code": Path(__file__)}.items()}
    return {"schema_version": 1, "status": "validated", "energy_tev": 40.,
        "validation_scope": "source mapping, artifacts, normalization and selection replay; generator quality flags require a separate physics assessment",
        "frozen_model_sha256": model_sha, "frozen_threshold": metadata["best_threshold"],
        "inputs": inputs, "classifier_replay": replay, "generation_domain_validation": domain_proof,
        "partition": {"variable": "original LHE Born diphoton mass", "focused_low_gev": 100.,
            "focused_high_gev": 150., "high_edge_exclusive": True,
            "mapping": "unique ordered hard-photon four-vector match with independently audited shower skips"},
        "cuts": broad[0]["metadata"]["cuts"],
        "shower_quality": {name: {key: item[2]["shower_completion"].get(key)
            for key in ("attempted_events", "saved_events", "discarded_events", "termination",
                        "max_momentum_violation_mev", "exception_counts", "exception_classes")}
            for name, item in (("broad", broad), ("focused", focused))},
        "selections": {name: combine_statistics(outer[name], inner[name]) for name in outer},
        "normalization": "sum native sigma_i * selected response sum in its disjoint domain / all attempted source events_i",
        "uncertainty_scope": "finite independent source sampling; native generator rates fixed for MC bounds; no theory, detector or luminosity uncertainty",
        "failed_shower_treatment": "all consumed attempts stay in the denominator; failed showers have zero simulated response; bounds do not cover their unknown physical response",
        "classifier_selection": "unchanged frozen score threshold, with no added reconstructed-mass window",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("broad-report", "focused-report", "broad-mapping", "focused-mapping", "frozen-classifier", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("use a fresh combination artifact")
    result = combine(args.broad_report, args.focused_report, args.broad_mapping, args.focused_mapping, args.frozen_classifier)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print(json.dumps({"status": result["status"], "selections": result["selections"]}, indent=2))


if __name__ == "__main__":
    main()
