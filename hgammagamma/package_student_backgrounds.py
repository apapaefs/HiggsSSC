#!/usr/bin/env python3
"""Publish portable event-level tables from the completed 40 TeV background analyses."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np


FEATURES = ["m_gg", "pt_gamma1", "eta_gamma1", "pt_gamma2", "eta_gamma2",
            "deltaR_gg", "deltaPhi_gg", "pt_gg", "y_gg", "n_selected_photons"]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verified_path(record):
    path = Path(record["path"])
    if sha256(path) != record["sha256"]:
        raise ValueError(f"source differs from its completed analysis: {path}")
    return path


def check_rate(actual, expected, context):
    if not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-13):
        raise ValueError(f"{context}: {actual} pb differs from verified {expected} pb")


def selected_rates(variables, weights, model, threshold):
    eligible = variables[:, -1] >= 2
    pair = variables[eligible]
    weights = weights[eligible]
    score = model.predict_proba(pair[:, :-1].astype(np.float32))[:, 1] if len(pair) else np.array([])
    if not np.all(np.isfinite(score)):
        raise ValueError("invalid frozen diphoton score")
    window = (pair[:, 0] >= 120) & (pair[:, 0] <= 130)
    return dict(pair_rate_pb=float(np.sum(weights)),
                mass_window_rate_pb=float(np.sum(weights[window])),
                frozen_classifier_rate_pb=float(np.sum(weights[score >= threshold])))


def loop_rows(combination):
    import uproot

    pieces = []
    provenance = []
    offset = 0
    for component, region in (("broad", "broad_outside"), ("focused", "focused_inside")):
        record = combination["inputs"][component + "_mapping"]
        mapping = json.loads(verified_path(record).read_text())
        response = verified_path(mapping["inputs"]["response"])
        source = combination["selections"]["cuts"]["components"][region]
        n = source["attempted_sources"]
        if (mapping["status"] != "validated" or mapping["attempted_events"] != n
                or n != source["normalization_denominator"] or mapping["run_tag"] != source["run_tag"]):
            raise ValueError(f"invalid source exposure in {region}")
        with uproot.open(response) as root:
            tree = root["Data2"]
            variables = tree["variables"].array(library="np")
            weights = tree["eventweight"].array(library="np").ravel()
            row_source = tree["sourceevent"].array(library="np").ravel()
        with np.load(verified_path(mapping["mapping_artifact"]), allow_pickle=False) as mapped:
            if not np.array_equal(row_source, mapped["response_sourceevent"]):
                raise ValueError(f"source rows differ from Born mapping for {component}")
            take = mapped["response_inside_window"]
            if component == "broad":
                take = ~take
            selected_sources = mapped["response_lhe_event_index"][take]
        if np.any(selected_sources < 0) or np.any(selected_sources >= n):
            raise ValueError("invalid hard-event source ID")
        pieces.append((variables[take], selected_sources + offset,
                       weights[take] * source["native_cross_section_pb"] / n))
        provenance.append(dict(region=region, run_tag=source["run_tag"], generated=n,
                               native_cross_section_pb=source["native_cross_section_pb"],
                               response_root_sha256=mapping["inputs"]["response"]["sha256"],
                               source_lhe_sha256=mapping["inputs"]["lhe"]["sha256"],
                               born_mapping_sha256=mapping["mapping_artifact"]["sha256"]))
        offset += n
    columns = dict(zip(("variables", "sourceevent", "weight_pb"),
                       (np.concatenate(column) for column in zip(*pieces))))
    return columns, provenance, offset


def jet_shard(index, shard, directory, analyzer, shards, model, threshold):
    run = directory / shard["directory"]
    inputs = list((run / "events").glob("*.root"))
    if len(inputs) != 1:
        raise ValueError(f"expected one HwSim ROOT file in {run}")
    root, = inputs
    n = shard["generated"]
    with tempfile.TemporaryDirectory() as tmp:
        moments, features = Path(tmp) / "moments.csv", Path(tmp) / "pairs.csv"
        subprocess.run([analyzer, "--dijet-bound", root, moments, str(shard["seed"] + 100), features], check=True)
        totals = np.genfromtxt(moments, names=True, delimiter=",", dtype=None, encoding="utf8", ndmin=1)
        pairs = np.genfromtxt(features, names=True, delimiter=",", ndmin=1)
    if not np.all(totals["events_read"] == n):
        raise ValueError(f"incomplete ROOT population for {run}")
    columns = np.column_stack([*(pairs[key] for key in FEATURES[:-1]), np.full(len(pairs), 2)])
    ids = pairs["source"].astype(np.int64)
    if np.any(ids != pairs["source"]) or np.any(ids < 0) or np.any(ids >= n):
        raise ValueError(f"invalid source-event index for {run}")
    # Each accepted pair carries a union-bound weight. If the cap at one were
    # active anywhere, a static table would need an explicit per-cut reducer.
    upper = pairs["fq"]**2
    lower = pairs["fg"]**2 * (pairs["eligible_jets"] == 2)
    if np.any(lower > upper) or np.any(np.bincount(ids, weights=upper, minlength=n) > 1):
        raise ValueError(f"double-fake probability exceeds the additive regime in {run}")
    factor_pb = 1.e3 * shard["generated_cross_section_nb"] / n
    for weight, scenario in ((upper, "all_quark_independent"), (lower, "all_gluon_two_eligible")):
        checks = selected_rates(columns, factor_pb * weight, model, threshold)
        for selection, prefix in (("mass_window_rate_pb", ""), ("frozen_classifier_rate_pb", "frozen_")):
            expected_pb = shard["estimates"][prefix + scenario + "@1"]["rate_fb"] / 1000
            check_rate(checks[selection], expected_pb, f"{shard['directory']} {scenario}/{selection}")
    # Equal-size independent samples are averaged, each with its own native
    # cross section. IDs cover the complete generated range, including gaps.
    columns = (columns, ids + index * n, factor_pb / shards * upper, factor_pb / shards * lower)
    receipt = dict(run=shard["directory"], seed=shard["seed"], generated=n,
                   generated_cross_section_nb=shard["generated_cross_section_nb"],
                   hwsim_sha256=sha256(root))
    return columns, receipt


def build(args, destination):
    from xgboost import XGBClassifier

    combination = json.loads(args.loop_combination.read_text())
    mb = json.loads(args.minbias_report.read_text())
    model_sha = sha256(args.model)
    if (combination["status"] != "validated" or combination["energy_tev"] != 40
            or combination["frozen_model_sha256"] != model_sha or mb["frozen_model"]["model_sha256"] != model_sha
            or mb["frozen_model"]["features"] != FEATURES[:-1] or sha256(args.analyzer) != mb["analyzer_sha256"]):
        raise ValueError("inputs differ from the completed 40 TeV analyses")
    threshold = combination["frozen_threshold"]
    if threshold != mb["frozen_model"]["threshold"]:
        raise ValueError("different classifier thresholds in the two analyses")
    classifier = XGBClassifier(n_jobs=1)
    classifier.load_model(args.model)
    classifier.set_params(n_jobs=1)

    manifest = dict(schema_version=1, energy_tev=40, detector_response="ssc", feature_names=FEATURES,
                    frozen_model_sha256=model_sha, frozen_threshold=threshold, samples={})

    def publish(name, columns, provenance, generated, mode, normalization):
        variables, weights, ids = (columns[k] for k in ("variables", "weight_pb", "sourceevent"))
        if (variables.shape != (len(weights), len(FEATURES)) or ids.shape != weights.shape
                or not np.all(np.isfinite(variables)) or not np.all(np.isfinite(weights))
                or np.any(weights < 0) or np.any(ids < 0) or np.any(ids >= generated)):
            raise ValueError(f"invalid event table for {name}")
        path = destination / (name + ".npz")
        np.savez_compressed(path, **columns)
        rows = selected_rates(variables, weights, classifier, threshold)
        manifest["samples"][name] = dict(file=path.name, sha256=sha256(path), rows=len(weights),
            generated_sources=generated, represented_sources=int(len(np.unique(ids))),
            sum_weight_pb=float(np.sum(weights)), response_mode=mode, normalization=normalization,
            frozen_threshold=threshold, inputs=provenance, **rows)
        return manifest["samples"][name]

    arrays, sources, generated = loop_rows(combination)
    loop = publish("bkg_loop_gg_aa", arrays, sources, generated, "genuine",
                   "Physical pb per SSC response hypothesis. Broad Born masses outside [100,150) GeV and "
                   "focused masses inside this interval use their own native rate and all attempted hard events.")
    for key, sel in (("mass_window_rate_pb", "cuts"), ("frozen_classifier_rate_pb", "frozen_classifier")):
        check_rate(loop[key], combination["selections"][sel]["selected_cross_section_pb"], "disjoint loop " + sel)
    loop["input_report_sha256"] = sha256(args.loop_combination)
    loop["shower_quality"] = {key: {"discarded_events": value.get("discarded_events"),
                              "max_momentum_violation_mev": value.get("max_momentum_violation_mev")}
                              for key, value in combination["shower_quality"].items()}

    shards = mb["shards"]
    if (len({shard["seed"] for shard in shards}) != len(shards)
            or len({shard["directory"] for shard in shards}) != len(shards)
            or len({shard["generated"] for shard in shards}) != 1
            or sum(shard["generated"] for shard in shards) != mb["generated_events"]):
        raise ValueError("expected disjoint, equally sized minimum-bias samples")
    def extract(args_for_shard):
        i, shard = args_for_shard
        return jet_shard(i, shard, args.minbias_dir, args.analyzer, len(shards), classifier, threshold)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        partial = list(pool.map(extract, enumerate(shards)))
    arrays = dict(zip(("variables", "sourceevent", "weight_pb", "lower_weight_pb"),
                      (np.concatenate(pieces) for pieces in zip(*(p[0] for p in partial)))))
    jet = publish("bkg_dijet_double_fake", arrays, [p[1] for p in partial], mb["generated_events"],
                  "dijet_double_fake", "Physical pb per eligible two-jet fake-pair hypothesis. The default "
                  "all-quark independent GEM upper endpoint is averaged over 132 equal-size collision samples, "
                  "each with its generated inelastic cross section. The lower alternative uses two eligible "
                  "jets with the gluon fake rate. Each source has total pair probability at most one.")
    jet["lower_endpoint"] = selected_rates(arrays["variables"], arrays["lower_weight_pb"], classifier, threshold)
    jet["input_report_sha256"] = sha256(args.minbias_report)
    jet["analyzer_sha256"] = mb["analyzer_sha256"]
    for endpoint, scenario in ((jet, "all_quark_independent"), (jet["lower_endpoint"], "all_gluon_two_eligible")):
        for key, prefix in (("mass_window_rate_pb", ""), ("frozen_classifier_rate_pb", "frozen_")):
            check_rate(endpoint[key], mb["totals"][prefix + scenario + "@1"]["rate_fb"] / 1000, "dijet " + key)
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for option in ("loop-combination", "minbias-report", "minbias-dir", "analyzer", "model", "output"):
        p.add_argument("--" + option, required=True, type=Path)
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    if args.output.exists() or args.workers < 1:
        p.error("supply a new output directory and a positive worker count")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.output.parent) as tmp:
        release = Path(tmp) / "release"
        release.mkdir()
        build(args, release)
        release.rename(args.output)
