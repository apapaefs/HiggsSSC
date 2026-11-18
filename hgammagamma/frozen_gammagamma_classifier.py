#!/usr/bin/env python3
"""Package and replay the saved diphoton classifier without training or tuning."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ho_xgboost_analysis import FEATURE_NAMES


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_metadata(metadata):
    if (metadata.get("mode") != "xgboost" or metadata.get("status") != "complete"
            or metadata.get("detector_response") != "ssc"
            or metadata.get("feature_names") != FEATURE_NAMES
            or metadata.get("evaluation_partition") != "test"):
        raise ValueError("expected the completed SSC diphoton classifier and exact saved feature order")
    threshold = metadata.get("best_threshold")
    if threshold is None or not np.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("missing or invalid frozen score threshold")


def build_bundle(report_dir, input_dir, output_dir):
    from hgammagamma.plot_ho_mass_diagnostics import load_events
    report_dir, output_dir = Path(report_dir).resolve(), Path(output_dir).resolve()
    if output_dir.exists():
        raise ValueError("classifier bundle exists; choose a new output directory")
    metadata, data, sources = load_events(report_dir, Path(input_dir).resolve() if input_dir else None)
    validate_metadata(metadata)
    summary = json.loads((report_dir / "summary.json").read_text())
    for key in ("feature_names", "best_threshold", "luminosity_fb", "detector_response", "run_tag"):
        if summary["metadata"].get(key) != metadata.get(key):
            raise ValueError("frozen model metadata and summary differ")
    output_dir.mkdir(parents=True)
    for name in ("signal_background_xgboost.json", "model_metadata.json", "summary.json"):
        shutil.copyfile(report_dir / name, output_dir / name)
    np.savez_compressed(output_dir / "score_reference.npz",
        features=np.column_stack([data[name] for name in FEATURE_NAMES]),
        score=data["score"], sample=data["sample"], partition=data["partition"],
        sourceevent=data["sourceevent"], tree_entry=data["tree_entry"],
        projected_physical_weight=data["projected_physical_weight"])
    names = ("signal_background_xgboost.json", "model_metadata.json", "summary.json", "score_reference.npz")
    proof = {"schema_version": 1, "source_report": str(report_dir), "source_provenance": sources,
        "feature_names": FEATURE_NAMES, "candidate_selection": "n_selected_photons >= 2",
        "score_selection": "score >= best_threshold; no additional mass window",
        "reference_rows": len(data["score"]), "files": {name: sha256(output_dir / name) for name in names}}
    (output_dir / "bundle.json").write_text(json.dumps(proof, indent=2, sort_keys=True) + "\n")
    return proof


def load_frozen_classifier(directory):
    import xgboost
    directory = Path(directory).resolve()
    proof = json.loads((directory / "bundle.json").read_text())
    names = ("signal_background_xgboost.json", "model_metadata.json", "summary.json", "score_reference.npz")
    if (proof.get("schema_version") != 1 or proof.get("feature_names") != FEATURE_NAMES
            or set(proof.get("files", {})) != set(names)):
        raise ValueError("invalid frozen diphoton classifier bundle")
    for name in names:
        if sha256(directory / name) != proof["files"][name]:
            raise ValueError(f"frozen classifier artifact changed: {name}")
    metadata = json.loads((directory / "model_metadata.json").read_text())
    validate_metadata(metadata)
    summary = json.loads((directory / "summary.json").read_text())
    model = xgboost.XGBClassifier(n_jobs=1)
    model.load_model(directory / "signal_background_xgboost.json")
    model.set_params(n_jobs=1)
    if model.get_booster().num_features() != len(FEATURE_NAMES):
        raise ValueError("frozen model feature count differs from its metadata")
    with np.load(directory / "score_reference.npz", allow_pickle=False) as reference:
        features, saved = reference["features"], reference["score"]
        if (features.shape != (proof["reference_rows"], len(FEATURE_NAMES))
                or saved.shape != (proof["reference_rows"],) or not np.all(np.isfinite(features))):
            raise ValueError("invalid saved-score reference dimensions or features")
        replay = model.predict_proba(features)[:, 1]
        threshold = metadata["best_threshold"]
        differences = np.abs(replay - saved)
        decisions = int(np.sum((replay >= threshold) != (saved >= threshold)))
        maximum = float(np.max(differences))
        if not np.all(np.isfinite(replay)) or maximum > 5e-7 or decisions:
            raise ValueError("frozen model does not reproduce saved scores and threshold decisions")
        selected = (reference["partition"] == "test") & (saved >= threshold)
        for row in summary["samples"]:
            weights = reference["projected_physical_weight"][selected & (reference["sample"] == row["sample"])]
            if not np.isclose(float(np.sum(weights)), float(row["expected_events"]), rtol=1e-9, atol=1e-10):
                raise ValueError("saved-score test yields differ from the frozen summary")
    diagnostic = {"bundle": str(directory), "bundle_sha256": sha256(directory / "bundle.json"),
        "files": proof["files"], "reference_rows": proof["reference_rows"],
        "maximum_absolute_score_difference": maximum, "threshold_decision_differences": decisions,
        "saved_test_yields_reproduced": True, "xgboost_version": xgboost.__version__}
    return model, metadata, summary, diagnostic


def score_rows(model, rows):
    eligible = np.asarray([row["n_selected_photons"] >= 2 for row in rows], dtype=bool)
    features = np.asarray([[row[name] for name in FEATURE_NAMES] for row in rows], dtype=float)
    if np.any(~np.isfinite(features[eligible])) or np.any(features[eligible] <= -900):
        raise ValueError("invalid eligible diphoton classifier features")
    scores = np.full(len(rows), np.nan)
    if np.any(eligible):
        scores[eligible] = model.predict_proba(features[eligible])[:, 1]
    if np.any(~np.isfinite(scores[eligible])):
        raise ValueError("nonfinite frozen classifier score")
    return eligible, scores


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", required=True, type=Path)
    parser.add_argument("--input-dir", type=Path, help="original SAMPLE.root copies for reports without embedded features")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(build_bundle(args.report_dir, args.input_dir, args.output_dir), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as error:
        raise SystemExit(str(error))
