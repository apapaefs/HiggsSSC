"""Project the frozen four-lepton test scores in 0.5 GeV bins across the 105–140 GeV preselection."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project(reference):
    files = {
        "original_scores": reference / "xgboost/scores.csv",
        "original_summary": reference / "xgboost/summary.json",
        "additional_scores": reference / "extended_test/additional_scores.csv",
        "extended_summary": reference / "extended_test/summary.json",
        "extension_plan": reference / "extended_test/plan.json",
    }
    hashes = {key: digest(path) for key, path in files.items()}
    original = json.loads(files["original_summary"].read_text())
    extended = json.loads(files["extended_summary"].read_text())
    plan = json.loads(files["extension_plan"].read_text())
    fixed_path = Path(__file__).resolve().parents[1] / "data/fourlepton-spectrum.json"
    fixed = json.loads(fixed_path.read_text())
    threshold = fixed["threshold"]
    if (hashes["original_scores"] != fixed["source_sha256"]["original"]
            or hashes["additional_scores"] != fixed["source_sha256"]["additional"]
            or hashes["extended_summary"] != fixed["extended_summary_sha256"]
            or hashes["original_summary"] != plan["frozen_files"]["summary.json"]["sha256"]
            or threshold != extended["threshold"] or threshold != original["threshold"]["threshold"]):
        raise ValueError("The score release differs from the verified fixed-mass inputs")
    if (original["card"]["preselection"] != {"m4l": {"min": 105., "max": 140.}}
            or original["card"]["cuts"] != {"m4l": {"min": 120., "max": 130.}}):
        raise ValueError("The stored baseline flag must represent only the fixed mass window")

    edges = np.linspace(105., 140., 71)
    by_sample = {name: np.zeros(70) for name in fixed["samples"]}
    for group in ("original", "additional"):
        with files[group + "_scores"].open(newline="") as stream:
            for row in csv.DictReader(stream):
                if group == "original" and row["partition"] != "test":
                    continue
                if float(row["score"]) < threshold:
                    continue
                mass = float(row["m4l"])
                if not edges[0] <= mass <= edges[-1] or row.get("eligible", "True") != "True":
                    raise ValueError("Score outside the frozen preselection")
                if (row["baseline"] == "True") != (120 <= mass <= 130):
                    raise ValueError("Unexpected cut other than the 125 GeV mass window")
                sample = row["sample"]
                rate = float(row["weight"]) / original["luminosity_fb"]
                if group == "original":
                    rate /= original["split_fractions"][row["source_group"]]["test"]
                if sample in extended["coefficients"]:
                    rate *= extended["coefficients"][sample][group + "_fraction"]
                index = min(int((mass - edges[0]) * 2), len(edges) - 2)
                by_sample[sample][index] += rate

    groups = {"signal_fb": "signal_", "qq_continuum_fb": "bkg_qq_", "gg_continuum_fb": "bkg_gg_"}
    output = dict(edges_gev=edges.tolist(), threshold=threshold, source_sha256=hashes,
                  fixed_mass_spectrum_sha256=digest(fixed_path))
    for name, bins in by_sample.items():
        central = bins[30:50].reshape(10, 2).sum(axis=1)
        np.testing.assert_allclose(central, fixed["samples"][name]["rate_fb"], rtol=1e-10, atol=1e-12,
                                   err_msg=f"Frozen fixed-mass bins differ for {name}")
    for key, prefix in groups.items():
        output[key] = np.sum([bins for name, bins in by_sample.items() if name.startswith(prefix)], axis=0).tolist()
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="Shared release directory with xgboost/ and extended_test/")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(project(args.reference), indent=2, allow_nan=False) + "\n")
