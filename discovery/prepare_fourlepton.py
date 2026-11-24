"""Project the frozen, independent four-lepton score files into selected mass bins."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project(reference):
    old, extra = reference / "xgboost", reference / "extended_test"
    original = json.loads((old / "summary.json").read_text())
    extended = json.loads((extra / "summary.json").read_text())
    plan = json.loads((extra / "plan.json").read_text())
    threshold = extended["threshold"]
    lumi = original["luminosity_fb"]
    sample_rows = extended["extended_test"]["samples"]
    edges = np.arange(120., 131.)
    cells = {row["sample"]: defaultdict(lambda: np.zeros(10)) for row in sample_rows}
    sources = {"original": old / "scores.csv", "additional": extra / "additional_scores.csv"}
    hashes = {key: digest(path) for key, path in sources.items()}
    if hashes["original"] != plan["frozen_files"]["scores.csv"]["sha256"] or original["threshold"]["threshold"] != threshold:
        raise ValueError("Frozen classifier inputs differ")
    for group, path in sources.items():
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                if group == "original" and row["partition"] != "test":
                    continue
                if row["baseline"] != "True" or float(row["score"]) < threshold:
                    continue
                sample = row["sample"]
                weight = float(row["weight"]) / lumi
                if group == "original":
                    weight /= original["split_fractions"][row["source_group"]]["test"]
                if sample in extended["coefficients"]:
                    weight *= extended["coefficients"][sample][group + "_fraction"]
                mass = float(row["m4l"])
                if not edges[0] <= mass <= edges[-1]:
                    raise ValueError(f"Selected event outside mass window for {sample}")
                index = min(int(mass - edges[0]), len(edges) - 2)
                cells[sample][(group, row["source"])][index] += weight
    hist = {}
    for row in sample_rows:
        sample = row["sample"]
        vectors = np.asarray(list(cells[sample].values())).reshape(-1, 10)
        rates = vectors.sum(axis=0)
        cov = vectors.T @ vectors
        if (not np.isclose(rates.sum(), row["sum_weight"] / lumi, atol=1.e-10)
                or not np.isclose(cov.sum(), row["variance"] / lumi**2, atol=1.e-10)):
            raise ValueError(f"Four-lepton test score replay failed for {sample}")
        hist[sample] = {"rate_fb": rates.tolist(), "covariance_fb2": cov.tolist()}
    return {"threshold": threshold, "edges_gev": edges.tolist(), "source_sha256": hashes,
            "extended_summary_sha256": digest(extra / "summary.json"), "samples": hist}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="Directory with xgboost/ and extended_test/")
    json.dump(project(parser.parse_args().reference), sys.stdout, separators=(",", ":"), allow_nan=False)
    print()
