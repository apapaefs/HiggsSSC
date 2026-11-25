"""Replay the native PYTHIA 5.6 electron-fake sample over the GEM tagged mass range."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

SELECTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SELECTION))
from selection_rates import model

WINDOWS = {
    "80_90": (80., 90.), "90_100": (90., 100.), "100_120": (100., 120.),
    "120_140": (120., 140.), "higgs": (123.3, 126.7),
}
FRACTIONS = np.arange(61) / 400.
EDGES = np.arange(80., 141.)


def generate(executable, reference, previous):
    seed, n = reference["seed"], reference["events"]
    rows = []
    sigma_pb = matched = None
    with subprocess.Popen([str(executable), str(seed), str(n), "1"],
                          stdout=subprocess.PIPE, text=True) as proc:
        for line in proc.stdout:
            if line.startswith('{"event":'):
                rows.append(json.loads(line))
            elif line.startswith("RATE "):
                sigma_pb = float(line.split()[-1])
            elif line.startswith("MATCH "):
                matched = list(map(int, line.split()[1:]))
        if proc.wait():
            raise RuntimeError(f"PYTHIA failed for seed {seed}")
    if matched != [reference["matched"], n]:
        raise ValueError(f"Native generator ancestry differs for seed {seed}")
    np.testing.assert_allclose(sigma_pb, reference["sigma_pb"], rtol=2.e-7)

    # Particle-level type probabilities already include 1/33 tracker fake rate.
    # Use GEM's exact 3% and its gluon-fusion multiplier on the qqbar cross section.
    norm = 1000. * sigma_pb / n * 1.20 * (0.03 * 33.)
    mass = np.array([row["mass"] for row in rows])
    variants, with_leak = model(rows)
    masks = [(lo <= mass) & (mass < hi) for lo, hi in WINDOWS.values()]
    if not np.all((80. <= mass) & (mass < 140.)):
        raise ValueError("Unexpected mass outside the tagged range")
    # Validate the generator on the full range *and* validate every previously
    # frozen 125 GeV response variation in the same seed before scanning.
    np.testing.assert_allclose(variants["nominal"].sum(), reference["probability"]
                               * reference["counts"]["80_140"]["pair_pt"][0], rtol=1.e-11)
    for key, per_row in variants.items():
        np.testing.assert_allclose(per_row[masks[-1]].sum() * norm,
                                   previous["rates_fb"][key], rtol=1.e-10)
    if masks[-1].sum() != previous["pre_isolation_assignments"]:
        raise ValueError(f"125 GeV assignment count differs for seed {seed}")

    rate, spectrum = [], []
    bin_number = np.searchsorted(EDGES, mass, side="right") - 1
    for fraction in FRACTIONS:
        weight = norm * with_leak(fraction)
        rate.append([weight[mask].sum() for mask in masks])
        spectrum.append(np.bincount(bin_number, weights=weight, minlength=len(EDGES)-1))
    rate = np.array(rate).T
    spectrum = np.array(spectrum).T
    for i, (lo, hi) in enumerate(list(WINDOWS.values())[:-1]):
        np.testing.assert_allclose(rate[i], spectrum[int(lo-80):int(hi-80)].sum(axis=0))
    print(f"Seed {seed} finished, {len(rows)} assignments", flush=True)
    return {"seed": seed, "generated": n, "sigma_pb": sigma_pb,
            "pre_isolation_assignments": [int(mask.sum()) for mask in masks],
            "rate_fb": rate.tolist(), "spectrum_rate_fb": spectrum.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    refs = [r for r in json.loads((SELECTION.parent / "vintage.json").read_text())["runs"]
            if r["flavours"] == 1]
    previous = {r["seed"]: r for r in json.loads((SELECTION / "native.json").read_text())["runs"]}
    with ThreadPoolExecutor(args.workers) as pool:
        runs = list(pool.map(lambda r: generate(args.executable.resolve(), r, previous[r["seed"]]), refs))
    report = {"window_order": list(WINDOWS), "windows_gev": WINDOWS,
              "spectrum_edges_gev": EDGES.tolist(), "leakage_fractions": FRACTIONS.tolist(),
              "conventions": "Each run contains independent cross-section estimates in fb, with axes [window or 1 GeV bin][leakage fraction]. The first four windows partition [80,140) GeV. The Higgs window lies inside the fourth. The 1 GeV bins partition the same range. Average the equal-exposure seeds. Use paired seed resampling across all columns.",
              "normalization": {"tracker_fake_probability": .03, "gluon_multiplier": 1.2,
                                "flavour": "forced Z to ee", "mass_unit": "GeV"},
              "runs": runs}
    raw = (json.dumps(report, separators=(",", ":"), allow_nan=False) + "\n").encode()
    args.output.write_bytes(gzip.compress(raw, mtime=0))


if __name__ == "__main__":
    main()
