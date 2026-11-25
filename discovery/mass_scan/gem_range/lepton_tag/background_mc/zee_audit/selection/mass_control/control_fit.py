"""Calibrate the electron-fake component at low mass and test its mass dependence."""

import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WINDOWS = {
    "80_90": (80, 90), "90_100": (90, 100), "100_120": (100, 120),
    "120_140": (120, 140), "higgs": (123.3, 126.7),
}


def integrate_anchors(masses, densities, lo, hi, log_linear=False):
    if lo < masses[0] or hi > masses[-1] or hi <= lo:
        raise ValueError(f"Window {lo}, {hi} falls outside the tabulation")
    x = np.r_[lo, masses[(masses > lo) & (masses < hi)], hi]
    if log_linear:
        log_y = np.interp(x, masses, np.log(densities))
        y = np.exp(log_y)
        mean = np.divide(np.diff(y), np.diff(log_y), out=y[:-1].copy(),
                         where=np.diff(log_y) != 0)
        return float(np.sum(np.diff(x) * mean))
    y = np.interp(x, masses, densities)
    return float(np.sum(np.diff(x) * (y[1:] + y[:-1]) / 2))


def references(masses, densities, log_linear=False):
    rates = {key: integrate_anchors(masses, densities, *window, log_linear=log_linear)
             for key, window in WINDOWS.items()}
    rates["training"] = rates["80_90"] + rates["90_100"]
    rates["low_ratio"] = rates["90_100"] / rates["80_90"]
    return rates


def summaries(rates, ref):
    """Rates have axes (full/leave-one-out, mass interval) at one leakage value."""
    low, upper_low, middle, high, higgs = rates.T
    training = low + upper_low
    scale = ref["training"] / training
    return {
        "raw_80_100_fb": training, "raw_125_fb": higgs,
        "scale": scale, "ratio_90_100_to_80_90": upper_low / low,
        "predicted_100_120_fb": scale * middle,
        "predicted_120_140_fb": scale * high, "predicted_125_fb": scale * higgs,
        "to_gem_120_140": scale * high / ref["120_140"],
        "to_gem_125": scale * higgs / ref["higgs"],
    }


def jackknife(values):
    replicates = values[1:]
    n = len(replicates)
    se = np.sqrt((n - 1) / n * np.sum((replicates - replicates.mean()) ** 2))
    return {"value": float(values[0]), "seed_mc_se": float(se)}


def match_low_mass(grid, samples, ref, match_shape=False):
    """Interpolate the rates linearly in leakage to match the low-mass yield or ratio."""
    fits, fitted_rates = [], []
    for sample in samples:
        if match_shape:
            residual = sample[1] - ref["low_ratio"] * sample[0]
        else:
            residual = sample[0] + sample[1] - ref["training"]
            if np.any(np.diff(residual) > 0):
                raise ValueError("Low-mass rate increases with added isolation leakage")
        at_knot = np.flatnonzero(residual == 0)
        crossings = np.flatnonzero(residual[:-1] * residual[1:] < 0)
        if len(at_knot) + len(crossings) != 1:
            raise ValueError("Expected one solution in the leakage grid")
        if len(at_knot):
            i = at_knot[0]
            fits.append(grid[i])
            fitted_rates.append(sample[:, i])
            continue
        i = crossings[0]
        weight = residual[i] / (residual[i] - residual[i + 1])
        fits.append(grid[i] + weight * (grid[i + 1] - grid[i]))
        fitted_rates.append(sample[:, i] + weight * (sample[:, i + 1] - sample[:, i]))
    return np.array(fits), np.array(fitted_rates)


def bounds(values):
    values = list(values)
    return [float(min(values)), float(max(values))]


def fit(native, legacy):
    grid = np.array(native["leakage_fractions"])
    order = [native["window_order"].index(key) for key in WINDOWS]
    if any(tuple(native["windows_gev"][key]) != edges for key, edges in WINDOWS.items()):
        raise ValueError("Native and GEM mass definitions differ")
    runs = native["runs"]
    if (len(runs) < 2 or len({run["seed"] for run in runs}) != len(runs)
            or len({run["generated"] for run in runs}) != 1):
        raise ValueError("The jackknife requires distinct seeds with equal generated counts")
    rates = np.array([run["rate_fb"] for run in runs])[:, order, :]
    if (rates.shape != (len(runs), len(WINDOWS), len(grid)) or grid[0] != 0
            or np.any(np.diff(grid) <= 0) or np.any(~np.isfinite(rates)) or np.any(rates < 0)):
        raise ValueError("Invalid native scan")
    mean = rates.mean(axis=0)
    samples = np.concatenate([mean[None], (rates.sum(axis=0) - rates) / (len(runs) - 1)])
    masses = np.array(legacy["masses_gev"])
    densities = np.array(legacy["selected_density_fb_per_gev"])
    gem = references(masses, densities)
    log_linear_gem = references(masses, densities, log_linear=True)
    corners = [references(masses, densities + shifts)
               for shifts in itertools.product([-.0005, .0005], repeat=len(masses))]
    scan = []
    for i, leakage in enumerate(grid):
        summary = summaries(samples[:, :, i], gem)
        scan.append({"leakage": float(leakage), **{k: jackknife(v) for k, v in summary.items()}})

    f_only, rates_only = match_low_mass(grid, samples, gem)
    pure = summaries(rates_only, gem)
    f_split, rates_split = match_low_mass(grid, samples, gem, match_shape=True)
    split = summaries(rates_split, gem)
    at_zero = summaries(samples[:, :, 0], gem)
    comparison = {k: jackknife(pure[k] - at_zero[k])
                  for k in ["ratio_90_100_to_80_90", "predicted_120_140_fb", "predicted_125_fb"]}
    rounding_pure, rounding_split = [], []
    for ref in corners:
        for match_shape, results in [(False, rounding_pure), (True, rounding_split)]:
            f_corner, rates_corner = match_low_mass(grid, samples[:1], ref, match_shape)
            estimates = summaries(rates_corner, ref)
            results.append({"leakage": f_corner[0], **{k: v[0] for k, v in estimates.items()}})
    rounding_keys = ["leakage", "scale", "predicted_125_fb", "predicted_120_140_fb",
                     "to_gem_125", "to_gem_120_140"]
    return {
        "gem_source": legacy["source"], "native_seed_count": len(runs),
        "events_per_seed": runs[0]["generated"], "windows_gev": WINDOWS,
        "procedure": "Piecewise-linear GEM rates. The 80, 90 and 100 GeV anchors alone determine the 80–100 GeV calibration. High-mass anchors enter only evaluation. Each leave-one-seed-out fit repeats the calibration across all bands and leakage values.",
        "historical_uncertainty": "Historical Monte Carlo errors and correlations are unreported. Rounding bounds vary each published anchor by +/-0.0005 fb/GeV and assume linear interpolation.",
        "gem_fb_except_low_ratio": gem,
        "gem_log_linear_interpolation_check": log_linear_gem,
        "gem_rounding_bounds": {k: bounds(c[k] for c in corners) for k in gem},
        "fixed_models_ratio_to_gem_rounding_multipliers": {
            k: bounds((c["training"] / gem["training"]) * (gem[k] / c[k]) for c in corners)
            for k in ["higgs", "120_140"]},
        "fixed_leakage_scan": scan,
        "unit_scale_fit": {"leakage": jackknife(f_only), **{k: jackknife(v) for k, v in pure.items()}},
        "unit_scale_fit_rounding_bounds": {
            k: bounds(c[k] for c in rounding_pure) for k in rounding_keys},
        "two_low_bins_fit": {"leakage": jackknife(f_split), **{k: jackknife(v) for k, v in split.items()}},
        "two_low_bins_fit_rounding_corner_range": {
            k: bounds(c[k] for c in rounding_split) for k in rounding_keys},
        "unit_scale_minus_zero_leakage": comparison,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, default=HERE / "native_scan.json.gz")
    parser.add_argument("--legacy", type=Path, default=HERE.parents[1] / "legacy.json")
    parser.add_argument("--output", type=Path, default=HERE / "control_fit.json")
    args = parser.parse_args()
    with gzip.open(args.native, "rt") as f:
        native = json.load(f)
    legacy = json.loads(args.legacy.read_text())["gem_1993"]
    result = fit(native, legacy)
    result["input_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [args.native, args.legacy]}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
