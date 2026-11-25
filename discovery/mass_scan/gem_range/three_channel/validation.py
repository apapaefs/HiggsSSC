"""Independent true-mass power and sufficient discovery exposures for the finite grid."""
import argparse
import json
import os
from pathlib import Path
import sys

if __name__ == "__main__":
    for key in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[key] = "1"

import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm, poisson

HERE = Path(__file__).resolve().parent
TAG_MC = HERE.parent / "lepton_tag/background_mc"
sys.path.insert(0, str(TAG_MC))
from refit import discovery_model
from likelihood import ALPHA


def reference_inputs():
    # The exclusive untagged score includes the untagged associated Higgs contribution.
    _, inputs = discovery_model()
    tag = json.loads((TAG_MC / "final/combined-mc.json").read_text())["scenarios"]["nominal"]
    four = json.loads((HERE.parents[2] / "published/results.json").read_text())["fourlepton"]
    rates = {
        name: (np.array([tag["signal_in_true_window_fb"], values["signal_fb"]]),
               np.array([tag["background_in_true_window_fb"], values["background_fb"]]))
        for name, values in four.items()}
    return inputs.fixed.q, rates


def true_mass_power(qgamma, luminosity, signal, background, alpha):
    """Average over both independent window counts, integrating the photon score exactly."""
    mu = np.sqrt(qgamma(luminosity))
    log_weights = np.log1p(signal / background)
    means = luminosity * (signal + background)
    # The neglected probability under either hypothesis is < 2e-17.
    maxima = poisson.ppf(1. - 1.e-13, means).astype(int) + 12
    if np.any(poisson.sf(maxima, means) >= 1.e-17):
        raise ValueError("Poisson grid does not cover the required tail")
    counts = [np.arange(last + 1) for last in maxima]
    weighted_counts = np.add.outer(log_weights[0] * counts[0], log_weights[1] * counts[1])

    def weights(rate):
        return np.multiply.outer(poisson.pmf(counts[0], luminosity * rate[0]),
                                 poisson.pmf(counts[1], luminosity * rate[1]))

    p0 = weights(background)
    def tail(z):
        return np.sum(p0 * ndtr(weighted_counts / mu - z)) - alpha

    threshold = brentq(tail, -12., weighted_counts.max() / mu + 12.)
    return float(np.sum(weights(signal + background) * ndtr(weighted_counts / mu + mu - threshold)))


def calculate():
    qgamma, rates = reference_inputs()
    n_masses = 161
    rows = []
    for scenario, (signal, background) in rates.items():
        for power in (.5, .9):
            row = {"fourlepton_background": scenario, "power": power}
            for name, alpha in [("fixed_mass", ALPHA), ("bonferroni_sufficient", ALPHA / n_masses)]:
                lumi = brentq(lambda l: true_mass_power(qgamma, l, signal, background, alpha) - power,
                              .2, 9., xtol=1.e-8)
                row[name] = {"luminosity_fb_inverse": lumi, "design_months": 1.2 * lumi}
            rows.append(row)
    return dict(mass_hypotheses=n_masses, bonferroni_local_z=float(norm.isf(ALPHA / n_masses)),
                tagged_window_signal_fb=float(rates["central"][0][0]),
                tagged_window_background_fb=float(rates["central"][1][0]),
                fourlepton_signal_fb=float(rates["central"][0][1]),
                fourlepton_background_fb={k: float(v[1][1]) for k, v in rates.items()},
                results=rows,
                interpretation="Bonferroni limits the grid-wide false-positive rate if every local test has valid size. Power here uses the true mass alone. The upper four-lepton background is an integrated control bound and has no measured spectrum for the other masses.")


def direct_power_check(scan, toys=60000, seed=972413):
    """Draw unconditioned signal experiments and check the power from conditional integration."""
    from scan_three import ANCHOR, JointScan, load_model

    model = load_model(scan["tag_background"])
    rows = []
    for index, row in enumerate(scan["crossings"]):
        lumi = row["luminosity_fb_inverse"]
        experiment = JointScan(model, lumi)
        alpha = norm.sf(row["local_z"])
        cuts = experiment.local(alpha)[0]
        rng = np.random.default_rng(seed + index)
        passed = passed_true = 0
        for i in range(0, toys, 512):
            n = min(toys - i, 512)
            z = rng.standard_normal((n, experiment.field.shape[1])) @ experiment.field.T
            scores = experiment.strength * (z + experiment.strength * experiment.response)
            for channel in experiment.channels:
                means = lumi * (channel.signal + channel.background)
                cell_counts = rng.poisson(means, size=(n, len(means)))
                prefix = np.column_stack([np.zeros(n, dtype=np.int64), cell_counts.cumsum(axis=1)])
                scores += channel.weights * (prefix[:, channel.stop] - prefix[:, channel.start])
            discoveries = scores >= cuts
            passed += discoveries.any(axis=1).sum()
            passed_true += discoveries[:, ANCHOR].sum()
        observed = passed / toys
        error = np.sqrt(observed * (1. - observed) / toys)
        true = passed_true / toys
        true_error = np.sqrt(true * (1. - true) / toys)
        if abs(observed - row["power"]) > 6 * np.hypot(error, row["power_mc_se"]):
            raise ValueError("Direct signal toys and fitted scan power disagree")
        if abs(true - row["true_mass_power_analytic"]) > 6 * true_error:
            raise ValueError("Direct true-mass toys and analytic probability disagree")
        rows.append(dict(target_power=row["target_power"], luminosity_fb_inverse=lumi,
                         scan_power=float(observed), scan_power_mc_se=float(error),
                         true_mass_power=float(true), true_mass_power_mc_se=float(true_error)))
    return dict(toys=toys, base_seed=seed, crossings=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "validation.json")
    parser.add_argument("--scan", type=Path, help="Check a completed scan with direct signal toys")
    args = parser.parse_args()
    result = calculate()
    if args.scan:
        result["independent_signal_toys"] = direct_power_check(json.loads(args.scan.read_text()))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
