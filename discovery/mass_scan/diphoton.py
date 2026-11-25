"""Approximate the global discovery reach for one SSC diphoton search.

Translate the reconstructed 125 GeV peak and profile the continuum in a
Gaussian score field. Use the Rice upcrossing approximation for the one-sided
global tail, and check it independently by conditional importance sampling.
Signal-injected Gaussian toys give the probability of discovery anywhere in
the scan, with the local strength matched to the calibrated 125 GeV likelihood.
"""
import argparse
import json
import os
from pathlib import Path
import sys

if __name__ == "__main__":
    for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[variable] = "1"

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, least_squares
from scipy.special import ndtr
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from input_data import BACKGROUNDS, SIGNAL, load_archive, pool
from likelihood import ALPHA, BLIND, EDGES, Calibration, Shape
from plots import BLUE, INK, ORANGE, save
import matplotlib.pyplot as plt

STEP_GEV = .05
TRUE_MASS = 125.
REFERENCE_LUMI = 5.   # Used for the across-mass covariance and the signal response.
NULL_TOYS = 40000
SIGNAL_TOYS = 100000
NULL_SEED = 67125
SIGNAL_SEED = 51625


def load_inputs():
    records, packets = load_archive()
    rates, covariances, _, _ = pool(records, packets)
    signal = rates[SIGNAL]
    observed = sum(rates[key] for key in BACKGROUNDS)
    covariance = sum(covariances[key] for key in BACKGROUNDS)
    shape = Shape("exponential", 2)
    theta, background, derivatives, _ = shape.fit(observed, ~BLIND)
    model = json.loads((HERE.parent / "published/results.json").read_text())["models"]["exponential2"]
    fixed = Calibration(shape, theta, background, signal, model["spurious_mu"], model["sigma_mu"])
    # Profile the three continuum coefficients in the Gaussian Poisson metric.
    continuum_basis = np.linalg.qr(derivatives / np.sqrt(background[:, None]))[0]
    noise_cov = covariance / np.sqrt(np.outer(background, background))
    return signal, background, continuum_basis, noise_cov, fixed


def fit_peak(signal):
    """Fit the integral of a common-mean Gaussian mixture over each actual mass bin."""
    def template(masses, parameters):
        mean, core, tail, tail_fraction = parameters
        offset = EDGES[None, :] - (masses[:, None] - TRUE_MASS) - mean
        cdf = (1 - tail_fraction) * ndtr(offset / core) + tail_fraction * ndtr(offset / tail)
        return signal.sum() * np.diff(cdf, axis=1)

    def residuals(parameters):
        prediction = template(np.array([TRUE_MASS]), parameters)[0]
        return (prediction - signal) / np.sqrt(np.maximum(signal, .2))

    best = least_squares(residuals, [125., .85, 1.7, .10],
                         bounds=([124., .5, 1.1, 0.], [126., 1.2, 5., .7]))
    return best.x, lambda masses: template(masses, best.x)


def score_field(templates, signal, background, continuum_basis, noise_cov, luminosity, true_index):
    """Each row of 'field' is a unit vector acting on the same Gaussian noise draw."""
    weighted = templates / np.sqrt(background)
    projected = weighted - (weighted @ continuum_basis) @ continuum_basis.T
    # A shared fluctuation of the background simulation perturbs each mass's
    # fitted spurious signal. Propagate the source covariance into those scores.
    transform = np.linalg.cholesky(np.eye(len(background)) + luminosity * noise_cov)
    vectors = projected @ transform
    norms = np.linalg.norm(vectors, axis=1)
    response = (projected @ (signal / np.sqrt(background))) / norms
    return vectors / norms[:, None], response / response[true_index]


def critical_value(field):
    # The expected number of Gaussian upcrossings is arc_length/(2π) exp(-z²/2)
    # for a smooth one-dimensional field of normalized templates. Count the
    # possibility of starting above threshold at the left boundary once.
    cosines = np.sum(field[1:] * field[:-1], axis=1)
    length = np.arccos(np.clip(cosines, -1., 1.)).sum()
    def upper_tail(z):
        return norm.sf(z) + (length / (2 * np.pi)) * np.exp(-.5 * z**2)
    return brentq(lambda z: upper_tail(z) - ALPHA, 5., 7.), norm.isf(upper_tail(5.))


def tail_check(field, z, *, seed=NULL_SEED, toys=NULL_TOYS):
    """Unbiased finite-grid trials factor at z using a conditional mixture."""
    rng = np.random.default_rng(seed)
    a = field.astype(np.float32)
    factors = []
    for start in range(0, toys, 500):
        size = min(500, toys - start)
        # Pick a mass uniformly; force this mass's score above threshold, with
        # the correct truncated-normal law. Count all excursions in the field.
        choice = rng.integers(len(a), size=size)
        gaussian = rng.standard_normal((size, a.shape[1]), dtype=np.float32)
        height = norm.isf(norm.sf(z) * rng.uniform(size=size))
        correction = height - np.sum(gaussian * a[choice], axis=1)
        gaussian += (correction[:, None] * a[choice]).astype(np.float32)
        over = np.count_nonzero(gaussian @ a.T > z, axis=1)
        if np.any(over == 0):
            raise ValueError("The conditional toy failed to cross at its sampled mass")
        factors.extend(len(a) / over)
    return float(np.mean(factors)), float(np.std(factors, ddof=1) / np.sqrt(len(factors)))


def signal_strengths(field, response, z):
    """Invert the maximum statistic as the common 125 GeV signal strength grows."""
    rng = np.random.default_rng(SIGNAL_SEED)
    a = field.astype(np.float32)
    positive, negative = response > 0, response < 0
    rising, falling, constant = [], [], []
    for start in range(0, SIGNAL_TOYS, 1000):
        size = min(1000, SIGNAL_TOYS - start)
        noise = rng.standard_normal((size, a.shape[1]), dtype=np.float32) @ a.T
        # A trial with positive signal response crosses above a lower limit on
        # the strength; a negative response crosses below an upper limit.
        rising.extend(np.min((z - noise[:, positive]) / response[positive], axis=1))
        falling.extend(np.max((z - noise[:, negative]) / response[negative],
                              axis=1, initial=-np.inf))
        constant.extend(np.any(noise[:, response == 0] >= z, axis=1))
    rising, falling, constant = np.array(rising), np.array(falling), np.array(constant)
    def power(strength):
        return np.mean((rising <= strength) | (falling >= strength) | constant)
    return [brentq(lambda strength: power(strength) - target, 0., 12.) for target in (.5, .9)]


def plot_expected(out, masses, response, critical_z, luminosity, local_strength):
    # This is the one-sided p-value of the expected (unfluctuated) score, in
    # a detector with a true 125 GeV Higgs. The scan power includes fluctuations.
    p = norm.sf(np.maximum(local_strength * response, 0.))
    fig, ax = plt.subplots(figsize=(7.7, 4.8), layout="constrained")
    ax.semilogy(masses, p, lw=2.2, color=BLUE, label="Expected local $p_0$")
    ax.axhline(ALPHA, color="#7a8592", ls=":", label=r"Local $5\sigma$")
    ax.axhline(norm.sf(critical_z), color=ORANGE, ls="--",
               label=rf"Global $5\sigma$ threshold (local {critical_z:.2f}$\sigma$)")
    true_p = p[np.argmin(abs(masses - TRUE_MASS))]
    ax.scatter([TRUE_MASS], [true_p], color=BLUE, s=27, zorder=3)
    ax.annotate("125 GeV signal", (TRUE_MASS, true_p), (TRUE_MASS + 2.0, 9.e-6),
                color=INK, arrowprops={"arrowstyle": "-", "color": INK, "lw": .8})
    ax.set(title="Diphoton mass scan at the median global discovery exposure",
           xlabel=r"Tested Higgs mass [GeV]", ylabel=r"Local background-only $p_0$",
           xlim=(masses[0], masses[-1]), ylim=(1.e-9, 2.0))
    ax.text(.025, .85, rf"SSC projection, 40 TeV, {luminosity:.2f} fb$^{{-1}}$",
            transform=ax.transAxes, va="top", color=INK)
    ax.legend(loc="center left", bbox_to_anchor=(.018, .45), fontsize=9, frameon=False)
    save(fig, out, "diphoton-local-p")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE)
    out = parser.parse_args().output
    out.mkdir(parents=True, exist_ok=True)
    signal, b, basis, mc, fixed = load_inputs()
    parameters, smooth = fit_peak(signal)
    empirical = PchipInterpolator(EDGES, np.r_[0., np.cumsum(signal)])
    raw_score = signal / np.sqrt(b)
    score = raw_score - basis @ (basis.T @ raw_score)
    linear_width = np.sqrt(score @ mc @ score) / (score @ score)
    def exposure(strength):
        return brentq(lambda lumi: fixed.q(lumi) - strength**2, .01, 25.)
    if abs(linear_width - fixed.width) > .001 or abs(exposure(5.) - 3.610034) > 1.e-5:
        raise ValueError("The local calibration differs from the published diphoton fit")
    result = dict(generator_mass_gev=TRUE_MASS, template_mean_core_tail_fraction=parameters.tolist(),
                  template_L1_fraction=float(np.abs(smooth(np.array([TRUE_MASS]))[0] - signal).sum() / signal.sum()),
                  linear_mc_spurious_sigma=linear_width, calibrated_mc_spurious_sigma=fixed.width,
                  fixed_median_fb_inverse=exposure(5.), fixed_90pct_fb_inverse=exposure(5 + norm.ppf(.9)),
                  null_importance_toys=NULL_TOYS, signal_toys=SIGNAL_TOYS,
                  null_seed=NULL_SEED, signal_seed=SIGNAL_SEED, grid_step_gev=STEP_GEV,
                  covariance_reference_fb_inverse=REFERENCE_LUMI, scans=[])
    for low, high in [(115, 130), (110, 135)]:
        masses = np.linspace(low, high, round((high - low) / STEP_GEV) + 1)
        index = np.argmin(abs(masses - TRUE_MASS))
        if masses[index] != TRUE_MASS:
            raise ValueError("The true mass must be on the scan grid")
        def field(templates, luminosity):
            return score_field(templates, signal, b, basis, mc, luminosity, index)
        templates = smooth(masses)
        a, response = field(templates, REFERENCE_LUMI)
        z, global_z_at_5 = critical_value(a)
        trials, error = tail_check(a, z)
        median, high_power = [exposure(x) for x in signal_strengths(a, response, z)]
        translated_edges = np.clip(EDGES[None, :] - masses[:, None] + TRUE_MASS, EDGES[0], EDGES[-1])
        alternative = np.diff(empirical(translated_edges), axis=1)
        z_alt = critical_value(field(alternative, REFERENCE_LUMI)[0])[0]
        expected_trials = ALPHA / norm.sf(z)
        if abs(trials - expected_trials) > 5 * error:
            raise ValueError("Upcrossing and importance-sampling tails differ")
        result["scans"].append(dict(
            mass_range_gev=[low, high], local_z_for_global_five=z, trials_factor=expected_trials,
            importance_trials_factor=trials, importance_mc_sigma=error,
            global_z_at_local_five=global_z_at_5, cdf_template_local_z=z_alt,
            local_z_without_mc=critical_value(field(templates, 0.)[0])[0],
            local_z_at_90pct_lumi=critical_value(field(templates, high_power)[0])[0],
            median_fb_inverse=median, power90_fb_inverse=high_power,
            true_mass_threshold_only_median_fb_inverse=exposure(z)))
        if (low, high) == (110, 135):
            plot_expected(out, masses, response, z, median, np.sqrt(fixed.q(median)))
    (out / "diphoton-scan.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
