"""Expected inclusive diphoton discovery on the finite 80–160 GeV GEM-range grid."""
import json
import os
from pathlib import Path
import sys

if __name__ == "__main__":
    for name in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = "1"

import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from diphoton import fit_peak, plot_expected, signal_strengths, tail_check
from grid_bounds import calibration_at_125
from input_data import BACKGROUNDS, SIGNAL
from likelihood import ALPHA
from plots import COLORS, INK, NAMES, ORANGE, save, style
import matplotlib.pyplot as plt


MASSES = np.arange(80., 160.01, .5)
RESOLUTION_MASS = np.arange(80., 161., 10.)
# June GEM Higgs study, Table 6. Only the relative change in width is used.
GEM_RESOLUTION = np.array([.66, .69, .77, .79, .84, .86, .91, .92, 1.])
Z_LO, Z_HI = 91.1876 - 10., 91.1876 + 10.
NULL_SEED, CHECK_SEED = 806125, 160625


def templates(edges, parameters, masses, constant_width=False):
    mean, core, tail, fraction = parameters
    widths = np.interp(masses, RESOLUTION_MASS, GEM_RESOLUTION)
    scale = np.ones_like(widths) if constant_width else widths / np.interp(125., RESOLUTION_MASS, GEM_RESOLUTION)
    centers = masses[:, None] + mean - 125.

    def cdf(points):
        u = (points - centers) / scale[:, None]
        return (1 - fraction) * ndtr(u / core) + fraction * ndtr(u / tail)

    area = np.diff(cdf(edges), axis=1)
    inside_z = np.diff(cdf(np.clip(edges, Z_LO, Z_HI)), axis=1)
    return area - (1 - (.86 / .96)**2) * inside_z


def score_weights(edges, counts, peaks, masses):
    """Profile smooth continuum variations locally, including the detector's Z band."""
    centers = .5 * (edges[1:] + edges[:-1])
    background = sum(counts[k] for k in BACKGROUNDS)
    fraction_z = np.diff(np.clip(edges, Z_LO, Z_HI)) / np.diff(edges)
    # A common variation in genuine-photon efficiency affects prompt and box
    # twice, and the one-photon + jet process once. Other central shapes
    # (including Drell-Yan and fakes) are already in the simulated spectrum.
    veto_shape = fraction_z * (2 * (counts[BACKGROUNDS[0]] + counts[BACKGROUNDS[1]])
                               + counts["bkg_gamma_j"])
    weights = np.zeros_like(peaks)
    for i, mass in enumerate(masses):
        window = (centers >= mass - 20.) & (centers < mass + 35.)
        b = background[window]
        x = (centers[window] - mass) / 25.
        design = np.column_stack((b, b*x, b*x*x, veto_shape[window]))
        whitened = design / np.sqrt(b[:, None])
        target = peaks[i, window] / np.sqrt(b)
        # lstsq handles the identically zero Z-band column in a high-mass window.
        beta = np.linalg.lstsq(whitened, target, rcond=1e-11)[0]
        weights[i, window] = (target - whitened @ beta) / np.sqrt(b)
    return weights


def score_field(weights, background, covariance, signal, luminosity, true_index):
    # Scores act on data minus the central MC expectation. Data Poisson noise
    # is Gaussian at these counts. Independent MC noise uses the source-event
    # covariance of the full simulated background and is shared across masses.
    transform = np.linalg.cholesky(np.diag(background) + luminosity * covariance)
    vectors = weights @ transform
    widths = np.linalg.norm(vectors, axis=1)
    response = weights @ signal / widths
    return vectors / widths[:, None], response / response[true_index], float(response[true_index])


def calibrate(field, seed, toys=40000):
    z = 5.87
    for _ in range(4):
        factor, error = tail_check(field, z, seed=seed, toys=toys)
        new_z = norm.isf(ALPHA / factor)
        if abs(new_z - z) < .0001:
            return float(new_z)
        z = new_z
    return float(z)


def plot_spectrum(edges, counts, covariance):
    fig, (top, low) = plt.subplots(2, 1, sharex=True, figsize=(8.7, 5.8),
                                   height_ratios=(2.3, 1), layout="constrained")
    for axis in (top, low):
        style(axis)
        axis.axvspan(edges[0], 80, color="#ececec", zorder=0)
        axis.axvspan(160, edges[-1], color="#ececec", zorder=0)
    bottom = np.zeros(len(edges) - 1)
    for sample, color, label in zip(BACKGROUNDS, COLORS, NAMES):
        top.stairs(bottom + counts[sample], edges, baseline=bottom, color=color,
                   fill=True, lw=0, label=label)
        bottom += counts[sample]
    top.errorbar((edges[:-1] + edges[1:])/2, bottom, np.sqrt(covariance.diagonal()),
                 fmt=".", ms=2, color=INK, elinewidth=.7, label="Total with MC error")
    top.set_ylabel(r"Background [fb / GeV]")
    top.legend(loc="upper right", frameon=False, fontsize=8.5, ncol=2)
    low.stairs(counts[SIGNAL], edges, color=ORANGE, lw=1.8)
    low.set(xlabel=r"Reconstructed $m_{\gamma\gamma}$ [GeV]", ylabel="125 GeV Higgs\n[fb / GeV]",
            xlim=(edges[0], edges[-1]), ylim=(0, max(counts[SIGNAL])*1.15))
    fig.suptitle("Inclusive diphotons at the 40 TeV SSC", x=.02, ha="left",
                 color=INK, weight="bold", fontsize=17)
    save(fig, HERE, "mass-spectrum")


def run():
    with np.load(HERE / "projection/pooled.npz", allow_pickle=False) as data:
        edges = data["edges_gev"]
        counts = dict(zip(data["sample"], data["rate_fb"]))
        covariance = sum(c for sample, c in zip(data["sample"], data["covariance_fb2"])
                         if sample in BACKGROUNDS)
    background = sum(counts[k] for k in BACKGROUNDS)
    plot_spectrum(edges, counts, covariance)
    signal = counts[SIGNAL]
    centers = (edges[1:] + edges[:-1]) / 2
    peak_params, _ = fit_peak(signal[(centers >= 105) & (centers < 160)])
    def make_weights(masses, constant_width=False):
        return score_weights(edges, counts, templates(edges, peak_params, masses, constant_width), masses)

    weights = make_weights(MASSES)
    _, fixed = calibration_at_125()

    def exposure(strength):
        return brentq(lambda l: fixed.q(l) - strength**2, .01, 30.)

    def discovery(power, initial, selected_masses, w, make_plot=False):
        true_index = np.flatnonzero(selected_masses == 125.).item()
        luminosity = initial
        for _ in range(5):
            field, response, abs_response = score_field(w, background, covariance, signal, luminosity, true_index)
            z = calibrate(field, NULL_SEED)
            strength = signal_strengths(field, response, z)[(.5, .9).index(power)]
            new_luminosity = exposure(strength)
            if abs(new_luminosity - luminosity) < .005:
                break
            luminosity = new_luminosity
        else:
            raise ValueError("Exposure and global score calibration failed to converge")
        check, error = tail_check(field, z, seed=CHECK_SEED, toys=80000)
        trials = ALPHA / norm.sf(z)
        if abs(check - trials) > 5 * error or new_luminosity > exposure(norm.isf(ALPHA / len(w)) + norm.ppf(power)):
            raise ValueError("Independent global tail or sufficient exposure check failed")
        variance = (w**2) @ background + luminosity * np.sum((w @ covariance) * w, axis=1)
        skew = (w**3) @ background / (np.sqrt(luminosity) * variance**1.5)
        most_skewed = np.argmax(abs(skew))
        if make_plot:
            plot_expected(HERE, selected_masses, response, z, new_luminosity, strength)
        return dict(power=power, luminosity_fb_inverse=new_luminosity,
            months=1.2*new_luminosity, local_z_for_global_five=z, effective_trials=trials,
            check_trials=check, check_trials_standard_error=error, field_luminosity_fb_inverse=luminosity,
            raw_true_mass_score=float(np.sqrt(luminosity)*abs_response),
            calibrated_true_mass_score=float(np.sqrt(fixed.q(luminosity))),
            largest_abs_poisson_skew=float(abs(skew[most_skewed])),
            mass_of_largest_abs_skew_gev=float(selected_masses[most_skewed]),
            max_first_order_local_poisson_threshold_shift=float(max(abs(skew))*(z*z - 1)/6))

    records = [discovery(.5, 5., MASSES, weights, make_plot=True), discovery(.9, 9., MASSES, weights)]
    inner = (MASSES >= 110) & (MASSES <= 135)
    narrow = discovery(.5, 4.8, MASSES[inner], weights[inner])
    fine_masses = np.arange(80., 160.01, .25)
    fine = discovery(.5, 5.3, fine_masses, make_weights(fine_masses))
    flat = make_weights(MASSES, constant_width=True)
    field = score_field(flat, background, covariance, signal, records[0]["luminosity_fb_inverse"],
                        np.flatnonzero(MASSES == 125.).item())[0]
    result = dict(channel="inclusive diphoton", true_mass_gev=125., scan_gev=[80,160],
        step_gev=.5, mass_hypotheses=len(MASSES), design_fb_inverse_per_year=10.,
        model="Gaussian local scores centred on the pooled background, normalized to the calibrated 125 GeV likelihood",
        relative_resolution_source="GEM-TN-93-373 Table 6",
        mean_core_tail_fraction_at_125=peak_params.tolist(), calibration_seed=NULL_SEED,
        calibration_experiments=40000, independent_null_seed=CHECK_SEED,
        independent_null_experiments=80000, signal_seed=51625, signal_experiments=100000,
        constant_width_median_local_z_for_global_five=calibrate(field, NULL_SEED),
        fixed_median_fb_inverse=exposure(5.), discoveries=records,
        same_model_110_135_median=narrow, quarter_gev_grid_80_160_median=fine)
    (HERE / "diphoton-scan.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    run()
