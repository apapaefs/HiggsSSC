"""Calibrate the combined fixed-grid SSC mass search and measure its power at 125 GeV.

The shared diphoton score acts on correlated Gaussian bins. Four-lepton events
are Poisson draws in disjoint mass cells used by every overlapping window.
"""
import argparse
import json
import os
from pathlib import Path

if __name__ == "__main__":
    for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[variable] = "1"

import diphoton
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm, poisson

from likelihood import ALPHA, joint_power

HERE = Path(__file__).resolve().parent
BATCH = 512


def physical_cells(raw, edges):
    """Remove negative signed-MC cells and conserve each sideband and the reference window."""
    positive = np.maximum(raw, 0.)
    for low, high in ((edges[0], 120.), (120., 130.), (130., edges[-1])):
        mask = (edges[:-1] >= low) & (edges[1:] <= high)
        total = raw[mask].sum()
        if total < 0 or positive[mask].sum() <= 0:
            raise ValueError("Mass region has no physical integrated rate")
        positive[mask] *= total / positive[mask].sum()
    return positive


class Experiment:
    def __init__(self, masses, photon_inputs, peak, four, exposure):
        self.lumi = exposure
        self.masses = masses
        anchor = self.anchor = int(np.flatnonzero(masses == 125.)[0])
        signal, continuum, basis, covariance, fixed = photon_inputs
        self.strength = np.sqrt(fixed.q(exposure))
        self.field, self.response = diphoton.score_field(peak(masses), signal, continuum, basis,
                                                        covariance, exposure, anchor)
        edges = np.asarray(four["edges_gev"])
        self.start = np.searchsorted(edges, masses - 5)
        self.stop = np.searchsorted(edges, masses + 5)
        np.testing.assert_allclose(edges[self.start], masses - 5)
        np.testing.assert_allclose(edges[self.stop], masses + 5)
        signal_cells = np.asarray(four["signal_fb"])
        background_cells = np.asarray(four["qq_continuum_fb"]) + np.asarray(four["gg_continuum_fb"])
        self.sig = physical_cells(signal_cells, edges)
        self.bkg = physical_cells(background_cells, edges)
        self.clipped_signal_fb = float(-signal_cells[signal_cells < 0].sum())
        self.clipped_background_fb = float(-background_cells[background_cells < 0].sum())
        self.b_window = np.array([self.bkg[i:j].sum() for i, j in zip(self.start, self.stop)])
        if np.any(self.b_window <= 0):
            raise ValueError("Each test needs a positive measured background in its four-lepton window")
        self.true_signal = self.sig[self.start[anchor]:self.stop[anchor]].sum()
        # Use the same positive reference signal yield at each tested mass.
        self.log_weight = np.log1p(self.true_signal / self.b_window)
        self.n = np.arange(1 + int(poisson.isf(1e-16, self.b_window.max() * exposure)))
        self.count_pmfs = poisson.pmf(self.n, self.b_window[:, None] * exposure)
        self.count_scores = self.log_weight[:, None] * self.n

    def thresholds(self, alpha):
        """A common one-sided local p-value with a different numerical score cut per mass."""
        cuts = []
        for offsets, pmf in zip(self.count_scores, self.count_pmfs):
            def null_tail(cut):
                return pmf @ ndtr((offsets - cut) / self.strength)
            cuts.append(brentq(lambda cut: null_tail(cut) - alpha, -10 * self.strength,
                               offsets[-1] + 12 * self.strength))
        return np.asarray(cuts)

    def window_counts(self, cells):
        sums = np.pad(cells.cumsum(axis=1), ((0, 0), (1, 0)))
        return sums[:, self.stop] - sums[:, self.start]

    def null_trials(self, alpha, n_toys, seed):
        """Conditional importance sampling of the union of all local exceedances.

        Choose j uniformly and draw an experiment conditional on passing cut j.
        The likelihood-ratio weight is M*alpha/K for K actual exceedances.
        """
        rng = np.random.default_rng(seed)
        cuts = self.thresholds(alpha)
        x_min = (cuts[:, None] - self.count_scores) / self.strength
        tails = ndtr(-x_min)
        conditional_counts = np.cumsum(self.count_pmfs * tails / alpha, axis=1)
        conditional_counts[:, -1] = 1.
        n_mass = len(cuts)
        trial_factors = []
        for first in range(0, n_toys, BATCH):
            n = min(BATCH, n_toys - first)
            j = rng.integers(n_mass, size=n)
            rows = np.arange(n)
            k = (rng.random((n, 1)) > conditional_counts[j]).sum(axis=1)
            # Truncated unit normal at j given N_j = k. All orthogonal noise is untouched.
            target = norm.isf(tails[j, k] * (1. - rng.random(n)))
            gaussian = rng.standard_normal((n, self.field.shape[1]))
            direction = self.field[j]
            gaussian += (target - np.sum(gaussian * direction, axis=1))[:, None] * direction
            photon = gaussian @ self.field.T
            # Start with independent Poisson cells. Overwrite just the chosen window
            # with k events distributed in proportion to their null rates therein.
            cells = rng.poisson(self.lumi * self.bkg, size=(n, len(self.bkg)))
            for mass_index in np.unique(j):
                group = rows[j == mass_index]
                lo, hi = self.start[mass_index], self.stop[mass_index]
                cells[group, lo:hi] = 0
                event_rows = np.repeat(group, k[group])
                probabilities = np.cumsum(self.bkg[lo:hi]) / self.b_window[mass_index]
                probabilities[-1] = 1.
                bins = lo + np.searchsorted(probabilities, rng.random(len(event_rows)))
                np.add.at(cells, (event_rows, bins), 1)
            values = self.strength * photon + self.log_weight * self.window_counts(cells)
            if np.any(values[rows, j] < cuts[j] - 1e-9):
                raise ValueError("Conditional simulation missed the selected local threshold")
            exceeds = values >= cuts
            exceeds[rows, j] = True
            trial_factors.extend(n_mass / exceeds.sum(axis=1))
        return float(np.mean(trial_factors)), float(np.std(trial_factors, ddof=1) / np.sqrt(n_toys))

    def ordinary_null(self, alpha, n_toys, seed):
        """Independent direct background ensembles for a moderate-tail IS validation."""
        rng = np.random.default_rng(seed)
        cuts = self.thresholds(alpha)
        passed = 0
        for first in range(0, n_toys, BATCH):
            n = min(BATCH, n_toys - first)
            noise = rng.standard_normal((n, self.field.shape[1])) @ self.field.T
            cells = rng.poisson(self.lumi * self.bkg, size=(n, len(self.bkg)))
            values = self.strength * noise + self.log_weight * self.window_counts(cells)
            passed += np.any(values >= cuts, axis=1).sum()
        p = passed / n_toys
        return p, np.sqrt(p * (1. - p) / n_toys)

    def signal_power(self, alpha, n_toys, seed):
        """Average over signal-injected counts and Gaussian noise orthogonal to the 125 score.

        Integrate the remaining one-dimensional Gaussian analytically. A positive
        or negative correlation with 125 gives an upper or lower tail, respectively.
        """
        rng = np.random.default_rng(seed)
        cuts = self.thresholds(alpha)
        corr = self.field @ self.field[self.anchor]
        positive = corr > 1e-12
        negative = corr < -1e-12
        residual_field = self.field - corr[:, None] * self.field[self.anchor]
        drift = self.strength * self.response
        scan_power, local_power = [], []
        for first in range(0, n_toys, BATCH):
            n = min(BATCH, n_toys - first)
            cells = rng.poisson(self.lumi * (self.bkg + self.sig), size=(n, len(self.bkg)))
            counts = self.window_counts(cells)
            noise = rng.standard_normal((n, self.field.shape[1])) @ residual_field.T
            excess = (cuts - self.log_weight * counts) / self.strength - drift - noise
            right = np.min(excess[:, positive] / corr[positive], axis=1)
            left = np.max(excess[:, negative] / corr[negative], axis=1) if negative.any() else -np.inf
            power = np.minimum(1., ndtr(-right) + ndtr(left))
            neutral = ~(positive | negative)
            if neutral.any():
                power[np.any(excess[:, neutral] <= 0., axis=1)] = 1.
            scan_power.extend(power)
            local_power.extend(ndtr(-excess[:, self.anchor]))
        def summary(values):
            return float(np.mean(values)), float(np.std(values, ddof=1) / np.sqrt(n_toys))
        return (*summary(scan_power), *summary(local_power))


def evaluate(masses, inputs, peak, four, lumi, toys, seed):
    exp = Experiment(masses, inputs, peak, four, lumi)
    # The trials factor changes slowly over this narrow range of local sizes.
    alpha = ALPHA / (.6 * len(masses))
    for _ in range(2):
        factor, _ = exp.null_trials(alpha, max(3000, toys // 3), seed)
        alpha = ALPHA / factor
    factor, error = exp.null_trials(alpha, toys, seed + 1)
    alpha = ALPHA / factor
    validation, validation_se = exp.null_trials(alpha, toys, seed + 2)
    p, dp, p125, dp125 = exp.signal_power(alpha, toys, seed + 3)
    exact_local = joint_power(inputs[-1].q, lumi, exp.true_signal, exp.b_window[exp.anchor], alpha=alpha)
    if abs(p125 - exact_local) > 6 * dp125 + 1e-7:
        raise ValueError("True-mass signal power differs from the published local-test construction")
    if abs(validation - factor) > 6 * np.hypot(validation_se, error) + .02:
        raise ValueError("Independent global-tail check failed")
    return dict(luminosity_fb_inverse=lumi, local_alpha=alpha, local_z=norm.isf(alpha),
                local_trials_factor=factor, trials_factor_mc_se=error, null_validation_p=alpha*validation,
                null_validation_mc_se=alpha*validation_se, power=p, power_mc_se=dp,
                true_mass_power=p125, true_mass_power_mc_se=dp125, true_mass_power_analytic=exact_local)


def plot_local_probabilities(exp, local_alpha, out):
    """Pointwise median local probabilities under injection of the actual 125 GeV signal."""
    import matplotlib.pyplot as plt
    from plots import BLUE, INK, ORANGE, save

    b = exp.lumi * exp.b_window
    s = exp.lumi * np.array([exp.sig[i:j].sum() for i, j in zip(exp.start, exp.stop)])
    photon_mean = exp.strength * exp.response
    counts = np.arange(1 + int(poisson.isf(1e-14, (s+b).max())))
    combined = []
    for i, weight in enumerate(exp.log_weight):
        alt_pmf = poisson.pmf(counts, (b+s)[i])
        alt_score = weight * counts + exp.strength * photon_mean[i]
        def alt_tail(score):
            return alt_pmf @ ndtr((alt_score - score) / exp.strength)
        median = brentq(lambda t: alt_tail(t) - .5,
                        alt_score.min() - 10 * exp.strength, alt_score.max() + 10 * exp.strength)
        combined.append(poisson.pmf(counts, b[i]) @ ndtr((weight * counts - median) / exp.strength))

    fig, ax = plt.subplots(figsize=(7.8, 4.9), layout="constrained")
    masses = exp.masses
    ax.semilogy(masses, combined, color=INK, lw=2.3, label=r"$\gamma\gamma+4\ell$")
    ax.plot(masses, norm.sf(photon_mean), color=BLUE, lw=1.9, ls="--", label=r"$\gamma\gamma$")
    ax.plot(masses, poisson.sf(poisson.ppf(.5, s+b)-1, b), color=ORANGE, lw=1.9,
            ls="-.", label=r"$4\ell$ (count)")
    for level, label, color in [
            (ALPHA, r"Local $5\sigma$", "#65707d"),
            (local_alpha, r"Combined global $5\sigma$", "#72334e")]:
        ax.axhline(level, color=color, lw=1.05, ls=":")
        ax.annotate(label, xy=(.015, level), xycoords=("axes fraction", "data"),
                    xytext=(0, 4), textcoords="offset points", color=color, fontsize=9.5)
    ax.set(title=rf"Single SSC experiment  ·  40 TeV  ·  {exp.lumi:.2f} fb$^{{-1}}$",
           xlabel=r"Tested Higgs mass [GeV]", ylabel=r"Median local $p$-value (true $m_H=125$ GeV)",
           xlim=(masses[0], masses[-1]), ylim=(1.5e-9, 1.4))
    ax.legend(loc="center left", frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    save(fig, out, "combined-local-p")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--toys", type=int, default=40000)
    parser.add_argument("--seed", type=int, default=265125)
    parser.add_argument("--range", nargs=2, type=float, default=(110., 135.))
    parser.add_argument("--output", type=Path, default=HERE / "combined-scan.json")
    args = parser.parse_args()
    lo, hi = args.range
    if not (110 <= lo <= 125 <= hi <= 135 and lo < hi
            and float(2*lo).is_integer() and float(2*hi).is_integer()):
        parser.error("--range requires half-GeV endpoints in [110,135], including 125 GeV")
    if args.toys < 3000 or args.seed < 0:
        parser.error("Use at least 3,000 toys and a nonnegative seed")
    inputs = diphoton.load_inputs()
    _, peak = diphoton.fit_peak(inputs[0])
    four = json.loads((HERE / "fourlepton.json").read_text())
    masses = np.arange(args.range[0], args.range[1] + .01, .5)
    published = json.loads((HERE.parent / "published/results.json").read_text())
    reference = published["fourlepton"]["central"]
    local = Experiment(np.array([125.]), inputs, peak, four, 1.83)
    np.testing.assert_allclose([local.true_signal, local.b_window[0]],
                               [reference["signal_fb"], reference["background_fb"]], atol=5e-9)
    results = dict(mass_range_gev=list(args.range), spacing_gev=.5, true_mass_gev=125., seed=args.seed,
                   toys=args.toys, central_fourlepton=True, removed_negative_cell_mass_fb={
                       "signal": local.clipped_signal_fb, "background": local.clipped_background_fb}, crossings=[])
    for power, endpoints in [(.5, [1.95, 2.45]), (.9, [3.35, 4.15])]:
        def test(lumi, n, seed):
            row = evaluate(masses, inputs, peak, four, lumi, n, seed)
            print(f"L={lumi:.5f}, P={row['power']:.5f} +/- {row['power_mc_se']:.5f}, "
                  f"trials factor={row['local_trials_factor']:.2f}", flush=True)
            return row
        a, b = [test(x, max(3000, args.toys // 4), args.seed) for x in endpoints]
        slope = (b["power"] - a["power"]) / (endpoints[1] - endpoints[0])
        estimate = endpoints[0] + (power - a["power"]) / slope
        for shift in (10, 20):
            near = test(estimate, args.toys, args.seed + shift)
            estimate += (power - near["power"]) / slope
        results["crossings"].append(dict(target_power=power, exposure_fb_inverse=estimate,
            design_months=1.2*estimate, power_slope_per_fb_inverse=slope,
            exposure_mc_se_from_power=near["power_mc_se"]/slope, observations=[a, b, near]))
    anchor_lumi = results["crossings"][0]["exposure_fb_inverse"]
    exp = Experiment(masses, inputs, peak, four, anchor_lumi)
    if (lo, hi) == (110., 135.):
        critical_alpha = results["crossings"][0]["observations"][-1]["local_alpha"]
        plot_local_probabilities(exp, critical_alpha, args.output.parent)
    factor, error = exp.null_trials(ALPHA, args.toys, args.seed + 91)
    p_scan = factor * ALPHA
    results["global_at_local_five"] = dict(luminosity_fb_inverse=anchor_lumi, trials_factor=factor,
        trials_factor_mc_se=error, global_p=p_scan, global_z=norm.isf(p_scan))

    # Validate the conditional proposal with ordinary null toys in an accessible tail.
    test_alpha = .002
    direct, direct_se = exp.ordinary_null(test_alpha, args.toys, args.seed + 92)
    factor, factor_se = exp.null_trials(test_alpha, args.toys, args.seed + 93)
    inferred, inferred_se = test_alpha * factor, test_alpha * factor_se
    if abs(direct-inferred) > 6 * np.hypot(direct_se, inferred_se):
        raise ValueError("Ordinary background ensembles disagree with conditional sampling")
    fixed_lumi = published["single_experiment"]["combined_median_fb_inverse"]
    single = Experiment(np.array([125.]), inputs, peak, four, fixed_lumi)
    p, dp, _, _ = single.signal_power(ALPHA, args.toys, args.seed + 94)
    if abs(p - .5) > 6*dp:
        raise ValueError("Single-mass benchmark differs from the original median reach")
    results["checks"] = dict(single_mass_power_at_original_median=p, single_mass_power_mc_se=dp,
        local_alpha_for_direct_comparison=test_alpha, direct_global_p=direct, direct_mc_se=direct_se,
        importance_global_p=inferred, importance_mc_se=inferred_se)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
