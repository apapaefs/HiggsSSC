"""Global discovery probability for the common GEM inclusive and lepton-tagged diphoton mass scan."""

import argparse
from copy import copy
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
sys.path.insert(0, str(HERE.parent))
import scan
from combined import Experiment
from gem_rates import background_density, resolution, selected_signal
from likelihood import ALL, ALPHA, BLIND, Calibration, Shape, joint_power


MASSES = scan.MASSES
ANCHOR = np.flatnonzero(MASSES == 125.).item()
TAG = MASSES <= 140.


def histogram_at_edges(edges, rates, target_edges):
    """Integrate piecewise uniform rates; the input values are fb per bin."""
    if not (len(edges) == len(rates) + 1 and np.all(np.diff(edges) > 0)
            and np.isfinite(rates).all() and np.isfinite(edges).all()):
        raise ValueError("Invalid associated diphoton histogram")
    return np.diff(np.interp(target_edges, edges, np.r_[0., np.cumsum(rates)]))


class Inputs:
    def __init__(self, ggf_veto=None, remove_tagged_background=False, associated=None):
        with np.load(HERE.parent / "projection/pooled.npz", allow_pickle=False) as data:
            edges = data["edges_gev"]
            counts = dict(zip(data["sample"], data["rate_fb"]))
            self.cov = sum(c for name, c in zip(data["sample"], data["covariance_fb2"])
                           if name in scan.BACKGROUNDS)
        self.background = sum(counts[k] for k in scan.BACKGROUNDS)
        # Exclusive-veto sensitivity: even if every tagged background had appeared in
        # the inclusive category, this is the fractional bin subtraction required.
        lower, upper = edges[:-1], edges[1:]
        density = (background_density(lower) + background_density(upper)) / 2
        self.overlap_fraction_bound = float(np.max(density / self.background))
        centers = (lower + upper) / 2
        fit_range = (centers >= 105) & (centers < 160)
        self.signal = counts[scan.SIGNAL]
        self.ggf_veto_fraction = 0.
        if ggf_veto is not None:
            # Count all ggF lepton candidates satisfying the pre-ID guard
            # as vetoed; give no credit for a ggF signal in the tag.
            guard = json.loads((HERE / "ggf_overlap.json").read_text())
            self.ggf_veto_fraction = guard["results"][ggf_veto]["fraction"]
            self.signal = self.signal * (1 - self.ggf_veto_fraction)
        if remove_tagged_background:
            # Assume the entire tagged density had also been in the inclusive
            # background. A binwise common rescaling preserves its mixture.
            for k in scan.BACKGROUNDS:
                counts[k] = counts[k] * (1 - density / self.background)
            self.background = self.background - density
        if associated is not None:
            self.signal = self.signal + histogram_at_edges(associated[0], associated[1], edges)
        if ggf_veto is not None or remove_tagged_background or associated is not None:
            shape = Shape("exponential", 2)
            b, s = self.background[fit_range], self.signal[fit_range]
            theta, continuum, _, _ = shape.fit(b, ~BLIND)
            pars, _, _, derivative = shape.fit(b, ALL, s)
            covariance = self.cov[np.ix_(fit_range, fit_range)]
            width = np.sqrt(derivative[-1] @ covariance @ derivative[-1])
            self.fixed = Calibration(shape, theta, continuum, s, pars[-1], width)
        else:
            _, self.fixed = scan.calibration_at_125()
        pars, _ = scan.fit_peak(self.signal[fit_range])
        self.weights = scan.score_weights(edges, counts, scan.templates(edges, pars, MASSES), MASSES)
        mass = MASSES[TAG]
        lo, hi = mass - 2 * resolution(mass), mass + 2 * resolution(mass)
        # Every window endpoint and every density knot is a cell edge. Trapezoid
        # integration is exact for the clamped, piecewise-linear GEM density.
        self.edges = np.unique(np.r_[lo, hi, [80., 90., 100., 120., 140.]])
        self.start = np.zeros(len(MASSES), dtype=int)
        self.stop = np.zeros(len(MASSES), dtype=int)
        self.start[TAG] = np.searchsorted(self.edges, lo)
        self.stop[TAG] = np.searchsorted(self.edges, hi)
        density = background_density(self.edges)
        self.tag_background = np.diff(self.edges) * (density[1:] + density[:-1]) / 2
        self.tag_signal_shape = np.diff(ndtr((self.edges - 125.) / resolution(125.)))
        self.clipped_tagged_fb = 0.

    def measured_tag(self, edges, spectrum, selected_fb):
        model = copy(self)
        raw = histogram_at_edges(edges, spectrum, self.edges)
        if selected_fb <= 0 or not np.isfinite(selected_fb):
            raise ValueError("A positive tagged signal rate is required")
        # Negative signed cells can occur in very sparse tails. Enforce
        # nonnegative Poisson means while keeping the true window yield and
        # the spectrum integral across all tested windows exactly.
        anchor_mask = np.zeros(len(raw), dtype=bool)
        anchor_mask[self.start[ANCHOR]:self.stop[ANCHOR]] = True
        positive = np.maximum(raw, 0.)
        for mask in (anchor_mask, ~anchor_mask):
            if raw[mask].sum() < 0:
                raise ValueError("Negative integrated tagged rate inside or outside the 125 window")
            if positive[mask].sum() > 0:
                positive[mask] *= raw[mask].sum() / positive[mask].sum()
        model.tag_signal_shape = positive / selected_fb
        model.clipped_tagged_fb = float(-np.minimum(raw, 0.).sum())
        return model


class TaggedExperiment(Experiment):
    """Use the shared rare-tail and signal-power engines for a Gaussian score plus Poisson windows."""

    def __init__(self, inputs, exposure, signal_fb, background_scale):
        self.lumi, self.masses, self.anchor = exposure, MASSES, ANCHOR
        self.strength = np.sqrt(inputs.fixed.q(exposure))
        self.field, self.response, _ = scan.score_field(inputs.weights, inputs.background, inputs.cov,
                                                      inputs.signal, exposure, ANCHOR)
        self.start, self.stop = inputs.start, inputs.stop
        self.bkg = background_scale * inputs.tag_background
        self.sig = signal_fb * inputs.tag_signal_shape
        self.b_window = self.window_counts(self.bkg[None, :])[0]
        self.true_signal = self.sig[self.start[ANCHOR]:self.stop[ANCHOR]].sum()
        # The signal at another mass enters as a test weight only. The
        # signal-injected experiment always has a 125 GeV Higgs.
        self.log_weight = np.zeros(len(MASSES))
        reference_window = selected_signal(MASSES[TAG]) / selected_signal(125.) * self.true_signal
        self.log_weight[TAG] = np.log1p(reference_window / self.b_window[TAG])
        self.n = np.arange(1 + int(poisson.isf(1e-16, self.b_window.max() * exposure)))
        self.count_pmfs = poisson.pmf(self.n, self.b_window[:, None] * exposure)
        self.count_scores = self.log_weight[:, None] * self.n


def fixed_reach(inputs, signal, background, power, alpha):
    def gap(lumi):
        return joint_power(inputs.fixed.q, lumi, signal, background, alpha=alpha) - power
    return brentq(gap, .05, 35.)


def search(inputs, signal_fb, bscale, toys, seed, max_iterations=4):
    def experiment(lumi):
        return TaggedExperiment(inputs, lumi, signal_fb, bscale)

    setup = experiment(5.)
    sig, bkg = setup.true_signal, setup.b_window[ANCHOR]
    result = dict(accepted_tagged_signal_fb=signal_fb, tagged_background_scale=bscale,
                  removed_ggf_fraction=inputs.ggf_veto_fraction,
                  signal_in_true_window_fb=sig, background_in_true_window_fb=bkg,
                  negative_tagged_cell_magnitude_fb=inputs.clipped_tagged_fb, crossings=[])
    for index, target in enumerate((.5, .9)):
        # A 125-only test with alpha/M is sufficient for a global discovery.
        # It supplies an upper bracket and a model-independent correlation check.
        local = fixed_reach(inputs, sig, bkg, target, ALPHA)
        upper = fixed_reach(inputs, sig, bkg, target, ALPHA / len(MASSES))
        lumi = .965 * upper
        alpha = ALPHA / 135.
        step_seed = seed + 100 * index
        for step in range(max_iterations):
            exp = experiment(lumi)
            n = max(4000, toys // 4) if step == 0 else toys
            trials, trials_se = exp.null_trials(alpha, n, step_seed + 1)
            alpha = ALPHA / trials
            power, error, p125, e125 = exp.signal_power(alpha, n, step_seed + 2)
            single = joint_power(inputs.fixed.q, lumi, sig, bkg, alpha=alpha)
            if abs(p125 - single) > 6 * e125 + 1.e-6:
                raise ValueError("Tagged true-mass toy power differs from its analytic Poisson mixture")
            def local_power(l):
                return joint_power(inputs.fixed.q, l, sig, bkg, alpha=alpha)
            slope = (local_power(lumi + .05) - local_power(lumi - .05)) / .1
            correction = (target - power) / slope
            if step > 0 and abs(correction) < max(.006, 1.4 * error / slope):
                break
            lumi += correction
        else:
            raise RuntimeError(f"Exposure failed to converge for target {target}")
        if not local <= lumi <= upper + .01:
            raise ValueError("Scan falls outside the fixed-mass and Bonferroni exposures")
        check, check_se = exp.null_trials(alpha, toys, step_seed + 3)
        if abs(check - trials) > 6 * np.hypot(check_se, trials_se):
            raise ValueError("Independent seeded global-tail calibration disagrees")
        result["crossings"].append(dict(target_power=target, luminosity_fb_inverse=lumi,
            design_months=1.2*lumi, local_z_for_global_five=norm.isf(alpha),
            effective_trials=trials, trials_mc_se=trials_se,
            independent_global_p=alpha*check, independent_global_p_mc_se=alpha*check_se,
            power=power, power_mc_se=error, exposure_mc_se_from_power=error/slope,
            true_mass_power_analytic=single, true_mass_power_toys=p125,
            fixed_125_local_luminosity_fb_inverse=local, bonferroni_sufficient_fb_inverse=upper))
    return result


def plot(inputs, nominal, output):
    import matplotlib.pyplot as plt
    from plots import BLUE, INK, ORANGE, save

    row = nominal["crossings"][0]
    exp = TaggedExperiment(inputs, row["luminosity_fb_inverse"], nominal["accepted_tagged_signal_fb"], 1.)
    photon = exp.strength * exp.response
    signal = exp.window_counts(exp.sig[None, :])[0]
    null = exp.lumi * exp.b_window
    alt = null + signal * exp.lumi
    k = np.arange(2 + int(poisson.isf(1.e-15, alt.max())))
    combined = []
    for weight, b, sb, shift in zip(exp.log_weight, null, alt, photon):
        def median_tail(t):
            return poisson.pmf(k, sb) @ ndtr((weight*k - t)/exp.strength + shift) - .5
        median = brentq(median_tail, -20*exp.strength, weight*k[-1] + 20*exp.strength)
        combined.append(poisson.pmf(k, b) @ ndtr((weight*k - median)/exp.strength))
    tagged_p = poisson.sf(poisson.ppf(.5, alt[TAG]) - 1., null[TAG])
    fig, ax = plt.subplots(figsize=(8.2, 4.7), layout="constrained")
    ax.semilogy(MASSES, combined, color=INK, lw=2, label="Combined")
    ax.semilogy(MASSES, norm.sf(photon), color=BLUE, ls="--", label="Untagged photons")
    ax.semilogy(MASSES[TAG], tagged_p, color=ORANGE, ls="-.", label="Lepton-tagged photons")
    ax.axhline(ALPHA, color="#777777", lw=1, ls=":", label=r"Local $5\sigma$")
    ax.axhline(norm.sf(row["local_z_for_global_five"]), color="#66334e", ls=":",
               label=r"Combined global $5\sigma$ threshold")
    ax.set(title=rf"One SSC experiment, 40 TeV, {exp.lumi:.2f} fb$^{{-1}}$",
           xlabel="Tested Higgs mass [GeV]", ylabel=r"Median local $p$-value for a 125 GeV Higgs",
           xlim=(80, 160), ylim=(1.e-10, 1.8))
    ax.legend(loc="center left", fontsize=8.4, frameon=False)
    save(fig, output, "combined-local-p")


def direct_models(data, transfer_fb, original_inputs):
    edges = np.asarray(data["mass_edges_gev"])
    tagged = np.asarray(data["tagged_spectrum_fb_per_bin"])
    inclusive = np.asarray(data["inclusive_spectrum_fb_per_bin"])
    vetoed = np.asarray(data["untagged_associated_spectrum_fb_per_bin"])
    rate = float(data["selected_fb"])
    rate_se = float(data["selected_mc_se_fb"])

    def with_tag(inputs):
        return inputs.measured_tag(edges, tagged, rate)

    nominal = with_tag(Inputs(ggf_veto="pt20_eta2p5", associated=(edges, vetoed)))
    guard = with_tag(Inputs(ggf_veto="pt15_eta3", remove_tagged_background=True, associated=(edges, vetoed)))
    return [
        ("nominal", rate, 1., nominal),
        ("same_production_inclusive_only", 0., 1., Inputs(associated=(edges, inclusive))),
        ("ggf_plus_tag_only", rate, 1., with_tag(Inputs(ggf_veto="pt20_eta2p5"))),
        ("double_tagged_background", rate, 2., nominal),
        ("exclusive_overlap_check", rate, 1., guard),
        ("tagged_signal_minus_one_mc_se", rate - rate_se, 1., nominal),
        ("tagged_signal_plus_one_mc_se", rate + rate_se, 1., nominal),
        ("modern_rates_gem_acceptance", transfer_fb, 1., original_inputs),
        ("original_gem_tagged_signal", float(selected_signal(125.)), 1., original_inputs),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--transfer-only", action="store_true",
                        help="Use historical tagged acceptance for WH/ttH with the previous ggF-only inclusive model")
    parser.add_argument("--signal-rate", type=float, help="Override the selected signal in fb for --transfer-only")
    parser.add_argument("--toys", type=int, default=40000)
    parser.add_argument("--seed", type=int, default=725081)
    args = parser.parse_args()
    if args.signal_rate is not None and not args.transfer_only:
        parser.error("--signal-rate requires --transfer-only")
    transfer_fb = json.loads((HERE / "signal/rates.json").read_text())["selected_fb"]
    if args.signal_rate is not None:
        transfer_fb = args.signal_rate
    if not np.isfinite(transfer_fb) or transfer_fb <= 0 or args.toys < 4000:
        parser.error("A positive finite signal and at least 4,000 experiments are required")
    inputs = Inputs()
    if args.transfer_only:
        models = [("nominal", transfer_fb, 1., inputs), ("double_tagged_background", transfer_fb, 2., inputs),
                  ("original_gem_tagged_signal", float(selected_signal(125.)), 1., inputs),
                  ("exclusive_overlap_check", transfer_fb, 1., Inputs("pt15_eta3", True))]
    else:
        models = direct_models(json.loads((HERE / "signal/direct.json").read_text()), transfer_fb, inputs)
    scenarios = {}
    for name, s, b, model in models:
        scenarios[name] = search(model, s, b, args.toys, args.seed)
        print(name, [(c["luminosity_fb_inverse"], c["design_months"]) for c in scenarios[name]["crossings"]],
              flush=True)
    nominal = scenarios["nominal"]
    nominal_inputs = models[0][-1]
    exp = TaggedExperiment(nominal_inputs, nominal["crossings"][0]["luminosity_fb_inverse"],
                           nominal["accepted_tagged_signal_fb"], 1.)
    alpha_test = .002
    direct, direct_se = exp.ordinary_null(alpha_test, args.toys, args.seed+402)
    moderate_factor, moderate_error = exp.null_trials(alpha_test, args.toys, args.seed+403)
    if abs(direct - alpha_test * moderate_factor) > 6 * np.hypot(direct_se, alpha_test*moderate_error):
        raise ValueError("Direct null toys disagree with conditional importance sampling")
    published = json.loads((HERE.parent / "diphoton-scan.json").read_text())
    control = []
    for row in published["discoveries"]:
        exp = TaggedExperiment(inputs, row["luminosity_fb_inverse"], 0., 1.)
        alpha = norm.sf(row["local_z_for_global_five"])
        power, power_se, _, _ = exp.signal_power(alpha, args.toys, args.seed+404)
        factor, error = exp.null_trials(alpha, args.toys, args.seed+405)
        if (abs(power - row["power"]) > 6*power_se + .004
                or abs(factor-row["effective_trials"]) > 6*np.hypot(error, row["check_trials_standard_error"])):
            raise ValueError("Disabling the tagged signal fails the original diphoton benchmark")
        control.append(dict(luminosity_fb_inverse=exp.lumi, power=power, power_mc_se=power_se,
                            trials_factor=factor, trials_factor_mc_se=error))
    args.output.mkdir(parents=True, exist_ok=True)
    plot(nominal_inputs, nominal, args.output)
    report = dict(inclusive_mass_range_gev=[80., 160.], tagged_mass_range_gev=[80., 140.],
        mass_step_gev=.5, true_mass_gev=125., design_fb_inverse_per_year=10.,
        signal_model="historical_acceptance_transfer" if args.transfer_only else "direct_exclusive_ggf_wh_tth",
        toys=args.toys, seed=args.seed, maximum_fraction_of_inclusive_background_in_tagged_spectrum=
        inputs.overlap_fraction_bound, scenarios=scenarios, disabled_tag_checks=control,
        moderate_tail_check=dict(local_alpha=alpha_test, direct_probability=direct,
            direct_standard_error=direct_se, importance_probability=alpha_test*moderate_factor,
            importance_standard_error=alpha_test*moderate_error))
    (args.output / "combined.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
