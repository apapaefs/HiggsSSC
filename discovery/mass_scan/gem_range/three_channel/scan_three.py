"""Global SSC Higgs discovery scan using an inclusive photon score and two Poisson counts."""
import argparse
from dataclasses import dataclass
from functools import lru_cache
import json
import os
from pathlib import Path
import sys

if __name__ == "__main__":
    for key in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[key] = "1"

import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri
from scipy.stats import norm, poisson

HERE = Path(__file__).resolve().parent
BACKGROUND = HERE.parent / "lepton_tag/background_mc"
sys.path.insert(0, str(BACKGROUND))
from refit import (ANCHOR, WINDOW, TaggedExperiment, add_zee_complement, bin_fractions,
                   discovery_model, load_spectrum, project_background)
from tagged_scan import MASSES
from combined import physical_cells
from likelihood import ALPHA

BATCH = 512
FOUR = (MASSES >= 110) & (MASSES <= 135)


@dataclass
class Channel:
    start: np.ndarray
    stop: np.ndarray
    background: np.ndarray
    signal: np.ndarray
    weights: np.ndarray

    def window_sums(self, cells):
        cumulative = np.pad(cells.cumsum(axis=-1), [(0, 0)]*(cells.ndim-1) + [(1, 0)])
        return cumulative[..., self.stop] - cumulative[..., self.start]

    def random_windows(self, rng, size, lumi, injected=False, indices=None, target=None):
        cells = rng.poisson(lumi*(self.background + self.signal if injected else self.background),
                            (size, len(self.background)))
        if indices is not None:
            # Given the count at the chosen trial mass, other cells remain independent.
            lo, hi = self.start[indices], self.stop[indices]
            cell_id = np.arange(len(self.background))
            cells[(cell_id >= lo[:, None]) & (cell_id < hi[:, None])] = 0
            cdf = np.r_[0., self.background.cumsum()]
            row = np.repeat(np.arange(size), target)
            u = cdf[lo[row]] + rng.random(len(row)) * (cdf[hi[row]] - cdf[lo[row]])
            col = np.searchsorted(cdf, u, side="right") - 1
            np.add.at(cells, (row, col), 1)
        counts = self.window_sums(cells)
        if indices is not None:
            np.testing.assert_array_equal(counts[np.arange(size), indices], target)
        return counts


def load_model(tag_background):
    signal, model = discovery_model()
    if tag_background != "gem-total":
        core = load_spectrum(BACKGROUND / "final/background.json.gz")
        add_zee_complement(core, load_spectrum(BACKGROUND / "zee-confirmation-spectrum.json.gz"))
        edges = np.asarray(core["edges_gev"])
        rates = np.asarray(core["rates_fb"])
        if tag_background == "gem-electron":
            electron = np.asarray(core["sources"]["z_ee_a"]["rates_fb"])
            groups = json.loads((BACKGROUND / "final/combined-mc.json").read_text())["window_components"]
            reference = groups["neutral_current_electron_fake"]["gem_fb"]
            rates = rates + (reference/(bin_fractions(edges, *WINDOW) @ electron) - 1.) * electron
        model.tag_background = project_background(model, edges, rates)
    four = json.loads((HERE.parents[1] / "fourlepton.json").read_text())
    edges = np.array(four["edges_gev"])
    lo = np.zeros(len(MASSES), np.int64)
    hi = lo.copy()
    lo[FOUR] = np.searchsorted(edges, MASSES[FOUR] - 5.)
    hi[FOUR] = np.searchsorted(edges, MASSES[FOUR] + 5.)
    np.testing.assert_allclose(edges[lo[FOUR]], MASSES[FOUR] - 5.)
    np.testing.assert_allclose(edges[hi[FOUR]], MASSES[FOUR] + 5.)
    four = Channel(lo, hi, physical_cells(np.array(four["qq_continuum_fb"]) + four["gg_continuum_fb"], edges),
                   physical_cells(np.array(four["signal_fb"]), edges), np.zeros(len(MASSES)))
    four.weights[FOUR] = np.log1p(four.window_sums(four.signal)[ANCHOR] /
                                  four.window_sums(four.background)[FOUR])
    return model, signal["selected_fb"], four


def poisson_mixture(means, weights):
    """The two independent window counts and their local log likelihood ratio offsets."""
    counts = [np.arange(1 + int(poisson.isf(1.e-16, mean))) for mean in means]
    first, second = counts
    pmf = np.outer(poisson.pmf(first, means[0]), poisson.pmf(second, means[1])).ravel()
    scores = (first[:, None]*weights[0] + second*weights[1]).ravel()
    return pmf, scores, len(second)


def threshold(strength, pmf, scores, alpha):
    return brentq(lambda cut: pmf @ ndtr((scores-cut)/strength) - alpha,
                  -10*strength, scores.max()+12*strength)


def fixed_power(strength, lumi, signals, backgrounds, alpha):
    weights = np.log1p(signals/backgrounds)
    pmf, scores, _ = poisson_mixture(lumi*backgrounds, weights)
    cut = threshold(strength, pmf, scores, alpha)
    pmf, scores, _ = poisson_mixture(lumi*(backgrounds + signals), weights)
    return pmf @ ndtr((scores-cut)/strength + strength)


def mean_se(values):
    return float(np.mean(values)), float(np.std(values, ddof=1) / np.sqrt(len(values)))


class JointScan:
    def __init__(self, model, lumi):
        inputs, signal_fb, four = model
        tagged = TaggedExperiment(inputs, lumi, signal_fb, 1.)
        self.lumi, self.field, self.response = lumi, tagged.field, tagged.response
        self.strength = tagged.strength
        self.channels = [Channel(tagged.start, tagged.stop, tagged.bkg, tagged.sig, tagged.log_weight), four]
        rates = np.column_stack([c.window_sums(c.background) for c in self.channels])
        self.mixtures = [poisson_mixture(b*lumi, w) for b,w in zip(rates,
                          np.column_stack([c.weights for c in self.channels]))]

    @lru_cache(maxsize=8)
    def local(self, alpha):
        cuts, targets, tails = [], [], []
        for pmf, scores, _ in self.mixtures:
            cut = threshold(self.strength, pmf, scores, alpha)
            tail = ndtr((scores-cut)/self.strength)
            cdf = np.cumsum(pmf*tail)
            if abs(cdf[-1]/alpha - 1.) > 1.e-7:
                raise RuntimeError("Failed to normalize the local conditional distribution")
            cdf /= cdf[-1]
            cuts.append(cut)
            targets.append(cdf)
            tails.append(tail)
        return np.array(cuts), targets, tails

    def null(self, alpha, toys, seed, conditional=True):
        rng = np.random.default_rng(seed)
        cuts, targets, tails = self.local(alpha)
        outcomes = []
        for first in range(0, toys, BATCH):
            n = min(BATCH, toys-first)
            noise = rng.standard_normal((n, self.field.shape[1]))
            if conditional:
                j = rng.integers(len(cuts), size=n)
                u = rng.random(n)
                truncated = rng.random(n)
                ks = np.zeros((2, n), int)
                for mass in np.unique(j):
                    mask = j == mass
                    cells = np.searchsorted(targets[mass], u[mask], side="right")
                    ks[:, mask] = np.divmod(cells, self.mixtures[mass][2])
                    truncated[mask] = -ndtri(tails[mass][cells] * (1.-truncated[mask]))
                direction = self.field[j]
                noise += (truncated - (noise*direction).sum(axis=1))[:, None] * direction
            scores = self.strength*(noise @ self.field.T)
            for i, channel in enumerate(self.channels):
                windows = channel.random_windows(rng, n, self.lumi, indices=j if conditional else None,
                                                 target=ks[i] if conditional else None)
                scores += windows * channel.weights
            exceeds = scores >= cuts
            if conditional:
                rows = np.arange(n)
                if np.any(scores[rows, j] < cuts[j] - 1.e-8):
                    raise RuntimeError("Toy misses its required conditional crossing")
                exceeds[rows, j] = True
                outcomes.extend(len(cuts)/exceeds.sum(axis=1))
            else:
                outcomes.extend(exceeds.any(axis=1))
        return mean_se(outcomes)

    def power(self, alpha, toys, seed):
        rng = np.random.default_rng(seed)
        cuts, _, _ = self.local(alpha)
        corr = self.field @ self.field[ANCHOR]
        positive, negative = corr > 1.e-12, corr < -1.e-12
        neutral = ~(positive | negative)
        residual = self.field - corr[:,None]*self.field[ANCHOR]
        all_power, true_power = [], []
        for first in range(0, toys, BATCH):
            n = min(BATCH, toys-first)
            scores = self.strength * self.response + rng.standard_normal((n, residual.shape[1])) @ residual.T
            for channel in self.channels:
                scores += channel.random_windows(rng, n, self.lumi, injected=True)*channel.weights/self.strength
            needed = cuts/self.strength - scores
            right = np.min(needed[:, positive]/corr[positive], axis=1)
            left = np.max(needed[:, negative]/corr[negative], axis=1) if np.any(negative) else -np.inf
            probability = np.minimum(1., ndtr(-right) + ndtr(left))
            if neutral.any():
                probability[np.any(needed[:,neutral] <= 0., axis=1)] = 1.
            all_power.extend(probability)
            true_power.extend(ndtr(-needed[:,ANCHOR]))
        return (*mean_se(all_power), *mean_se(true_power))


def run(toys, seed, variant):
    model = load_model(variant)
    first = JointScan(model, 3.)
    signals = np.array([c.window_sums(c.signal)[ANCHOR] for c in first.channels])
    backgrounds = np.array([c.window_sums(c.background)[ANCHOR] for c in first.channels])

    def fixed(lumi, alpha):
        return fixed_power(np.sqrt(model[0].fixed.q(lumi)), lumi, signals, backgrounds, alpha)

    report = dict(toys=toys, seed=seed, tag_background=variant, true_mass_gev=125.,
                  inclusive_range_gev=[80,160], tagged_range_gev=[80,140], fourlepton_range_gev=[110,135],
                  step_gev=.5, design_fb_inverse_per_year=10, fourlepton_background="central qq+gg",
                  tagged_signal_fb=signals[0], tagged_background_fb=backgrounds[0], four_signal_fb=signals[1],
                  four_background_fb=backgrounds[1], crossings=[])
    for index, target in enumerate((.5,.9)):
        local = brentq(lambda l: fixed(l, ALPHA) - target, .1, 12.)
        upper = brentq(lambda l: fixed(l, ALPHA/len(MASSES)) - target, .1, 12.)
        lumi, alpha = .95*upper, ALPHA/133.
        base = seed + index*100
        for step in range(12):
            experiment = JointScan(model, lumi)
            n = max(4000, toys//4) if step == 0 else toys
            trials, trials_se = experiment.null(alpha, n, base)
            updated = ALPHA/trials
            power, error, p125, e125 = experiment.power(updated, n, base+1)
            slope = (fixed(lumi+.03, updated)-fixed(lumi-.03, updated))/.06
            correction = (target - power)/slope
            print(f"{target:.0%}: L={lumi:.5f}, power {power:.5f} +/- {error:.5f}, trials {trials:.2f}",
                  flush=True)
            if step and abs(correction) < max(.004, 1.25*error/slope):
                alpha = updated
                break
            lumi += correction
            alpha = updated
        else:
            raise RuntimeError("Exposure failed to converge")
        if abs(p125 - fixed(lumi, alpha)) > 6*e125+1e-6 or lumi > upper+.01:
            raise RuntimeError("Anchored true-mass power or exposure bound failed")
        check, check_se = experiment.null(alpha, toys, base+2)
        if abs(check-trials) > 6*np.hypot(check_se,trials_se):
            raise RuntimeError("Rare-tail null validation failed")
        report["crossings"].append(dict(target_power=target, luminosity_fb_inverse=lumi, design_months=1.2*lumi,
            power=power, power_mc_se=error, exposure_mc_se_from_power=error/slope, local_z=norm.isf(alpha),
            effective_trials=trials, effective_trials_mc_se=trials_se, independent_global_p=alpha*check,
            independent_global_p_mc_se=alpha*check_se, true_mass_power_toys=p125, true_mass_power_mc_se=e125,
            true_mass_power_analytic=fixed(lumi, alpha), fixed_mass_local_luminosity_fb_inverse=local,
            bonferroni_sufficient_fb_inverse=upper))
    exp = JointScan(model, report["crossings"][0]["luminosity_fb_inverse"])
    alpha = .002
    direct, de = exp.null(alpha, toys, seed+500, conditional=False)
    imp, ie = exp.null(alpha, toys, seed+501)
    if abs(direct-alpha*imp) > 6*np.hypot(de, alpha*ie):
        raise RuntimeError("Direct and sampled moderate null tails disagree")
    report["moderate_null_check"] = dict(local_alpha=alpha, direct_p=direct, direct_mc_se=de,
                                         importance_p=alpha*imp, importance_mc_se=alpha*ie)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--toys", type=int, default=40000)
    parser.add_argument("--seed", type=int, default=923125)
    parser.add_argument("--tag-background", choices=("nominal","gem-total","gem-electron"), default="nominal")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.toys < 4000 or args.seed < 0:
        parser.error("Use at least 4000 toys and a nonnegative seed")
    report = run(args.toys, args.seed, args.tag_background)
    path = args.output or HERE/f"{args.tag_background}.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
