"""Sufficient exposures from Bonferroni testing on a predeclared finite mass grid.

Only the on-grid 125 GeV test contributes to the power calculation below.
The global size bound assumes each individual mass test has valid local size.

Generate the inclusive diphoton result from the repository root:
    python3 discovery/mass_scan/grid_bounds.py --gem-diphoton > discovery/mass_scan/gem_range/diphoton-grid-bound.json
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import norm

DISCOVERY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DISCOVERY))
from input_data import SIGNAL
from likelihood import ALPHA, BLIND, Calibration, Shape, crossing, fourlepton_counting_reach, joint_power


SPACING_GEV = 0.5


def calibration_at_125():
    published = DISCOVERY / "published"
    inputs = json.loads((published / "results.json").read_text())
    spectrum = np.genfromtxt(published / "mass_spectrum.csv", delimiter=",", names=True)
    shape = Shape("exponential", 2)
    theta, continuum, _, _ = shape.fit(spectrum["background_fb"], ~BLIND)
    model = inputs["models"]["exponential2"]
    cal = Calibration(shape, theta, continuum, spectrum[SIGNAL], model["spurious_mu"], model["sigma_mu"])
    return inputs, cal


def finite_grid_threshold(lo, hi):
    count = int((hi - lo) / SPACING_GEV) + 1
    alpha = ALPHA / count
    return count, alpha, norm.isf(alpha)


def diphoton_exposure(cal, critical, power):
    target = (critical + norm.ppf(power))**2
    return crossing(lambda luminosity: cal.q(luminosity) - target)


def evaluate():
    inputs, cal = calibration_at_125()
    rows = []
    for lo, hi in [(125, 125), (115, 130), (110, 135)]:
        count, alpha, critical = finite_grid_threshold(lo, hi)
        for power in (.5, .9):
            exposures = {"diphoton": diphoton_exposure(cal, critical, power)}
            for scenario, four in inputs["fourlepton"].items():
                exposures["fourlepton_" + scenario], _ = fourlepton_counting_reach(
                    four["signal_fb"], four["background_fb"], power, alpha)
                def power_at_lumi(l):
                    return joint_power(cal.q, l, four["signal_fb"], four["background_fb"], alpha=alpha)
                exposures["combined_" + scenario] = crossing(lambda l: power_at_lumi(l) - power)
            rows.append(dict(scan_gev=[lo, hi], spacing_gev=SPACING_GEV, mass_hypotheses=count,
                             true_mass_gev=125., local_critical_z=critical, power=power,
                             sufficient_luminosity_fb_inverse=exposures))
    return rows


def gem_diphoton_bound():
    inputs, cal = calibration_at_125()
    lo, hi = 80, 160
    count, alpha, critical = finite_grid_threshold(lo, hi)
    rate = inputs["design_fb_inverse_per_year"]
    rows = []
    for power in (.5, .9):
        lumi = diphoton_exposure(cal, critical, power)
        rows.append(dict(power=power, sufficient_luminosity_fb_inverse=lumi,
                         design_months=12 * lumi / rate))
    return dict(channel="inclusive diphoton", scan_gev=[lo, hi], spacing_gev=SPACING_GEV,
                mass_hypotheses=count, true_mass_gev=125.,
                grid_definition="Inclusive endpoints span the original GEM Table 2-3 diphoton benchmarks. "
                                "The 0.5 GeV spacing is a modern analysis choice.",
                global_one_sided_size_upper_bound=ALPHA, local_one_sided_size=alpha,
                local_critical_z=critical,
                size_condition="Bonferroni bounds the false-positive probability at a fixed exposure "
                               "provided each of the 161 local tests has valid size under the null. "
                               "The off-mass tests have not been calibrated here.",
                power_condition="Only the on-grid 125 GeV test contributes to this power calculation, "
                                "using the published calibrated sideband likelihood and Gaussian score "
                                "approximation. Under this model, the exposures suffice for the grid to "
                                "reach at least the stated power.",
                design_fb_inverse_per_year=rate, results=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gem-diphoton", action="store_true",
                        help="evaluate the inclusive-photon 80–160 GeV finite-grid bound")
    args = parser.parse_args()
    print(json.dumps(gem_diphoton_bound() if args.gem_diphoton else evaluate(),
                     indent=2, allow_nan=False))
