"""Sufficient exposures from Bonferroni testing on a predeclared finite mass grid.

Only the on-grid 125 GeV test contributes to the power calculation below.
The global size bound assumes each individual mass test has valid local size.
"""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import norm

DISCOVERY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DISCOVERY))
from input_data import SIGNAL
from likelihood import ALPHA, BLIND, Calibration, Shape, crossing, fourlepton_counting_reach, joint_power


def evaluate():
    published = DISCOVERY / "published"
    inputs = json.loads((published / "results.json").read_text())
    spectrum = np.genfromtxt(published / "mass_spectrum.csv", delimiter=",", names=True)
    shape = Shape("exponential", 2)
    theta, continuum, _, _ = shape.fit(spectrum["background_fb"], ~BLIND)
    model = inputs["models"]["exponential2"]
    cal = Calibration(shape, theta, continuum, spectrum[SIGNAL], model["spurious_mu"], model["sigma_mu"])
    rows = []
    for lo, hi in [(125, 125), (115, 130), (110, 135)]:
        count = 2 * (hi - lo) + 1
        alpha = ALPHA / count
        critical = norm.isf(alpha)
        for power in (.5, .9):
            exposures = {"diphoton": crossing(lambda l: cal.q(l) - (critical + norm.ppf(power))**2)}
            for scenario, four in inputs["fourlepton"].items():
                exposures["fourlepton_" + scenario], _ = fourlepton_counting_reach(
                    four["signal_fb"], four["background_fb"], power, alpha)
                def power_at_lumi(l):
                    return joint_power(cal.q, l, four["signal_fb"], four["background_fb"], alpha=alpha)
                exposures["combined_" + scenario] = crossing(lambda l: power_at_lumi(l) - power)
            rows.append(dict(scan_gev=[lo, hi], spacing_gev=.5, mass_hypotheses=count,
                             true_mass_gev=125., local_critical_z=critical, power=power,
                             sufficient_luminosity_fb_inverse=exposures))
    return rows


if __name__ == "__main__":
    print(json.dumps(evaluate(), indent=2, allow_nan=False))
