"""Predict the high-mass eeγ rates using the nominal modern spectrum and GEM at 80–100 GeV."""
import gzip
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import legacy


def calculate():
    reports = []
    for name in ["final/background.json.gz", "zee-confirmation-spectrum.json.gz"]:
        with gzip.open(legacy.BACKGROUND / name, "rt") as stream:
            reports.append(json.load(stream))
    if [r["generation_region"] for r in reports] != ["nominal", "outside_nominal"]:
        raise ValueError("Expected complementary generation regions")
    edges = np.array(reports[0]["edges_gev"])
    np.testing.assert_array_equal(edges, reports[1]["edges_gev"])
    sources = [r["sources"]["z_ee_a"] for r in reports]
    rates = sum(np.asarray(s["rates_fb"]) for s in sources)
    cov = sum(np.asarray(s["covariance_fb2"]) for s in sources)
    gem = json.loads((legacy.BACKGROUND.parent / "gem_inputs.json").read_text())
    masses = np.array(gem["tagged_mass_gev"])
    densities = np.array(gem["background_components_fb_per_gev"]["Z_gamma_electron_fake"])

    low = (80., 100.)
    low_mask = legacy.bin_mask(edges, low)
    gem_low = legacy.interpolate_integral(masses, densities, low)
    mc_low = low_mask @ rates
    factor = gem_low / mc_low
    results = []
    for window in [(100., 120.), (120., 140.), legacy.WINDOW]:
        mask = legacy.bin_mask(edges, window)
        mc_window = mask @ rates
        predicted = factor * mc_window
        # The same events contribute to the calibration and prediction.
        gradient = factor * (mask - mc_window / mc_low * low_mask)
        results.append(dict(interval_gev=window, unscaled_fb=mc_window, predicted_fb=predicted,
                            predicted_mc_se_fb=np.sqrt(gradient @ cov @ gradient),
                            gem_fb=legacy.interpolate_integral(masses, densities, window)))
    return dict(calibration_interval_gev=low, gem_calibration_fb=gem_low, mc_calibration_fb=mc_low,
                scale=factor, scale_mc_se=factor * np.sqrt(low_mask @ cov @ low_mask) / mc_low,
                comparison=results)


if __name__ == "__main__":
    print(json.dumps(calculate(), allow_nan=False, indent=2))
