"""Reconstruct the published lepton-tagged diphoton rates in GEM Tables 2-9/2-10."""
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm, poisson


HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "gem_inputs.json").read_text())
TAG_MASSES = np.asarray(DATA["tagged_mass_gev"])


def background_density(mass):
    """Accepted continuum in fb/GeV, constant beyond the outer anchors for window edges."""
    return np.interp(mass, TAG_MASSES, DATA["background_density_total_fb_per_gev"])


def resolution(mass):
    """Gaussian photon-pair mass width in GeV from GEM Table 6."""
    return np.interp(mass, DATA["resolution_mass_gev"], DATA["resolution_sigma_gev"])


def selected_signal(mass):
    """Historical WH + ttH signal in fb after all cuts, integrated over the entire peak."""
    return np.interp(mass, TAG_MASSES, DATA["accepted_signal_before_mass_window_fb"])


def historical_generated_signal(mass):
    """Historical cross section in fb into the lepton + diphoton final state before cuts."""
    return np.interp(mass, TAG_MASSES, DATA["sigma_times_decay_fb"])


def historical_table_check():
    """Interpolate the exact Poisson tail in the expected count as specified in GEM Eq. 10."""
    half_width = DATA["mass_window_half_width_sigma"]
    signal_fb = selected_signal(TAG_MASSES) * (2 * norm.cdf(half_width) - 1)
    background_fb = 2 * half_width * resolution(TAG_MASSES) * background_density(TAG_MASSES)
    rows = []
    for lumi, published in DATA["tabulated_peak_significance"].items():
        background = float(lumi) * background_fb
        count = float(lumi) * signal_fb + background
        low = np.floor(count)
        lower_tail = poisson.sf(low - 1, background)
        upper_tail = poisson.sf(low, background)
        local_z = norm.isf(lower_tail + (count - low) * (upper_tail - lower_tail))
        np.testing.assert_allclose(local_z, published, rtol=0., atol=.08)
        rows.append(dict(luminosity_fb_inverse=float(lumi), mass_gev=TAG_MASSES.tolist(),
                         reconstructed_z=local_z.tolist(), reported_z=published))
    return rows


if __name__ == "__main__":
    print(json.dumps(historical_table_check(), indent=2))
