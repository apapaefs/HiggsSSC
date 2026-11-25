"""Compare the selected eeγ spectrum with GEM's historical rates and cut study."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BACKGROUND = HERE.parent
BROAD = (80., 140.)
WINDOW = (123.3, 126.7)


def bin_mask(edges, interval):
    for endpoint in interval:
        if not np.any(np.isclose(edges, endpoint, atol=1.e-8, rtol=0)):
            raise ValueError(f"Missing boundary {endpoint} in simulated spectrum")
    return ((edges[:-1] >= interval[0]) & (edges[1:] <= interval[1])).astype(float)


def interpolate_integral(masses, densities, interval):
    lo, hi = interval
    assert masses[0] <= lo < hi <= masses[-1]
    grid = np.r_[lo, masses[(masses > lo) & (masses < hi)], hi]
    y = np.interp(grid, masses, densities)
    return np.diff(grid) @ ((y[1:] + y[:-1]) / 2.)


def run():
    core = BACKGROUND / "final/background.json.gz"
    guard = BACKGROUND / "zee-confirmation-spectrum.json.gz"
    gem_file = BACKGROUND.parent / "gem_inputs.json"
    reports = []
    for path in (core, guard):
        with gzip.open(path, "rt") as stream:
            reports.append(json.load(stream))
    edges = np.array(reports[0]["edges_gev"])
    np.testing.assert_equal(edges, reports[1]["edges_gev"])
    if [d["generation_region"] for d in reports] != ["nominal", "outside_nominal"]:
        raise ValueError("Expected nominal and complementary generation regions")
    sources = [d["sources"]["z_ee_a"] for d in reports]
    rate = sum(np.asarray(s["rates_fb"]) for s in sources)
    covariance = sum(np.asarray(s["covariance_fb2"]) for s in sources)
    gem = json.loads(gem_file.read_text())
    masses = np.array(gem["tagged_mass_gev"])
    densities = np.array(gem["background_components_fb_per_gev"]["Z_gamma_electron_fake"])

    def compare(interval):
        mask = bin_mask(edges, interval)
        measured = mask @ rate
        error = np.sqrt(mask @ covariance @ mask)
        historical = interpolate_integral(masses, densities, interval)
        return dict(interval_gev=interval, simulated_fb=measured, simulated_mc_se_fb=error,
                    gem_interpolated_fb=historical, ratio=measured/historical, ratio_mc_se=error/historical)

    broad = compare(BROAD)
    window = compare(WINDOW)
    fraction = window["simulated_fb"]/broad["simulated_fb"]
    gradient = (bin_mask(edges, WINDOW) - fraction*bin_mask(edges, BROAD))/broad["simulated_fb"]
    se = np.sqrt(gradient @ covariance @ gradient)
    gem_fraction = window["gem_interpolated_fb"]/broad["gem_interpolated_fb"]
    np.testing.assert_allclose(window["ratio"], broad["ratio"] * fraction / gem_fraction)

    # Roundoff bounds assume independently rounded decimal entries and linear
    # interpolation. They carry no interpretation as statistical errors.
    rounded = []
    for shift in itertools.product((-.0005, .0005), repeat=len(masses)):
        y = densities + shift
        b = interpolate_integral(masses, y, BROAD)
        w = interpolate_integral(masses, y, WINDOW)
        rounded.append([b, w, w/b])
    ranges = np.array([np.min(rounded, axis=0), np.max(rounded, axis=0)]).T
    note92 = "https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-92-126.pdf"
    note93 = "https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf"
    return dict(
        sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (core, guard, gem_file)},
        gem_1993=dict(source=note93+"#page=61", table=19, masses_gev=masses.tolist(),
                      selected_density_fb_per_gev=densities.tolist(),
                      interpolation="Piecewise linear in reconstructed diphoton mass",
                      pdf_source=note93+"#page=31", pdf="EHLQ-1",
                      generator_source=note93+"#page=60", enhancement_on_qq_to_z_gamma=1.2,
                      electron_fake_probability_in_tag=.03, photon_track_acceptance_at_working_point=.96,
                      isolation_definition=note93+"#page=17",
                      historical_component_mc_error="Unreported in Table 19",
                      historical_preselection="Table 19 provides intermediate rates for the Higgs signal and post-selection rates for the backgrounds."),
        broad=broad, window=window,
        shape=dict(simulated_window_fraction=fraction, simulated_fraction_mc_se=se,
                   gem_window_fraction=gem_fraction, ratio=fraction/gem_fraction, ratio_mc_se=se/gem_fraction,
                   relation="window ratio = broad selected-rate ratio times window-fraction ratio"),
        gem_rounding_only=dict(broad_fb=ranges[0].tolist(), window_fb=ranges[1].tolist(),
                               window_fraction=ranges[2].tolist(),
                               assumption="The five published densities are rounded to the nearest 0.001 fb/GeV. The bounds exclude modelling and historical MC uncertainty."),
        adjoining_mass_intervals=[compare((float(a), float(b))) for a,b in zip(masses[:-1],masses[1:])],
        gem_1992=dict(
            generator_source=note92+"#page=31", tabulated_hard_rate_pb=14.4,
            process="qqbar/gg -> Z gamma -> l+ l- gamma", lepton_definition=note92+"#page=2",
            flavour_caveat="The report defines l as e or mu and labels the total llgamma. An electron-only weight is unspecified for this entry. Halving the quoted rate assumes it includes both flavors with equal rates.",
            hard_cut="pT(hard photon), pT(Z) >= 10 GeV", qqbar_enhancement=1.2,
            electron_fake_probability=.03, fake_source=note92+"#page=33",
            cut_comparison_source=note92+"#page=34", table=15, interval_gev=[75.,165.],
            photon_isolation_radius=.45, lepton_isolation_radius=.30, lumi_fb_inverse=10.,
            fake_events_before_pair_pt_40=51, fake_events_after_pair_pt_40=32,
            pair_pt_efficiency=32/51, ideal_calo_selected_fb=3.2, approximate_lar_acceptance=.56,
            lar_corrected_selected_fb=3.2*.56,
            geometry="Table 15 gives the ideal calorimeter yields and the preceding paragraph gives an overall LAr acceptance.",
            comparison="The earlier study used a wider counting interval and a different detector simulation. Its numbers are a separate historical cross-check."))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "legacy.json")
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
