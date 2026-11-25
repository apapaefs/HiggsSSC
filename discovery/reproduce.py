#!/usr/bin/env python3
"""Refit the archived 40 TeV data, verify the campaign and regenerate the discovery figures."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import tempfile

# The fit solves small dense systems. One BLAS thread avoids oversubscribing shared hosts.
if __name__ == "__main__":
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[variable] = "1"

import numpy as np
from scipy.stats import chi2, norm

from input_data import BACKGROUNDS, DATA, SCENARIOS, SIGNAL, fourlepton_rates, fourlepton_spectrum, load_archive, pool
from likelihood import ALL, BLIND, EDGES, MODELS, Calibration, Shape, crossing, fourlepton_median, joint_power


def check(actual, expected, label):
    if not np.allclose(actual, expected, atol=3.e-6, rtol=2.e-6):
        raise ValueError(f"Campaign replay differs for {label}")


def analyze():
    records, packets = load_archive()
    fit, reference = records["records/final_fit.json"], records["records/final_reach.json"]
    hist, cov, total, added = pool(records, packets)
    signal = hist[SIGNAL]
    background = sum(hist[key] for key in BACKGROUNDS)
    background_cov = sum(cov[key] for key in BACKGROUNDS)
    check(background, fit["background_per_bin_fb"], "diphoton background")
    check(signal, fit["signal_per_bin_fb"], "diphoton signal")
    if (added != reference["new_generated_events"]
            or any(total[k] != n for k, n in reference["total_generated_by_stratum"].items())):
        raise ValueError("Event exposure differs from campaign")
    s4, b4 = fourlepton_rates(records["records/extended_fourlepton.json"], records["records/reducible_bound.json"])
    rows, models, spectra = [], {}, {}
    def append(model, treatment, channel, scenario, luminosity):
        rows.append(dict(model=model, treatment=treatment, channel=channel, heavy_flavour=scenario,
                         luminosity_fb_inverse=luminosity, design_months=12 * luminosity / 10))
    four = {}
    for scenario, tag in SCENARIOS.items():
        source = fit["fourlepton_rates_fb"][tag]
        check([s4, b4[scenario]], [source["signal_fb"], source["background_fb"]], scenario)
        luminosity, count = fourlepton_median(s4, b4[scenario])
        key = "irreducible_only" if scenario == "central" else "plus_old_joint_95_bound"
        exact = records["records/reducible_bound.json"]["discovery_scenarios"][key]["exact_median_local_poisson"]
        check([luminosity, count], [exact["luminosity_fb_inverse"], exact["critical_count"]], f"4l {scenario}")
        four[scenario] = dict(signal_fb=s4, background_fb=b4[scenario], critical_count=count)
        append("exact_poisson", "count", "fourlepton", scenario, luminosity)
    for family, degree in MODELS:
        name = f"{family}{degree}"
        shape = Shape(family, degree)
        theta, continuum, jac, sideband_grad = shape.fit(background, ~BLIND)
        truth = continuum + signal
        _, null, _, _ = shape.fit(truth, ALL)
        information = float(2 * np.sum(null - truth + truth * np.log(truth / null)))
        pars, _, _, full_grad = shape.fit(background, ALL, signal)
        g = full_grad[-1]
        variance = {key: float(g @ cov[key] @ g) for key in BACKGROUNDS}
        width = math.sqrt(sum(variance.values()))
        center = float(pars[-1])
        calibration = Calibration(shape, theta, continuum, signal, center, width)
        ref = reference["continuum_models"][name]
        check([center, width], [ref["spurious_mu"], ref["sigma_mu"]], f"{name} spurious signal")
        smooth = 25 / information
        calibrated = crossing(lambda l: calibration.q(l) - 25) if width < .2 else None
        if calibrated is None:
            raise ValueError(f"No finite 5-sigma diphoton crossing for {name}")
        check([smooth, calibrated], [ref["smooth_diphoton_fb_inverse"], ref["calibrated_diphoton_fb_inverse"]], name)
        models[name] = dict(spurious_mu=center, sigma_mu=width, variance_by_sample=variance,
                            luminosity_penalty_percent=100 * (calibrated / smooth - 1))
        for treatment, lumi, q in (("smooth", smooth, lambda l: information * l),
                                   ("calibrated", calibrated, calibration.q)):
            append(name, treatment, "diphoton", "both", lumi)
            for scenario, tag in SCENARIOS.items():
                joint = crossing(lambda l: joint_power(q, l, s4, b4[scenario]) - .5)
                published = ref[treatment + "_combined"][tag]["gaussian_score_exact_fourlepton_fb_inverse"]
                check(joint, published["luminosity_fb_inverse"], f"{name} {treatment} {scenario}")
                append(name, treatment, "combined", scenario, joint)
        if name == "exponential2":
            normalization = records["records/ihixs_normalization.json"]
            native_rate = normalization["native_cross_section_pb"]
            ihixs_rate = normalization["cross_section_pb"]
            scale = native_rate / ihixs_rate
            # Express the same physical background bias and covariance in native signal units.
            native = Calibration(shape, theta, continuum, scale * signal, center / scale, width / scale)
            native_lumi = crossing(lambda l: native.q(l) - 25)
            signal_check = dict(ihixs_pb=ihixs_rate, native_hjminnlo_pb=native_rate,
                                branching_fraction=normalization["weight_scale"],
                                native_diphoton_fb_inverse=native_lumi,
                                native_diphoton_design_months=1.2 * native_lumi,
                                native_combined_fb_inverse={k: crossing(lambda l: joint_power(
                                    native.q, l, s4, rate) - .5) for k, rate in b4.items()})
            response = np.eye(len(background)) - jac @ sideband_grad
            residual_cov = response @ background_cov @ response.T
            difference = background - continuum
            center_cov = residual_cov[np.ix_(BLIND, BLIND)]
            rank = int(np.linalg.matrix_rank(center_cov, hermitian=True))
            statistic = float(difference[BLIND] @ np.linalg.pinv(center_cov, hermitian=True) @ difference[BLIND])
            closure = dict(window_chi2=statistic, window_dof=rank, window_p=float(chi2.sf(statistic, rank)),
                           spurious_two_sided_p=float(2 * norm.sf(abs(center / width))))
            check(statistic, fit["fixed_125_GeV_background_only_closure"]["window_120_130"]["chi2"], "closure")
            spectra = dict(low_gev=EDGES[:-1], high_gev=EDGES[1:], **hist,
                           background_fb=background, background_mc_sigma_fb=np.sqrt(background_cov.diagonal()),
                           sideband_fit_fb=continuum, residual_fb=difference,
                           residual_mc_sigma_fb=np.sqrt(residual_cov.diagonal()))
            nominal_lumi = next(r["luminosity_fb_inverse"] for r in rows
                                if (r["model"], r["treatment"], r["channel"], r["heavy_flavour"])
                                == (name, "calibrated", "combined", "central"))
            fixed_exposure = dict(
                combined_90pct_fb_inverse={scenario: crossing(
                    lambda l: joint_power(calibration.q, l, s4, rate) - .9) for scenario, rate in b4.items()},
                combined_median_fb_inverse=nominal_lumi,
                yields_120_130_gev_at_combined_median={
                    "diphoton_signal": nominal_lumi * signal[BLIND].sum(),
                    "diphoton_calibrated_background": nominal_lumi * (continuum + center * signal)[BLIND].sum(),
                    "fourlepton_signal": nominal_lumi * s4,
                    "fourlepton_central_background": nominal_lumi * b4["central"]})
    with (DATA / "gem-table18.csv").open() as file:
        gem = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(file)]
    for row in gem:
        lumi = row["reference_luminosity_fb_inverse"] * (5 / row["diphoton_significance"])**2
        row.update(estimated_five_sigma_fb_inverse=lumi, design_months=1.2 * lumi)
    result = dict(energy_tev=40, higgs_mass_gev=125, design_fb_inverse_per_year=10,
                  target_met=all(m["luminosity_penalty_percent"] <= 20 for m in models.values()),
                  campaign_fit_sha256=reference["fit_sha256"], final_campaign_utc=reference["updated_utc"],
                  prior_packets=reference["precampaign_receipts"], final_packets=reference["extra_receipts"],
                  total_generated_by_stratum=total, final_campaign_generated_by_stratum=added,
                  closure=closure, models=models, fourlepton=four, reach=rows,
                  signal_normalization_check=signal_check, gem_table18=gem, single_experiment=fixed_exposure)
    if result["target_met"] != reference["precision_target"]["achieved"]:
        raise ValueError("Campaign precision decision changed")
    four_spectrum = fourlepton_spectrum(records["records/extended_fourlepton.json"], records["records/reducible_bound.json"])
    return result, spectra, four_spectrum


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/discovery-ssc40"))
    out = parser.parse_args().output
    from plots import plot_fourlepton, plot_reach, plot_spectrum
    from discovery_spectra import plot_discovery_spectra
    from discovery_power import plot_discovery_power
    results, spectra, four_spectrum = analyze()
    fiducial = np.genfromtxt(DATA / "fourlepton-fiducial.csv", delimiter=",", names=True)
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out) as directory:
        temp = Path(directory)
        (temp / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
        write_csv(temp / "reach.csv", results["reach"])
        write_csv(temp / "mass_spectrum.csv", [dict(zip(spectra, row)) for row in zip(*spectra.values())])
        write_csv(temp / "fourlepton_mass_spectrum.csv", [dict(zip(four_spectrum, row)) for row in zip(*four_spectrum.values())])
        plot_fourlepton(four_spectrum, temp)
        plot_spectrum(spectra, temp)
        plot_reach(results["reach"], results["gem_table18"], temp)
        plot_discovery_spectra(spectra, four_spectrum, results["single_experiment"]["combined_median_fb_inverse"],
                               temp, spurious_mu=results["models"]["exponential2"]["spurious_mu"],
                               fiducial_fourlepton=fiducial)
        plot_discovery_power(spectra, results, temp)
        for path in temp.iterdir():
            path.replace(out / path.name)
    for row in results["reach"]:
        if (row["model"] == "exponential2" and row["treatment"] == "calibrated") or row["channel"] == "fourlepton":
            print(f"{row['channel']:12} {row['heavy_flavour']:8} {row['luminosity_fb_inverse']:6.3f} fb^-1  {row['design_months']:6.3f} months")
    print(f"Verified {results['prior_packets'] + results['final_packets']} packets against the final Timur fit. Outputs in {out}")


if __name__ == "__main__":
    main()
