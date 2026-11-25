"""Repeat the one-experiment GEM-range diphoton search with a measured tagged continuum."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys

if __name__ == "__main__":
    for var in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"

import numpy as np

HERE = Path(__file__).resolve().parent
TAG_DIR = HERE.parent
sys.path.insert(0, str(TAG_DIR))
from tagged_scan import ANCHOR, TAG, Inputs, TaggedExperiment, histogram_at_edges, plot, search
from gem_rates import DATA, TAG_MASSES, background_density
from tagged_spectrum import draw_spectrum
from plots import BLUE, ORANGE, save, style

WINDOW = (123.3, 126.7)
GEM_WINDOW_FB = 0.408
PRECISION_FB = 0.0408

PLOT_GROUPS = {
    "charged_current": (r"$\ell\nu\gamma\gamma$", "#9dc8d2"),
    "neutral_current_two_real_photons": (r"$\ell^+\ell^-\gamma\gamma$", "#b8d4c2"),
    "neutral_current_electron_fake": (r"$ee\gamma$, electron fake", "#e5c783"),
    "all_top": ("All top backgrounds", "#c2abc9"),
}


def bin_fractions(edges, lower, upper):
    left = np.maximum(edges[:-1], lower)
    right = np.minimum(edges[1:], upper)
    return np.maximum(right - left, 0.) / np.diff(edges)


def project_background(model, edges, rates):
    if edges[0] > model.edges[0] or edges[-1] < model.edges[-1]:
        raise ValueError("Tagged background must cover every mass window, including the endpoints")
    if (rates < 0).any():
        raise ValueError("Negative background bin in an unweighted LO sample")
    for endpoint in WINDOW:
        if not np.isclose(edges, endpoint, rtol=0, atol=1e-8).any():
            raise ValueError("Source spectrum must have both 125 GeV window boundaries as bin edges")
    cells = histogram_at_edges(edges, rates, model.edges)
    # Interpolation can produce roundoff in otherwise empty bins.
    if np.any(cells < -1.e-12):
        raise ValueError("Negative Poisson mean after rebinning")
    cells = np.maximum(cells, 0.)
    sums = np.r_[0., np.cumsum(cells)]
    windows = sums[model.stop] - sums[model.start]
    if np.any(windows[TAG] <= 0.):
        raise ValueError("An entire tagged mass window has zero simulated continuum")
    np.testing.assert_allclose(windows[ANCHOR], bin_fractions(edges, *WINDOW) @ rates, rtol=1e-9)
    return cells


def spectrum_plot(edges, rate, covariance, sources, groups, output):
    import matplotlib.pyplot as plt

    bins = np.arange(80., 141., 5.)
    transform = np.array([bin_fractions(edges, lo, hi) / (hi-lo) for lo, hi in zip(bins, bins[1:])])
    density = transform @ rate
    sigma = np.sqrt(np.maximum(0., np.diag(transform @ covariance @ transform.T)))
    fig, ax = plt.subplots(figsize=(8.3, 4.7), layout="constrained")
    style(ax)
    bottom = np.zeros(len(bins)-1)
    for key, group in groups.items():
        contribution = sum(np.asarray(sources[name]["rates_fb"]) for name in group["samples"])
        height = bottom + transform @ contribution
        label, color = PLOT_GROUPS[key]
        ax.stairs(height, bins, baseline=bottom, fill=True, color=color, label=label)
        bottom = height
    ax.axvspan(*WINDOW, color="#444444", alpha=.09, label="125 GeV window")
    ax.plot(bins, background_density(bins), lw=2, color=ORANGE, label="GEM Table 19")
    ax.errorbar((bins[:-1]+bins[1:])/2, density, xerr=2.5, yerr=sigma, fmt="o", ms=4,
                color=BLUE, label="Simulated background, MC uncertainty")
    ax.set(xlabel=r"Reconstructed $m_{\gamma\gamma}$ [GeV]", ylabel="Selected cross section [fb/GeV]",
           xlim=(80., 140.), ylim=(0, None), title="Lepton-tagged diphotons at the 40 TeV SSC")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.17), frameon=False, ncol=3, fontsize=8.5)
    save(fig, output, "tagged-background")


def component_comparison(sources):
    gem = DATA["background_components_fb_per_gev"]
    groupings = {
        "charged_current": (["w_lnu_aa"], ["W_gamma_gamma", "W_gamma_lepton_FSR"]),
        "neutral_current_two_real_photons": (["z_ll_aa"], ["Z_gamma_lepton_FSR"]),
        "neutral_current_electron_fake": (["z_ee_a"], ["Z_gamma_electron_fake"]),
        "all_top": ([key for key in sources if key.startswith("tt_")], ["ttbar_gamma_gamma", "ttbar_fakes"]),
    }
    result = {}
    for key, (samples, historical_rows) in groupings.items():
        reference = (WINDOW[1]-WINDOW[0]) * sum(np.interp(125., TAG_MASSES, gem[name]) for name in historical_rows)
        windows = [sources[name]["window_125"] for name in samples]
        result[key] = dict(samples=samples, selected_fb=sum(w["rate_fb"] for w in windows),
                           mc_se_fb=np.linalg.norm([w["error_fb"] for w in windows]), gem_fb=reference,
                           attempted_events=sum(sources[name]["attempted_events"] for name in samples),
                           contributing_window_events=sum(w["contributing_events"] for w in windows))
    return result


def load_spectrum(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return json.load(stream)


def discovery_model():
    signal = json.loads((TAG_DIR / "signal/direct.json").read_text())
    edges = np.asarray(signal["mass_edges_gev"])
    model = Inputs(ggf_veto="pt20_eta2p5", associated=(edges, signal["untagged_associated_spectrum_fb_per_bin"]))
    return signal, model.measured_tag(edges, signal["tagged_spectrum_fb_per_bin"], signal["selected_fb"])


def add_zee_complement(summary, guard):
    """Add independently normalized, disjoint eeγ events to the fit histograms and covariance."""
    if guard["generation_region"] != "outside_nominal" or set(guard["sources"]) != {"z_ee_a"}:
        raise ValueError("The acceptance control must contain only the complementary eeγ phase space")
    np.testing.assert_allclose(guard["edges_gev"], summary["edges_gev"])
    if {s["source_lhe_sha256"] for s in guard["shards"]} & {s["source_lhe_sha256"] for s in summary["shards"]}:
        raise ValueError("An eeγ control integration also appears in the nominal inputs")
    for base, extra in ((summary, guard), (summary["sources"]["z_ee_a"], guard["sources"]["z_ee_a"])):
        for key in ("rates_fb", "gluon_jet_rates_fb", "covariance_fb2", "mc_covariance_fb2",
                    "integration_covariance_fb2"):
            base[key] = np.asarray(base[key])+np.asarray(extra[key])
        window, more = base["window_125"], extra["window_125"]
        for key in ("rate_fb", "contributing_events"):
            window[key] += more[key]
        for key in ("error_fb", "mc_error_fb", "integration_error_fb"):
            window[key] = float(np.hypot(window[key], more[key]))
        window["effective_events"] = (window["rate_fb"]/window["mc_error_fb"])**2
    for key in ("attempted_events", "saved_events", "shards"):
        summary["sources"]["z_ee_a"][key] += guard["sources"]["z_ee_a"][key]


def run(path, output, toys, seed, pilot=False, migration_path=None):
    summary = load_spectrum(path)
    if not summary["complete_sources"]:
        raise ValueError("Incomplete background campaign")
    migration = None
    if migration_path is not None:
        guard = load_spectrum(migration_path)
        restricted_rates = np.asarray(summary["rates_fb"])
        migration = dict(spectrum_sha256=hashlib.sha256(migration_path.read_bytes()).hexdigest(),
                         restricted_window_fb=summary["window_125"]["rate_fb"],
                         restricted_window_mc_se_fb=summary["window_125"]["error_fb"],
                         additional_window_fb=guard["window_125"]["rate_fb"],
                         additional_window_mc_se_fb=guard["window_125"]["error_fb"])
        add_zee_complement(summary, guard)
    edges = np.asarray(summary["edges_gev"], dtype=float)
    rates = np.asarray(summary["rates_fb"], dtype=float)
    cov = np.asarray(summary["covariance_fb2"], dtype=float)
    if (cov.shape != (len(rates), len(rates)) or not np.isfinite(cov).all()
            or len(edges) != len(rates)+1 or not np.all(np.diff(edges) > 0)):
        raise ValueError("Malformed tagged-background spectrum or covariance")
    np.testing.assert_allclose(cov, cov.T, atol=1.e-12)
    if np.linalg.eigvalsh(cov)[0] < -1.e-10*max(1., cov.diagonal().max()):
        raise ValueError("Indefinite tagged-background covariance")
    fraction = bin_fractions(edges, *WINDOW)
    measured = float(fraction @ rates)
    error = float(np.sqrt(max(0., fraction @ cov @ fraction)))
    unseen = summary["window_125"].get("zero_strata_upper_fb_95", 0.)
    np.testing.assert_allclose([measured, error],
                               [summary["window_125"]["rate_fb"], summary["window_125"]["error_fb"]])
    if error <= 0. or measured <= 0. or (error + unseen > PRECISION_FB and not pilot):
        raise ValueError(f"Window {measured:.5g} ± {error:.5g} fb, upper limit for empty sources {unseen:.5g} fb")
    signal, model = discovery_model()
    se = np.asarray(signal["mass_edges_gev"])
    model.tag_background = project_background(model, edges, rates)
    # Each variant uses its own coherent background for both null calibration
    # and signal-injected experiments. This measures sensitivity to MC normalization.
    results = {}
    scales = {"nominal": 1., "window_mc_high": 1.+error/measured}
    if measured > error:
        scales["window_mc_low"] = 1.-error/measured
    if unseen:
        scales["empty_sources_upper_with_total_shape"] = 1.+unseen/measured
    for name, scale in scales.items():
        results[name] = search(model, signal["selected_fb"], scale, toys, seed, max_iterations=12)
    # GEM added 20% to its quark-initiated Zγ estimate to account for
    # gluon fusion. Applying that historical prescription to the modern Z
    # samples is a sensitivity test, independent of their MC precision.
    z_rates = sum(np.asarray(summary["sources"][key]["rates_fb"]) for key in ("z_ll_aa", "z_ee_a"))
    nominal_cells = model.tag_background
    model.tag_background = project_background(model, edges, rates + .20*z_rates)
    results["gem_z_enhancement"] = search(model, signal["selected_fb"], 1., toys, seed, max_iterations=12)
    model.tag_background = project_background(model, edges, np.asarray(summary["gluon_jet_rates_fb"]))
    results["unmatched_jets_gluon_response"] = search(model, signal["selected_fb"], 1., toys, seed, max_iterations=12)
    if migration_path is not None:
        model.tag_background = project_background(model, edges, restricted_rates)
        results["restricted_zee_acceptance"] = search(model, signal["selected_fb"], 1., toys, seed, max_iterations=12)
        # Change only the electron-fake normalization to GEM's Table 19
        # value at 125 GeV. Other sources and the simulated shape stay fixed.
        electron_rates = np.asarray(summary["sources"]["z_ee_a"]["rates_fb"])
        historical = DATA["background_components_fb_per_gev"]["Z_gamma_electron_fake"]
        gem_electron = (WINDOW[1]-WINDOW[0]) * np.interp(125., TAG_MASSES, historical)
        ratio = gem_electron / (fraction @ electron_rates)
        model.tag_background = project_background(model, edges, rates + (ratio-1.)*electron_rates)
        results["gem_electron_fake_normalization"] = search(
            model, signal["selected_fb"], 1., toys, seed, max_iterations=12)
    model.tag_background = nominal_cells
    median = results["nominal"]["crossings"][0]
    exp = TaggedExperiment(model, median["luminosity_fb_inverse"], signal["selected_fb"], 1.)
    alpha = .002
    ordinary, ordinary_se = exp.ordinary_null(alpha, toys, seed+402)
    factor, factor_se = exp.null_trials(alpha, toys, seed+403)
    if abs(ordinary-alpha*factor) > 6.*np.hypot(ordinary_se, alpha*factor_se):
        raise ValueError("Direct null toys disagree with the importance-sampled global tail")
    history = json.loads((TAG_DIR / "combined.json").read_text())
    groups = component_comparison(summary["sources"])
    np.testing.assert_allclose(sum(g["selected_fb"] for g in groups.values()), measured)
    report = dict(toys=toys, seed=seed, pilot=pilot, target_mc_se_fb=PRECISION_FB,
                  measured_tagged_background_fb=measured,
                  measured_tagged_background_mc_se_fb=error, gem_tagged_background_fb=GEM_WINDOW_FB,
                  empty_sources_upper_fb_95=unseen,
                  measured_over_gem=measured/GEM_WINDOW_FB, window_gev=WINDOW,
                  source_summary_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  historical_scan_sha256=hashlib.sha256((TAG_DIR / "combined.json").read_bytes()).hexdigest(),
                  design_fb_inverse_per_year=10., inclusive_mass_range_gev=[80,160], tagged_mass_range_gev=[80,140],
                  mass_step_gev=.5, true_mass_gev=125., scenarios=results, window_components=groups,
                  wider_zee_acceptance=migration,
                  comparison={name:history["scenarios"][name]["crossings"]
                              for name in ("nominal", "same_production_inclusive_only")},
                  global_tail_control=dict(local_alpha=alpha, direct_probability=ordinary,
                      direct_mc_se=ordinary_se, importance_probability=alpha*factor, importance_mc_se=alpha*factor_se))
    output.mkdir(parents=True, exist_ok=True)
    plot(model, results["nominal"], output)
    spectrum_plot(edges, rates, cov, summary["sources"], groups, output)
    event_edges = np.arange(80., 141.)
    signal_bins = histogram_at_edges(se, np.asarray(signal["tagged_spectrum_fb_per_bin"]), event_edges)
    names = list(groups)
    component_bins = []
    for name in names:
        spectrum = sum(np.asarray(summary["sources"][k]["rates_fb"]) for k in groups[name]["samples"])
        component_bins.append(histogram_at_edges(edges, spectrum, event_edges))
    labels, colors = zip(*(PLOT_GROUPS[name] for name in names))
    draw_spectrum(results["nominal"], event_edges, signal_bins, component_bins, names, colors, labels, output)
    (output / "combined-mc.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    for row in results["nominal"]["crossings"]:
        print(f"{row['target_power']:.0%} global discovery power at {row['luminosity_fb_inverse']:.3f} fb^-1 "
              f"({row['design_months']:.3f} design months)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--toys", type=int, default=40000)
    parser.add_argument("--seed", type=int, default=725081)
    parser.add_argument("--pilot", action="store_true",
                        help="Refit an interim spectrum and record its preliminary status")
    parser.add_argument("--zee-complement", type=Path,
                        help="Include disjoint eeγ acceptance migrations in the nominal fit and compare restricted cuts")
    args = parser.parse_args()
    if args.toys < 4000:
        parser.error("Use at least 4,000 experiments")
    if args.pilot and args.output.resolve() == HERE:
        parser.error("Choose a separate --output directory for an interim fit")
    run(args.spectrum, args.output, args.toys, args.seed, args.pilot, args.zee_complement)
