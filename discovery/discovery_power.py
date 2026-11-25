"""Probability that one SSC experiment passes the local 5σ test at a given exposure."""
import csv

import numpy as np
from scipy.special import ndtr

from input_data import SIGNAL
from likelihood import BLIND, Calibration, Shape, joint_power
from plots import BLUE, INK, ORANGE, save
import matplotlib.pyplot as plt


def plot_discovery_power(spectra, results, out):
    nominal = results["models"]["exponential2"]
    shape = Shape("exponential", 2)
    theta, continuum, _, _ = shape.fit(spectra["background_fb"], ~BLIND)
    np.testing.assert_allclose(continuum, spectra["sideband_fit_fb"])
    calibrated = Calibration(shape, theta, continuum, spectra[SIGNAL],
                             nominal["spurious_mu"], nominal["sigma_mu"])
    rows = [r for r in results["reach"]
            if r["model"] == "exponential2" and r["treatment"] == "calibrated"]
    medians = {(r["channel"], r["heavy_flavour"]): r["luminosity_fb_inverse"] for r in rows}
    # Evaluate the profile and exact four-lepton sum at every point, including the crossings.
    lum = np.unique(np.r_[np.linspace(.05, 6., 120), list(medians.values())])
    q = np.array([calibrated.q(x) for x in lum])
    curves = {"diphoton": ndtr(np.sqrt(q) - 5.)}
    for tag in ("central", "upper"):
        four = results["fourlepton"][tag]
        curves[tag] = np.array([joint_power(calibrated.q, x, four["signal_fb"], four["background_fb"])
                                for x in lum])
    settings = [
        ("central", "combined", "central", BLUE, "-", r"$\gamma\gamma+4\ell$, central background"),
        ("upper", "combined", "upper", ORANGE, "--", r"$\gamma\gamma+4\ell$, heavy-flavour upper bound"),
        ("diphoton", "diphoton", "both", "#545c67", "-.", r"$\gamma\gamma$ alone"),
    ]
    fig, ax = plt.subplots(figsize=(7.4, 4.8), layout="constrained")
    ax.axhline(.5, ls=":", lw=1, color="#84939d")
    for key, channel, scenario, color, line, label in settings:
        median = medians[channel, scenario]
        idx = np.searchsorted(lum, median)
        if not np.isclose(curves[key][idx], .5, atol=1.e-7, rtol=0):
            raise ValueError(f"Discovery curve misses the published median for {key}")
        ax.plot(lum, curves[key], lw=2.2, color=color, ls=line,
                label=label + f"  ({median:.2f} " + r"fb$^{-1}$)")
        ax.scatter(median, .5, s=35, color=color, zorder=3)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.text(.035, .965, "SSC projection", transform=ax.transAxes, ha="left", va="top",
            color=INK, fontweight="bold", fontsize=13)
    ax.text(.035, .887, r"$pp$, 40 TeV   ·   $m_H=125$ GeV", transform=ax.transAxes,
            ha="left", va="top", color=INK, fontsize=10)
    ax.set(xlim=(0, 6), ylim=(0, 1.01), xlabel=r"Integrated luminosity [fb$^{-1}$]",
           ylabel=r"Probability of a local $5\sigma$ discovery")
    ax.set_yticks(np.arange(5) / 4, labels=["0%", "25%", "50%", "75%", "100%"])
    ax.tick_params(direction="in", right=True)
    design_months = 12 / results["design_fb_inverse_per_year"]
    top = ax.secondary_xaxis("top", functions=(lambda x: x * design_months,
                                               lambda x: x / design_months))
    rate = results["design_fb_inverse_per_year"]
    top.set_xlabel(rf"Months at {rate:g} fb$^{{-1}}$ per year")
    top.tick_params(direction="in")
    save(fig, out, "discovery-probability")

    with (out / "discovery_probability.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["luminosity_fb_inverse", "design_months", "combined_central_probability",
                         "combined_upper_probability", "diphoton_probability"])
        writer.writerows(zip(lum, lum * design_months, curves["central"], curves["upper"],
                             curves["diphoton"]))
