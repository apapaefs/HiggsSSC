"""Plot the exclusive tagged category at the combined median discovery exposure."""
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import ndtr
from scipy.stats import poisson

from gem_rates import DATA, TAG_MASSES, background_density, resolution

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from plots import INK, ORANGE, save, style
import matplotlib.pyplot as plt


def run():
    report = json.loads((HERE / "combined.json").read_text())
    result = report["scenarios"]["nominal"]
    lumi = result["crossings"][0]["luminosity_fb_inverse"]
    edges = np.arange(80., 141.)
    centers = .5 * (edges[1:] + edges[:-1])
    # Every kink in the linearly interpolated rates is at an integer bin edge.
    # Midpoints integrate each line exactly; reconcile the rounded component
    # entries to the independently tabulated total in each exported bin.
    density = background_density(centers)
    rows = DATA["background_components_fb_per_gev"]
    components = np.array([np.interp(centers, TAG_MASSES, row) for row in rows.values()])
    components *= density / components.sum(axis=0)
    if report["signal_model"] == "direct_exclusive_ggf_wh_tth":
        direct = json.loads((HERE / "signal/direct.json").read_text())
        cdf = np.r_[0., np.cumsum(direct["tagged_spectrum_fb_per_bin"])]
        signal = np.diff(np.interp(edges, direct["mass_edges_gev"], cdf))
    else:
        signal = result["accepted_tagged_signal_fb"] * np.diff(ndtr((edges - 125.)/resolution(125.)))
    with (HERE / "tagged-spectrum.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["low_gev", "high_gev", "signal_fb", *rows, "background_fb"])
        writer.writerows(zip(edges[:-1], edges[1:], signal, *components, density))

    colors = ("#739ab7", "#a1c9c3", "#d4a466", "#a77895", "#9aab76", "#70808a")
    labels = [r"$t\bar t\gamma\gamma$", r"$W\gamma\gamma$", r"$W\gamma$ + lepton radiation",
              r"$Z\gamma$ + lepton radiation", r"$Z\gamma$ + electron fake", r"$t\bar t$ fakes"]
    fig, (spec, count) = plt.subplots(1, 2, figsize=(11, 4.4), width_ratios=(1.55, 1),
                                      layout="constrained")
    for ax in (spec, count):
        style(ax)
    bottom = np.zeros(len(centers))
    for contribution, color, label in zip(components, colors, labels):
        upper = bottom + lumi * contribution
        spec.stairs(upper, edges, baseline=bottom, color=color, lw=0, fill=True, label=label)
        bottom = upper
    spec.stairs(bottom + lumi * signal, edges, color=ORANGE, lw=2, label="Signal + background")
    spec.set(xlabel=r"Reconstructed $m_{\gamma\gamma}$ [GeV]", ylabel="Expected events / GeV",
             xlim=(80, 140), ylim=(0, 1.48 * max(bottom + lumi * signal)))
    spec.legend(loc="upper left", frameon=False, fontsize=8, ncol=2)

    b = lumi * result["background_in_true_window_fb"]
    s = lumi * result["signal_in_true_window_fb"]
    n = np.arange(int(poisson.ppf(.998, b+s)) + 1)
    count.bar(n-.18, poisson.pmf(n, b), width=.36, color=INK,
              label=f"Background ({b:.2f} expected)")
    count.bar(n+.18, poisson.pmf(n, s+b), width=.36, color=ORANGE,
              label=f"Signal + background ({s+b:.2f})")
    count.set(xlabel="Tagged events in 123.3–126.7 GeV", ylabel="Probability",
              xticks=n[::2], ylim=(0, .42))
    count.legend(loc="upper right", frameon=False, fontsize=8)
    fig.suptitle(rf"Lepton-tagged diphotons at the 40 TeV SSC, {lumi:.2f} fb$^{{-1}}$",
                 x=.015, ha="left", fontweight="bold", color=INK)
    save(fig, HERE, "tagged-spectrum")


if __name__ == "__main__":
    run()
