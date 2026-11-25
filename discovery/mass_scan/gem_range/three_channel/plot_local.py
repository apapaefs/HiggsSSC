"""Plot the individual and joint median local p-values at the global discovery exposure."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtr
from scipy.stats import norm, poisson

from scan_three import ALPHA, FOUR, MASSES, JointScan, load_model, poisson_mixture, threshold
from plots import save, style

HERE = Path(__file__).resolve().parent


def plot(report, output):
    row = next(r for r in report["crossings"] if r["target_power"] == .5)
    lumi = row["luminosity_fb_inverse"]
    exp = JointScan(load_model(report["tag_background"]), lumi)
    alt_means = np.column_stack([lumi*c.window_sums(c.background + c.signal) for c in exp.channels])
    null_means = np.column_stack([lumi*c.window_sums(c.background) for c in exp.channels])
    weights = np.column_stack([c.weights for c in exp.channels])
    drift = exp.strength * exp.response
    combined = []
    for i, (pmf, offsets, _) in enumerate(exp.mixtures):
        alt_pmf, alt_offsets, _ = poisson_mixture(alt_means[i], weights[i])
        # Translation by the mean of the injected Gaussian commutes with the median.
        median_score = threshold(exp.strength, alt_pmf, alt_offsets, .5) + exp.strength * drift[i]
        combined.append(pmf @ ndtr((offsets - median_score) / exp.strength))
    fig, ax = plt.subplots(figsize=(8.9, 5.1), layout="constrained")
    style(ax)
    ax.semilogy(MASSES, combined, color="#213d50", lw=2.3, label="Three-channel search")
    ax.semilogy(MASSES, norm.sf(drift), color="#347c9f", ls="--", label="Untagged diphotons")
    for col, mask, color, label in [(0, MASSES <= 140, "#ae622f", "Lepton-tagged diphotons"),
                                    (1, FOUR, "#6d4b80", "Four leptons")]:
        median_count = poisson.ppf(.5, alt_means[mask, col])
        ax.semilogy(MASSES[mask], poisson.sf(median_count-1, null_means[mask, col]),
                    color=color, lw=1.4, label=label)
    local_alpha = norm.sf(row["local_z"])
    for p, label in [(ALPHA, "Local 5σ"), (local_alpha, "Global 5σ threshold")]:
        ax.axhline(p, color="#62686c", ls=":", lw=1)
        ax.text(159.6, p*1.32, label, ha="right", va="bottom", fontsize=9)
    ax.set(xlabel="Tested Higgs mass [GeV]", ylabel="Median local background p-value",
           ylim=(min(local_alpha / 3, min(combined)/2), 2.3), xlim=(80., 160.),
           title=rf"One 40 TeV SSC experiment, $m_H=125$ GeV, {lumi:.2f} fb$^{{-1}}$")
    ax.legend(frameon=False, fontsize=9, loc="center left")
    save(fig, output, "local-p")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=HERE / "nominal.json")
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    plot(json.loads(args.input.read_text()), args.output)
