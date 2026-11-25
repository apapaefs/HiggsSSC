"""Plot the selected eeγ rates and the dependence on visible shower leakage."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(input_file, output):
    report = json.loads(input_file.read_text())
    results = report["results"]
    ref = report["gem_fb"]
    dark, blue, rust = "#13485b", "#236d87", "#99441f"
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.2, 4.1),
                                      gridspec_kw={"width_ratios": [1.35, 1.]},
                                      layout="constrained")
    names = ["nominal", "five_by_five", "cell_threshold", "noise",
             "compactness_hcal", "with_tracks"]
    labels = ["Previous selection", "5×5 cluster subtraction", "Cell threshold",
              "Cone noise", "Cluster shape + HCAL", "Additional tracks"]
    for y, name in enumerate(names):
        obs = results[name]
        left.errorbar(obs["rate_fb"], y, xerr=obs["mc_se_fb"], fmt="o", capsize=3, color=blue)
    left.set_yticks(range(len(names)), labels)
    left.set_ylim(len(names) - .5, -.5)
    left.axvline(ref, label="GEM Table 19", color=rust, linestyle="--", linewidth=1.4)
    left.set_xlabel("Selected background [fb]")
    left.set_title("Successive selection changes", loc="left", fontsize=11)
    left.legend(loc="lower left", bbox_to_anchor=(0,1.1), frameon=False, fontsize=9)
    left.grid(axis="x", alpha=.2)
    left.set_xlim(0., max(results[k]["rate_fb"] + results[k]["mc_se_fb"] for k in names) * 1.15)

    x = [0., .5, 1., 2., 3., 4., 5., 6., 7., 8.]
    keys = ["with_tracks"] + [f"leak_{f:g}pct" for f in x[1:]]
    y = np.array([results[k]["rate_fb"] for k in keys])
    se = np.array([results[k]["mc_se_fb"] for k in keys])
    right.axhline(ref, color=rust, linestyle="--", linewidth=1.4)
    right.fill_between(x, y-se, y+se, color=blue, alpha=.15)
    right.plot(x[:7], y[:7], "o-", color=dark, markersize=4)
    right.plot(x[6:], y[6:], "o:", color=dark, markersize=4)
    right.set(xlabel="Primary transverse energy counted\nin isolation outside the cluster [%]",
              ylabel="Selected background [fb]", ylim=(0., max(y+se)*1.12), xticks=range(0,9,2))
    right.set_title("Sensitivity to shower leakage", loc="left", fontsize=11)
    right.grid(axis="y", alpha=.2)
    fig.suptitle(r"$ee\gamma$ in the tagged 125 GeV Higgs window at 40 TeV",
                 fontsize=15)
    for axis in (left, right):
        axis.spines[["right", "top"]].set_visible(False)
        axis.set_axisbelow(True)
    for ext in ["pdf","png"]:
        fig.savefig(output/f"gem-selection.{ext}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plot(args.input, args.output)
