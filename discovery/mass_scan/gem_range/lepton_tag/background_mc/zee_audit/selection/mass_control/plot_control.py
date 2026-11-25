"""Plot how a low-mass calibration changes the eeγ predictions at higher masses."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent


def plot(report, out):
    pure = report["unit_scale_fit"]
    f_only = 100 * pure["leakage"]["value"]
    xmax = np.ceil(f_only + 0.5)
    scan = [s for s in report["fixed_leakage_scan"] if 100 * s["leakage"] <= xmax]
    x = [100 * s["leakage"] for s in scan]
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.4, 4.2), layout="constrained")
    for axis in (left, right):
        axis.spines[["top", "right"]].set_visible(False)
        axis.axhline(1, color="#555555", ls="--", lw=1)
        axis.set(xlim=(0, xmax), xlabel="Primary transverse energy added to isolation [%]")
        axis.grid(axis="y", alpha=.18)

    def curve(axis, key, color, label):
        y = np.array([s[key]["value"] for s in scan])
        err = np.array([s[key]["seed_mc_se"] for s in scan])
        axis.plot(x, y, color=color, label=label, lw=2)
        axis.fill_between(x, y - err, y + err, color=color, alpha=.14)
        axis.plot(f_only, pure[key]["value"], "o", color=color, mec="white", zorder=4)

    curve(left, "scale", "#26637a", None)
    left.set(ylabel="Required rate multiplier", title="Calibration in 80–100 GeV")
    curve(right, "to_gem_125", "#933f26", "123.3–126.7 GeV")
    curve(right, "to_gem_120_140", "#26637a", "120–140 GeV")
    right.set(ylabel="Predicted rate / GEM rate", title="Prediction at higher mass")
    right.legend(loc="best", frameon=False, title="Diphoton mass")
    fig.suptitle("Electron fakes in the lepton-tagged diphoton search", fontsize=14)
    for ext in ["png", "pdf"]:
        fig.savefig(out / f"mass-control.{ext}", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=HERE / "control_fit.json")
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    plot(json.loads(args.input.read_text()), args.output)
