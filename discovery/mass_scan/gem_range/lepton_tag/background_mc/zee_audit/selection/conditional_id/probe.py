"""Test the fake electron's calorimeter selection after photon tracking in the existing samples."""
import argparse
import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
SELECTION = HERE.parent
sys.path.insert(0, str(SELECTION))
from mass_control.control_fit import jackknife
from paired_fit import evaluate, load
from selection_rates import columns, model


def herwig_probe():
    packet = load(SELECTION / "herwig_candidates.json.gz", SELECTION / "paired_weights.json.gz")
    rows = packet["assignments"]
    selected, _ = model(rows)
    primary = columns(rows, "primary_et")
    inner = columns(rows, "core3_em_et")
    outer = columns(rows, "core_em_et")
    hadronic = columns(rows, "core_had_et")
    cluster = ((.95 * primary + inner) / (primary + outer) > .90) & (hadronic < .10 * (primary + outer))
    # The saved object probability already includes the clean-shower efficiency
    # for each object. Only the fake electron's 90% factor is divided out.
    others = cluster[:, 0] & cluster[:, 2] & (columns(rows, "extra_track_pt") <= 1.).all(axis=1)
    denominator = selected["noise"] * others / .90
    local_shape = denominator * cluster[:, 1]
    observed = .90 * local_shape
    np.testing.assert_allclose(observed, selected["with_tracks"])
    fit = evaluate(packet, {"before_fake_calo": denominator, "surrounding_deposits": local_shape,
                            "including_clean_shower_efficiency": observed}, baseline="before_fake_calo")
    keys = ["rate_fb", "mc_se_fb", "relative_to_baseline", "ratio_paired_se"]
    return {key: {field: obs[field] for field in keys} for key, obs in fit["results"].items()}


def effective_efficiency():
    mass = SELECTION / "mass_control"
    with gzip.open(mass / "native_scan.json.gz", "rt") as stream:
        native = json.load(stream)
    reference = json.loads((mass / "control_fit.json").read_text())["gem_fb_except_low_ratio"]
    assert native["leakage_fractions"][0] == 0.
    rate = np.array([run["rate_fb"] for run in native["runs"]])[:, :, 0]
    n = len(rate)
    samples = np.concatenate([rate.mean(axis=0)[None], (rate.sum(axis=0) - rate) / (n - 1)])
    index = native["window_order"].index
    low = .90 * reference["training"] / (samples[:, index("80_90")] + samples[:, index("90_100")])
    high = .90 * reference["higgs"] / samples[:, index("higgs")]
    return {"80_100": jackknife(low), "higgs_window": jackknife(high),
            "difference_high_minus_low": jackknife(high - low)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "probe.json")
    args = parser.parse_args()
    out = {"herwig_fake_calo_given_other_cuts": herwig_probe(),
           "native_effective_clean_shower_efficiency_at_zero_leakage": effective_efficiency(),
           "assumptions": "The Herwig denominator includes the 3% tracking weight and every other event cut. "
                          "The native inversion assigns the full residual to the fake electron's 90% clean-shower "
                          "probability. The uncertainties include only simulation statistics with GEM fixed."}
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
