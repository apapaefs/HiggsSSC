"""Measure isolation and cluster-selection changes on the same tagged eeγ collisions."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from noise import pass_probability
from paired_fit import combine, evaluate, load


HERE = Path(__file__).resolve().parent
ROLES = ("real", "fake", "tag")
RADII = np.array([.45, .45, .30])
GEM_FB = .0374


def columns(rows, name, grid=None):
    return np.asarray([[(row[role][grid] if grid else row[role])[name]
                        for role in ROLES] for row in rows], dtype=float)


def model(rows, grid=None):
    """Kinematics and the object ID factors are frozen across the variations."""
    base = np.array([r["probability"] for r in rows])
    primary = columns(rows, "primary_et", grid)
    activity = columns(rows, "cell_iso_et", grid)
    ring = columns(rows, "core_em_et", grid)
    center = columns(rows, "core3_em_et", grid)
    had = columns(rows, "core_had_et", grid)
    tracks = columns(rows, "extra_track_pt", grid)
    track_pass = (tracks <= 1.).all(axis=1)

    def shape(share=.95):
        # An isolated primary shower already carries the 90% plateau. This
        # mask measures the extent to which other particles spoil an otherwise
        # passing cluster. The clean E3/E5 share is varied independently.
        return ((share * primary + center) / (primary + ring) > .90).all(axis=1)

    def isolation(leak=0., power=1., field=activity, calibrated=True):
        return pass_probability(field + leak * primary, RADII, scale_power=power,
                                calibrated=calibrated)

    had_pass = (had < .10 * (primary + ring)).all(axis=1)
    cluster_pass = shape() & had_pass
    selected = cluster_pass & track_pass
    passed = isolation()
    iso = base * passed.prod(axis=1)
    with_cluster = iso * cluster_pass
    variants = {
        "nominal": np.array([r["nominal_probability"] for r in rows]),
        "five_by_five": base * (columns(rows, "core_iso", grid) < 5.).all(axis=1),
        "cell_threshold": base * (activity < 5.).all(axis=1),
        "noise": iso,
        "compactness": iso * shape(),
        "compactness_hcal": with_cluster,
        "with_tracks": with_cluster * track_pass,
        "compactness_pointlike": iso * shape(1.) * had_pass * track_pass,
        "compactness_near_cut": iso * shape(.905) * had_pass * track_pass,
        "all_neighbour_em_above_1gev_veto": iso * selected * (columns(rows, "extra_em_pt", grid) <= 1.).all(axis=1),
        "hcal_tag_1gev": iso * selected * (columns(rows, "core_had_e", grid)[:, 2] < 1.),
        "tag_track_central_three": with_cluster * (tracks[:, :2] <= 1.).all(axis=1)
                                     * (columns(rows, "extra_track3_pt", grid)[:, 2] <= 1.),
        "tag_noise_same_as_photon": base * selected * isolation(power=0.).prod(axis=1),
        "tag_noise_pileup_scaling": base * selected * isolation(power=1.88).prod(axis=1),
        "noise_uncalibrated": base * selected * isolation(calibrated=False).prod(axis=1),
        "noise_joint_lower": base * selected * np.maximum(0., passed.sum(axis=1) - 2.),
        "noise_joint_upper": base * selected * passed.min(axis=1),
    }

    def with_leak(fraction):
        # Effective additional transverse energy *counted* by the isolation
        # after thresholding. Forward/back loss and subthreshold cells do not
        # automatically count. The scan includes all three electromagnetic objects.
        return base * selected * isolation(leak=fraction).prod(axis=1)

    for percent in (0.5, 1., 2., 3., 4., 5., 6., 7., 8.):
        variants[f"leak_{percent:g}pct"] = with_leak(percent / 100.)
    for field, label in (("mean92_em_angle_rad", "gem92_angular_cut"),
                         ("mean_em_angle_rad", "angular_cut_in_1993_cluster")):
        # The precursor applied this multiphoton discriminator to the two
        # photon candidates. The electron promoted to photon counts as an EM shower.
        if not grid:
            angles = columns(rows, field)
            angular_pass = (angles[:, :2] < .005).all(axis=1)
            variants[label] = iso * selected * angular_pass
            variants[label + "_after_isolation"] = iso * angular_pass
    return variants, with_leak


def run(features, hard):
    packet = load(features, hard)
    rows = packet["assignments"]
    variants, with_leak = model(rows)
    baseline = variants["nominal"]
    cone_check = (columns(rows, "cone") < 5.).all(axis=1) * [r["probability"] for r in rows]
    np.testing.assert_allclose(baseline, cone_check, atol=1.e-12)
    uniform, _ = model(rows, "uniform")
    variants["uniform_barrel_grid"] = uniform["with_tracks"]
    report = evaluate(packet, variants, reference_fb=GEM_FB)
    results = report["results"]
    np.testing.assert_allclose(results["nominal"]["rate_fb"], .13406788141483655, rtol=1.e-9)
    np.testing.assert_allclose(results["nominal"]["modern_hard_rate_fb"], .25745677493885, rtol=1.e-9)
    np.testing.assert_allclose(results["nominal"]["mc_se_fb"], .01579847293359152, rtol=1.e-9)
    obs = report["observable_order"]
    covariance = np.array(report["covariance_fb2"])

    def paired_change(first, last):
        i, j = obs.index(first), obs.index(last)
        return {"rate_change_fb": results[last]["rate_fb"] - results[first]["rate_fb"],
                "paired_mc_se_fb": np.sqrt(max(0., covariance[j,j] + covariance[i,i] - 2*covariance[i,j]))}

    report["changes"] = {
        "noise_after_cell_threshold": paired_change("cell_threshold", "noise"),
        "concentration_after_noise": paired_change("noise", "compactness"),
        "hadronic_after_concentration": paired_change("compactness", "compactness_hcal"),
        "extra_track_after_cluster": paired_change("compactness_hcal", "with_tracks"),
        "shape_and_hadronic_after_noise": paired_change("noise", "compactness_hcal"),
        "precursor_angular_after_tracks": paired_change("with_tracks", "gem92_angular_cut"),
    }
    central = results["with_tracks"]
    eff = .90 * GEM_FB / central["rate_fb"]
    report["required_fake_conditional_efficiency"] = {
        "value": eff, "mc_se": eff * central["mc_se_fb"] / central["rate_fb"]}

    def leak_rate(leak):
        prob = np.bincount(packet["assignment_events"], weights=with_leak(leak),
                           minlength=len(packet["events"]))
        rates, cov = combine(prob[:,None], packet, True)
        return rates[0], np.sqrt(cov[0,0])

    crossing = brentq(lambda leak: leak_rate(leak)[0] - GEM_FB, 0., .5)
    step = 1.e-5
    slope = (leak_rate(crossing + step)[0] - leak_rate(crossing - step)[0]) / (2.*step)
    report["required_visible_leakage_fraction_in_scan"] = crossing
    report["required_visible_leakage_fraction_mc_se"] = leak_rate(crossing)[1] / abs(slope)
    report["noise_only_pass_R045_at_3_4_5"] = pass_probability(0., .45, margin=np.array([3.,4.,5.])).tolist()
    report["noise_only_pass_R030_at_5_assumed"] = float(pass_probability(0., .30))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("features", type=Path)
    parser.add_argument("weights", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = run(args.features, args.weights)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
