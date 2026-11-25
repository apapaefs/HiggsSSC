"""Measure the separate Zee migration sample at the prespecified confirmation exposure."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import tempfile

import numpy as np
from scipy.stats import beta

from response_summary import WINDOW, read_campaign, summarize


def measure(groups, max_response, alpha):
    count = sum(m.events for _, _, m in groups)
    selected = sum(m.window_count for _, _, m in groups)
    weighted_sum = 0.
    variance = 0.
    integration_variance = 0.
    for sigma, error, m in groups:
        if m.events < 2 or m.events != m.sumw or m.events != m.sumw2:
            raise ValueError("The prespecified binomial limit requires independent positive unit-weight events")
        if m.max_window > max_response + 1.e-9:
            raise ValueError("Measured event probability exceeds the two-assignment Zee response bound")
        first = float(m.s1[WINDOW].sum())
        second = float(m.s2[np.ix_(WINDOW, WINDOW)].sum())
        weighted_sum += sigma * first
        variance += sigma**2 * max(0., second - first**2 / m.events) * m.events / (m.events - 1.)
        integration_variance += (error * first)**2
    probability_upper = float(beta.ppf(1.-alpha, selected+1, count-selected)) if selected<count else 1.
    raw_upper = max(sigma for sigma, _, _ in groups) * probability_upper
    return dict(attempted=count, selected=selected, rate_fb=weighted_sum/count,
                mc_error_fb=math.sqrt(variance)/count,
                integration_error_fb=math.sqrt(integration_variance)/count,
                unweighted_upper_fb_95=raw_upper, detector_upper_fb_95=raw_upper*max_response)


def result(root, plan):
    groups, _ = read_campaign(root, {"z_ee_a"}, False, "outside_nominal", False)
    rounds = {}
    for hard in sorted(root.glob("z_ee_a-edge-*/hard.json")):
        data = json.loads(hard.read_text())
        group = groups.get(("z_ee_a", data["integration_id"]))
        if group and data["generation_region"] == "outside_nominal":
            rounds[data["hard_seed"]] = group
    ids = plan["confirmation"]["hard_seeds"]
    pilot_id = plan["pilot"]["hard_seed"]
    counts = {k:v[2].events for k,v in rounds.items()}
    if any(counts.get(k)!=plan["confirmation"]["events_per_round"] for k in ids):
        raise ValueError(f"Confirmation incomplete: attempted events by hard seed {counts}")
    cap = plan["maximum_response_per_collision"]
    alpha = plan["alpha_one_sided"]
    confirm = [rounds[k] for k in ids]
    output = dict(primary_confirmation=measure(confirm, cap, alpha),
                  max_response=cap, alpha=alpha, hard_seeds=ids,
                  hard_cross_sections_pb=[s/1000. for s,_,_ in confirm],
                  n_confirmation_shards=sum(g[2].shards for g in confirm))
    if pilot_id in rounds:
        if counts[pilot_id] != plan["pilot"]["events"]:
            raise ValueError("The accompanying pilot sample is incomplete")
        pilot = rounds[pilot_id]
        output.update(pilot_events=pilot[2].events, pilot_selected=pilot[2].window_count,
                      descriptive_full_sample={k:v for k,v in measure([pilot]+confirm,cap,alpha).items()
                                               if "upper" not in k})
    return output


def confirmation_spectrum(root, plan, fixed):
    seeds = set(plan["confirmation"]["hard_seeds"])
    root = root.resolve()
    # Keep exactly the same relative receipt/packet paths as the physical archive.
    # The standard reader validates hashes, source intervals, normalization and covariance.
    with tempfile.TemporaryDirectory(prefix="ssc-zee-confirmation-") as name:
        view = Path(name)
        for receipt_path in sorted((root / "showers").glob("*/receipt.json")):
            receipt = json.loads(receipt_path.read_text())
            if receipt["sample"] != "z_ee_a" or receipt["hard_seed"] not in seeds:
                continue
            if receipt["generation_region"] != "outside_nominal":
                raise ValueError("invalid confirmation region")
            packet = receipt_path.with_name("packet.csv")
            folder = view / receipt_path.parent.relative_to(root)
            folder.mkdir(parents=True)
            for path in (receipt_path, packet):
                (folder / path.name).symlink_to(path)
        spectrum = summarize(view, region="outside_nominal")
    component = spectrum["sources"]["z_ee_a"]
    if component["attempted_events"] != plan["confirmation"]["events"]:
        raise ValueError("binned confirmation exposure and plan differ")
    measured = fixed["primary_confirmation"]
    for key in ("rate_fb", "mc_error_fb", "integration_error_fb"):
        if not math.isclose(spectrum["window_125"][key], measured[key], rel_tol=1.e-9, abs_tol=1.e-12):
            raise ValueError(f"binned and direct confirmation disagree for {key}")
    spectrum["confirmation"] = dict(plan=plan, generation_region="outside_nominal",
                                    primary_confirmation=measured,
                                    reference_root=str(root), hard_seeds=sorted(seeds))
    return spectrum


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--plan", type=Path, default=Path(__file__).with_name("ew_cards")/"zee_edge_plan.json")
    p.add_argument("--spectrum", type=Path, help="Full standard binned rates, covariance and shard provenance for confirmation seeds")
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    fixed = result(a.root, plan)
    a.output.write_text(json.dumps(fixed, indent=2)+"\n")
    if a.spectrum:
        model = confirmation_spectrum(a.root, plan, fixed)
        model["confirmation"]["plan_sha256"] = hashlib.sha256(a.plan.read_bytes()).hexdigest()
        opener = gzip.open if a.spectrum.suffix == ".gz" else open
        with opener(a.spectrum, "wt") as stream:
            json.dump(model, stream, separators=(",",":"),allow_nan=False)
            stream.write("\n")
