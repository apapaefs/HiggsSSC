"""Evaluate SSC dijet double-fake envelopes for the mass window and frozen classifier."""
import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hgammagamma import ho_shower_completion as completion
from hgammagamma.dijet_response import load_model, response_rows


def evaluate(bin_dirs, analyzer, model=None):
    bins = []
    totals = {}
    previous_high = None
    sources_seen = set()
    roots_seen = set()
    for directory in map(Path, bin_dirs):
        source = json.loads((directory / "source.json").read_text())
        low, high = source["bin_GeV"]
        if low <= 0 or (high != -1 and high <= low) or (previous_high is not None and previous_high != low):
            raise ValueError("supply consecutive, non-overlapping Born bins in increasing order")
        if source["lhe_sha256"] in sources_seen:
            raise ValueError("repeated dijet source")
        previous_high = high
        sources_seen.add(source["lhe_sha256"])

        proof = json.loads((directory / "shower.json").read_text())
        n = source["events"]
        completion.validate_completion(
            proof, dict(sample="dijet", nevents_requested=n, lhe=dict(events=n)), directory=directory)
        if (proof["source_weights"]["artifact"]["sha256"] != source["lhe_sha256"] or
                proof["source_weights"]["sum_weight"] != n or proof["attempted_events"] != n or n < 2):
            raise ValueError("inconsistent hard source or insufficient population")

        groups = {}
        for root in proof["artifacts"]["roots"]:
            if root["sha256"] in roots_seen:
                raise ValueError("reused HwSim ROOT")
            roots_seen.add(root["sha256"])
            for row in response_rows(analyzer, root["path"], source["seed"] + 2, model):
                key = row["scenario"] + "@" + row["fake_scale"]
                moments = groups.setdefault(key, [0., 0., 0, 0])
                moments[0] += float(row["sum_probability"])
                moments[1] += float(row["sum_probability_squared"])
                moments[2] += int(row["events_read"])
                moments[3] += int(row["contributing_events"])

        generated_fb = 1000 * source["cross_section_pb"]
        estimates = {}
        for key, (sum_p, sum_p2, read, count) in groups.items():
            if read != proof["saved_events"]:
                raise ValueError("analyzer source count differs from finite shower")
            rate = generated_fb * sum_p / n
            variance = generated_fb**2 * max(0, sum_p2 - sum_p**2 / n) / (n * (n-1))
            integration_error = 1000 * source["integration_error_pb"] * sum_p / n
            estimates[key] = dict(rate_fb=rate, mc_error_fb=math.sqrt(variance),
                                  integration_error_fb=integration_error, contributing_sources=count)
            total = totals.setdefault(key, [0., 0., 0.])
            total[0] += rate
            total[1] += variance
            total[2] += integration_error**2
        bins.append(dict(bin_GeV=[low, high], generated=n, saved=proof["saved_events"], seed=source["seed"],
                         source_sha256=source["lhe_sha256"], cross_section_pb=source["cross_section_pb"],
                         estimates=estimates))

    return dict(
        window_GeV=[120, 130], bins=bins, frozen_model=model[1] if model else None,
        analyzer_sha256=completion._artifact(analyzer)["sha256"],
        interpretation="Pair-probability envelopes in the generated Born domain. MC errors describe sampled events. "
                       "Empty bins can still contribute. The calculation covers only the generated range.",
        totals={key: dict(rate_fb=rate, mc_error_fb=math.sqrt(variance), integration_error_fb=math.sqrt(mg_variance))
                for key, (rate, variance, mg_variance) in totals.items()},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bins", type=Path, nargs="+", required=True)
    parser.add_argument("--analyzer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    args = parser.parse_args()
    result = evaluate(args.bins, args.analyzer.resolve(), load_model(args.frozen_dir))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
