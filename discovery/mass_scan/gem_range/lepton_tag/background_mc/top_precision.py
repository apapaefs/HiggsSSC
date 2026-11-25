"""Finish a resolved top sample and extend it until its window precision is adequate."""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

TWO_PHOTON = ("tt_1a_1dec", "tt_0a_2dec_same", "tt_0a_2dec_split")
ONE_DECAY = ("tt_0a_1dec", "tt_0a_1dec_ll")


def sparse_upper(rate_fb, attempted, selected, sample):
    if selected >= 5:
        return None
    from scipy.stats import beta
    max_response = .988 * (.90 * .96)**(2 if sample in TWO_PHOTON else 1)
    return float(beta.ppf(.95, selected + 1, attempted - selected)) * rate_fb * max_response


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("out", "generator", "worker", "repo", "mg5", "plugin", "response", "photon-matcher", "summary"):
        parser.add_argument("--" + key, required=True, type=Path)
    parser.add_argument("--sample", required=True, choices=(*TWO_PHOTON, *ONE_DECAY))
    parser.add_argument("--cores", required=True, type=int)
    parser.add_argument("--target", required=True, type=float, help="window MC error budget in fb")
    parser.add_argument("--block", default=500, type=int)
    parser.add_argument("--wait-complete", action="store_true", help="Wait for all initial detector shards from another worker")
    parser.add_argument("--ready-marker", type=Path, help="Wait for a prerequisite file before reducing this source")
    a = parser.parse_args()
    if a.cores < 1 or a.block < 1 or a.target <= 0.:
        parser.error("invalid allocation or precision goal")
    out = a.out.resolve()
    source = Path(__file__).resolve().parent
    packets = out / "packets" / a.sample
    # One independent process owns all adaptive rounds of this source.
    import fcntl
    lock = (out / a.sample / "precision.lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    sys.path.insert(0, str(a.summary.resolve().parent))
    from response_summary import summarize
    if a.ready_marker:
        while not a.ready_marker.is_file():
            time.sleep(10)
    round_number = 0
    next_events = None
    status = out / a.sample / "precision.json"
    if status.is_file():
        previous = json.loads(status.read_text())
        if previous["target_fb"] != a.target or previous["sample"] != a.sample:
            raise ValueError("precision settings differ from the existing run")
        if previous["complete"]:
            available = sum(json.loads(p.read_text())["attempted_events"] for p in packets.glob("*/*/receipt.json"))
            if available == previous["attempted_events"]:
                return
            round_number = previous["round"]
        else:
            round_number = previous["round"] + 1
            next_events = previous["next_events"]
    while True:
        hard = out / a.sample / f"round-{round_number:03}"
        if round_number == 0:
            # The pilot dispatcher checks the same 500 source events, with an
            # absolute four-momentum comparison of each input photon to its history.
            while not (packets / "round-000/0000000/receipt.json").is_file():
                time.sleep(10)
            if a.wait_complete:
                total = json.loads((hard / "hard.json").read_text())["generated_events"]
                while sum(json.loads(p.read_text())["attempted_events"]
                          for p in (packets / "round-000").glob("*/receipt.json")) < total:
                    time.sleep(10)
        elif not (hard / "hard.json").is_file():
            if shutil.disk_usage(out).free < 8 * 1024**3:
                raise RuntimeError("Less than 8 GiB free for the next independent radiative integration")
            subprocess.run([sys.executable, source / "top_radiative.py", "--out", out, "--generator", a.generator,
                            "--mg5", a.mg5, "--cores", str(a.cores), "--events", str(next_events),
                            "--round", str(round_number), "--samples", a.sample], check=True)
        subprocess.run([sys.executable, source / "top_shower.py", "--out", out, "--worker", a.worker,
                        "--repo", a.repo, "--plugin", a.plugin, "--response", a.response,
                        "--photon-matcher", a.photon_matcher, "--cores", str(a.cores), "--radiative",
                        "--samples", a.sample, "--round", str(round_number), "--block", str(a.block)], check=True)
        report = summarize(packets, radiative=True, dileptonic=a.sample in ONE_DECAY)["sources"][a.sample]
        window = report["window_125"]
        n = report["attempted_events"]
        k = window["contributing_events"]
        se = window["mc_error_fb"]
        desired = n * (se / a.target)**2
        # A handful of selected collisions can give a small empirical SE by chance.
        # Use a 95% one-sided binomial bound for this sparse source. The detector
        # probability is bounded by the matched photon IDs and the trigger efficiency.
        sigma = max(group["cross_section_fb"] for group in report["integration_groups"])
        upper = sparse_upper(sigma, n, k, a.sample)
        if upper is not None and upper > 2 * a.target:
            desired = max(desired, n * upper / (2 * a.target))
        complete = se <= a.target and (upper is None or upper <= 2 * a.target)
        next_events = max(10000, min(200000, math.ceil(1.25 * (desired - n) / a.block) * a.block))
        current = dict(sample=a.sample, round=round_number, complete=complete, target_fb=a.target,
                       attempted_events=n, contributing_events=k, window=window, sparse_upper_95_fb=upper,
                       next_events=None if complete else next_events)
        temp = status.with_suffix(".partial")
        temp.write_text(json.dumps(current, indent=2) + "\n")
        os.replace(temp, status)
        print(json.dumps(current), flush=True)
        if complete:
            return
        round_number += 1


if __name__ == "__main__":
    main()
