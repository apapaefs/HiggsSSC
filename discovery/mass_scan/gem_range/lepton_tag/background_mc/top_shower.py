"""Reduce disjoint blocks of the three top integrations, with bounded parallelism."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

from top_campaign import SAMPLES
from top_radiative import SAMPLES as RADIATIVE


def shower(job):
    packet, command = job
    if (packet / "receipt.json").exists():
        return
    packet.mkdir(parents=True, exist_ok=True)
    with (packet / "worker.log").open("a") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("out", "worker", "plugin", "response", "repo"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--cores", type=int, default=36)
    parser.add_argument("--round", type=int, default=0)
    parser.add_argument("--block", type=int, default=500)
    parser.add_argument("--radiative", action="store_true")
    parser.add_argument("--samples", nargs="+", choices=(*SAMPLES, *RADIATIVE))
    parser.add_argument("--photon-matcher", type=Path)
    parser.add_argument("--limit-events", type=int, help="Only shower this initial interval of each integration")
    parser.add_argument("--packets", type=Path, help="Separate packet root for a shower variation")
    parser.add_argument("--shower-card", type=Path)
    parser.add_argument("--lanes", type=int, default=1, help="Partition block indices into this many disjoint lanes")
    parser.add_argument("--lane", type=int, default=0)
    args = parser.parse_args()
    if args.cores < 1 or args.block < 1 or args.round < 0 or not 0 <= args.lane < args.lanes:
        parser.error("invalid parallelism, block size or round")
    if args.limit_events is not None and (args.limit_events <= 0 or args.limit_events % args.block):
        parser.error("limit-events must cover a positive integer number of full blocks")
    jobs = []
    round_id = f"round-{args.round:03}"
    samples = RADIATIVE if args.radiative else SAMPLES
    if args.radiative and args.photon_matcher is None:
        parser.error("radiative decay samples require LHE photon matching")
    if args.samples and not set(args.samples).issubset(samples):
        parser.error("requested samples do not match the production family")
    packets = args.packets or args.out / "packets"
    selected = args.samples or (RADIATIVE[:4] if args.radiative else SAMPLES)
    for i, sample in enumerate(samples):
        if sample not in selected:
            continue
        hard = args.out / sample / round_id
        info = json.loads((hard / "hard.json").read_text())
        total = info["generated_events"]
        if args.limit_events is not None:
            total = min(total, args.limit_events)
        for j, start in enumerate(range(0, total, args.block)):
            if j % args.lanes != args.lane:
                continue
            count = min(args.block, total - start)
            packet = packets / sample / round_id / f"{start:07}"
            seed = (82100000 if args.radiative else 81100000) + 100000 * args.round + 10000 * i + j
            if j >= 10000:
                raise ValueError("too many blocks for this seed range")
            command = [sys.executable, args.worker, "--hard", hard, "--out", packet, "--sample", sample,
                       "--seed", str(seed), "--start", str(start), "--events", str(count),
                       "--plugin", args.plugin, "--response", args.response, "--repo", args.repo]
            if j == 0:
                command.append("--keep-root")
            if args.photon_matcher:
                command += ["--photon-matcher", args.photon_matcher]
            if args.shower_card:
                command += ["--card", args.shower_card]
            jobs.append((j, i, packet, command))
    # Interleave samples so every stratum contributes to early closure checks.
    with ThreadPoolExecutor(max_workers=args.cores) as pool:
        for _ in pool.map(shower, [(p, c) for _, _, p, c in sorted(jobs)]):
            pass


if __name__ == "__main__":
    main()
