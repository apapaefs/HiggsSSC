"""Generate the three independent top-background matrix elements on Timur."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys


SAMPLES = ("tt_0a", "tt_1a", "tt_2a")


def hard_job(sample, args):
    index = SAMPLES.index(sample)
    work = args.out / sample / f"round-{args.round:03}"
    work.mkdir(parents=True, exist_ok=True)
    if (work / "hard.json").exists():
        return
    source = Path(__file__).resolve().parent
    options = []
    if index:
        options += ["--set", "pta=14.", "--set", "etaa=3.5"]
    if index == 2:
        options += ["--set", "draa=0.", "--set", "mmaa=60.", "--set", "mmaamax=170."]
    seed = 721000 + 100 * args.round + index
    command = [sys.executable, source / "ew_generate.py", "--mg5", args.mg5, "--process", sample,
               "--work", work, "--events", str(args.events[index]), "--cores", str(args.cores // len(args.samples)),
               "--seed", str(seed), "--pdf", "335900", "--card", source / "top_cards" / (sample + ".mg5"),
               *options]
    if (work / "madgraph/Events/hard/unweighted_events.lhe.gz").exists():
        command.append("--resume-completed")
    with (work / "driver.log").open("a") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--mg5", required=True, type=Path)
    parser.add_argument("--cores", default=36, type=int)
    parser.add_argument("--events", nargs=3, default=[30000, 30000, 30000], type=int,
                        metavar=("TT", "TTG", "TTGG"))
    parser.add_argument("--round", default=0, type=int)
    parser.add_argument("--samples", nargs="+", choices=SAMPLES, default=list(SAMPLES))
    args = parser.parse_args()
    if len(set(args.samples)) != len(args.samples) or args.cores < len(args.samples) or min(args.events) <= 0 or args.round < 0:
        parser.error("invalid resource count, event count or round")
    with ThreadPoolExecutor(max_workers=len(args.samples)) as pool:
        jobs = [pool.submit(hard_job, sample, args) for sample in args.samples]
        for job in jobs:
            job.result()


if __name__ == "__main__":
    main()
