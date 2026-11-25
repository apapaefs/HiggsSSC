"""Run independent seeds of GEM's generator and apply the project's tagged response."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import statistics
import subprocess

WINDOWS = ("75_165", "80_140", "123p3_126p7")
SEEDS = tuple(range(827101, 827113))
SOURCES = {"pythia56.f":"0470ac77f0cf53ae30bd02a8fcff76156096d6057e5b61562dc3b501275c82b4",
           "jetset73.f":"d33a6c72fcb1bf5797828a34cfa76169119078419ebf9899e0a0378f69c1bf6f"}


def run_seed(exe, cfg):
    seed, events, flavours = cfg
    text = subprocess.check_output([exe, str(seed), str(events), str(flavours)], text=True)
    result = {"seed":seed, "events":events, "flavours":flavours}
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "RATE":
            if list(map(int, fields[1:4])) != list(cfg):
                raise ValueError("Generator seed, exposure or decay convention changed")
            result["sigma_pb"] = float(fields[4])
        if fields[0] == "MATCH":
            result["matched"], generated = map(int, fields[1:])
            if result["matched"]/events < .999 or generated != events:
                raise ValueError(f"Incomplete native identification of Z daughters: {line}")
        if fields[0] == "PROB":
            result["probability"] = float(fields[1])
        if fields[0] == "BIN":
            window, stage, total, squares = fields[1:]
            result.setdefault("counts", {}).setdefault(WINDOWS[int(window)], {})[stage] = (int(total), int(squares))
    if "sigma_pb" not in result or (flavours == 1 and "counts" not in result):
        raise ValueError("Incomplete native generator output")
    return result


def estimate(values):
    return {"mean":statistics.mean(values), "se":statistics.stdev(values)/math.sqrt(len(values))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("executable", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--events", type=int, default=80000)
    a = p.parse_args()
    if a.workers < 1 or a.events < 1000:
        p.error("Use a positive worker count and at least 1000 events per seed")
    configs = [(s, a.events, 1) for s in SEEDS] + [(s+100, 10000, 2) for s in SEEDS]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        runs = list(pool.map(lambda c: run_seed(a.executable.resolve(), c), configs))
    electron = [r for r in runs if r["flavours"]==1]
    fits = {}
    for window in WINDOWS:
        rates = {}
        for stage in electron[0]["counts"][window]:
            rate = [r["sigma_pb"]*1000*r["probability"]*r["counts"][window][stage][0]/r["events"] for r in electron]
            rates[stage] = estimate(rate)
        fits[window] = rates
    report = {"process":"PYTHIA 5.6 subprocess 19, Z decays to ee or ee + mumu",
              "pp_tev":40, "pdf":"EHLQ1", "partonic_pt_min_gev":10, "official_fortran_sha256":SOURCES,
              "reported_qq_only":True, "historical_multiplier":1.20,
              "uncertainties":"Standard error across the independent generator seeds, including integration and selection",
              "sigma_ee_pb":estimate([r["sigma_pb"] for r in electron]),
              "sigma_ee_mumu_pb":estimate([r["sigma_pb"] for r in runs if r["flavours"]==2]),
              "rates_with_object_probabilities_fb":fits, "runs":runs}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")

if __name__ == "__main__":
    main()
