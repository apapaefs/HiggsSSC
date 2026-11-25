"""Repeat the detector variations on independent, native PYTHIA 5.6 showers."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from selection_rates import model

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vintage_runs import SEEDS


def generate(executable, reference):
    seed, n = reference["seed"], reference["events"]
    output = subprocess.check_output([str(executable), str(seed), str(n), "1"], text=True)
    lines = output.splitlines()
    entries = [json.loads(line) for line in lines if line.startswith('{"event":')]
    rate, = [float(line.split()[-1]) for line in lines if line.startswith("RATE ")]
    matched, = [list(map(int,line.split()[1:])) for line in lines if line.startswith("MATCH ")]
    if matched != [reference["matched"], n]:
        raise ValueError(f"Native generator ancestry differs for seed {seed}")
    np.testing.assert_allclose(rate, reference["sigma_pb"], rtol=2.e-7)
    variants, _ = model(entries)
    rates = {k: float(v.sum()) * 1000 * rate / n * 1.2 * .99 for k,v in variants.items()}
    expected = 1000 * rate / n * reference["probability"] * reference["counts"]["123p3_126p7"]["pair_pt"][0] * 1.2 * .99
    np.testing.assert_allclose(rates["nominal"], expected, rtol=1.e-10)
    return {"seed":seed, "events":n, "pre_isolation_assignments":len(entries), "rates_fb":rates}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    path = Path(__file__).resolve().parents[1] / "vintage.json"
    refs = [r for r in json.loads(path.read_text())["runs"] if r["flavours"] == 1]
    assert tuple(r["seed"] for r in refs) == SEEDS
    with ThreadPoolExecutor(args.workers) as pool:
        runs = list(pool.map(lambda r: generate(args.executable.resolve(), r), refs))
    names = list(runs[0]["rates_fb"])
    rates = np.array([[run["rates_fb"][key] for key in names] for run in runs])
    cov = np.cov(rates, rowvar=False, ddof=1) / len(runs)
    values = rates.mean(axis=0)
    result = {}
    for i, key in enumerate(names):
        gradient = np.zeros(len(names))
        gradient[i] += 1.
        gradient[names.index("nominal")] -= 1.
        result[key] = {"fb": values[i], "mc_se_fb": np.sqrt(cov[i,i]),
                       "change_from_nominal_fb": gradient @ values,
                       "paired_mc_se_fb": np.sqrt(max(0., gradient @ cov @ gradient))}
    report = {"generated":sum(r["events"] for r in runs), "observable_order":names,
              "results":result, "covariance_fb2":cov.tolist(), "runs":runs}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
