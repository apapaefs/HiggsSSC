"""Repeat the tagged eeγ selection with individual native PYTHIA 5.6 switches."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

SELECTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SELECTION))
from selection_rates import model

# The first three entries switch off one ingredient. The scale factors
# multiply kT squared, and correspond to halving or doubling the ISR scale.
SETTINGS = {
    "no_isr": (0, 2, 1, 1.), "no_timelike_qed": (1, 1, 1, 1.),
    "no_mpi": (1, 2, 0, 1.), "isr_half_scale": (1, 2, 1, .25),
    "isr_double_scale": (1, 2, 1, 4.), "nominal": (1, 2, 1, 1.),
}
GEM_FB = .0374


def measure(output, mode, seed, n, baseline):
    rows = []
    metadata = {}
    for line in output.splitlines():
        if line.startswith('{"event":'):
            rows.append(json.loads(line))
        else:
            items = line.split()
            if items and items[0] in ("RATE", "MATCH", "PROB", "SWITCHES", "ISR_SCALE"):
                metadata[items[0]] = items[1:]
            if items[:3] == ["BIN", "2", "pair_pt"]:
                nominal_count = int(items[3])
    if list(map(int, metadata["RATE"][:3])) != [seed, n, 1]:
        raise ValueError("Wrong seed, exposure or Z decay")
    matched, total = map(int, metadata["MATCH"])
    if total != n or matched / n < .999:
        raise ValueError("Unexpected hard-object ancestry loss")
    if list(map(int, metadata["SWITCHES"])) != list(SETTINGS[mode][:3]):
        raise ValueError("Unexpected ISR, final-state QED or MPI switch")
    # The first off-switch runs used a driver that printed only the switches.
    if mode.startswith("isr_") or "ISR_SCALE" in metadata:
        np.testing.assert_equal(list(map(float, metadata["ISR_SCALE"])), [2, SETTINGS[mode][3]])
    sigma_pb = float(metadata["RATE"][3])
    norm = 1000. * sigma_pb / n * 1.2 * (.03 * 33.)
    weights = np.array([r["probability"] for r in rows])
    variants, _ = model(rows)
    selected = variants["with_tracks"]
    np.testing.assert_allclose(variants["nominal"].sum(), float(metadata["PROB"][0]) * nominal_count)
    rate = norm * selected.sum()
    if mode == "nominal":
        ref = baseline[seed]
        np.testing.assert_allclose(rate, ref["rates_fb"]["with_tracks"], rtol=1e-11)
        if ref["pre_isolation_assignments"] != len(rows):
            raise ValueError("Historical baseline assignments changed")
    pt = np.array([[r[k]["true_pt"] for k in ("real", "fake", "tag")] for r in rows])
    phi = np.array([[r[k]["true_phi"] for k in ("real", "fake", "tag")] for r in rows])
    recoil = np.hypot((pt * np.cos(phi)).sum(axis=1), (pt * np.sin(phi)).sum(axis=1))
    cone = np.array([r["fake"]["cell_iso_et"] for r in rows])
    return dict(seed=seed, generated=n, matched=matched, hard_sigma_pb=sigma_pb,
                jetset_error_messages=output.count("Error type"),
                assignments=len(rows), before_iso_fb=norm * weights.sum(), selected_fb=rate,
                recoil_numerator_gev_fb=norm * (weights @ recoil),
                fake_cone_numerator_gev_fb=norm * (weights @ cone))


def summary(runs):
    def estimate(blocks):
        observables = ["hard_sigma_pb", "before_iso_fb", "selected_fb",
                       "recoil_numerator_gev_fb", "fake_cone_numerator_gev_fb"]
        sums = {key: np.mean([b[key] for b in blocks]) for key in observables}
        pre = sums["before_iso_fb"]
        return dict(hard_sigma_pb=sums["hard_sigma_pb"], selected_fb=sums["selected_fb"],
                    before_iso_fb=pre, pass_fraction=sums["selected_fb"] / pre,
                    mean_bare_recoil_gev=sums["recoil_numerator_gev_fb"] / pre,
                    mean_fake_cone_et_gev=sums["fake_cone_numerator_gev_fb"] / pre)

    full = estimate(runs)
    n = len(runs)
    leave = [estimate(runs[:i] + runs[i+1:]) for i in range(n)]
    result = {}
    for key, value in full.items():
        v = np.array([rep[key] for rep in leave])
        result[key] = dict(value=float(value), mc_se=float(np.sqrt((n-1)/n * ((v-v.mean())**2).sum())))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--read-logs", type=Path, help="Read an existing complete set of generator logs")
    args = parser.parse_args()
    refs = {r["seed"]: r for r in json.loads((SELECTION / "native.json").read_text())["runs"]}
    jobs = []
    for j, mode in enumerate(SETTINGS):
        for i, reference in enumerate(refs.values()):
            seed = reference["seed"] if mode == "nominal" else 7200001 + 10000*j + 97*i
            jobs.append((mode, seed, reference["events"]))
    if len({job[1] for job in jobs}) != len(jobs):
        raise ValueError("Seed reused between settings")

    def execute(job):
        mode, seed, n = job
        if args.read_logs:
            output = (args.read_logs / f"{mode}-{seed}.txt").read_text()
        else:
            output = subprocess.check_output([str(args.executable.resolve()), str(seed),
                                              str(n), "1", mode], text=True)
        return mode, measure(output, mode, seed, n, refs)

    runs = {mode: [] for mode in SETTINGS}
    with ThreadPoolExecutor(args.workers) as pool:
        for mode, block in pool.map(execute, jobs):
            runs[mode].append(block)
    results = {mode: summary(blocks) for mode, blocks in runs.items()}
    ref = results["nominal"]["selected_fb"]
    for mode, result in results.items():
        rate = result["selected_fb"]
        ratio = rate["value"] / ref["value"]
        result["ratio_to_gem"] = rate["value"] / GEM_FB
        result["ratio_to_nominal"] = dict(value=ratio, mc_se=0. if mode == "nominal" else
            ratio * float(np.hypot(rate["mc_se"] / rate["value"], ref["mc_se"] / ref["value"])))
    args.output.write_text(json.dumps(dict(window_gev=[123.3,126.7], gem_fb=GEM_FB,
        seed_blocks_per_setting=len(refs), switches_order=["MSTP(61)", "MSTJ(41)", "MSTP(81)", "PARP(64)"],
        switches=SETTINGS, runs=runs, results=results), indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
