"""Run the combined reach on Timur once the requested tagged samples and precision are available."""
import argparse
from collections import Counter
import gzip
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from top_precision import sparse_upper

HERE = Path(__file__).resolve().parent


def wait_for_fit(root, output, expected, toys, interval, zee_complement):
    output.mkdir(parents=True, exist_ok=True)
    spectrum = output / "background.json"
    if (output / "combined-mc.json").is_file():
        raise FileExistsError(output / "combined-mc.json")
    extra_se = 0.
    options = []
    if zee_complement is not None:
        opener = gzip.open if zee_complement.suffix == ".gz" else open
        with opener(zee_complement, "rt") as stream:
            extra_se = json.load(stream)["window_125"]["error_fb"]
        options = ["--zee-complement", zee_complement]
    previous = None
    while True:
        counts = Counter()
        for packet_root in (root/"ew/showers", root/"top/packets"):
            for path in packet_root.rglob("receipt.json"):
                receipt = json.loads(path.read_text())
                counts[receipt["sample"]] += receipt["attempted_events"]
        if counts != previous:
            print(time.strftime("%Y-%m-%d %H:%M:%S %Z"), dict(sorted(counts.items())), flush=True)
            previous = counts
            if all(counts[name] >= minimum for name, minimum in expected.items()):
                mode = "--top-dileptonic" if "tt_0a_1dec_ll" in expected else "--top-radiative"
                subprocess.run([sys.executable, HERE/"response_summary.py", root, spectrum, mode],
                               check=True)
                model = json.loads(spectrum.read_text())
                window = model["window_125"]
                sample = "tt_0a_2dec_split"
                source = model["sources"][sample]
                split = source["window_125"]
                upper = sparse_upper(max(g["cross_section_fb"] for g in source["integration_groups"]),
                                     source["attempted_events"], split["contributing_events"], sample)
                if (model["complete_sources"] and split["error_fb"] <= .011
                        and (upper is None or upper <= .022)
                        and math.hypot(window["error_fb"], extra_se) + window["zero_strata_upper_fb_95"] <= .0408):
                    subprocess.run([sys.executable, HERE/"refit.py", spectrum, "--output", output,
                                    "--toys", str(toys), *options], check=True)
                    return
        time.sleep(interval)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("minimum_counts", type=Path, help="JSON object giving attempted hard events by source")
    parser.add_argument("output", type=Path)
    parser.add_argument("--toys", type=int, default=40000)
    parser.add_argument("--interval", type=int, default=60)
    parser.add_argument("--zee-complement", type=Path)
    a = parser.parse_args()
    required = json.loads(a.minimum_counts.read_text())
    if not required or any(type(n) is not int or n <= 0 for n in required.values()) or a.interval < 1:
        parser.error("Use positive event counts and poll interval")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[key] = "1"
    wait_for_fit(a.campaign.resolve(), a.output.resolve(), required, a.toys, a.interval, a.zee_complement)
