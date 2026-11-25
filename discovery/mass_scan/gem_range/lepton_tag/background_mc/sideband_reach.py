"""Repeat the mass scan with the tagged continuum normalized to the blinded window estimates."""
import argparse
import hashlib
import json
import os
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np

from refit import add_zee_complement, discovery_model, load_spectrum, project_background, search


def calculate(core, guard, sidebands, toys, seed):
    report = load_spectrum(core)
    add_zee_complement(report, load_spectrum(guard))
    check = json.loads(sidebands.read_text())
    for path, field in ((core, "core_sha256"), (guard, "zee_complement_sha256")):
        if hashlib.sha256(path.read_bytes()).hexdigest() != check["inputs"][field]:
            raise ValueError("Sideband estimates were calculated from different input spectra")
    direct = report["window_125"]["rate_fb"]
    np.testing.assert_allclose(direct, check["nominal_direct_window"]["rate_fb"])
    signal, model = discovery_model()
    model.tag_background = project_background(model, np.asarray(report["edges_gev"]), np.asarray(report["rates_fb"]))
    fits = {}
    for label, rate in zip(("lower_rate", "upper_rate"), check["fit_rate_envelope_fb"]):
        fits[label] = search(model, signal["selected_fb"], rate/direct, toys, seed, max_iterations=12)
        np.testing.assert_allclose(fits[label]["background_in_true_window_fb"], rate)
    return dict(method="Full simulated tagged spectrum rescaled to the two endpoints of the blinded 125 GeV window estimates",
                direct_window_fb=direct, sideband_range_fb=check["fit_rate_envelope_fb"],
                sidebands_sha256=hashlib.sha256(sidebands.read_bytes()).hexdigest(),
                toys=toys, seed=seed, scenarios=fits)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("core", "guard", "sidebands", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--toys", type=int, default=40000)
    parser.add_argument("--seed", type=int, default=725081)
    args = parser.parse_args()
    if args.toys < 4000:
        parser.error("Use at least 4,000 experiments")
    report = calculate(args.core, args.guard, args.sidebands, args.toys, args.seed)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
