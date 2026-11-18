"""Read the normalized SSC 40 TeV loop and double-fake response tables."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

RELEASE = Path(__file__).resolve().parent / "HOAnalysis/student_backgrounds/ssc40_v1"
NAMES = ("bkg_loop_gg_aa", "bkg_dijet_double_fake")
DIJET_SCENARIO = ("The dijet double-fake yield uses the upper estimate from independent "
                  "quark-to-photon fake probabilities for both jets.")
FEATURES = ("m_gg", "pt_gamma1", "eta_gamma1", "pt_gamma2", "eta_gamma2",
            "deltaR_gg", "deltaPhi_gg", "pt_gg", "y_gg", "n_selected_photons")


def inventory(requested, directory=None):
    """Return metadata for the requested components in the published release."""
    directory = RELEASE if directory is None else Path(directory)
    manifest = directory / "manifest.json"
    index = json.loads(manifest.read_text())
    if (index["schema_version"] != 1 or index["energy_tev"] != 40 or index["detector_response"] != "ssc"
            or tuple(index["feature_names"]) != FEATURES or set(index["samples"]) != set(NAMES)):
        raise ValueError(f"invalid 40 TeV diphoton background inventory: {manifest}")
    records = {}
    for name in NAMES:
        if requested is not None and name not in requested:
            continue
        record = index["samples"][name]
        if record["file"] != name + ".npz":
            raise ValueError(f"unexpected diphoton background filename: {name}")
        path = directory / record["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"student background checksum differs: {path}")
        records[name] = (record, path, manifest)
    return records


def arrays(record, path):
    """Return all hypotheses with their differential cross sections in pb."""
    import numpy as np
    with np.load(path, allow_pickle=False) as table:
        variables, weights, sources = (table[key] for key in ("variables", "weight_pb", "sourceevent"))
    n = record["rows"]
    if (variables.shape != (n, len(FEATURES)) or weights.shape != (n,) or sources.shape != (n,)
            or sources.dtype.kind not in "iu" or np.any(sources < 0)
            or not np.isfinite(variables).all() or not np.isfinite(weights).all() or np.any(weights < 0)
            or len(np.unique(sources)) != record["represented_sources"]):
        raise ValueError(f"invalid SSC response population: {path}")
    for mask, key in ((slice(None), "sum_weight_pb"), (variables[:, -1] >= 2, "pair_rate_pb")):
        if not math.isclose(float(weights[mask].sum()), record[key], rel_tol=1.e-10, abs_tol=1.e-18):
            raise ValueError(f"{key} differs from the published response table: {path}")
    return variables, weights, sources
