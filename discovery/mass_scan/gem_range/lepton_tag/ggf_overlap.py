"""Measure the ggF signal exposed to a lepton veto in the original diphoton sample."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import uproot


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def lepton_mask(arrays, pt_min, eta_max):
    passed = np.zeros(len(arrays["evweight"]), dtype=bool)
    for name in ("Electrons", "Positrons", "Muons", "antiMuons"):
        _, px, py, pz = arrays["the" + name].transpose(1, 0, 2)
        present = np.arange(px.shape[1]) < arrays["num" + name][:, None]
        pt = np.hypot(px, py)
        passed |= np.any(present & (pt > pt_min) & (abs(pz) < np.sinh(eta_max) * pt), axis=1)
    return passed


def measure(raw, release):
    manifest_path = release / "diphoton_samples.json"
    row = json.loads(manifest_path.read_text())["samples"]["signal_gg_h_aa"]
    response = release / row["var_file"]
    checked = json.loads((release / "dataset.json").read_text())["files"]
    for path in (manifest_path, response):
        if file_hash(path) != checked[str(path.relative_to(release))]["sha256"]:
            raise ValueError(f"Student data hash differs for {path.name}")
    with uproot.open(response) as f:
        data = f["Data2"].arrays(["variables", "sourceevent", "eventweight"], library="np", how=dict)
    ids = data["sourceevent"].ravel()
    variables = data["variables"]
    source_count = row["events_read"]
    if not (ids.min() == 0 and ids.max() == source_count - 1):
        raise ValueError("Diphoton response IDs are inconsistent with the event manifest")
    keep = (variables[:, -1] >= 2) & (variables[:, 0] >= 120) & (variables[:, 0] < 130)
    weights = data["eventweight"].ravel() * (1000 * row["cross_section_pb"] * row["weight_scale"]
                                             / row["normalization_sum_weight"])
    size = int(source_count)
    reference = np.bincount(ids[keep], weights[keep], minlength=size)
    possible = keep & (variables[:, 7] > 40.)
    pair_pt = np.bincount(ids[possible], weights[possible], minlength=size)
    if not np.isclose(reference.sum(), 168.7822888236, rtol=1e-8):
        raise ValueError("Selected 125 GeV signal differs from published diphoton rate")
    names = ["evweight"] + [p + n for n in ("Electrons", "Positrons", "Muons", "antiMuons")
                             for p in ("num", "the")]
    nominal, loose, signed = [], [], []
    with uproot.open(raw) as f:
        tree = f["Data"]
        if tree.num_entries != size:
            raise ValueError("Shower source count differs from diphoton response file")
        for a in tree.iterate(names, step_size=1024, library="np", how=dict):
            nominal.append(lepton_mask(a, 20., 2.5))
            loose.append(lepton_mask(a, 15., 3.))
            signed.extend(a["evweight"])
    # The weighted response enumerates all photon outcomes for each source event.
    summed_response = np.bincount(ids, data["eventweight"].ravel(), minlength=size)
    np.testing.assert_allclose(summed_response, np.asarray(signed) * row["weight_scale"], rtol=1e-9, atol=1e-13)
    results = {}
    for name, parts in (("pt20_eta2p5", nominal), ("pt15_eta3", loose)):
        subset = pair_pt[np.concatenate(parts)]
        results[name] = dict(eligible_sources=int(np.count_nonzero(subset)), signed_fb=float(subset.sum()),
                             absolute_weight_fb=float(np.abs(subset).sum()), mc_error_fb=float(np.linalg.norm(subset)),
                             fraction=float(subset.sum()/reference.sum()))
    return dict(shower_file_sha256=file_hash(raw), response_file_sha256=file_hash(response),
                generated_events=row["shower_quality"]["attempted_events"], shower_events=size,
                mass_window_gev=[120, 130], inclusive_ggf_fb=float(reference.sum()),
                pair_pt_min_gev=40, isolation_and_id_applied_to_leptons=False,
                measure="Candidates with any shower electron or muon passing the stated kinematic cuts. "
                        "The signed rates precede lepton isolation and identification.",
                results=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path, help="Original ggF shower HwSim ROOT file")
    parser.add_argument("release", type=Path, help="Student diphoton release directory")
    args = parser.parse_args()
    print(json.dumps(measure(args.raw, args.release), indent=2, allow_nan=False))
