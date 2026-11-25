"""Measure lepton and jet activity alongside baseline prompt and photon+jet candidates.

Run on Timur with uproot and NumPy. The native HwSim trees supply the leptons
and jets. Selected photon pairs and their signed weights come from the student
detector release. This check does not implement jet-to-lepton misidentification.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


BASELINE = Path("/home/shared/HiggsSSC/diphoton/ssc40-20261118")
NATIVE = Path("/home/apapaefs/Projects/HiggsSSC/hgammagamma/HOAnalysis/runs/ho_bg_1m_01/Backgrounds/events")
SAMPLES = ("bkg_prompt_aa", "bkg_gamma_j")
LEPTONS = ("Electrons", "Positrons", "Muons", "antiMuons")
JETS = ("Jets", "bJets")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1024**2), b""):
            digest.update(part)
    return digest.hexdigest()


def objects(arrays, key, row, pt_cut):
    n = min(arrays["num"+key][row], 100)
    _, px, py, pz = arrays["the"+key][row, :4, :n]
    pt = np.hypot(px, py)
    keep = pt > pt_cut
    eta = np.arcsinh(pz[keep] / pt[keep])
    phi = np.arctan2(py[keep], px[keep])
    return np.column_stack((eta, phi))


def activity(arrays, row):
    leptons = np.concatenate([objects(arrays, k, row, 20.) for k in LEPTONS])
    leptons = leptons[abs(leptons[:, 0]) < 2.5]
    jets = np.concatenate([objects(arrays, k, row, 5.) for k in JETS])
    hard_jets = np.concatenate([objects(arrays, k, row, 20.) for k in JETS])
    isolated = 0
    for eta, phi in leptons:
        dphi = (jets[:, 1] - phi + np.pi) % (2*np.pi) - np.pi
        isolated += int(np.all(np.hypot(jets[:, 0] - eta, dphi) > .3))
    return len(leptons), isolated, np.count_nonzero(abs(hard_jets[:, 0]) < 2.5)


def check(sample, baseline, native):
    import uproot

    release = json.loads((baseline / "diphoton_samples.json").read_text())["samples"][sample]
    detector = baseline / release["var_file"]
    archived_sha = json.loads((baseline / "dataset.json").read_text())["files"][release["var_file"]]["sha256"]
    if sha256(detector) != archived_sha:
        raise ValueError("Unrecognized detector file")
    with uproot.open(detector) as r:
        t = r["Data2"]
        v = t["variables"].array(library="np")
        cut = (v[:, 9] == 2) & (v[:, 7] > 40) & (v[:, 0] >= 100) & (v[:, 0] < 150)
        ids = t["sourceevent"].array(library="np").ravel()[cut]
        generator = t["generatorweight"].array(library="np").ravel()[cut]
        scale = 1000 * release["cross_section_pb"] / release["normalization_sum_weight"]
        weights = t["eventweight"].array(library="np").ravel()[cut] * scale
        masses = v[cut, 0]
    if np.any(np.diff(ids) <= 0) or release["weight_scale"] != 1:
        raise ValueError("Expected one selected hypothesis per source and unit generator scale")
    tree_path = native / sample / "herwig/events" / (sample + ".root")
    n_lepton, n_isolated, n_jet = np.zeros((3, len(ids)), dtype=int)
    with uproot.open(tree_path) as r:
        t = r["Data"]
        if t.num_entries != release["events_read"]:
            raise ValueError("Native and response exposure differ")
        np.testing.assert_array_equal(t["evweight"].array(library="np").ravel()[ids], generator)
        branches = [p+k for k in (*JETS, *LEPTONS) for p in ("num", "the")]
        for start in range(0, t.num_entries, 3000):
            stop = min(start + 3000, t.num_entries)
            positions = np.flatnonzero((ids >= start) & (ids < stop))
            if not len(positions):
                continue
            arrays = t.arrays(branches, entry_start=start, entry_stop=stop, library="np", how=dict)
            for pos in positions:
                row = ids[pos] - start
                n_lepton[pos], n_isolated[pos], n_jet[pos] = activity(arrays, row)
    windows = []
    for lo, hi in ((100., 150.), (120., 130.), (123.3, 126.7)):
        keep = (masses >= lo) & (masses < hi)
        selection = {}
        for label, values in (("pair_only", np.ones(len(ids))), ("fiducial_lepton", n_lepton > 0),
                              ("lepton_jet_isolation_proxy", n_isolated > 0),
                              ("jet_candidates_sum_including_fake_photon_jet", n_jet)):
            w = weights[keep] * values[keep]
            selection[label] = dict(source_count=int(np.count_nonzero(w)), rate_fb=float(w.sum()),
                                    mc_error_fb=float(np.linalg.norm(w)))
        windows.append(dict(window_gev=[lo, hi], results=selection))
    return dict(sample=sample, native_path=str(tree_path), native_sha256=sha256(tree_path),
                detector_sha256=archived_sha, maximum_abs_selected_source_weight_fb=float(abs(weights).max()),
                native_entries=release["events_read"], window_checks=windows)


def check_dijet(tables):
    """Count unused jet candidates from the complete list of unordered fake-photon pairs."""
    manifest = json.loads((tables / "manifest.json").read_text())["samples"]["bkg_dijet_double_fake"]
    if sha256(tables / manifest["file"]) != manifest["sha256"]:
        raise ValueError("Unrecognized dijet pair table")
    with np.load(tables / manifest["file"]) as data:
        pairs = data["variables"]
        keep = ((pairs[:, 9] == 2) & (pairs[:, 7] > 40)
                & (pairs[:, 0] >= 123.3) & (pairs[:, 0] < 126.7))
        _, source, multiplicity = np.unique(data["sourceevent"], return_inverse=True, return_counts=True)
        # There is one row for each unordered pair of distinct, eligible jets.
        jets = np.rint((1 + np.sqrt(1 + 8*multiplicity)) / 2).astype(int)
        np.testing.assert_array_equal(jets * (jets - 1) // 2, multiplicity)
        pairs_fb = 1000 * data["weight_pb"][keep]
        unused_jets = jets[source[keep]] - 2
    results = {}
    for label, w in (("pair_only", pairs_fb), ("unused_fiducial_jet_candidates", pairs_fb * unused_jets)):
        source_weights = np.bincount(source[keep], weights=w)
        results[label] = dict(source_count=int(np.count_nonzero(source_weights)),
                              rate_fb=float(source_weights.sum()),
                              mc_error_fb=float(np.linalg.norm(source_weights)))
    return dict(sample="bkg_dijet_double_fake", table_sha256=manifest["sha256"],
                generated_events=manifest["generated_sources"], selected_pairs=int(keep.sum()),
                window_gev=[123.3, 126.7], results=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=BASELINE)
    parser.add_argument("--native", type=Path, default=NATIVE)
    parser.add_argument("--dijet", action="store_true", help="Run the selected 125 GeV dijet control alone")
    parser.add_argument("--tables", type=Path, default=Path("hgammagamma/HOAnalysis/student_backgrounds/ssc40_v1"))
    args = parser.parse_args()
    result = check_dijet(args.tables) if args.dijet else [check(s, args.baseline, args.native) for s in SAMPLES]
    print(json.dumps(result, indent=2))
