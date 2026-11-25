"""Reproject the November 18 student diphoton samples and check the published overlap."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import tarfile

import numpy as np
from scipy.sparse import coo_matrix
import uproot


EDGES = np.arange(70., 171.)


def verify_hash(path, expected):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected:
        raise ValueError(f"Source hash mismatch for {path}")
    return digest


def project(variables, weights, ids, count):
    masses = variables[:, 0]
    if ids.dtype.kind not in "iu" or ids.min() < 0 or ids.max() >= count:
        raise ValueError("Invalid source-event IDs")
    if not np.isfinite(weights).all() or not np.isfinite(masses).all():
        raise ValueError("Nonfinite source weights or masses")
    keep = (variables[:, -1] >= 2) & (masses >= EDGES[0]) & (masses < EDGES[-1])
    bins = np.searchsorted(EDGES, masses[keep], side="right") - 1
    contributions = coo_matrix((weights[keep], (ids[keep], bins)),
                               shape=(count, len(EDGES) - 1)).tocsr()
    rate = np.asarray(contributions.sum(axis=0)).ravel()
    covariance = (contributions.T @ contributions).toarray()
    return rate, covariance, int(np.unique(ids[keep]).size)


def run(dataset, tables, archive, output):
    dataset_files = json.loads((dataset / "dataset.json").read_text())["files"]
    sample_file = dataset / "diphoton_samples.json"
    verify_hash(sample_file, dataset_files[sample_file.name]["sha256"])
    native = json.loads(sample_file.read_text())["samples"]
    packaged = json.loads((tables / "manifest.json").read_text())["samples"]
    with tarfile.open(archive) as tar:
        baseline = json.load(tar.extractfile("baseline.json"))["datasets"]
    rates, covariances, metadata = [], [], []
    for source in baseline:
        sample, stratum, count = (source[k] for k in ("sample", "stratum", "generated_events"))
        if sample in native and sample not in packaged:
            item = native[sample]
            path = dataset / item["var_file"]
            digest = verify_hash(path, dataset_files[item["var_file"]]["sha256"])
            if (item["shower_quality"]["attempted_events"], stratum) != (count, "full"):
                raise ValueError(f"Generator count or stratum mismatch for {sample}")
            with uproot.open(path) as root:
                tree = root["Data2"]
                variables = tree["variables"].array(library="np")
                ids = tree["sourceevent"].array(library="np").reshape(-1)
                raw_weights = tree["eventweight"].array(library="np").reshape(-1)
            factor = 1000 * item["cross_section_pb"] * item["weight_scale"] / item["normalization_sum_weight"]
            rate, covariance, selected = project(variables, raw_weights * factor, ids, count)
        else:
            item = packaged[sample]
            path = tables / item["file"]
            digest = verify_hash(path, item["sha256"])
            with np.load(path, allow_pickle=False) as data:
                ids = data["sourceevent"].reshape(-1)
                if stratum == "full":
                    start, stop = 0, item["generated_sources"]
                else:
                    start = 0
                    for part in item["inputs"]:
                        stop = start + part["generated"]
                        if part["region"] == stratum:
                            break
                        start = stop
                    else:
                        raise ValueError(f"Unknown packaged stratum {stratum}")
                if stop - start != count:
                    raise ValueError(f"Generator count mismatch for {sample}/{stratum}")
                domain = (ids >= start) & (ids < stop)
                rate, covariance, selected = project(data["variables"][domain],
                    1000 * data["weight_pb"][domain].reshape(-1), ids[domain] - start, count)
        old_edges = np.asarray(source["edges_gev"])
        low = int(np.searchsorted(EDGES, old_edges[0]))
        high = low + len(old_edges) - 1
        np.testing.assert_array_equal(EDGES[low:high+1], old_edges)
        rate_diff = rate[low:high] - source["rate_fb"]
        cov_diff = covariance[low:high, low:high] - source["covariance_fb2"]
        np.testing.assert_allclose(rate[low:high], source["rate_fb"], atol=1e-9, rtol=1e-9)
        np.testing.assert_allclose(covariance[low:high, low:high], source["covariance_fb2"], atol=1e-9, rtol=1e-9)
        metadata.append(dict(sample=sample, stratum=stratum, independent_id=source["independent_id"],
            source_file=path.name, sha256=digest, selected_sources=selected,
            max_overlap_rate_difference_fb=float(np.abs(rate_diff).max()),
            max_overlap_covariance_difference_fb2=float(np.abs(cov_diff).max())))
        rates.append(rate)
        covariances.append(covariance)
    np.savez_compressed(output, edges_gev=EDGES, rate_fb=rates, covariance_fb2=covariances,
        sample=[x["sample"] for x in baseline], stratum=[x["stratum"] for x in baseline],
        generated_events=[x["generated_events"] for x in baseline], metadata_json=json.dumps(metadata))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True, help="ssc40-20261118 native ROOT release")
    parser.add_argument("--tables", type=Path, required=True, help="ssc40_v1 response table directory")
    parser.add_argument("--archive", type=Path, required=True, help="discovery/data/spectra.tar.gz")
    parser.add_argument("--output", type=Path, help="Write NPZ here (defaults to stdout)")
    args = parser.parse_args()
    run(args.dataset, args.tables, args.archive, args.output or sys.stdout.buffer)
