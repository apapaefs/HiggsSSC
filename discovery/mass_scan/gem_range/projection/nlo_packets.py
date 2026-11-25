"""Project the prior and final prompt-diphoton and photon+jet response trees.

Run on Timur with --release pointing at the November 23 student discovery release.
The numeric packet archive and all source trees are read-only inputs.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import tarfile

import numpy as np
from scipy.sparse import coo_matrix
import uproot


EDGES = np.arange(70., 171.)
OVERLAP = slice(35, 90)
SAMPLES = ("bkg_prompt_aa", "bkg_gamma_j")


def project_tree(path, cross_section_pb, parent_events, denominator, shard_events):
    ids, masses, weights = [], [], []
    with uproot.open(path) as source:
        for arrays in source["Data2"].iterate(["variables", "sourceevent", "eventweight"],
                                              step_size="64 MB", library="np", how=dict):
            v = arrays["variables"]
            eligible = v[:, 9] >= 2
            if not (np.isfinite(v[eligible, 0]).all() and np.isfinite(arrays["eventweight"]).all()):
                raise ValueError(f"Nonfinite mass or weight in {path}")
            keep = eligible & (v[:, 0] >= EDGES[0]) & (v[:, 0] < EDGES[-1])
            masses.append(v[keep, 0])
            ids.append(arrays["sourceevent"].reshape(-1)[keep])
            weights.append(arrays["eventweight"].reshape(-1)[keep])
    if not ids:
        return np.zeros(100), np.zeros((100, 100))
    sources, rows = np.unique(np.concatenate(ids), return_inverse=True)
    bins = np.searchsorted(EDGES, np.concatenate(masses), side="right") - 1
    raw = np.concatenate(weights)
    if shard_events is None:
        scaled = raw * 1000 * cross_section_pb / denominator
    else:
        scaled = raw * 1000 * cross_section_pb * parent_events / (denominator * shard_events)
    matrix = coo_matrix((scaled, (rows, bins)), shape=(len(sources), len(EDGES) - 1)).tocsr()
    return np.asarray(matrix.sum(axis=0)).ravel(), (matrix.T @ matrix).toarray()


def project_packet(task):
    sources, name, packet = task
    sample, n = packet["sample"], packet["generated_events"]
    origin = packet["provenance"]
    if "shard_directory" not in origin:
        run = Path(origin["path"]).parts[-4]
        roots = [(Path(path).name, digest) for path, digest in origin["sha256"].items()
                 if path.endswith("_var.root")]
        ((filename, digest),) = roots
        tree = sources / "earlier_ho" / run / sample / filename
        hasher = hashlib.sha256()
        with tree.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                hasher.update(chunk)
        if hasher.hexdigest() != digest:
            raise ValueError(f"Tree hash failed for {name}")
        rate, cov = project_tree(tree, packet["cross_section_pb"], n,
                                 packet["normalization_sum_weight"], None)
        ntrees = 1
    else:
        suffix = Path(origin["shard_directory"]).name
        folder = (sources / "sideband_20pct_20261120/nlo/parents" if name.startswith("packets/final/")
                  else sources / "sideband-20261119/nlo") / suffix
        receipts = sorted(folder.glob("shard_*/receipt.json"))
        ntrees = len(receipts)
        if ntrees != origin["complete_shards"] or not ntrees:
            raise ValueError(f"Incomplete shower receipts in {folder}")
        covered = signed = 0
        rate, cov = np.zeros(100), np.zeros((100, 100))
        denom = origin["source_signed_denominator"]
        for receipt_path in receipts:
            receipt = json.loads(receipt_path.read_bytes())
            size = receipt["generated_events"]
            if (receipt["first_source"] != covered or receipt["last_source_exclusive"] != covered + size
                    or receipt["parent_events"] != n or receipt["global_denominator"] != denom):
                raise ValueError(f"Inconsistent source coverage in {receipt_path}")
            tree = receipt_path.with_name("ssc-side_var.root")
            shard_rate, shard_cov = project_tree(tree, origin["native_cross_section_pb"], n, denom, size)
            rate += shard_rate * (size / n)
            cov += shard_cov * (size / n) ** 2
            covered += size
            signed += receipt["signed_source_sum"]
        if covered != n or signed != denom:
            raise ValueError(f"Wrong generated or signed event count for {name}")
    np.testing.assert_array_equal(EDGES[35:91], packet["edges_gev"])
    dr = rate[OVERLAP] - packet["rate_fb"]
    dc = cov[OVERLAP, OVERLAP] - packet["covariance_fb2"]
    for derived, published in ((rate[OVERLAP], packet["rate_fb"]),
                               (cov[OVERLAP, OVERLAP], packet["covariance_fb2"])):
        np.testing.assert_allclose(derived, published, atol=1e-11, rtol=1e-11, err_msg=name)
    return sample, n, rate, cov, dict(packet=name, trees=ntrees, generated_events=n,
        independent_id=packet["independent_id"], max_overlap_rate_difference_fb=float(abs(dr).max()),
        max_overlap_covariance_difference_fb2=float(abs(dc).max()))


def run(release, output, workers, benchmark):
    archive = release / "analysis/discovery/data/spectra.tar.gz"
    tasks = []
    with tarfile.open(archive) as bundle:
        for entry in bundle:
            if entry.name.startswith(("packets/prior/", "packets/final/")):
                packet = json.load(bundle.extractfile(entry))
                if packet["sample"] in SAMPLES:
                    tasks.append((release / "sources", entry.name, packet))
    if benchmark:
        tasks = [t for t in tasks if t[1].startswith("packets/final/")][::40][:12]
    identities = [p["independent_id"] for _, _, p in tasks]
    if len(set(identities)) != len(tasks):
        raise ValueError("Duplicate NLO packet")
    totals = {s: sum(p["generated_events"] for _, _, p in tasks if p["sample"] == s) for s in SAMPLES}
    rates = np.zeros((len(SAMPLES), 100))
    covariances = np.zeros((len(SAMPLES), 100, 100))
    details = {s: [] for s in SAMPLES}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for j, (sample, n, rate, covariance, meta) in enumerate(pool.map(project_packet, tasks), 1):
            index = SAMPLES.index(sample)
            fraction = n / totals[sample]
            rates[index] += rate * fraction
            covariances[index] += covariance * fraction ** 2
            details[sample].append(meta)
            if j % 50 == 0:
                print(f"Verified {j}/{len(tasks)} NLO packets", file=sys.stderr, flush=True)
    for s, packets in details.items():
        print(f"{s}: {len(packets)} packets, {sum(x['trees'] for x in packets)} trees, "
              f"{totals[s]} generated events", file=sys.stderr)
    if not benchmark:
        metadata = [dict(sample=s, packets=details[s], stratum="full",
                    archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest()) for s in SAMPLES]
        np.savez_compressed(output, edges_gev=EDGES, sample=SAMPLES, stratum=["full"] * len(SAMPLES),
            generated_events=[totals[s] for s in SAMPLES], rate_fb=rates, covariance_fb2=covariances,
            metadata_json=json.dumps(metadata))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, required=True)
    parser.add_argument("--output", type=Path, help="Destination NPZ, defaults to stdout")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--benchmark", action="store_true", help="Time twelve packets without writing output")
    args = parser.parse_args()
    run(args.release, args.output or sys.stdout.buffer, args.workers, args.benchmark)
