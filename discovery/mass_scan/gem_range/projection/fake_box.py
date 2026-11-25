"""Project every prior and final dijet and box packet in the frozen SSC campaign."""

import argparse
from collections import defaultdict
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

import numpy as np
from scipy.sparse import coo_matrix
import uproot


EDGES = np.arange(70., 171.)
NBINS = len(EDGES) - 1


def digest(data):
    return hashlib.sha256(data).hexdigest()


def verified(path, expected):
    data = path.read_bytes()
    if digest(data) != expected:
        raise ValueError(f"Source hash mismatch for {path}")
    return data


def project(masses, ids, weights):
    if not (np.isfinite(masses).all() and np.isfinite(weights).all()):
        raise ValueError("Nonfinite response")
    keep = (masses >= EDGES[0]) & (masses < EDGES[-1])
    bins = np.searchsorted(EDGES, masses[keep], side="right") - 1
    sources, rows = np.unique(ids[keep], return_inverse=True)
    response = coo_matrix((weights[keep], (rows, bins)), shape=(len(sources), NBINS)).tocsr()
    return np.asarray(response.sum(axis=0)).ravel(), (response.T @ response).toarray()


def project_dijet(packet, directory, checksums):
    rate = np.zeros(NBINS)
    covariance = np.zeros((NBINS, NBINS))
    count = packet["generated_events"]
    if sum(shard["events"] for shard in packet["shards"]) != count:
        raise ValueError("Dijet shard counts do not sum to the packet denominator")
    for shard in packet["shards"]:
        spectrum = directory / shard["path"]
        receipt = json.loads(verified(spectrum, shard["sha256"]))
        if receipt["generated_events"] != shard["events"] or receipt["generator_seed"] != shard["seed"]:
            raise ValueError(f"Inconsistent dijet shard identity at {spectrum}")
        segments = [(spectrum.parent, receipt)]
        if "components" in receipt:
            segments = []
            for component in receipt["components"].values():
                path = directory / component["path"]
                segment = json.loads(verified(path, component["sha256"]))
                if segment["generated_events"] != component["size"]:
                    raise ValueError(f"Inconsistent repaired segment at {path}")
                segments.append((path.parent, segment))
            n_segments = sum(segment["generated_events"] for _, segment in segments)
            if n_segments + receipt["zero_response_native_aborts"] != receipt["generated_events"]:
                raise ValueError(f"Incomplete set of repaired segments at {spectrum}")
        for folder, segment in segments:
            path = folder / "pairs.csv"
            data = verified(path, checksums[path])
            columns = np.loadtxt(io.BytesIO(data), delimiter=",", skiprows=1, usecols=(0, 4, 6), ndmin=2)
            ids, fq, masses = columns.T
            if np.any(ids != np.floor(ids)) or np.any(ids < 0) or np.any(ids >= segment["generated_events"]):
                raise ValueError(f"Invalid pair source IDs at {path}")
            # Each independent segment has its own inclusive rate in nb. The pooling
            # weight n_segment/n_packet cancels the segment denominator.
            weights = fq**2 * segment["sigma_nb"] * 1e6 / count
            partial_rate, partial_cov = project(masses, ids, weights)
            rate += partial_rate
            covariance += partial_cov
    return rate, covariance


def broad_mask(packet, summary, ids, by_digest, repository):
    # The public source archive includes the validated skip list and LHE, while the
    # original working directory also contains the larger per-response mapping.
    if (summary["status"] != "validated" or not summary["synthetic_skip_recovery_passed"]
            or summary["inputs"]["response"]["sha256"] != packet["response_root_sha256"]
            or summary["inputs"]["lhe"]["sha256"] != packet["source_lhe_sha256"]):
        raise ValueError("Born association belongs to different events or has not passed validation")
    sys.path.insert(0, str(repository))
    try:
        from hgammagamma.born_mass_mapping import lhe_photons, region
    finally:
        sys.path.pop(0)
    path = by_digest[packet["source_lhe_sha256"]][0]
    verified(path, packet["source_lhe_sha256"])
    _, masses, _, intervals = lhe_photons(path)
    partition = summary["partition"]
    low, high, guard = (partition[k] for k in ("low_gev", "high_gev", "serialization_boundary_guard_gev"))
    if (low, high, partition["high_edge_exclusive"]) != (100., 150., True):
        raise ValueError("Unknown broad/focused partition")
    for edge in (low, high):
        if np.any((intervals[:, 0] - guard <= edge) & (edge <= intervals[:, 1] + guard)):
            raise ValueError("Born mass is ambiguous at the partition boundary")
    inside = region(masses, low, high, guard)
    # The archived matching audit established that the shower preserves LHE order,
    # with exactly these skipped hard collisions. Detector IDs index the saved rows.
    saved = np.delete(inside, summary["discarded_lhe_indices"])
    if (len(inside), inside.sum(), len(saved), saved.sum()) != (
            summary["attempted_events"], partition["attempted_inside"],
            summary["saved_events"], partition["saved_inside"]):
        raise ValueError("Reconstructed LHE partition disagrees with the original mapping")
    np.testing.assert_array_equal(np.unique(ids), np.arange(len(saved)))
    return ~saved[ids]


def project_box(packet, by_digest, repository):
    response_sha = packet["response_root_sha256"]
    path = by_digest[response_sha][0]
    with uproot.open(io.BytesIO(verified(path, response_sha))) as root:
        tree = root["Data2"]
        values = tree["variables"].array(library="np")
        ids = tree["sourceevent"].array(library="np").reshape(-1)
        weights = tree["eventweight"].array(library="np").reshape(-1)
    keep = values[:, -1] >= 2
    if packet["stratum"] == "broad_outside":
        summary_sha = packet["born_map_sha256"]
        summary_path = by_digest[summary_sha][0]
        summary = json.loads(verified(summary_path, summary_sha))
        keep &= broad_mask(packet, summary, ids, by_digest, repository)
    elif packet["stratum"] != "focused_inside":
        raise ValueError("Unexpected box phase-space stratum")
    weights = weights * (packet["native_cross_section_pb"] * 1000 / packet["generated_events"])
    return project(values[keep, 0], ids[keep], weights[keep])


def run(release, archive, output, samples):
    checksums = {}
    by_digest = defaultdict(list)
    listing = (release / "source-files.sha256").read_bytes()
    for line in listing.decode().splitlines():
        sha, name = line.split("  ", 1)
        path = release / name
        checksums[path] = sha
        by_digest[sha].append(path)
    rates, covariances, metadata = [], [], []
    with tarfile.open(archive) as tar:
        for member in tar:
            if not member.isfile() or not member.name.startswith(("packets/prior/", "packets/final/")):
                continue
            data = tar.extractfile(member).read()
            packet = json.loads(data)
            sample = packet["sample"]
            if sample not in samples:
                continue
            if sample == "bkg_dijet_double_fake":
                campaign = ("sideband_campaign_20261119/dijet" if member.name.startswith("packets/prior/")
                            else "sideband_20pct_20261120/fake_loop/dijet")
                rate, covariance = project_dijet(packet, release / "sources" / campaign, checksums)
            elif sample == "bkg_loop_gg_aa":
                rate, covariance = project_box(packet, by_digest, release / "analysis")
            else:
                continue
            edges = packet.get("edges_gev", packet.get("mass_bin_edges_GeV"))
            low = int(np.searchsorted(EDGES, edges[0]))
            high = low + len(edges) - 1
            np.testing.assert_array_equal(EDGES[low:high+1], edges)
            old_rate = packet.get("rate_fb", packet.get("rates_fb"))
            old_cov = packet["covariance_fb2"]
            np.testing.assert_allclose(rate[low:high], old_rate, rtol=1e-9, atol=1e-9, err_msg=member.name)
            np.testing.assert_allclose(covariance[low:high, low:high], old_cov, rtol=1e-9, atol=1e-9, err_msg=member.name)
            metadata.append(dict(sample=sample, stratum=packet.get("stratum", "full"), packet=member.name,
                packet_sha256=digest(data), independent_id=packet["independent_id"],
                generated_events=packet["generated_events"],
                max_overlap_rate_difference_fb=float(np.max(np.abs(rate[low:high] - old_rate))),
                max_overlap_covariance_difference_fb2=float(np.max(np.abs(covariance[low:high,low:high] - old_cov)))))
            rates.append(rate)
            covariances.append(covariance)
    if not rates:
        raise ValueError("No matching packets in the archive")
    np.savez_compressed(output, edges_gev=EDGES, rate_fb=rates, covariance_fb2=covariances,
        sample=[x["sample"] for x in metadata], stratum=[x["stratum"] for x in metadata],
        generated_events=[x["generated_events"] for x in metadata], metadata_json=json.dumps(metadata),
        source_index_sha256=digest(listing))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, required=True, help="Timur student discovery source release")
    parser.add_argument("--archive", type=Path, required=True, help="discovery/data/spectra.tar.gz")
    parser.add_argument("--samples", nargs="+", default=["bkg_dijet_double_fake", "bkg_loop_gg_aa"],
                        choices=["bkg_dijet_double_fake", "bkg_loop_gg_aa"])
    parser.add_argument("--output", type=Path, help="Write NPZ here (defaults to stdout)")
    args = parser.parse_args()
    run(args.release, args.archive, args.output or sys.stdout.buffer, args.samples)
