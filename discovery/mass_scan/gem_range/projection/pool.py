"""Pool the wide diphoton inputs and verify the 105–160 GeV archive overlap.

Run from any directory with the discovery dependencies. The three local NPZ
inputs are the original dataset projection, the fake and box packet projections,
and the two grouped NLO packet projections. The signed-off archive defines
which physical exposures belong to each stratum.
"""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

HERE = Path(__file__).resolve().parent
DISCOVERY = HERE.parents[2]
sys.path.insert(0, str(DISCOVERY))
from input_data import BACKGROUNDS, DATA, SIGNAL, load_archive, pool
from likelihood import EDGES as ARCHIVE_EDGES

SAMPLES = (SIGNAL, *BACKGROUNDS)
FILES = ("baseline.npz", "fake_box.npz", "nlo_packets.npz")
EDGES = np.arange(70., 171.)
OVERLAP = slice(int(ARCHIVE_EDGES[0] - EDGES[0]), int(ARCHIVE_EDGES[-1] - EDGES[0]))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(derived, archived, label):
    np.testing.assert_allclose(derived, archived, rtol=1e-9, atol=1e-9, err_msg=label)
    return float(np.abs(derived - archived).max())


def archive_overlap(members, count):
    rate = np.zeros(len(ARCHIVE_EDGES) - 1)
    cov = np.zeros((rate.size, rate.size))
    for item in members:
        weight = item["generated_events"] / count
        rate += weight * np.asarray(item.get("rate_fb", item.get("rates_fb")))
        cov += weight**2 * np.asarray(item["covariance_fb2"])
    return rate, cov


def run(input_dir, output_dir):
    absent = [name for name in FILES if not (input_dir / name).is_file()]
    if absent:
        raise SystemExit("Complete projections still required: " + ", ".join(absent))

    objects, packets = load_archive()
    ref_rates, ref_covs, ref_counts, _ = pool(objects, packets)
    archive_hash = sha256(DATA / "spectra.tar.gz")
    canonical = {p: objects[p] for p in packets}
    packet_hashes = {p: src["sha256"] for p, src in zip(
        packets, objects["records/final_fit.json"]["independent_new_inputs"])}
    for item in objects["baseline.json"]["datasets"]:
        canonical["baseline/" + item["independent_id"]] = item
    if len(canonical) != len(packets) + len(objects["baseline.json"]["datasets"]):
        raise ValueError("Repeated baseline independent ID")

    strata = defaultdict(list)
    used = set()
    stage_counts = defaultdict(lambda: defaultdict(lambda: {"inputs": 0, "generated_events": 0}))
    projection_audit = {}
    for filename in FILES:
        path = input_dir / filename
        row_checks, packet_checks = [], []
        with np.load(path, allow_pickle=False) as f:
            np.testing.assert_array_equal(f["edges_gev"], EDGES)
            sample, stratum, counts = f["sample"], f["stratum"], f["generated_events"]
            rates, covs = f["rate_fb"], f["covariance_fb2"]
            metadata = json.loads(f["metadata_json"].item())
            if (counts.dtype.kind not in "iu" or counts.shape != (len(metadata),)
                    or sample.shape != counts.shape or stratum.shape != counts.shape
                    or np.any(counts <= 0) or rates.shape != (len(metadata), 100)
                    or covs.shape != (len(metadata), 100, 100)
                    or not np.isfinite(rates).all() or not np.isfinite(covs).all()):
                raise ValueError(f"Invalid event counts or spectra in {filename}")
            np.testing.assert_allclose(covs, covs.transpose(0, 2, 1), atol=1e-10, rtol=1e-12)
            eig = np.linalg.eigvalsh(covs)
            if np.any(eig[:, 0] < -1e-9 * np.maximum(1., eig[:, -1])):
                raise ValueError(f"Invalid covariance in {filename}")
            for s, t, n, rate, cov, record in zip(sample, stratum, counts, rates, covs, metadata):
                if (record.get("sample", s), record.get("stratum", t)) != (s, t):
                    raise ValueError(f"Row labels differ from metadata in {filename}")
                if record.get("archive_sha256", archive_hash) != archive_hash:
                    raise ValueError(f"Archive digest mismatch in {filename}")
                if filename == "baseline.npz":
                    members = [("baseline/" + record["independent_id"], record)]
                elif filename == "nlo_packets.npz":
                    members = [(p["packet"], p) for p in record["packets"]]
                else:
                    members = [(record["packet"], record)]
                rows = []
                for key, provenance in members:
                    if key in used or key not in canonical:
                        raise ValueError(f"Unknown or duplicate archive exposure {key}")
                    used.add(key)
                    entry = canonical[key]
                    if ((entry["sample"], entry.get("stratum", "full")) != (s, t)
                            or entry["independent_id"] != provenance["independent_id"]
                            or entry["generated_events"] != provenance.get("generated_events", n)
                            or ("packet_sha256" in provenance and provenance["packet_sha256"] != packet_hashes[key])):
                        raise ValueError(f"Exposure identity or generated count differs for {key}")
                    check = [provenance["max_overlap_rate_difference_fb"],
                             provenance["max_overlap_covariance_difference_fb2"]]
                    archived_arrays = (entry.get("rate_fb", entry.get("rates_fb")), entry["covariance_fb2"])
                    for difference, archived in zip(check, archived_arrays):
                        if not 0 <= difference <= 1e-9 * (1 + np.max(np.abs(archived))):
                            raise ValueError(f"Invalid individual overlap check for {key}")
                    packet_checks.append(check)
                    rows.append(entry)
                    stage = "baseline" if filename == "baseline.npz" else key.split("/")[1]
                    tally = stage_counts[f"{s}/{t}"][stage]
                    tally["inputs"] += 1
                    tally["generated_events"] += entry["generated_events"]
                if n != sum(r["generated_events"] for r in rows):
                    raise ValueError(f"Generated event total differs for {s}/{t} in {filename}")
                ref_rate, ref_cov = archive_overlap(rows, n)
                row_checks.append((compare(rate[OVERLAP], ref_rate, filename),
                                   compare(cov[OVERLAP, OVERLAP], ref_cov, filename)))
                strata[(s, t)].append((int(n), rate, cov))
        max_rate, max_cov = np.max(row_checks, axis=0)
        max_packet_rate, max_packet_cov = np.max(packet_checks, axis=0)
        projection_audit[filename] = dict(sha256=sha256(path), npz_rows=len(row_checks),
            archive_inputs=len(packet_checks), max_row_overlap_rate_difference_fb=max_rate,
            max_row_overlap_covariance_difference_fb2=max_cov,
            max_individual_overlap_rate_difference_fb=max_packet_rate,
            max_individual_overlap_covariance_difference_fb2=max_packet_cov)

    if missing := canonical.keys() - used:
        raise ValueError(f"Missing {len(missing)} archive inputs, including {sorted(missing)[:3]}")

    rates, covs = np.zeros((len(SAMPLES), 100)), np.zeros((len(SAMPLES), 100, 100))
    totals = {}
    for (sample, stratum), members in strata.items():
        index = SAMPLES.index(sample)
        n = sum(row[0] for row in members)
        totals[f"{sample}/{stratum}"] = n
        for size, rate, cov in members:
            fraction = size / n
            rates[index] += fraction * rate
            covs[index] += fraction**2 * cov
    if totals != ref_counts:
        raise ValueError("Pooled exposures differ by sample or stratum")
    per_sample = {}
    for i, sample in enumerate(SAMPLES):
        per_sample[sample] = dict(
            max_rate_difference_fb=compare(rates[i, OVERLAP], ref_rates[sample], sample),
            max_covariance_difference_fb2=compare(covs[i, OVERLAP, OVERLAP], ref_covs[sample], sample))
    audit = dict(archive_sha256=archive_hash, archive_packets=len(packets),
                 baseline_strata=len(canonical) - len(packets), generated_events=sum(totals.values()),
                 exposures_by_stratum_and_stage=stage_counts, source_projections=projection_audit,
                 reference_overlap_gev=[int(ARCHIVE_EDGES[0]), int(ARCHIVE_EDGES[-1])],
                 pooled_overlap_by_sample=per_sample)
    output_dir.mkdir(exist_ok=True, parents=True)
    with TemporaryDirectory(dir=output_dir) as temp:
        tmp = Path(temp)
        np.savez_compressed(tmp / "pooled.npz", edges_gev=EDGES, sample=SAMPLES, rate_fb=rates,
                            covariance_fb2=covs, metadata_json=json.dumps(audit, sort_keys=True))
        with (tmp / "rates.csv").open("w", newline="") as stream:
            table = csv.writer(stream, lineterminator="\n")
            table.writerow(["low_gev", "high_gev", *SAMPLES, "background_fb", "background_mc_sigma_fb"])
            table.writerows(zip(EDGES[:-1], EDGES[1:], *rates, rates[1:].sum(axis=0),
                                np.sqrt(covs[1:].sum(axis=0).diagonal())))
        (tmp / "audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True, allow_nan=False) + "\n")
        for name in ("pooled.npz", "rates.csv", "audit.json"):
            (tmp / name).replace(output_dir / name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=HERE)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    run(args.input_dir, args.output_dir)
