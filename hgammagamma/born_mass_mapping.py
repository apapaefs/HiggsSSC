"""Audit an LO gg->gamma gamma LHE/raw-ROOT/response Born-mass association.

No generator or nominal output is modified. A unique ordered four-vector match
allows missing shower events, bounded by the independently audited skip count.
All normalization denominators remain the full consumed source populations.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import gzip
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np

if str(Path(__file__).resolve().parents[1]) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hgammagamma import ho_shower_completion as completion


def lhe_photons(path):
    momenta, masses, intervals, weights = [], [], [], []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as stream:
        root = None
        for action, node in ET.iterparse(stream, events=("start", "end")):
            if root is None:
                root = node
            if action != "end" or node.tag != "event":
                continue
            rows = [line.split() for line in (node.text or "").splitlines() if line.strip() and not line.lstrip().startswith("#")]
            particles = rows[1:1+int(rows[0][0])]
            incoming = [p for p in particles if int(p[1]) == -1]
            photons = [p for p in particles if int(p[0]) == 22 and int(p[1]) == 1]
            if len(particles) != 4 or len(photons) != 2 or [int(p[0]) for p in incoming] != [21, 21]:
                raise ValueError("mapping requires an exclusive LO gg->gamma gamma source")
            momenta.append([[float(p[i]) for i in (9, 6, 7, 8)] for p in photons])
            with localcontext() as ctx:
                ctx.prec = 50
                total = [sum(Decimal(p[i]) for p in photons) for i in (9, 6, 7, 8)]
                rounding = [sum(Decimal(1).scaleb(Decimal(p[i]).as_tuple().exponent) / 2
                                for p in photons) for i in (9, 6, 7, 8)]
                squared = total[0]**2 - sum(x*x for x in total[1:])
                if squared <= 0:
                    raise ValueError("nonpositive source diphoton invariant mass")
                masses.append(float(squared.sqrt()))
                ranges = [(q-d, q+d) for q,d in zip(total, rounding)]
                def square_range(bounds):
                    lo, hi = bounds
                    return (Decimal(0) if lo <= 0 <= hi else min(lo*lo, hi*hi), max(lo*lo, hi*hi))
                squares = [square_range(bounds) for bounds in ranges]
                lower = max(Decimal(0), squares[0][0]-sum(item[1] for item in squares[1:]))
                upper = max(Decimal(0), squares[0][1]-sum(item[0] for item in squares[1:]))
                intervals.append([float(lower.sqrt()), float(upper.sqrt())])
            weights.append(float(rows[0][2]))
            root.clear()
    if not momenta:
        raise ValueError("empty source")
    return np.asarray(momenta), np.asarray(masses), np.asarray(weights), np.asarray(intervals)


def match_subsequence(lhe, raw, discarded, *, atol=1.e-6, rtol=1.e-10):
    """Require a unique ordered match, including possible leading/trailing skips."""
    if len(lhe) - len(raw) != discarded or discarded < 0:
        raise ValueError("source/raw population differs from audited shower losses")
    cursor, used_skips, mapping, errors = 0, 0, [], []
    for vector in raw:
        stop = min(len(lhe), cursor + discarded - used_skips + 1)
        candidates = lhe[cursor:stop]
        tolerance = atol + rtol * np.abs(candidates)
        direct = np.all(np.abs(candidates-vector) <= tolerance, axis=(1, 2))
        swapped = np.all(np.abs(candidates-vector[::-1]) <= tolerance, axis=(1, 2))
        matches = np.flatnonzero(direct | swapped)
        if len(matches) != 1:
            raise ValueError(f"raw event {len(mapping)} has {len(matches)} possible LHE matches; no row join is safe")
        local = int(matches[0]); index = cursor + local
        error = min(np.max(np.abs(lhe[index]-vector)), np.max(np.abs(lhe[index]-vector[::-1])))
        errors.append(float(error)); mapping.append(index)
        used_skips += local; cursor = index + 1
    used_skips += len(lhe) - cursor
    if used_skips != discarded:
        raise ValueError("matched skips differ from completion proof")
    return np.asarray(mapping, dtype=np.int64), np.asarray(errors)


def region(mass, low, high, guard):
    if not low < high or guard < 0:
        raise ValueError("invalid partition window")
    ambiguous = np.flatnonzero(np.minimum(np.abs(mass-low), np.abs(mass-high)) <= guard)
    if len(ambiguous):
        raise ValueError(f"Born-window boundary is within serialization guard for source indices {ambiguous[:10].tolist()}")
    return (mass >= low) & (mass < high)


def audit(manifest_path, response_path, output, *, low=100., high=150., guard=1.e-4):
    import uproot
    if output.exists():
        raise ValueError("use a fresh diagnostic output directory")
    manifest_path, response_path = manifest_path.resolve(), response_path.resolve()
    manifest = json.loads(manifest_path.read_text())
    proof = completion.validate_completion(manifest["shower_completion"], manifest,
                                           directory=manifest_path.parent)
    source = Path(manifest["lhe_file"]).resolve()
    chain_path = Path(manifest["analysis"]["input"]).resolve()
    expected_response = chain_path.with_name(chain_path.stem + "-" + manifest["run_tag"] + "_var.root")
    if response_path != expected_response:
        raise ValueError("response path differs from the recorded analysis input and run tag")
    roots = [Path(line.strip()).resolve() for line in chain_path.read_text().splitlines() if line.strip()]
    if [str(path) for path in roots] != [row["path"] for row in proof["root_metadata"]]:
        raise ValueError("response raw-ROOT chain order differs from audited inventory")
    inputs = {"manifest": completion._artifact(manifest_path), "lhe": completion._artifact(source),
              "raw_chain": completion._artifact(chain_path), "response": completion._artifact(response_path),
              "code": completion._artifact(Path(__file__))}
    lhe, born_mass, weights, mass_intervals = lhe_photons(source)
    if len(lhe) != proof["attempted_events"] or proof["attempted_events"] != proof["lhe_events"]:
        raise ValueError("mapping requires the full finite LHE to have been consumed")
    unit = weights / np.max(np.abs(weights))
    if not np.allclose(unit, 1., rtol=1.e-12, atol=1.e-12):
        raise ValueError("this first mapping utility requires positive unit source weights")
    raw, raw_weights = [], []
    for path in roots:
        with uproot.open(path) as file:
            rows = file["Data"].arrays(["numoutgoing", "partons", "evweight"], library="np")
        if not np.all(rows["numoutgoing"] == 2) or not np.all(rows["partons"][:, 4, :2] == 22):
            raise ValueError("saved hard subprocess is not exactly two photons")
        raw.extend(rows["partons"][:, :4, :2].transpose(0, 2, 1))
        raw_weights.extend(rows["evweight"].reshape(-1))
    raw, raw_weights = np.asarray(raw), np.asarray(raw_weights)
    if len(raw) != proof["saved_events"] or not np.allclose(raw_weights, 1., rtol=1.e-10, atol=1.e-12):
        raise ValueError("raw ROOT count/weights differ from completion and unit-source proof")
    mapping, errors = match_subsequence(lhe, raw, proof["discarded_events"])
    lost = np.setdiff1d(np.arange(len(lhe)), mapping)
    for edge in (low, high):
        if np.any((mass_intervals[:, 0]-guard <= edge) & (edge <= mass_intervals[:, 1]+guard)):
            raise ValueError("Born-window boundary overlaps the LHE numeric-rounding interval")
    focus = region(born_mass, low, high, guard)
    with uproot.open(response_path) as file:
        rows = file["Data2"].arrays(["sourceevent", "generatorweight", "eventweight", "responseweight"], library="np")
    ids = rows["sourceevent"].reshape(-1)
    if ids.dtype.kind not in "iu" or not np.array_equal(np.unique(ids), np.arange(len(raw))):
        raise ValueError("response sourceevent does not cover each saved raw row exactly as an index")
    probability = rows["responseweight"].reshape(-1)
    scale = float(manifest["weight_scale"])
    if (np.any(probability < 0) or np.any(probability > 1 + 1.e-12)
            or not np.allclose(np.bincount(ids, weights=probability, minlength=len(raw)), 1., atol=1.e-10, rtol=1.e-10)
            or not np.allclose(rows["generatorweight"].reshape(-1), raw_weights[ids]*scale, atol=1.e-10, rtol=1.e-10)
            or not np.allclose(rows["eventweight"].reshape(-1), raw_weights[ids]*scale*probability, atol=1.e-10, rtol=1.e-10)):
        raise ValueError("response hypotheses do not close independently to their indexed raw source events")
    # A read-only in-memory deletion test exercises leading, internal and final
    # skips even when the physical smoke happened to have no shower failures.
    deleted = np.unique([0, len(raw)//2, len(raw)-1])
    keep = np.setdiff1d(np.arange(len(raw)), deleted)
    recovered, _ = match_subsequence(lhe, raw[keep], proof["discarded_events"]+len(deleted))
    if not np.array_equal(recovered, mapping[keep]):
        raise ValueError("synthetic shower-skip recovery failed")
    total = raw.sum(axis=1)
    raw_mass = np.sqrt(np.maximum(0., total[:, 0]**2 - (total[:, 1:]**2).sum(axis=1)))
    record = {"schema_version": 1, "status": "validated", "sample": manifest["sample"],
        "run_tag": manifest["run_tag"], "inputs": inputs,
        "raw_roots": proof["artifacts"]["roots"], "attempted_events": len(lhe), "saved_events": len(raw),
        "discarded_events": len(lost), "discarded_lhe_indices": lost.tolist(),
        "mapping": "unique ordered two-photon four-vector match; response sourceevent indexes saved raw ROOT rows",
        "momentum_match_atol_gev": 1.e-6, "momentum_match_rtol": 1.e-10,
        "maximum_absolute_fourvector_difference_gev": float(np.max(errors)),
        "maximum_raw_vs_lhe_mass_difference_gev": float(np.max(np.abs(raw_mass-born_mass[mapping]))),
        "mass_calculation": "50-digit Decimal arithmetic on the original LHE text, then converted to float",
        "maximum_lhe_mass_rounding_interval_width_gev": float(np.max(mass_intervals[:, 1]-mass_intervals[:, 0])),
        "partition": {"low_gev": low, "high_gev": high, "high_edge_exclusive": True,
                      "serialization_boundary_guard_gev": guard,
                      "attempted_inside": int(np.count_nonzero(focus)),
                      "saved_inside": int(np.count_nonzero(focus[mapping])),
                      "discarded_inside": int(np.count_nonzero(focus[lost]))},
        "response_rows": len(ids), "response_source_index_coverage": "all saved raw rows",
        "synthetic_deleted_raw_rows": deleted.tolist(), "synthetic_skip_recovery_passed": True,
        "normalization_denominator": proof["source_weights"]["sum_weight"],
        "combination_rule": "Use the original broad sigma/denominator with the outside-window response mask, plus the native focused sigma/denominator with the inside-window mask. Keep all consumed source records in both denominators; never normalize the broad outside subset as a full-rate sample.",
        "requirements_before_combining": ["Validate each actual completed sample independently, including real shower skips.",
            "Prove the focused generation differs only by the stated Born mass interval and uses an independent seed/source.",
            "Apply the same Born-window convention to both samples; keep reconstructed-mass and frozen-score selections unchanged.",
            "Combine independent per-source response variances; do not add the full broad prediction and the focused prediction."]}
    for item in inputs.values():
        completion._validate_artifact(item, verify_file=True)
    output.mkdir(parents=True)
    np.savez_compressed(output / "mapping.npz", raw_sourceevent=np.arange(len(raw)),
                        lhe_event_index=mapping, born_mgg_gev=born_mass[mapping],
                        born_mgg_rounding_interval_gev=mass_intervals[mapping],
                        inside_window=focus[mapping], response_sourceevent=ids,
                        response_lhe_event_index=mapping[ids], response_born_mgg_gev=born_mass[mapping[ids]],
                        response_inside_window=focus[mapping[ids]], discarded_lhe_event_index=lost)
    record["mapping_artifact"] = completion._artifact(output / "mapping.npz")
    (output / "summary.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--response-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.manifest, args.response_root, args.output_dir.resolve())
    print(json.dumps({key: result[key] for key in ("status", "attempted_events", "saved_events", "discarded_events", "partition")}))
