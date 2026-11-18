"""Diagnostic native-decay control/importance-sampling pilot on disjoint LHE ranges.

Reuse a completed inclusive-heavy-flavour-v2 hard campaign. The two arms never
share source records; both load exactly the same native Herwig B/D decay table.
The weighted arm multiplies direct B->hadron e/mu modes by four, with complete
proposal support and independently audited physical p/q weights. Nothing here
promotes the pilot into a nominal background estimate.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, replace
import fcntl
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

if str(Path(__file__).resolve().parents[1]) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import run_reducible_controls as controls
from hfourlepton import hf_decay_audit
from hfourlepton.minimal_runtime import _decode
from hgammagamma import ho_shower_completion as completion


def code_identity():
    paths = [Path(__file__), Path(four.__file__), Path(controls.__file__),
             Path(hf_decay_audit.__file__), Path(completion.__file__),
             Path(__file__).parent / "Herwig/LHEBranchingRatioReweighter.cc"]
    return {str(path.resolve()): four.sha256_file(path) for path in paths}


def slice_lhe(source, destination, start, stop):
    """Write an exact half-open slice, preserving the native inclusive init rate."""
    source, destination = Path(source), Path(destination)
    if start < 0 or stop <= start or destination.exists():
        raise ValueError("invalid or already existing slice")
    destination.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(destination, "wt") as output:
        for kind, data in controls.lhe_parts(source):
            if kind == "event":
                if start <= count < stop:
                    output.write(data)
                count += 1
            else:
                output.write(data)
    if count < stop:
        raise ValueError("requested slice exceeds the source shard")
    return destination


def select_range(record, start, events, destination):
    stop = start + events
    if start < 0 or events < 1 or stop > record["source_events"]:
        raise ValueError("pilot range is outside its hard parent")
    result, cursor = [], start
    for shard in sorted(record["shards"], key=lambda row: row["start"]):
        lo, hi = max(start, shard["start"]), min(stop, shard["stop"])
        if lo >= hi:
            continue
        if lo != cursor or four.sha256_file(Path(shard["lhe"])) != shard["sha256"]:
            raise ValueError("parent shards are incomplete, overlapping, or changed")
        path = Path(shard["lhe"])
        if lo != shard["start"] or hi != shard["stop"]:
            path = slice_lhe(path, Path(destination) / f"{lo:09d}-{hi:09d}.lhe.gz",
                             lo - shard["start"], hi - shard["start"])
        result.append({"index": len(result), "start": lo, "stop": hi, "events": hi-lo,
                       "lhe": str(path), "sha256": four.sha256_file(path)})
        cursor = hi
    if cursor != stop:
        raise ValueError("parent inventory does not cover the complete requested range")
    return result


def worker(job):
    payload = json.loads(Path(job).read_text(), object_hook=_decode)
    if payload["code_identity"] != code_identity():
        raise ValueError("pilot code changed after job preparation")
    cfg = four.Config(**payload["config"])
    sample = next(sample for sample in four.selected_samples(cfg) if sample.name == payload["sample"])
    four.require_inputs(cfg, [sample])
    if four.sha256_file(cfg.herwig_b_decays) != payload["decay_backend"]["heavy_flavour_decay_table_sha256"]:
        raise ValueError("native decay table changed after preparation")
    planned = four.build_manifest(cfg, [sample])
    existing = four.load_resumable_manifest(cfg)
    if existing and existing["configuration_fingerprint"] != planned["configuration_fingerprint"]:
        raise ValueError("worker configuration changed")
    manifest = four.merge_manifest_samples(existing, planned)
    manifest.update(status="running", reducible_control_partition=payload["partition"],
                    weighted_pilot_diagnostic_only=True, weighted_pilot_code_identity=payload["code_identity"])
    for entry in manifest["samples"]:
        entry.update(payload["decay_backend"])
    four.write_manifest(manifest, cfg)
    original = four.run_herwig
    four.run_herwig = lambda sample, lhe, cfg: controls.run_finite_shower(sample, lhe, cfg, native_decay_table=True)
    started = time.monotonic()
    try:
        four.run_sample(sample, [sample], manifest, cfg)
        manifest["status"] = "complete"
    except BaseException as exc:
        manifest.update(status="failed", error=str(exc))
        raise
    finally:
        four.run_herwig = original
        manifest["last_attempt_wall_seconds"] = time.monotonic() - started
        four.write_manifest(manifest, cfg)


def campaign(args):
    if not 1 <= args.workers <= 40:
        raise ValueError("pilot reserves at most 40 logical workers")
    if max(args.weighted_start, args.control_start) < min(args.weighted_start + args.events_per_arm, args.control_start + args.events_per_arm):
        raise ValueError("weighted and control source ranges must be disjoint")
    parent = json.loads(args.parent_index.read_text(), object_hook=_decode)
    if parent.get("status") != "complete" or set(parent["samples"]) != set(controls.SAMPLES):
        raise ValueError("pilot requires a completed three-process hard campaign")
    cfg = replace(four.Config(**parent["identity"]["config"]), run_tag=args.run_tag,
                  nb_core=1, nevents=args.events_per_arm, herwig_b_decays=args.herwig_b_decays.resolve(),
                  heavy_flavour_bias="unbiased", postprocess_only=True)
    if not cfg.herwig_b_decays.is_file():
        raise ValueError("native Herwig decay table is missing")
    cfg.run_dir.mkdir(parents=True, exist_ok=True)
    index_path = cfg.run_dir / "reducible_controls.json"
    backend = {"heavy_flavour_decay_backend": "herwig_native_table",
               "heavy_flavour_decay_table": str(cfg.herwig_b_decays),
               "heavy_flavour_decay_table_sha256": four.sha256_file(cfg.herwig_b_decays)}
    identity = {"parent_index": completion._artifact(args.parent_index), "code": code_identity(),
                "weighted_start": args.weighted_start, "control_start": args.control_start,
                "events_per_arm": args.events_per_arm, "config": asdict(cfg), **backend}
    with (cfg.run_dir / ".weighted_pilot.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if index_path.exists():
            index = json.loads(index_path.read_text(), object_hook=_decode)
            if index["identity"] != identity:
                raise ValueError("weighted pilot identity changed; use a new tag")
        else:
            index = {"schema_version": 1, "identity": identity, "status": "prepared", "samples": {},
                     "diagnostic_only": True,
                     "normalization": "Raw pre-decay LHE denominator. Pool disjoint slices separately within each sampling arm; do not add arms or promote without closure."}
            for base in controls.SAMPLES:
                original = parent["samples"][base]
                if original.get("generation_policy", {}).get("name") != "inclusive_heavy_flavour_v2":
                    raise ValueError("weighted pilot requires the production v2 hard phase space")
                if four.sha256_file(Path(original["source_lhe"])) != original["source_sha256"]:
                    raise ValueError("hard parent LHE changed")
                for biased, start in ((True, args.weighted_start), (False, args.control_start)):
                    name = base + "_hfbiased" if biased else base
                    record = {key: value for key, value in original.items() if key != "shards"}
                    record.update(backend)
                    record.update(sampling_stratum="weighted_native" if biased else "unbiased_native", shards=[])
                    shards = select_range(original, start, args.events_per_arm, cfg.run_dir / "slices" / name)
                    for shard in shards:
                        tag = f"{args.run_tag}-{name}-{shard['index']:05d}"
                        worker_cfg = replace(cfg, run_tag=tag, nevents=shard["events"], bias_control_events=shard["events"],
                            sample_nevents={name: shard["events"]}, run_samples=base,
                            heavy_flavour_bias="semileptonic" if biased else "unbiased",
                            source_backend="external-lhe", external_lhe={base: Path(shard["lhe"])},
                            external_perturbative_order="LO", seed_base=1 + four.stable_offset(tag))
                        job = cfg.run_dir / "jobs" / f"{name}-{shard['index']:05d}.json"
                        partition = {key: record[key] for key in ("source_lhe", "source_sha256", "production_seed", "generation_policy", "sampling_stratum")}
                        partition.update(shard)
                        controls.write_json(job, {"config": asdict(worker_cfg), "sample": name,
                                                  "partition": partition, "code_identity": code_identity(),
                                                  "decay_backend": backend})
                        record["shards"].append({**shard, "run_tag": tag, "manifest": str(worker_cfg.manifest_path),
                                                "job": str(job), "log": str(job.with_suffix(".log")), "status": "prepared"})
                    index["samples"][name] = record
            controls.write_json(index_path, index)
        print(f"Pilot inventory: {index_path}", flush=True)
        if args.prepare_only:
            return
        # All mutations of shared build products finish before any worker starts.
        four.run_runtime_command(["make", "-C", cfg.analysis_build_dir, cfg.analysis_target], cfg)
        four.build_decay_reweighter(cfg)
        jobs = [shard for record in index["samples"].values() for shard in record["shards"]]
        stop = threading.Event()
        def launch(shard):
            if stop.is_set():
                return None
            with Path(shard["log"]).open("a") as out:
                result = subprocess.run([cfg.mg5_python or sys.executable, Path(__file__).resolve(), "--worker", shard["job"]],
                    stdout=out, stderr=subprocess.STDOUT,
                    env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONUNBUFFERED": "1"})
            if result.returncode:
                stop.set()
            return result.returncode
        index.update(status="showering", workers=args.workers)
        controls.write_json(index_path, index)
        failures = []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(launch, shard): shard for shard in jobs}
            for future in as_completed(futures):
                shard, code = futures[future], future.result()
                shard["status"] = "prepared" if code is None else "complete" if code == 0 else "failed"
                if code:
                    failures.append(shard["log"])
                controls.write_json(index_path, index)
                print(f"{shard['run_tag']}: {shard['status']}", flush=True)
        index["status"] = "failed" if failures else "complete"
        controls.write_json(index_path, index)
        if failures:
            raise RuntimeError(f"{len(failures)} pilot shards failed; first log: {failures[0]}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--parent-index", type=Path)
    parser.add_argument("--run-tag")
    parser.add_argument("--events-per-arm", type=int, default=100000)
    parser.add_argument("--weighted-start", type=int, default=0)
    parser.add_argument("--control-start", type=int, default=100000)
    parser.add_argument("--workers", type=int, default=40)
    parser.add_argument("--herwig-b-decays", type=Path, default=Path("/home/shared/Herwig/share/Herwig/defaults/HerwigBDecays.in"))
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args(argv)
    if args.worker:
        worker(args.worker)
    elif args.parent_index and args.run_tag:
        campaign(args)
    else:
        parser.error("--parent-index and --run-tag are required")


if __name__ == "__main__":
    main()
