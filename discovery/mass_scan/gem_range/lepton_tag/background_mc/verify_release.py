"""Replay the fixed nominal and wider-acceptance spectra from the two student archives."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return json.load(stream)


def extract_verified(archive, directory):
    hashes = {}
    with tarfile.open(archive, "r|gz") as source:
        for item in source:
            member = Path(item.name)
            if (not item.isfile() or member.is_absolute() or ".." in member.parts
                    or member in hashes):
                raise ValueError(f"Unexpected member in {archive.name}: {item.name}")
            output = directory / member
            output.parent.mkdir(parents=True, exist_ok=True)
            if shutil.disk_usage(directory).free <= item.size:
                raise OSError(f"Insufficient replay space for {member} in {directory}")
            digest = hashlib.sha256()
            with source.extractfile(item) as incoming, output.open("xb") as dest:
                for chunk in iter(lambda: incoming.read(1024 * 1024), b""):
                    digest.update(chunk)
                    dest.write(chunk)
            hashes[member] = digest.hexdigest()
    return hashes


def check_sources(directory, hashes, model):
    parents = {(row["integration_id"], row["source_lhe_sha256"]) for row in model["shards"]}
    found = set()
    used_software = set()
    for row in model["shards"]:
        receipt_path = Path(row["receipt"])
        receipt = load_json(directory / receipt_path)
        if (hashes[receipt_path] != row["receipt_sha256"]
                or hashes[receipt_path.with_name("packet.csv")] != row["packet_sha256"]
                or hashes[receipt_path.with_name("tagged.in")] != receipt["card_sha256"]):
            raise ValueError(f"Invalid archived packet, receipt or card at {receipt_path}")
        used_software.update([receipt["response_sha256"], receipt["plugin_sha256"]])
    binaries = {digest for path, digest in hashes.items()
                if len(path.parts) == 3 and path.parts[0] == "software" and path.parts[1] == digest}
    if not used_software.issubset(binaries):
        raise ValueError("Missing an original detector executable or shower plugin")
    raw_hashes = set(hashes.values())
    for member in hashes:
        if member.name != "hard.json":
            continue
        record = load_json(directory / member)
        parent = (record["integration_id"], record["lhe_sha256"])
        if (parent not in parents or parent in found
                or hashes[member.with_name("events.lhe.gz")] != parent[1]):
            raise ValueError(f"Invalid or unrelated hard-event source: {member}")
        if "cut_function_sha256" in record:
            if hashes[member.with_name("outside_nominal.f")] != record["cut_function_sha256"]:
                raise ValueError(f"Missing generation acceptance at {member}")
        if "gridpack_sha256" in record:
            if record["gridpack_sha256"] not in raw_hashes:
                raise ValueError(f"Missing native integration gridpack at {member}")
            if any(lot["lhe_sha256"] not in raw_hashes for lot in record["native_lots"]):
                raise ValueError(f"Missing native lot for {member}")
        found.add(parent)
    if parents != found:
        raise ValueError("Archive is missing a hard-event source used by the spectrum")
    return dict(hard_integrations=len(found), shower_shards=len(model["shards"]),
                software_sha256=sorted(used_software))


def replay(release, scratch, kind):
    core = kind == "core"
    archive = release / "raw" / ("core-inputs.tar.gz" if core else "zee-confirmation-inputs.tar.gz")
    frozen = release / ("background.json.gz" if core else "zee-confirmation-spectrum.json.gz")
    expected = load_json(frozen)
    with tempfile.TemporaryDirectory(dir=scratch, prefix=f"tagged-ssc-{kind}-") as temp:
        directory = Path(temp)
        hashes = extract_verified(archive, directory)
        provenance = check_sources(directory, hashes, expected)
        if core:
            output = directory / "replay.json"
            subprocess.run([sys.executable, release / "code/core/response_summary.py", directory, output,
                            "--top-dileptonic"], check=True, stdout=subprocess.DEVNULL)
            if load_json(output) != expected:
                raise ValueError("Core replay differs from the frozen spectrum or its source provenance")
        else:
            guard_code = release / "code/guard"
            fixed = directory / "replay-window.json"
            shape = directory / "replay-shape.json.gz"
            subprocess.run([sys.executable, guard_code / "ew_edge_result.py", directory, fixed,
                            "--plan", guard_code / "ew_cards/zee_edge_plan.json", "--spectrum", shape], check=True)
            if load_json(fixed)["primary_confirmation"] != load_json(release / "zee-confirmation.json")["primary_confirmation"]:
                raise ValueError("Changed confirmation-only migration measurement or confidence limits")
            rebuilt = load_json(shape)
            # This informational field gives the absolute input path on each machine.
            rebuilt["confirmation"]["reference_root"] = expected["confirmation"]["reference_root"]
            if rebuilt != expected:
                raise ValueError("Guard replay differs from the frozen spectrum or its source provenance")
        return dict(archive_sha256=sha256(archive), archive_bytes=archive.stat().st_size,
                    archive_members=len(hashes), spectrum_sha256=sha256(frozen),
                    reproduced_bins=len(expected["rates_fb"]), window_fb=expected["window_125"]["rate_fb"],
                    **provenance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release", type=Path)
    parser.add_argument("--scratch", type=Path, help="Writable temporary directory for the extracted archive")
    parser.add_argument("--output", type=Path, help="Optional JSON audit file")
    args = parser.parse_args()
    result = {kind: replay(args.release.resolve(), args.scratch, kind) for kind in ("core", "guard")}
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
