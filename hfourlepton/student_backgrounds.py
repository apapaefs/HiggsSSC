"""Package and read the two detector responses of the published heavy-flavour controls."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from hfourlepton import shared_events as shared

DIRECTORY = REPO / "hfourlepton/student_backgrounds/ssc40_v1"


def controls_for(raw, directory=DIRECTORY):
    """Choose only the response calculated for this 40 TeV SSC selection."""
    profile = raw.get("detector_profile")
    if (raw.get("order") != "ho" or raw.get("collider", {}).get("sqrt_s_tev") != 40
            or raw.get("selection_profile") != "ssc" or profile not in ("perfect", "ssc")
            or {s.get("selection_version") for s in raw["samples"]} != {"ssc_v1"}):
        raise ValueError("the published heavy-flavour controls require the 40 TeV HO campaign and SSC v1 selection")
    manifest = json.loads((directory / "manifest.json").read_text())
    response = manifest["profiles"][profile]
    if manifest["schema_version"] != 1 or raw["detector"] != response["detector"]:
        raise ValueError("the detector settings differ from the published heavy-flavour controls")
    return response["samples"]


def package(release, destination):
    """Copy the exact public analysis products after verifying their pinned digests."""
    from hfourlepton.pooled_controls import PROCESSES
    index_path = release / "dataset.json"
    index = json.loads(index_path.read_bytes())
    if index.get("kind") != "fourlepton_shared_events" or index.get("status") != "complete":
        raise ValueError("a completed shared four-lepton release is required")
    if destination.exists():
        raise ValueError(f"choose an unused destination: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        staging = Path(tmp) / "release"
        staging.mkdir()
        profiles = {}
        for profile in ("ssc", "perfect"):
            relative = index["manifests"][profile]
            manifest_path = shared.safe_path(release, relative)
            if shared.digest(manifest_path) != index["files"][relative]["sha256"]:
                raise ValueError(f"shared {profile} manifest checksum changed")
            raw = json.loads(manifest_path.read_bytes())
            refs = raw["pooled_backgrounds"]
            if set(refs) != set(PROCESSES):
                raise ValueError(f"incomplete {profile} control inventory")
            samples = {}
            for name, info in sorted(refs.items()):
                source = shared.safe_path(release, info["path"])
                if shared.digest(source) != info["sha256"] or info["sha256"] != index["files"][info["path"]]["sha256"]:
                    raise ValueError(f"shared control checksum changed: {source}")
                filename = f"{name}_{profile}.json.gz"
                target = staging / filename
                target.write_bytes(gzip.compress(source.read_bytes(), mtime=0))
                samples[name] = dict(path=filename, sha256=shared.digest(target), public_path=info["path"],
                                     public_sha256=info["sha256"])
            profiles[profile] = dict(detector=raw["detector"], samples=samples)
        shared.save(staging / "manifest.json", dict(schema_version=1, release=release.name,
                    dataset_sha256=shared.digest(index_path), profiles=profiles))
        staging.rename(destination)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    package(args.release, args.output)
