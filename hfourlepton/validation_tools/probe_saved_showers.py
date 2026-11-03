#!/usr/bin/env python3
"""Read-only CMS-v2 definition probes on existing HwSim files; no generators.

Outputs are diagnostic variations, not validated physics corrections. In
particular, angular matching of pre-shower hard leptons is not ancestry.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()
    if not 1 <= args.jobs <= 40:
        p.error("jobs must be between 1 and 40")
    manifest = args.manifest.resolve()
    raw = json.loads(manifest.read_text())
    if raw.get("selection_profile") != "cms" or raw.get("detector_profile") != "perfect":
        p.error("this probe requires a CMS/perfect manifest")
    samples = raw["samples"]
    if not samples or any(s.get("status") != "complete" or s.get("selection_version") != "cms_fiducial_v2" for s in samples):
        p.error("this diagnostic is defined for completed cms_fiducial_v2 samples")
    recorded = raw["tracked_inputs"]["analyzer_source"]
    analyzer = Path(recorded["path"])
    if digest(analyzer) != recorded["sha256"]:
        p.error("recorded analyzer source changed")
    for sample in samples:
        roots = sample["shower_completion"]["artifacts"]["roots"]
        if len(roots) != 1:
            p.error("this campaign diagnostic supports one saved HwSim file per sample")
        root = roots[0]
        if digest(Path(root["path"])) != root["sha256"]:
            p.error("saved shower changed: " + sample["name"])
    output = args.output_dir.resolve()
    if output.exists():
        p.error("choose a new diagnostic output directory")
    output.mkdir(parents=True)
    source = Path(__file__).with_name("acceptance_probe.cc")
    shutil.copyfile(source, output / source.name)
    shutil.copyfile(analyzer, output / analyzer.name)
    flags = shlex.split(subprocess.check_output(["root-config", "--cflags", "--libs"], text=True))
    executable = output / "acceptance_probe"
    with (output / "build.log").open("w") as log:
        subprocess.run(["c++", "-O2", "-std=c++17", str(output / source.name), "-o", str(executable), *flags],
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)
    started = time.monotonic()

    def inspect(sample):
        root = sample["shower_completion"]["artifacts"]["roots"][0]
        result = subprocess.run([str(executable), root["path"]], text=True,
                                capture_output=True, check=True, timeout=240)
        stats = json.loads(result.stdout)
        if stats["all"]["entries"] != sample["shower_completion"]["saved_events"]:
            raise ValueError("saved event count differs: " + sample["name"])
        return dict(name=sample["name"], mode=sample.get("production_mode", ""),
                    normalization=sample["normalization"], stats=stats)

    result = dict(manifest=str(manifest), manifest_sha256=digest(manifest),
                  probe_sha256=digest(source), analyzer_source_sha256=digest(analyzer), samples=[])
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(inspect, sample) for sample in samples]):
            sample = future.result()
            result["samples"].append(sample)
            print("Inspected", sample["name"], flush=True)
    result["seconds"] = time.monotonic() - started
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print("Diagnostic variations only:", output / "results.json")


if __name__ == "__main__":
    main()
