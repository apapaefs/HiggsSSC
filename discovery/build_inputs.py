"""Rebuild the portable input archive from the finished campaign and student records."""
import argparse
from contextlib import ExitStack
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def build(campaign, records, output):
    control = campaign / "control"
    fit_bytes = (control / "final_fit.json").read_bytes()
    reach_bytes = (control / "final_reach.json").read_bytes()
    fit, reach = json.loads(fit_bytes), json.loads(reach_bytes)
    if sha(fit_bytes) != reach["fit_sha256"]:
        raise ValueError("Final campaign fit checksum differs from the reach record")
    sources = []
    def add(name, path, expected=None):
        data = path.read_bytes()
        if expected is not None and expected != sha(data):
            raise ValueError(f"Source checksum mismatch for {name}")
        sources.append((name, data))
    add("baseline.json", campaign / "fit/prior/baseline.json", fit["reference_sha256"])
    for i, packet in enumerate(fit["independent_new_inputs"]):
        parts = Path(packet["path"]).parts
        start = parts.index("sideband_20pct_20261120") + 1
        local = Path(*parts[start:])
        phase = "prior" if i < reach["precampaign_receipts"] else "final"
        if phase == "final" and reach["new_inputs_sha256"][packet["path"]] != packet["sha256"]:
            raise ValueError("Final packet disagrees between the two signed-off outputs")
        add(f"packets/{phase}/{i:04d}-{local.name}", campaign / local, packet["sha256"])
    if len(sources) != 1 + reach["precampaign_receipts"] + reach["extra_receipts"]:
        raise ValueError("Missing campaign inputs")
    for name in ("extended_fourlepton.json", "reducible_bound.json", "ihixs_normalization.json"):
        add("records/" + name, records / name)
    add("records/final_fit.json", control / "final_fit.json", reach["fit_sha256"])
    add("records/final_reach.json", control / "final_reach.json")
    manifest = {"schema": 1, "files": {name: sha(data) for name, data in sources}}
    sources.insert(0, ("manifest.json", (json.dumps(manifest, indent=2) + "\n").encode()))
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        target = Path(tmp) / "spectra.tar.gz"
        with ExitStack() as stack:
            stream = stack.enter_context(target.open("wb"))
            zipped = stack.enter_context(gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0, compresslevel=9))
            tar = stack.enter_context(tarfile.open(fileobj=zipped, mode="w|", format=tarfile.USTAR_FORMAT))
            for name, data in sources:
                entry = tarfile.TarInfo(name)
                entry.size = len(data)
                tar.addfile(entry, io.BytesIO(data))
        target.replace(output)
    return sha(output.read_bytes())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("campaign", type=Path)
    p.add_argument("records", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    print(build(args.campaign, args.records, args.output))
