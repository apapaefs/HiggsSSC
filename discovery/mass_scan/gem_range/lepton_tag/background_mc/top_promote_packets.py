"""Atomically promote a reclassified, event-equivalent top-decay packet directory on Linux."""
import argparse
import csv
import ctypes
import itertools
import json
import os
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("production", type=Path)
    p.add_argument("staged", type=Path)
    p.add_argument("marker", type=Path)
    p.add_argument("--expected-events", required=True, type=int)
    a = p.parse_args()
    if a.marker.is_file():
        return
    packets = lambda root: {r.relative_to(root) for r in root.glob("*/*/receipt.json")}
    names = packets(a.production)
    if not names or names != packets(a.staged):
        raise ValueError("old and new packet intervals differ")
    count = 0
    for name in sorted(names):
        old, new = (json.loads((base / name).read_text()) for base in (a.production, a.staged))
        for key in ("sample", "integration_id", "source_first_event", "source_stop_event", "seed",
                    "event_block_sha256", "attempted_events", "saved_events", "plugin_sha256"):
            if old[key] != new[key]:
                raise ValueError(f"nonidentical generation or shower settings in {name}: {key}")
        count += new["attempted_events"]
        def rows(base):
            with (base / name).with_name("packet.csv").open() as f:
                reader = csv.DictReader(f)
                if base == a.staged and "direct_top_dilepton" not in reader.fieldnames:
                    raise ValueError("missing decay classification")
                for row in reader:
                    row.pop("direct_top_dilepton", None)
                    if row["component"] == "decay_electron_ll":
                        row["component"] = "decay_electron"
                    yield row
        for x, y in itertools.zip_longest(rows(a.production), rows(a.staged)):
            if x != y:
                raise ValueError(f"detector outcomes changed in {name}")
    if count != a.expected_events:
        raise ValueError(f"expected {a.expected_events} attempted events, found {count}")
    # Atomic exchange on the same Linux filesystem. The former production files
    # remain in the staging location for an independent cross-check.
    rename = ctypes.CDLL(None, use_errno=True).renameat2
    result = rename(-100, os.fsencode(a.production), -100, os.fsencode(a.staged), 2)
    if result:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    temporary = a.marker.with_suffix(".partial")
    temporary.write_text(json.dumps(dict(events=count, shards=len(names), previous_packets=str(a.staged),
                                         response_sha256=new["response_sha256"]), indent=2) + "\n")
    temporary.replace(a.marker)


if __name__ == "__main__":
    main()
