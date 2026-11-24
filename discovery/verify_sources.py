"""Check the final diphoton campaign against the independently accessible Timur source archive."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np

from input_data import load_archive


def verify(release):
    all_data, ordered = load_archive()
    campaign = release / "sources/sideband_20pct_20261120"
    manifest = {}
    by_digest = defaultdict(list)
    with (release / "source-files.sha256").open() as listing:
        for line in listing:
            digest, name = line.rstrip("\n").split("  ", 1)
            path = release / name
            manifest[path] = digest
            by_digest[digest].append(path)

    def require(path, digest=None):
        if not path.is_file() or path not in manifest or (digest and manifest[path] != digest):
            raise ValueError(f"Missing or inconsistent source {path}")

    def read_spectrum(path, expected_hash=None):
        require(path, expected_hash)
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != manifest[path]:
            raise ValueError(f"Changed source spectrum {path}")
        return json.loads(raw)

    counts = Counter()
    for name in ordered:
        if not name.startswith("packets/final/"):
            continue
        packet = all_data[name]
        sample, size = packet["sample"], packet["generated_events"]
        if sample == "bkg_loop_gg_aa":
            # A repaired shower can live under a different path from its first attempt.
            folder = campaign / "fake_loop/loop"
            for key in ("response_root_sha256", "source_lhe_sha256"):
                if not any(folder in path.parents and path.is_file() for path in by_digest[packet[key]]):
                    raise ValueError(f"Missing loop source {key} for {name}")
        else:
            rates = np.zeros(55)
            covariance = np.zeros((55, 55))
            covered = 0
            if sample == "bkg_dijet_double_fake":
                chunks = [read_spectrum(campaign / "fake_loop/dijet" / shard["path"], shard["sha256"])
                          for shard in packet["shards"]]
            else:
                origin = packet["provenance"]
                folder = campaign / "nlo/parents" / Path(origin["parent_lhe"]).parent.name
                require(folder / "events.lhe.gz", origin["sha256"])
                chunks = []
                for path in sorted(folder.glob("shard_*/receipt.json")):
                    require(path.with_name("ssc-side_var.root"))
                    chunks.append(read_spectrum(path))
                if len(chunks) != origin["complete_shards"]:
                    raise ValueError(f"Missing NLO shower shards for {name}")
            for chunk in chunks:
                n = chunk["generated_events"]
                if "first_source" in chunk and (chunk["first_source"] != covered
                        or chunk["last_source_exclusive"] != covered + n):
                    raise ValueError(f"Incomplete event coverage for {name}")
                covered += n
                rates += np.asarray(chunk["rate_fb"]) * (n / size)
                covariance += np.asarray(chunk["covariance_fb2"]) * (n / size)**2
            if (covered != size or not np.allclose(rates, packet["rate_fb"], rtol=1.e-9, atol=1.e-10)
                    or not np.allclose(covariance, packet["covariance_fb2"], rtol=1.e-9, atol=1.e-10)):
                raise ValueError(f"Source shards disagree with {name}")
        counts[sample] += 1
    print(f"Verified {sum(counts.values())} final campaign packets against the source archive")
    for sample, count in sorted(counts.items()):
        print(f"{sample}: {count}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path,
                        default=Path("/home/shared/HiggsSSC/student-results/ssc40-discovery-20261123"))
    verify(parser.parse_args().archive)
