"""Archive the exact hard events, cards, detector packets and receipts used in a tagged result."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import tarfile


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            sha.update(block)
    return sha.hexdigest()


def archive(root, report_path, target, software_root=None):
    if target.exists():
        raise FileExistsError(target)
    opener = gzip.open if report_path.suffix == ".gz" else open
    with opener(report_path, "rt") as stream:
        report = json.load(stream)
    chosen = set()
    parents = set()
    response_hashes = set()
    plugins = {}
    for row in report["shards"]:
        receipt = root / row["receipt"]
        packet = receipt.with_name("packet.csv")
        if digest(receipt) != row["receipt_sha256"] or digest(packet) != row["packet_sha256"]:
            raise ValueError(f"Detector packet or receipt changed in {receipt.parent}")
        parents.add((row["integration_id"], row["source_lhe_sha256"]))
        audit = json.loads(receipt.read_text())
        config = receipt.with_name("tagged.in")
        response_hashes.add(audit["response_sha256"])
        for line in config.read_text().splitlines():
            if line.startswith("library /"):
                plugin = Path(line.split(maxsplit=1)[1].strip())
                plugins[plugin] = audit["plugin_sha256"]
        chosen.update([receipt, packet, config])
    # Only include integrations which actually contributed to the published estimate.
    found = set()
    grid_hashes = set()
    hard_paths = itertools.chain(root.glob("top/tt_*/round-*/hard.json"), root.glob("ew/*/hard.json"),
                                 root.glob("*/hard.json"), root.glob("tt_*/round-*/hard.json"))
    grid_paths = set()
    for path in hard_paths:
        info = json.loads(path.read_text())
        grid = path.with_name("native-gridpack.tar.gz")
        if grid.is_file():
            grid_paths.add(grid)
        source = (info.get("integration_id"), info.get("lhe_sha256"))
        if source not in parents:
            continue
        lhe = path.with_name("events.lhe.gz")
        if source in found or digest(lhe) != source[1]:
            raise ValueError(f"Ambiguous or altered hard events for {source}")
        found.add(source)
        chosen.update([path, lhe, path.with_name("run_card.dat"),
                       path.with_name("param_card.dat"), path.with_name("proc_card_mg5.dat")])
        if "cut_function_sha256" in info:
            cut = path.with_name("outside_nominal.f")
            if digest(cut) != info["cut_function_sha256"]:
                raise ValueError(f"Generator acceptance function changed in {path.parent}")
            chosen.add(cut)
        if "gridpack_sha256" in info:
            grid_hashes.add(info["gridpack_sha256"])
            for lot in info["native_lots"]:
                native = Path(lot["lhe"]).resolve()
                if not native.is_relative_to(root) or digest(native) != lot["lhe_sha256"]:
                    raise ValueError(f"Unverified native gridpack events in {path.parent}")
                chosen.add(native)
        else:
            chosen.add(path.with_name("banner.txt"))
    if found != parents:
        raise ValueError(f"Missing hard-event inputs for {parents-found}")
    if grid_hashes:
        packs = {digest(p):p for p in grid_paths}
        if grid_hashes - packs.keys():
            raise ValueError(f"Missing adaptive grids for {grid_hashes-packs.keys()}")
        chosen.update(packs[sha] for sha in grid_hashes)
    executables = {}
    for path in (software_root or root).glob("response*"):
        if path.is_file() and path.stat().st_mode & 0o111:
            executables[digest(path)] = path
    if response_hashes - executables.keys():
        raise ValueError(f"Missing detector binaries for {response_hashes - executables.keys()}")
    software = {executables[sha]: f"software/{sha}/{executables[sha].name}" for sha in response_hashes}
    for path, sha in plugins.items():
        if digest(path) != sha:
            raise ValueError(f"Shower plugin hash mismatch at {path}")
        software[path] = f"software/{sha}/{path.name}"
        proof = path.with_name("build.json")
        if proof.is_file():
            software[proof] = f"software/{sha}/build.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(target.name+".partial")
    # A compact packet suffices to rerun the detector-level statistics.
    # The archived LHE events also allow the full shower to be regenerated.
    with tarfile.open(tmp, "x:gz", compresslevel=2, dereference=True) as out:
        for path in sorted(chosen):
            out.add(path, arcname=path.relative_to(root), recursive=False)
        for path, name in sorted(software.items()):
            out.add(path, arcname=name, recursive=False)
    with tarfile.open(tmp, "r:gz") as check:
        if len(check.getmembers()) != len(chosen) + len(software):
            raise ValueError("Incomplete tagged background archive")
    tmp.rename(target)
    print(f"{digest(target)}  {target.name}")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("campaign", type=Path)
    cli.add_argument("spectrum", type=Path)
    cli.add_argument("archive", type=Path)
    cli.add_argument("--software-root", type=Path, help="Directory holding the detector-response binaries")
    args = cli.parse_args()
    archive(args.campaign.resolve(), args.spectrum, args.archive, args.software_root)
