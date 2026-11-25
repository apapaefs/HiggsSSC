"""Re-create the original eeγ HwSim trees for the frozen shower shards."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ew_shower import isolate


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def spectrum(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return json.load(stream)


def records(snapshot, packet_root, hard_root):
    parents = {json.loads(p.read_text())["lhe_sha256"]: p.with_name("events.lhe.gz")
               for p in hard_root.glob("z_ee_a*/hard.json")}
    for row in spectrum(snapshot)["shards"]:
        if row["sample"] != "z_ee_a":
            continue
        yield packet_root / row["receipt"], parents[row["source_lhe_sha256"]], row["receipt_sha256"], packet_root


def replay(record, destination):
    receipt, parent, expected, packet_root = record
    out = destination / receipt.parent.name
    out.mkdir(parents=True, exist_ok=True)
    done = out / "tree.json"
    if done.exists():
        return
    if digest(receipt) != expected:
        raise ValueError(f"Changed input receipt {receipt}")
    ref = json.loads(receipt.read_text())
    old = receipt.parent / "events/tagged.root"
    card = receipt.with_name("tagged.in")
    if digest(card) != ref["card_sha256"]:
        raise ValueError(f"Changed shower card {card}")
    if old.exists():
        if digest(old) != ref["root_sha256"]:
            raise ValueError(f"Original ROOT file changed {old}")
        root = old
    else:
        (out / "events").mkdir(exist_ok=True)
        lhe = out / "source.lhe.gz"
        sha = isolate(parent, lhe, ref["source_first_event"], ref["attempted_events"])
        if sha != ref["event_block_sha256"]:
            raise ValueError(f"Changed hard event block {receipt}")
        text = card.read_text()
        prefix = "set LesHouchesReader:FileName "
        stored_lhe = next(line[len(prefix):] for line in text.splitlines() if line.startswith(prefix))
        original_plugin = Path(next(line[8:] for line in text.splitlines() if line.startswith("library /")))
        archived_plugin = packet_root / "software" / ref["plugin_sha256"] / original_plugin.name
        plugin = archived_plugin.resolve() if archived_plugin.is_file() else original_plugin
        text = text.replace(prefix + stored_lhe, prefix + str(lhe))
        (out / "tagged.in").write_text(text.replace("library " + str(original_plugin), "library " + str(plugin)))
        if digest(plugin) != ref["plugin_sha256"]:
            raise ValueError(f"Changed HwSim plugin {plugin}")
        for name,cmd in [("read", ["read","tagged.in"]),("run", ["run","tagged.run", f"-N{ref['attempted_events']+2}"])]:
            with (out / f"{name}.log").open("w") as stream:
                subprocess.run(["Herwig", *cmd, "-L", str(plugin.parent)], cwd=out,
                               stdout=stream, stderr=subprocess.STDOUT, check=(name=="read"))
        root = out / "events/tagged.root"
        if not root.is_file() or root.stat().st_size == 0:
            raise ValueError(f"Missing shower tree {out}")
    tmp = done.with_suffix(".tmp")
    tmp.write_text(json.dumps(dict(root=str(root), receipt=str(receipt), original_root=(old==root)))+"\n")
    os.replace(tmp, done)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("core","core_packets","core_hard","wide","wide_packets","wide_hard","out"):
        p.add_argument(key,type=Path)
    p.add_argument("--workers",type=int,default=8)
    args=p.parse_args()
    jobs=[*records(args.core,args.core_packets,args.core_hard),*records(args.wide,args.wide_packets,args.wide_hard)]
    # Later pilot seeds in the wider region have no entry in the frozen model.
    if len({j[0].parent.name for j in jobs}) != len(jobs):
        raise ValueError("Shower block directory names should be unique")
    with ThreadPoolExecutor(args.workers) as pool:
        for count, _ in enumerate(pool.map(lambda j:replay(j,args.out.resolve()),jobs),1):
            if count%25==0 or count==len(jobs): print(count,len(jobs),flush=True)
