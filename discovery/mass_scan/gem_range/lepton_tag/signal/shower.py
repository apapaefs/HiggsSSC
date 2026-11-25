"""Reshower independent 125 GeV WH and ttH hard events with H -> gamma gamma."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

MODES = ("wp_h", "wm_h", "tt_h")
FLAVOURS = ("4e", "4mu", "2e2mu")
BR_GAMMAGAMMA = .00227


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_job(events, plugin, out, name, parent, header, block, seed):
    directory = out / name
    directory.mkdir(parents=True, exist_ok=False)
    (directory / "events").mkdir()
    data = directory / "hard.lhe.gz"
    with gzip.open(data, "wt") as stream:
        stream.write(header + "".join(block) + "</LesHouchesEvents>\n")
    card_path = events / parent / "herwig" / (parent + ".in")
    card = card_path.read_text()
    card = re.sub(r"(?m)^set LesHouchesReader:FileName .+$", "set LesHouchesReader:FileName " + str(data), card)
    card = re.sub(r"(?m)^set EventGenerator:RandomNumberGenerator:Seed \d+$", f"set EventGenerator:RandomNumberGenerator:Seed {seed}", card)
    card = card.replace("set /Herwig/Particles/h0:NominalMass 125.0*GeV", """set /Herwig/Particles/h0:NominalMass 125.0*GeV
do /Herwig/Particles/h0:SelectDecayModes h0->gamma,gamma;
set /Herwig/Particles/h0/h0->gamma,gamma;:BranchingRatio 1.0""")
    card = card.replace("library HwSim.so", "library " + str(plugin))
    card = card.replace("set /Herwig/Analysis/HwSim:SaveObjects Yes", "set /Herwig/Analysis/HwSim:SaveObjects Yes\nset /Herwig/Analysis/HwSim:SaveAncestry Yes")
    card = card.replace("saverun " + parent, "saverun tagged")
    (directory / "tagged.in").write_text(card)
    audit = dict(seed=seed, events=len(block), parent=parent, parent_card_sha256=digest(card_path))
    for term, command in (("read", ["Herwig", "read", "tagged.in"]),
                          ("run", ["Herwig", "run", "tagged.run", "-N"+str(len(block)+2)])):
        with (directory/(term+".log")).open("w") as log:
            result = subprocess.run(command + ["-L", str(plugin.parent)], cwd=directory, stdout=log, stderr=subprocess.STDOUT)
        if term == "read" and result.returncode:
            raise RuntimeError("Herwig read failed in " + str(directory))
    # Check finite-source termination from the project’s existing audit.
    from hgammagamma import ho_shower_completion as completion
    import uproot
    roots = list((directory / "events").glob("*.root"))
    if len(roots) != 1:
        raise ValueError(f"Expected one ROOT file in {directory}")
    import numpy as np
    with uproot.open(roots[0]) as source:
        signed = source["Data"]["evweight"].array(library="np")
    audit.update(completion._parse_logs(directory, "tagged", "source_exhausted"))
    if audit["attempted_events"] != len(block) or audit["saved_events"] != len(signed):
        raise ValueError(f"Incomplete source in {directory}")
    norm = completion.source_weight_summary(data, len(block))
    audit.update(signed_attempted=norm["sum_weight"], signed_saved=float(signed.sum()),
                 negative_saved=int((signed < 0).sum()), lhe_sha256=digest(data))
    (directory / "receipt.json").write_text(json.dumps(audit, indent=2)+"\n")
    print("Completed", directory.name, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("events", type=Path)
    parser.add_argument("plugin", type=Path)
    parser.add_argument("repo", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pilot", type=int, default=0, help="Use this many collisions per mode from the first source only")
    parser.add_argument("--workers", type=int, default=18)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo.resolve()))
    events, plugin, out = [p.resolve() for p in (args.events, args.plugin, args.output)]
    out.mkdir(parents=True, exist_ok=False)
    jobs, inputs = [], []
    for i, mode in enumerate(MODES):
        for j, flavour in enumerate(FLAVOURS[:1] if args.pilot else FLAVOURS):
            parent = f"signal_{mode}_{flavour}"
            path = events / parent / f"hard-{flavour}.lhe.gz"
            with gzip.open(path, "rt") as stream:
                raw = stream.read()
            head = raw.index("<event>")
            blocks = re.findall(r"<event>.*?</event>\s*", raw[head:], flags=re.S)
            if len(blocks) != 5000:
                raise ValueError("Unexpected number of hard events at " + str(path))
            inputs.append(dict(relative_path=str(path.relative_to(events)), sha256=digest(path), generated=len(blocks)))
            if args.pilot:
                blocks = blocks[:args.pilot]
            for start in range(0, len(blocks), 1000):
                name = f"{mode}-{flavour}-{start:04d}"
                jobs.append((events,plugin,out,name,parent,raw[:head],blocks[start:start+1000],411000+i*1000+j*100+start//1000))
    (out/"inputs.json").write_text(json.dumps(dict(inputs=inputs, plugin_sha256=digest(plugin)), indent=2)+"\n")
    with ThreadPoolExecutor(max_workers=args.workers) as workers:
        list(workers.map(lambda a: write_job(*a), jobs))


if __name__ == "__main__":
    main()
