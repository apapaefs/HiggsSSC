"""Match saved HwSim events to LHE collisions and export prompt photon four-vectors."""
import argparse
import gzip
import math
from pathlib import Path

import numpy as np


def lhe_events(path):
    stream = gzip.open(path, "rt") if path.suffix == ".gz" else path.open()
    with stream:
        for line in stream:
            if line.strip() != "<event>":
                continue
            n = int(next(stream).split()[0])
            particles = [next(stream).split() for _ in range(n)]
            incoming = sorted(((int(p[0]), float(p[8]), float(p[9])) for p in particles if p[1]=="-1"),
                              key=lambda p:-p[1])
            if len(incoming) != 2:
                raise ValueError("Expected two incoming LHE partons")
            photons = []
            for i, p in enumerate(particles):
                if p[0] != "22" or p[1] != "1":
                    continue
                # Resolve common mother ranges in the original, one-based LHE graph.
                decay = False
                visited = set()
                stack = [i+1]
                while stack:
                    j = stack.pop()
                    if j in visited:
                        continue
                    visited.add(j)
                    parent = particles[j-1]
                    if abs(int(parent[0])) in (6,24):
                        decay = True
                    first, last = int(parent[2]), int(parent[3])
                    if first:
                        if last < first or last >= j:
                            raise ValueError("Malformed LHE parent indices")
                        stack.extend(range(first,last+1))
                photons.append((int(decay), *(float(x) for x in p[6:10])))
            yield incoming, photons


def export(lhe, incoming, output):
    source = iter(enumerate(lhe_events(lhe)))
    with output.open("w") as out:
        out.write("saved_entry,lhe_entry,decay,px,py,pz,energy\n")
        for i, arr in enumerate(incoming):
            beams = sorted(((int(arr[4,j]),arr[3,j],arr[0,j]) for j in range(2)), key=lambda p:-p[1])
            def matches(other):
                return all(p[0]==q[0] and all(math.isclose(x,y,rel_tol=2.e-8,abs_tol=1.e-7)
                           for x,y in zip(p[1:],q[1:])) for p,q in zip(beams,other))
            for source_index, (hard_beams, photons) in source:
                if matches(hard_beams):
                    break
            else:
                raise ValueError(f"Cannot identify saved entry {i} by both incoming flavours, pz and energies")
            if not photons:
                raise ValueError(f"No hard photon in LHE event {source_index}")
            for decay, *p in photons:
                out.write(",".join(map(str, (i,source_index,decay,*p)))+"\n")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("lhe","root","output"):
        parser.add_argument(name,type=Path)
    args = parser.parse_args()
    import uproot
    with uproot.open(args.root) as f:
        arr = np.asarray(f["Data"]["incoming"].array(library="np"))
    export(args.lhe, arr, args.output)
