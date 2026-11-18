"""Loose-isolation diagnostics on stored SSC four-lepton response hypotheses.

For rows failing nominal isolation, the analyzer stores the four leading
crack-accepted leptons. Use that same quartet in both decay samples. All
reported quantities precede the trigger. Masses for this diagnostic use the
stored pT, eta and phi and the massless approximation.
"""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np


def second_mass(arrays):
    pt, eta, phi = (arrays["lepton_" + name] for name in ("pt", "eta", "phi"))
    charge, flavour = (arrays["lepton_" + name] for name in ("charge", "flavour"))
    def pair(i, j):
        m2 = 2 * pt[:, i] * pt[:, j] * (np.cosh(eta[:, i]-eta[:, j]) - np.cos(phi[:, i]-phi[:, j]))
        ossf = (charge[:, i] != charge[:, j]) & (flavour[:, i] == flavour[:, j])
        return np.sqrt(np.maximum(0, m2)), ossf
    second = np.full(len(pt), np.nan)
    distance = np.full(len(pt), np.inf)
    for i,j,k,l in ((0,1,2,3), (0,2,1,3), (0,3,1,2)):
        a, ab = pair(i,j)
        b, cd = pair(k,l)
        da, db = abs(a-91.1876), abs(b-91.1876)
        better = ab & cd & (np.minimum(da,db) < distance)
        second[better] = np.where(da <= db, b, a)[better]
        distance[better] = np.minimum(da,db)[better]
    return second


def shard_moments(manifest_path):
    import uproot
    manifest = json.loads(Path(manifest_path).read_text())
    s, = manifest["samples"]
    if manifest["status"] != "complete" or s["status"] != "complete" or s["selection_version"] != "ssc_v1":
        raise ValueError(f"incomplete/non-SSC manifest: {manifest_path}")
    root = Path(s["output_root"])
    digest = hashlib.sha256()
    with root.open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20), b""):
            digest.update(chunk)
    if digest.hexdigest() != s["output_root_sha256"]:
        raise ValueError(f"changed input: {root}")
    branches = ["cut_mask","source_index","pretrigger_event_weight","lepton_isolation",
                "lepton_pt","lepton_eta","lepton_phi","lepton_charge","lepton_flavour"]
    bins = ("all", "unpaired", "m2_lt12", "m2_12to30", "m2_30to120", "m2_ge120")
    with uproot.open(root) as f:
        n = f["SourceEvents"].num_entries
        sums = {f"{tag}/{reg}/{b}":np.zeros(n) for tag in ("open","lt15")
                for reg in ("3P1F","2P2F") for b in bins}
        for a in f["FourLepton"].iterate(branches,step_size=40000,library="np"):
            selected = (a["cut_mask"] & 15) == 7
            if not selected.any():
                continue
            a = {k:v[selected] for k,v in a.items()}
            iso = a["lepton_isolation"]
            fail = (iso >= 5.).sum(axis=1)
            mass = second_mass(a)
            for tag, loose in (("open", np.ones(len(iso),dtype=bool)), ("lt15",iso.max(axis=1)<15)):
                for reg, mask in (("3P1F",fail==1),("2P2F",fail==2)):
                    for b, mb in zip(bins,(np.ones(len(iso),dtype=bool),np.isnan(mass),mass<12,
                                         (mass>=12)&(mass<30),(mass>=30)&(mass<120), mass>=120)):
                        take = mask & loose & mb
                        np.add.at(sums[f"{tag}/{reg}/{b}"], a["source_index"][take], a["pretrigger_event_weight"][take])
    moments = {k:[float(v.sum()),float(v@v),int(np.count_nonzero(v))] for k,v in sums.items()}
    return s["name"], s["normalization"]["denominator_event_count"], s["normalization"]["native_cross_section_pb"], moments


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--index",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--start",type=int,default=0,help="first parent event (must coincide with a shard boundary)")
    p.add_argument("--stop",type=int,help="exclusive end in each parent, aligned to a shard boundary")
    args=p.parse_args()
    index_bytes = args.index.read_bytes()
    index=json.loads(index_bytes)
    if index["status"]!="complete": raise ValueError("index is incomplete")
    events, rates, stats = {}, {}, {}
    manifests=[]
    for record in index["samples"].values():
        cursor=0
        end=args.stop if args.stop is not None else record["source_events"]
        if not 0 <= args.start < end <= record["source_events"]:
            raise ValueError("invalid parent range")
        for shard in sorted(record["shards"],key=lambda s:s["start"]):
            if shard["start"]!=cursor or shard["stop"]<=cursor or shard["status"]!="complete": raise ValueError("invalid source partition")
            if shard["start"] < args.start < shard["stop"] or shard["start"] < end < shard["stop"]:
                raise ValueError("range splits a source shard")
            if args.start <= shard["start"] < end:
                manifests.append(shard["manifest"])
            cursor=shard["stop"]
        if cursor!=record["source_events"]: raise ValueError("incomplete source coverage")
    if len(manifests)!=len(set(manifests)): raise ValueError("duplicate manifest")
    inputs = []
    for name in manifests:
        raw = Path(name).read_bytes()
        sample, = json.loads(raw)["samples"]
        inputs.append(dict(manifest=name, manifest_sha256=hashlib.sha256(raw).hexdigest(),
                           output_root_sha256=sample["output_root_sha256"]))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for name, n, xs, moments in pool.map(shard_moments,manifests):
            events[name]=events.get(name,0)+n
            if rates.setdefault(name,xs)!=xs: raise ValueError("unequal cross sections within a parent")
            dest=stats.setdefault(name,defaultdict(lambda:np.zeros(3)))
            for key, value in moments.items(): dest[key]+=value
    results={}
    for name, moments in stats.items():
        scale=1000*rates[name]/events[name]
        results[name]={key:dict(rate_fb=a*scale,mc_error_fb=np.sqrt(b)*scale,
                                sources=int(c),effective_sources=a*a/b if b else 0.)
                       for key,(a,b,c) in moments.items()}
    if any(hashlib.sha256(Path(item["manifest"]).read_bytes()).hexdigest() != item["manifest_sha256"]
           for item in inputs) or args.index.read_bytes() != index_bytes:
        raise ValueError("input index or manifest changed during sideband evaluation")
    args.output.write_text(json.dumps(dict(index=str(args.index),start=args.start,stop=args.stop,source_events=events,
        index_sha256=hashlib.sha256(index_bytes).hexdigest(),
        diagnostic_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), inputs=inputs,
        definition="Pretrigger, four leading crack-accepted leptons in failed-isolation hypotheses; massless nearest-Z pairing; lt15 requires all four isolation ET <15 GeV, open has no upper bound. Regions compare generator responses.",
        results=results),indent=2)+"\n")

if __name__=="__main__":
    main()
