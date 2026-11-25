"""Recover and reweight the hard collision of every candidate before the isolation variations."""
import argparse
from collections import defaultdict
import ctypes
import gzip
import json
from pathlib import Path
import sys

import numpy as np

AUDIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(AUDIT))
from hard_matrix import matrix_sources, evaluate
from pdf_reweight import PDF_NAME, SOURCES, densities, digest, read_events, recovered_indices


def read(path):
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt") as f:
        return json.load(f)


def write(path, value):
    data = (json.dumps(value,separators=(",",":"),allow_nan=False)+"\n").encode()
    if path.suffix == ".gz":
        with path.open("wb") as f, gzip.GzipFile(fileobj=f,mode="wb",filename="",mtime=0) as out:
            out.write(data)
    else:
        path.write_bytes(data)


def event_requests(features):
    requested = defaultdict(set)
    for row in features:
        requested[row["shard"]].add(int(row["event"]))
    return requested


def extract_momenta(path, events):
    with gzip.open(path,"rt") as stream:
        index = -1
        for line in stream:
            if line.strip() != "<event>":
                continue
            index += 1
            head = next(stream).split()
            rows = [next(stream).split() for _ in range(int(head[0]))]
            if index not in events:
                continue
            inc = sorted((r for r in rows if r[1] == "-1"), key=lambda r: -int(r[0]))
            out = {int(r[0]):r for r in rows if r[1] == "1"}
            flavour = int(inc[0][0])
            if len(inc)!=2 or flavour not in (1,2,3,4) or int(inc[1][0])!=-flavour or set(out)!={11,-11,22}:
                raise ValueError("Unknown hard final state")
            momenta = [[float(r[k]) for k in (9,6,7,8)] for r in inc+[out[k] for k in (-11,11,22)]]
            events[index].update(hard_flavour=flavour,p4=momenta,scale_gev=float(head[3]))
    if any("p4" not in row for row in events.values()):
        raise ValueError("Could not recover all selected LHE events")


def collect(a):
    import lhapdf
    for source, expected in SOURCES.items():
        if digest(a.library.resolve().with_name(source))!=expected:
            raise ValueError("Historical generator sources differ from original audit")
    lib=ctypes.CDLL(str(a.library.resolve()))
    lhapdf.setVerbosity(0)
    pdf=lhapdf.mkPDF(PDF_NAME,0)
    requested=event_requests(read(a.features))
    result=dict(groups={},events=[],source_spectra={})
    for name,root,snapshot in (("core",a.core,a.core_spectrum),("complement",a.complement,a.complement_spectrum)):
        model=read(snapshot)
        result["source_spectra"][name]=digest(snapshot)
        source=model["sources"]["z_ee_a"]
        records=defaultdict(list)
        for shard in model["shards"]:
            if shard["sample"]=="z_ee_a":
                records[shard["integration_id"]].append(shard)
        parents={read(p)["integration_id"]:p for par in (root,root/"ew") for p in par.glob("z_ee_a*/hard.json")}
        for group in source["integration_groups"]:
            key=group["integration_id"]
            result["groups"][key]={**group,"region":name,"region_generated":source["attempted_events"]}
            meta=parents[key]
            info=read(meta)
            path=meta.with_name("events.lhe.gz")
            if digest(path)!=info["lhe_sha256"] or info["lhe_init"][6:8]!=["335900","335900"]:
                raise ValueError("LHE parent or PDF changed")
            n=group["generated"]
            beam=float(info["lhe_init"][2])
            ids,x,q=read_events(path,n,beam)
            events={}
            for shard in records[key]:
                receipt=root/shard["receipt"]
                if digest(receipt)!=shard["receipt_sha256"] or shard["source_lhe_sha256"]!=info["lhe_sha256"]:
                    raise ValueError("Published shower receipt mismatch")
                tag=receipt.parent.name
                entries=sorted(requested.pop(tag,[]))
                r=read(receipt)
                if not entries:
                    continue
                if min(entries)<0 or max(entries)>=r["saved_events"]:
                    raise ValueError("Detector index outside saved shower range")
                if r["saved_events"]==r["attempted_events"]:
                    mapping={i:r["source_first_event"]+i for i in entries}
                else:
                    mapping=recovered_indices(entries,ids,x,beam,receipt,a.trees)
                for entry,i in mapping.items():
                    if i in events:
                        raise ValueError("Hard collision assigned more than once")
                    events[i]=dict(shard=tag,event=entry,hard_event=i,integration=key)
            if events:
                indices=np.array(sorted(events))
                pdf_ratios=densities(lib,pdf,ids[indices],x[indices],q[indices])
                for i,w in zip(indices,pdf_ratios):
                    events[i]["pdf_ratio"]=float(w)
                extract_momenta(path,events)
                result["events"].extend(events[i] for i in indices)
    if requested:
        raise ValueError(f"Unmatched detector blocks {tuple(requested)}")
    write(a.output,result)


def weight(a):
    data=read(a.input)
    a.work.mkdir(parents=True,exist_ok=True)
    work=a.work.resolve()
    matrix_sources(a.mg5.resolve(),work)
    for up in (True,False):
        rows=[r for r in data["events"] if (r["hard_flavour"]%2==0)==up]
        if not rows:
            continue
        full=evaluate(work,"full",up,rows,a.card.resolve())
        z=evaluate(work,"production_z",up,rows,a.card.resolve())
        if np.any(full<=0):
            raise ValueError("Vanishing full matrix element")
        for row,ratio in zip(rows,z/full):
            p=row["p4"][-1]
            row["production_ratio"]=float(ratio)
            row["hard_photon_pt_gev"]=float(np.hypot(p[1],p[2]))
    data["parameter_sha256"]=digest(a.card)
    data["mg5_version"]=(a.mg5.resolve().parents[1]/"VERSION").read_text().strip()
    write(a.output,data)


def seal(a):
    data=read(a.input)
    keys={(r["shard"],int(r["event"])) for r in read(a.features)}
    rows=[r for r in data["events"] if (r["shard"],r["event"]) in keys]
    if {(r["shard"],r["event"]) for r in rows}!=keys or len(rows)!=len(keys):
        raise ValueError("Incomplete hard-event support for detector candidates")
    data["events"]=rows
    data["features_sha256"]=digest(a.features)
    write(a.output,data)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(required=True)
    for task,names,fn in (("collect","features core complement core_spectrum complement_spectrum trees library output",collect),
                          ("weight","input mg5 card work output",weight),
                          ("seal","input features output",seal)):
        p=sub.add_parser(task)
        p.set_defaults(fn=fn)
        for name in names.split():
            p.add_argument(name,type=Path)
    args=parser.parse_args()
    args.fn(args)
