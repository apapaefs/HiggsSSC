"""Measure the eeγ detector cutflow with a collision-level covariance and a frozen-packet closure."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from replay_showers import records

STAGES = [
    "truth_matched_reconstructed_mass", "real_photon_kinematics",
    "fake_electron_kinematics", "tag_electron_kinematics",
    "all_em_fiducial", "real_photon_isolation", "fake_electron_isolation",
    "tag_electron_isolation", "pair_pt_above_40", "real_photon_id",
    "tag_electron_id", "electron_fake_id_and_exclusivity", "trigger",
    "all_ids_before_pair_pt", "final_fake_electron", "final_fake_positron",
    "final_both_charges_same_collision"
]
BANDS = {"higgs_window": (123.3,126.7), "wide_mass": (75.,165.)}


def measure(job, out, executable, executable_sha):
    path, _, sha, _ = job
    d = out / path.parent.name
    while not (d/"tree.json").is_file():
        time.sleep(2)
    tree = Path(json.loads((d/"tree.json").read_text())["root"])
    r = json.loads(path.read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest() != sha:
        raise ValueError(f"Changed frozen receipt {path}")
    packet=path.with_name("packet.csv")
    if hashlib.sha256(packet.read_bytes()).hexdigest()!=r["packet_sha256"]:
        raise ValueError(f"Changed published detector outcomes {packet}")
    table=d/"cutflow-models.csv"
    meta=table.with_suffix(".provenance.json")
    cache=dict(response_executable_sha256=executable_sha, frozen_receipt_sha256=sha)
    if table.is_file() and meta.is_file():
        cache["table_sha256"]=hashlib.sha256(table.read_bytes()).hexdigest()
    if not meta.is_file() or json.loads(meta.read_text())!=cache or not table.is_file():
        tmp=table.with_suffix(".partial")
        proc=subprocess.run([executable,tree,tmp,str(r["response_seed"])],check=True,
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        if f"{r['saved_events']} saved collisions" not in proc.stderr:
            raise ValueError(f"Replay lost showers in {path}")
        os.replace(tmp,table)
        cache["table_sha256"]=hashlib.sha256(table.read_bytes()).hexdigest()
        meta.write_text(json.dumps(cache)+"\n")
    n=r["attempted_events"]
    results = np.zeros((3,2,n,len(STAGES)))
    keyed = defaultdict(list)
    with table.open() as stream:
        for row in csv.DictReader(stream):
            ev=int(row["event"]); mode=int(row["mode"])
            mass=float(row["mass"]); ncut=int(row["passed"]); prob=float(row["probability"])
            pdg=int(row["pdg_fake"]); tags=int(row["eligible_tags"])
            if ncut==8 and prob>0 and mode==1:
                keyed[ev].append((mass,prob))
            for ib,(lo,hi) in enumerate(BANDS.values()):
                if not lo<=mass<hi:continue
                event=results[mode,ib,ev]
                event[:ncut+1]+=1
                if ncut>=7:
                    event[13]+=prob
                if ncut==8:
                    event[9]+=.90*.96
                    event[10]+=.90*.96*(1.-(1.-.90)**tags)
                    event[11]+=prob/.988
                    event[12]+=prob
                    event[14 if pdg==11 else 15]+=prob
    reference=defaultdict(list)
    with packet.open() as stream:
        for row in csv.DictReader(stream):
            m=float(row["mass_gev"])
            if row["component"]=="electron_hard" and 75.<=m<165.:
                reference[int(row["source_event"])].append((m,float(row["probability"])))
    if set(keyed)!=set(reference):
        raise ValueError(f"Replay and published MC differ in selected event identities {path}")
    for event,rows in keyed.items():
        np.testing.assert_allclose(sorted(rows),sorted(reference[event]),rtol=3.e-11,atol=1.e-14,
                                   err_msg=f"Different selected pair {event} in {path}")
    results[:,:,:,16]=np.where((results[:,:,:,14]>0)&(results[:,:,:,15]>0),results[:,:,:,12],0.)
    nwin=sum(any(123.3<=m<126.7 for m,p in h) for h in keyed.values())
    return (r["integration_id"],n,r["saved_events"],1000.*r["cross_section_pb"],1000.*r["cross_section_error_pb"],
            len(keyed),nwin,results.sum(axis=2),np.einsum("abij,abik->abjk",results,results))


def summarize(measurements):
    grouped=defaultdict(list)
    for m in measurements: grouped[m[0]].append(m)
    nall=sum(m[1] for m in measurements)
    rates=np.zeros((3,2,len(STAGES)))
    cov=np.zeros((3,2,len(STAGES),len(STAGES)))
    for source,group in grouped.items():
        n=sum(m[1] for m in group)
        vals=sum(m[7] for m in group)
        outer=sum(m[8] for m in group)
        sigma,dx=group[0][3:5]
        if any(tuple(m[3:5])!=(sigma,dx) for m in group):
            raise ValueError(f"Changed normalization for integration {source}")
        rates+=sigma/nall*vals
        product=np.einsum("abi,abj->abij",vals,vals)
        cov+=(sigma/nall)**2*n/(n-1)*(outer-product/n)+(dx/nall)**2*product
    return dict(generated=nall,saved=sum(m[2] for m in measurements),matched_events_75_165=sum(m[5] for m in measurements),
                matched_events_123p3_126p7=sum(m[6] for m in measurements),
                generated_cross_section_fb=sum(m[1]*m[3] for m in measurements)/nall,
                generated_cross_section_mc_error_fb=np.linalg.norm([
                    group[0][4]*sum(m[1] for m in group)/nall for group in grouped.values()])),rates,cov


def tables(rates,cov):
    result={}
    for ismeared,label in enumerate(("unsmeared","nominal","precursor_selection")):
        result[label]={}
        for ib,band in enumerate(BANDS):
            r,c=rates[ismeared,ib],cov[ismeared,ib]
            result[label][band]=dict(rate_fb=r.tolist(), mc_se_fb=np.sqrt(np.maximum(np.diag(c),0.)).tolist(),
                                      covariance_fb2=c.tolist())
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("core","core_packets","core_hard","wide","wide_packets","wide_hard","trees","executable","output"):
        p.add_argument(key,type=Path)
    p.add_argument("--workers",type=int,default=8)
    a=p.parse_args()
    full_rates=np.zeros((3,2,len(STAGES)))
    full_cov=np.zeros((3,2,len(STAGES),len(STAGES)))
    result=dict(stages=STAGES, mass_ranges_gev=BANDS,
                baseline="Candidate assignments in the indicated mass range with one truth-matched prompt photon and separate fake and tag electrons in the saved stable-particle record. Pre-ID rates count assignments.",
                core_sha256=hashlib.sha256(a.core.read_bytes()).hexdigest(),
                wider_sha256=hashlib.sha256(a.wide.read_bytes()).hexdigest(),regions={})
    executable_sha=hashlib.sha256(a.executable.read_bytes()).hexdigest()
    for label,snap,root,hard in (("core",a.core,a.core_packets,a.core_hard),("wider_complement",a.wide,a.wide_packets,a.wide_hard)):
        jobs=list(records(snap,root,hard))
        findings=[]
        with ThreadPoolExecutor(a.workers) as pool:
            for i,m in enumerate(pool.map(lambda j:measure(j,a.trees,a.executable,executable_sha),jobs),1):
                findings.append(m)
                if i%25==0 or i==len(jobs): print(label,i,len(jobs),flush=True)
        record,rates,cov=summarize(findings)
        source=dict(record,shower_shards=len(jobs),**tables(rates,cov))
        result["regions"][label]=source
        full_rates+=rates
        full_cov+=cov
        print(label,source["nominal"]["higgs_window"]["rate_fb"][12],flush=True)
    result["combined"]=tables(full_rates,full_cov)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
