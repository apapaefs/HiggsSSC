"""Reduce every disjoint LO event block in one EW integration using a bounded pool."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("hard","root","plugin","response","repo"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--cores",type=int,required=True)
    p.add_argument("--seed",type=int,required=True)
    p.add_argument("--block",type=int,default=500)
    args=p.parse_args()
    hard=args.hard.resolve()
    raw=json.loads((hard/"hard.json").read_text())
    root=args.root.resolve()
    root.mkdir(parents=True,exist_ok=True)
    ranges=[]
    # Verify and skip blocks successfully measured before this launch, including a short pilot.
    for existing in root.glob(hard.name+"-*/receipt.json"):
        rec=json.loads(existing.read_text())
        if rec["integration_id"] != raw["integration_id"]:
            raise ValueError("mismatched sample in "+str(existing))
        ranges.append((rec["source_first_event"],rec["source_stop_event"]))
    total=raw["generated_events"]
    missing=[]
    i=0
    for start,stop in sorted(ranges)+[(total,total)]:
        if not i<=start<=stop<=total:
            raise ValueError("overlapping or invalid event intervals")
        while i<start:
            n=min(args.block,start-i)
            missing.append((i,n))
            i+=n
        i=stop
    def shower(item):
        i,n=item
        out=root/f"{hard.name}-{i:06d}"
        out.mkdir(parents=True,exist_ok=True)
        cmd=[sys.executable,Path(__file__).with_name("ew_shower.py"),"--out",out,"--hard",hard,
             "--sample",raw["sample"],"--start",str(i),"--events",str(n),"--seed",str(args.seed+i),
             "--plugin",args.plugin,"--response",args.response,"--repo",args.repo]
        if i==0:
            cmd.append("--keep-root")
        with (out/"worker.log").open("w") as log:
            subprocess.run([str(x) for x in cmd],stdout=log,stderr=subprocess.STDOUT,check=True)
        print(out.name,"complete",flush=True)
    with ThreadPoolExecutor(max_workers=args.cores) as pool:
        list(pool.map(shower,missing))


if __name__=="__main__":
    main()
