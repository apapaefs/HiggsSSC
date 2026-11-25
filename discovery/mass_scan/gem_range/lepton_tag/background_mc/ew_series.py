"""Run fixed-size independent hard integrations and disjoint Herwig showers in one CPU allocation."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample",required=True)
    for key in ("work","plugin","response","repo","mg5"):
        p.add_argument("--"+key,type=Path,required=True)
    for key in ("cores","seed","shower-seed","first","last","events"):
        p.add_argument("--"+key,type=int,required=True)
    p.add_argument("--block",type=int,default=1250)
    p.add_argument("--after",help="Wait for an earlier integration's disjoint intervals to finish")
    a=p.parse_args()
    here=Path(__file__).resolve().parent
    work=a.work.resolve()
    packets=work/"showers"
    if a.after:
        while True:
            receipts=[json.loads(f.read_text()) for f in packets.glob(a.after+"-*/receipt.json")]
            completed=sum(r["attempted_events"] for r in receipts)
            meta=work/a.after/"hard.json"
            if meta.is_file() and completed==json.loads(meta.read_text())["generated_events"]:
                break
            time.sleep(5)
    for i in range(a.first,a.last+1):
        job=work/f"{a.sample}-{i:03d}"
        if not (job/"hard.json").exists():
            cmd=[sys.executable,here/"ew_generate.py","--mg5",a.mg5,"--process",a.sample,
                 "--work",job,"--events",a.events,"--seed",a.seed+i,"--cores",a.cores]
            if (job/"madgraph/Events/hard/unweighted_events.lhe.gz").is_file():
                cmd.append("--resume-completed")
            subprocess.run(list(map(str,cmd)),check=True)
        cmd=[sys.executable,here/"ew_campaign.py","--hard",job,"--root",packets,"--plugin",a.plugin,
             "--response",a.response,"--repo",a.repo,"--cores",a.cores,"--block",a.block,
             "--seed",a.shower_seed+i*100000]
        subprocess.run(list(map(str,cmd)),check=True)


if __name__=="__main__":
    main()
