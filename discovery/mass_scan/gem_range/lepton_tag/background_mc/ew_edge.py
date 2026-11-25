"""Run a complementary EW generation and shower check in a separate validation directory."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("nominal","work","plugin","response","repo","mg5"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--after",required=True,help="Final nominal integration whose CPU allocation is transferred")
    p.add_argument("--sample",required=True,choices=("w_lnu_aa","z_ll_aa","z_ee_a"))
    for key in ("seed","shower-seed","cores","events"):
        p.add_argument("--"+key,type=int,required=True)
    a=p.parse_args()
    nominal,work=a.nominal.resolve(),a.work.resolve()
    if nominal.parent==work or nominal.parent in work.parents:
        p.error("place validation outside the full nominal campaign tree")
    job=work/(a.sample+"-edge-001")
    here=Path(__file__).resolve().parent
    while True:
        base=nominal/a.after/"hard.json"
        if base.is_file():
            counts=sum(json.loads(f.read_text())["attempted_events"] for f in
                       (nominal/"showers").glob(a.after+"-*/receipt.json"))
            if counts==json.loads(base.read_text())["generated_events"]:
                break
        time.sleep(6)
    if not (job/"hard.json").exists():
        cmd=[sys.executable,here/"ew_generate.py","--mg5",a.mg5,"--process",a.sample,
             "--work",job,"--events",a.events,"--seed",a.seed,"--cores",a.cores,"--outside-nominal"]
        if (job/"madgraph/Events/hard/unweighted_events.lhe.gz").is_file():
            cmd.append("--resume-completed")
        subprocess.run(list(map(str,cmd)),check=True)
    cmd=[sys.executable,here/"ew_campaign.py","--hard",job,"--root",work/"showers","--plugin",a.plugin,
         "--response",a.response,"--repo",a.repo,"--cores",a.cores,"--block",1500,
         "--seed",a.shower_seed]
    subprocess.run(list(map(str,cmd)),check=True)


if __name__=="__main__":
    main()
