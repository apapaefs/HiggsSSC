"""Extract the 125 GeV pre-isolation assignments and reproduce the frozen nominal packets."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_shard(item):
    tree_record, shard, executable = item
    tree = json.loads(tree_record.read_text())
    receipt_path = Path(tree["receipt"])
    receipt = json.loads(receipt_path.read_text())
    packet = receipt_path.with_name("packet.csv")
    if (digest(receipt_path) != shard["receipt_sha256"] or
            digest(packet) != shard["packet_sha256"]):
        raise ValueError(f"Frozen receipt or packet changed for {tree_record}")
    with tempfile.NamedTemporaryFile(suffix=".jsonl") as out:
        run = subprocess.run([executable, tree["root"], out.name, str(receipt["response_seed"])],
                             capture_output=True, text=True, check=True)
        if run.stderr.strip().splitlines()[-1] != f"{receipt['saved_events']} saved collisions":
            raise ValueError(f"Missing replayed events in {tree_record}")
        rows = [json.loads(line) for line in out]
    model = sorted((r["event"], r["mass"], r["nominal_probability"]) for r in rows if r["nominal_probability"])
    with packet.open() as fd:
        frozen = sorted((int(r["source_event"]), float(r["mass_gev"]), float(r["probability"]))
                        for r in csv.DictReader(fd) if r["component"] == "electron_hard" and
                        123.3 <= float(r["mass_gev"]) < 126.7)
    if len(model) != len(frozen):
        raise ValueError(f"Changed nominal count in {tree_record}")
    for (i,m,p),(j,n,q) in zip(model, frozen):
        if i != j or not math.isclose(m,n,rel_tol=1.e-12) or not math.isclose(p,q,rel_tol=1.e-12):
            raise ValueError(f"Nominal pair mismatch for {tree_record} event {i}")
    for row in rows:
        row["shard"] = tree_record.parent.name
        if all(row[k]["cone"]<5. for k in ("real","fake","tag")) != bool(row["nominal_probability"]):
            raise ValueError("Continuous scalar-cone comparison differs from the frozen selection")
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("core_spectrum","complement_spectrum","trees","executable","output"):
        p.add_argument(name,type=Path)
    p.add_argument("--workers",type=int,default=16)
    a=p.parse_args()
    jobs=[]
    for file in (a.core_spectrum,a.complement_spectrum):
        with gzip.open(file,"rt") as stream:
            for shard in json.load(stream)["shards"]:
                if shard["sample"]=="z_ee_a":
                    tree=a.trees/Path(shard["receipt"]).parent.name/"tree.json"
                    jobs.append((tree,shard,a.executable.resolve()))
    rows=[]
    with ThreadPoolExecutor(a.workers) as pool:
        for i,items in enumerate(pool.map(run_shard,jobs),1):
            rows.extend(items)
            if i%50==0: print(i,len(jobs),len(rows),flush=True)
    ids=lambda rows: {(r["shard"],r["event"]) for r in rows}
    nominal=[r for r in rows if r["nominal_probability"]]
    if len(jobs)!=498 or len(ids(nominal))!=149:
        raise ValueError("Frozen eeγ campaign attempt coverage or nominal Higgs window changed")
    packed=json.dumps(rows,allow_nan=False,separators=(",",":")).encode()
    with a.output.open("wb") as file, gzip.GzipFile(filename="",mode="wb",fileobj=file,mtime=0) as out:
        out.write(packed)
    print(f"{len(rows)} possible pairs, {len(ids(rows))} independent collisions, "
          f"{len(nominal)} nominal pairs in {len(ids(nominal))} independent collisions")

if __name__=="__main__":
    main()
