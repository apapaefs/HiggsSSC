"""Measure the signal and validate inclusive rates event by event with the existing SSC executable."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import subprocess
from tempfile import TemporaryDirectory
import numpy as np
import uproot


def run(pair):
    folder, binary, inclusive, work = pair
    work = work / folder.name
    work.mkdir()
    receipt = json.loads((folder/"receipt.json").read_text())
    root = folder/"events/tagged.root"
    in_file = work/"inclusive.input"
    in_file.write_text(str(root) + "\n")
    seed = str(receipt["seed"] + 99901)
    with (work/"inclusive.log").open("w") as log:
        subprocess.run([inclusive,str(in_file),"--response-mode","genuine","--seed",seed],check=True,
                       stdout=log,stderr=subprocess.STDOUT)
    csv = work/"selected.csv"
    subprocess.run([binary,root,csv,seed],check=True)
    events = np.atleast_1d(np.genfromtxt(csv,delimiter=",",names=True))
    if len(events) != receipt["saved_events"] or not np.isclose(events["evweight"].sum(), receipt["signed_saved"]):
        raise ValueError("Incomplete shower event records")
    with uproot.open(work/"inclusive_var.root") as f:
        t = f["Data2"]
        var = t["variables"].array(library="np")
        mask = var[:,9]>=2
        source = t["sourceevent"].array(library="np").ravel()[mask]
        np.testing.assert_array_equal(source, np.flatnonzero(events["p_inclusive"]>0))
        np.testing.assert_allclose(var[mask,0],events["mass_inclusive"][source],rtol=1.e-7)
        np.testing.assert_allclose(t["eventweight"].array(library="np").ravel()[mask],
                                   events["evweight"][source]*events["p_inclusive"][source],rtol=1.e-7)
    return events, receipt


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("showers", type=Path)
    p.add_argument("measurement_exe", type=Path)
    p.add_argument("inclusive_exe", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--workers", type=int, default=12)
    args = p.parse_args()
    names = sorted(x.parent for x in args.showers.glob("*/receipt.json"))
    if len(names)!=45: raise ValueError("The release comprises 45 disjoint 1000-event jobs")
    with TemporaryDirectory(prefix="ssc-tag-measure-") as temporary:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            jobs = [(p, args.measurement_exe.resolve(), args.inclusive_exe.resolve(), Path(temporary)) for p in names]
            results = list(pool.map(run, jobs))
    inputs = json.loads((args.showers/"inputs.json").read_text())
    meta = []
    events = []
    for i,(folder,(rows,r)) in enumerate(zip(names,results)):
        meta.append(dict(name=folder.name, signed_attempted=r["signed_attempted"],
                    generated=r["attempted_events"], saved=r["saved_events"], signed_saved=r["signed_saved"],
                    seed=r["seed"], parent=r["parent"]))
        events.append(np.column_stack((np.full(len(rows),i),*[rows[k] for k in rows.dtype.names])))
    columns = ["shard", *results[0][0].dtype.names]
    np.savez_compressed(args.output, values=np.vstack(events), columns=columns,
                        metadata=json.dumps(dict(shards=meta,**inputs), sort_keys=True))


if __name__=="__main__":
    main()
