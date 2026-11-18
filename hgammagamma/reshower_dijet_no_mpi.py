"""Shower one verified dijet bin with multipartonic interactions disabled."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hgammagamma.run_dijet_bound import shower, completion, mg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("destination", type=Path)
    args=parser.parse_args()
    src=json.loads((args.original/"source.json").read_text())
    proof=json.loads((args.original/"shower.json").read_text())
    n=src["events"]
    completion.validate_completion(proof,dict(sample="dijet",nevents_requested=n,lhe=dict(events=n)),directory=args.original)
    lhe=Path(src["lhe"])
    if mg.sha256_file(lhe)!=src["lhe_sha256"] or src["lhe_sha256"]!=proof["source_weights"]["artifact"]["sha256"]:
        raise ValueError("hard events differ from the original shower")
    args.destination.mkdir(parents=True,exist_ok=False)
    src.update(shower_variation="MPI disabled",baseline_shower_fingerprint=proof["fingerprint"])
    (args.destination/"source.json").write_text(json.dumps(src,indent=2)+"\n")
    shower(args.destination.resolve(),lhe,proof["source_weights"],n,src["seed"]+1,"NNPDF40_lo_as_01180",disable_mpi=True)


if __name__=="__main__":
    main()
