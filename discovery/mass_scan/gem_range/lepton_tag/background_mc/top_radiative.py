"""Generate resolved radiative top-decay strata with loose lab-frame cuts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import shutil
import subprocess
import sys

SAMPLES = ("tt_0a_1dec", "tt_1a_1dec", "tt_0a_2dec_same", "tt_0a_2dec_split", "tt_0a_1dec_ll")


def massive_lepton_model(work, mg5):
    """Use the standard SM with massless initial b and physical charged leptons."""
    destination = work / "sm_massive_leptons"
    original = mg5.resolve().parent.parent / "models/sm"
    if not destination.is_dir():
        shutil.copytree(original, destination, ignore=shutil.ignore_patterns("*.pkl", "*.pyc", "__pycache__"))
    restriction = (original / "restrict_lepton_masses.dat").read_text()
    restriction, n = re.subn(r"(?im)^(\s*5\s+)\S+(\s*#\s*y?mb\b)", r"\g<1>0.000000e+00\2", restriction)
    if n != 2:
        raise ValueError("unexpected SM mass/Yukawa restriction")
    (destination / "restrict_default.dat").write_text(restriction)
    return destination


def generate(i, args):
    sample = SAMPLES[i]
    work = args.out / sample / f"round-{args.round:03}"
    work.mkdir(parents=True, exist_ok=True)
    if (work / "hard.json").is_file():
        return
    # The t and internal W poles decay through their leading-order SM amplitudes.
    # No cuts on their quark or lepton daughters beyond the QED collinear regulator.
    cuts = dict(pta=14., etaa=3.5, cut_decays="True", ptj=0., ptb=0., ptl=0.,
                etaj=-1., etab=-1., etal=-1., drjj=0., drbj=0., drbb=0., drjl=0.,
                drbl=0., drll=0., mmll=0., mmjj=0., draj=args.dr, drab=args.dr,
                dral=0., draa=0., ptgmin=0., maxjetflavor=5)
    if sample in ("tt_1a_1dec", "tt_0a_2dec_same", "tt_0a_2dec_split"):
        cuts.update(mmaa=60., mmaamax=170.)
    options = [v for key, value in cuts.items() for v in ("--set", f"{key}={value}")]
    code = Path(__file__).resolve().parent
    model = massive_lepton_model(work, args.mg5)
    card = work / "radiative.mg5"
    card.write_text((code / "top_cards" / (sample + ".mg5")).read_text().replace(
        "import model sm-no_b_mass", f"import model {model.resolve()}"))
    command = [sys.executable, args.generator, "--mg5", args.mg5, "--process", sample,
               "--work", work, "--events", str(args.events), "--cores", str(args.cores // len(args.samples)),
               "--seed", str(751000 + 100 * args.round + i), "--pdf", "335900", "--card",
               card, *options]
    if (work / "madgraph/Events/hard/unweighted_events.lhe.gz").is_file():
        command.append("--resume-completed")
    with (work / "driver.log").open("a") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("out", "generator", "mg5"):
        p.add_argument("--" + name, required=True, type=Path)
    p.add_argument("--cores", type=int, default=36)
    p.add_argument("--events", type=int, default=10000)
    p.add_argument("--dr", type=float, default=.05)
    p.add_argument("--round", type=int, default=0)
    p.add_argument("--samples", nargs="+", choices=SAMPLES, default=list(SAMPLES[:4]))
    a = p.parse_args()
    if a.cores < len(a.samples) or len(set(a.samples)) != len(a.samples) or a.round < 0 or not 0 < a.dr < 1. or a.events <= 0:
        p.error("invalid sample list, resource count or generation cuts")
    with ThreadPoolExecutor(max_workers=len(a.samples)) as pool:
        jobs = [pool.submit(generate, SAMPLES.index(s), a) for s in a.samples]
        for job in jobs:
            job.result()


if __name__ == "__main__":
    main()
