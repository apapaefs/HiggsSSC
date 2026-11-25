"""Generate resolved radiative top-decay strata with loose lab-frame cuts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import errno
import json
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


def check_collar(work, lower):
    from top_diagnostics import separations
    distances = separations(work / "events.lhe.gz")
    if min(distances) < lower-1.e-5 or max(distances) >= .05+1.e-5:
        raise ValueError("hard events violate the exclusive photon-quark collar")
    result = dict(events=len(distances), minimum=min(distances), maximum=max(distances))
    (work / "collar-check.json").write_text(json.dumps(result) + "\n")


def generate(i, args):
    sample = SAMPLES[i]
    work = args.out / sample / f"round-{args.round:03}"
    work.mkdir(parents=True, exist_ok=True)
    if (work / "hard.json").is_file():
        info = json.loads((work / "hard.json").read_text())
        settings = info["run_settings"]
        if (info["sample"] != sample or info["generated_events"] != args.events
                or info["hard_seed"] != 751000 + 100 * args.round + i
                or any(float(settings[k]) != args.dr for k in ("draj", "drab"))
                or args.quark_collar != ("quark_photon_collar.f" in str(settings.get("custom_fcts", "")))):
            raise ValueError(f"Existing integration has different settings in {work}")
        if args.quark_collar and not (work / "collar-check.json").is_file():
            check_collar(work, args.dr)
        return
    base = work.parent / "round-000"
    if args.round > 0 and not args.fresh_integral and (base / "hard.json").is_file():
        settings = json.loads((base / "hard.json").read_text())["run_settings"]
        if not all(float(settings[k]) == args.dr for k in ("draj", "drab")):
            raise ValueError("a changed generation cut needs --fresh-integral")
        from top_gridpack import generate as continue_grid
        try:
            continue_grid(base, work, args.events, 751000 + 100 * args.round + i, args.cores // len(args.samples),
                          args.native_scratch)
            return
        except (OSError, ValueError, subprocess.CalledProcessError) as error:
            if isinstance(error, OSError) and error.errno == errno.ENOSPC:
                raise
            # A failed grid or normalization check cannot contribute events.
            # Start an independent integration while preserving its diagnostic.
            (work / "gridpack-failure.txt").write_text(str(error) + "\n")
            if (work / "events.lhe.gz").is_file():
                (work / "events.lhe.gz").rename(work / "unvalidated-gridrun.lhe.gz")
    # The t and internal W poles decay through their leading-order SM amplitudes.
    # No cuts on their quark or lepton daughters beyond the QED collinear regulator.
    cuts = dict(pta=14., etaa=3.5, cut_decays="True", ptj=0., ptb=0., ptl=0.,
                etaj=-1., etab=-1., etal=-1., drjj=0., drbj=0., drbb=0., drjl=0.,
                drbl=0., drll=0., mmll=0., mmjj=0., draj=args.dr, drab=args.dr,
                dral=0., draa=0., ptgmin=0., maxjetflavor=5)
    if sample in ("tt_1a_1dec", "tt_0a_2dec_same", "tt_0a_2dec_split"):
        cuts.update(mmaa=60., mmaamax=170.)
    code = Path(__file__).resolve().parent
    if args.quark_collar:
        cuts["custom_fcts"] = str([str(code / "top_cards/quark_photon_collar.f")])
    options = [v for key, value in cuts.items() for v in ("--set", f"{key}={value}")]
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
    if args.quark_collar:
        check_collar(work, args.dr)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("out", "generator", "mg5"):
        p.add_argument("--" + name, required=True, type=Path)
    p.add_argument("--cores", type=int, default=36)
    p.add_argument("--events", type=int, default=10000)
    p.add_argument("--dr", type=float, default=.05)
    p.add_argument("--quark-collar", action="store_true", help="Exclusive lower-cut variation with min dR(gamma,quark) < 0.05")
    p.add_argument("--round", type=int, default=0)
    p.add_argument("--fresh-integral", action="store_true", help="Independently integrate this round")
    p.add_argument("--native-scratch", type=Path, help="Private temporary native grid directories on fast scratch")
    p.add_argument("--samples", nargs="+", choices=SAMPLES, default=list(SAMPLES[:4]))
    a = p.parse_args()
    if a.cores < len(a.samples) or len(set(a.samples)) != len(a.samples) or a.round < 0 or not 0 < a.dr < 1. or a.events <= 0:
        p.error("invalid sample list, resource count or generation cuts")
    if a.quark_collar and (not a.fresh_integral or a.dr >= .05):
        p.error("the quark collar needs --fresh-integral and a lower --dr than 0.05")
    with ThreadPoolExecutor(max_workers=len(a.samples)) as pool:
        jobs = [pool.submit(generate, SAMPLES.index(s), a) for s in a.samples]
        for job in jobs:
            job.result()


if __name__ == "__main__":
    main()
