"""Generate one independent 40 TeV MG5 LO integration and positive-weight LHE file."""
import argparse
import ast
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def set_card(path, params):
    data = path.read_text()
    for key, value in params.items():
        if not re.fullmatch(r"[a-z][a-z0-9_]*", key):
            raise ValueError(key)
        pattern = rf"(?im)^[^#!\n=]+(=\s*{key}\b[^\n]*)$"
        data, n = re.subn(pattern, lambda m: f" {value} {m[1]}", data)
        if n > 1:
            raise ValueError(f"Repeated {key}")
        if not n:
            data += f"\n {value} = {key}\n"
    path.write_text(data)


def run(command, path, log):
    with log.open("w") as stream:
        subprocess.run([str(v) for v in command], stdin=subprocess.DEVNULL, stdout=stream,
                       stderr=subprocess.STDOUT, cwd=path, check=True)


def lhe_info(path):
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.strip() == "<init>":
                break
        header = next(stream).split()
        if [int(h) for h in header[:2]] != [2212,2212] or [float(h) for h in header[2:4]] != [20000.,20000.]:
            raise ValueError("LHE beam mismatch")
        rows = [next(stream).split() for _ in range(int(header[-1]))]
        for line in stream:
            if line.strip() == "</init>":
                break
            if line.strip() and not line.lstrip().startswith("<generator "):
                raise ValueError(f"unexpected LHE initialization: {line}")
        else:
            raise ValueError("unterminated LHE initialization")
        n = sum(1 for line in stream if line.strip() == "<event>")
    return dict(cross_section_pb=sum(float(row[0]) for row in rows),
                cross_section_error_pb=math.sqrt(sum(float(row[1])**2 for row in rows)),
                lhe_init=header, generated_events=n)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mg5", required=True, type=Path)
    p.add_argument("--process", required=True)
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--seed", required=True, type=int)
    p.add_argument("--events", required=True, type=int)
    p.add_argument("--cores", required=True, type=int)
    p.add_argument("--pdf", default=335900, type=int)
    p.add_argument("--card", type=Path)
    p.add_argument("--resume-completed", action="store_true", help="Build receipt from an existing finished integration")
    p.add_argument("--outside-nominal", action="store_true", help="Generate the complementary looser-cut control for EW processes")
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    a = p.parse_args()
    work, mg5 = a.work.resolve(), a.mg5.resolve()
    if not re.fullmatch(r"[a-z0-9_]+", a.process) or min(a.seed,a.events,a.cores) <= 0:
        p.error("invalid sample, seed, event count or cores")
    card = a.card or Path(__file__).with_name("ew_cards")/(a.process+".dat")
    commands = card.read_text()
    work.mkdir(parents=True, exist_ok=True)
    target = work/"madgraph"
    params = dict(nevents=a.events, iseed=a.seed, ebeam1=20000., ebeam2=20000.,
                  pdlabel="lhapdf", lhaid=a.pdf, use_syst="False", systematics_program="none")
    if not a.card:
        params.update(pta=14., etaa=3., dral=0., ptl=0., etal=-1., drll=0., mmll=0.)
        if a.process == "z_ee_a":
            params.update(ptl=14., etal=3., dral=0.20)
        else:
            params.update(draa=0., mmaa=60., mmaamax=170.)
        if a.process == "w_lnu_aa":
            params.update(ptl=13., etal=3.5)
        if a.process == "z_ll_aa":
            params.update(xptl=13.)
    if a.outside_nominal:
        if a.card or a.set:
            p.error("complement is defined only for the standard EW cards and cuts")
        custom=Path(__file__).resolve().with_name("ew_cards")/"outside_nominal.f"
        params.update(pta=8.,etaa=4.,custom_fcts=str([str(custom)]))
        if a.process=="z_ee_a":
            params.update(ptl=8.,etal=4.,dral=.08)
        else:
            params.update(mmaa=50.,mmaamax=180.)
            if a.process=="w_lnu_aa":
                params.update(ptl=8.,etal=4.5)
            else:
                params.update(xptl=8.)
    params.update(dict(x.split("=",1) for x in a.set))
    if (work/"hard.json").exists():
        raise FileExistsError(work/"hard.json")
    if not a.resume_completed:
        if target.exists():
            raise FileExistsError(f"Existing process directory {target}; choose a new integration path")
        (work/"proc.mg5").write_text(commands.rstrip()+"\noutput "+str(target)+" -f\n")
        run([sys.executable,mg5,work/"proc.mg5"],work,work/"export.log")
        if not (target/"bin/generate_events").is_file():
            raise RuntimeError(f"MG5 export failed; see {work/'export.log'}")
        if path := shutil.which("lhapdf-config"):
            with (target/"Cards/me5_configuration.txt").open("a") as stream:
                stream.write(f"\nlhapdf_py3 = {path}\n")
        set_card(target/"Cards/run_card.dat",params)
        run([sys.executable,target/"bin/generate_events","hard","-f",f"--nb_core={a.cores}"],target,work/"generate.log")
    lhe = target/"Events/hard/unweighted_events.lhe.gz"
    info = lhe_info(lhe)
    if info["generated_events"] != a.events or not info["cross_section_pb"] > 0:
        raise RuntimeError(f"incomplete LHE {info}")
    final = work/"events.lhe.gz"
    if final.exists():
        if not final.samefile(lhe):
            raise ValueError("source identity changed")
    else:
        os.link(lhe,final)
    # Validate reported settings against the actual generated cards when recovering.
    banner=target/"Events/hard/hard_tag_1_banner.txt"
    actual=banner.read_text().split("<MGRunCard>",1)[1].split("</MGRunCard>",1)[0]
    for key,value in params.items():
        m=re.search(rf"(?im)^\s*([^#!\n=]+?)\s*=\s*{key}\b",actual)
        if m is None:
            raise ValueError(f"missing archived run-card field {key}")
        try:
            same = float(m[1]) == float(value)
        except ValueError:
            same = m[1].strip().lower() == str(value).lower()
        if key=="custom_fcts":
            expected=ast.literal_eval(value)
            same=m[1].strip()==expected[0] if len(expected)==1 else False
        if not same:
            raise ValueError(f"archived run-card mismatch for {key}={value}")
    (work/"run_card.dat").write_text(actual.strip()+"\n")
    shutil.copy2(banner,work/"banner.txt")
    for name in ("param_card.dat","proc_card_mg5.dat"):
        shutil.copy2(target/"Cards"/name,work/name)
    info.update(sample=a.process, hard_seed=a.seed, run_settings=params, lhe_sha256=digest(final),
                run_card_sha256=digest(work/"run_card.dat"), proc_card_sha256=digest(work/"proc_card_mg5.dat"),
                mg5_version=(mg5.parent.parent/"VERSION").read_text().strip())
    if a.outside_nominal:
        from ew_diagnostics import verify_outside
        if verify_outside(lhe,a.process)!=info["generated_events"]:
            raise ValueError("complement event count mismatch")
        shutil.copy2(custom,work/"outside_nominal.f")
        info.update(generation_region="outside_nominal",cut_function_sha256=digest(custom))
    info["integration_id"] = f"{a.process}-{a.seed}-{info['lhe_sha256'][:16]}"
    tmp=work/"hard.json.tmp"
    tmp.write_text(json.dumps(info,indent=2)+"\n")
    tmp.replace(work/"hard.json")
    print(json.dumps(info),flush=True)


if __name__ == "__main__":
    main()
