"""Produce one disjoint LO QCD dijet Born-pT bin in an activated MG5/Herwig runtime."""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from string import Template

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hfourlepton import run_four_lepton_campaign as mg
from hfourlepton.run_reducible_controls import lhe_parts
from hgammagamma import ho_shower_completion as completion
from hgammagamma.run_gammagamma_ho_campaign import ROOT_METADATA_SCRIPT


def set_cuts(card, values):
    text = card.read_text()
    for key, value in values.items():
        pattern = rf"(?im)^([ \t]*)[^\s=#]+([ \t]*=[ \t]*{key}\b)"
        text, count = re.subn(pattern, lambda m: f"{m[1]}{value}{m[2]}", text)
        if count == 0:
            text += f"\n{value} = {key}\n"
        elif count != 1:
            raise ValueError(f"duplicate {key} in run card")
    card.write_text(text)


def check_include(path, cuts):
    text = path.read_text()
    for key, val in cuts.items():
        match = re.search(rf"(?im)^\s*{key}\s*=\s*([^\s!]+)", text)
        if match is None or float(match[1].lower().replace("d","e")) != val:
            raise ValueError(f"effective {key} differs from {val}: {path}")


def audit_source(src, n, low, high):
    count = 0
    for kind, block in lhe_parts(src):
        if kind != "event":
            continue
        header, *rest = block.splitlines()[1:]
        particles = [x.split() for x in rest[:int(header.split()[0])]]
        final = [p for p in particles if p[1]=="1"]
        if len(final)!=2 or any(abs(int(p[0])) not in (1,2,3,4,5,21) for p in final):
            raise ValueError("expected two final QCD partons")
        for p in final:
            pt=math.hypot(float(p[6]),float(p[7]))
            if pt < low-1e-4 or (high>0 and pt>high+1e-4):
                raise ValueError("generated event outside its Born bin")
        count+=1
    init=mg.parse_lhe_init_records(src)[0].split()
    if count!=n or list(map(int,init[:2]))!=[2212,2212] or list(map(float,init[2:4]))!=[20000.,20000.] or list(map(int,init[6:8]))!=[331900,331900]:
        raise ValueError("incorrect beam, PDF or source population")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dir",type=Path,required=True)
    p.add_argument("--low",type=float,required=True)
    p.add_argument("--high",type=float,required=True,help="-1 for an open upper edge")
    p.add_argument("--events",type=int,required=True)
    p.add_argument("--seed",type=int,required=True)
    p.add_argument("--cores",type=int,default=2)
    p.add_argument("--mg5-dir",type=Path,required=True)
    a=p.parse_args()
    if not (a.low>0 and (a.high>a.low or a.high==-1) and a.events>0 and a.cores>0):
        p.error("invalid bin, event or CPU count")
    dest=a.dir.resolve()
    dest.mkdir(parents=True,exist_ok=False)
    cfg=mg.parse_config(["--order","lo","--higgs-mass","125.09","--ebeam","20000","--run-tag","dijet",
        "--nevents",str(a.events),"--nb-core",str(a.cores),"--mg5-dir",str(a.mg5_dir),
        "--mg5-python",sys.executable,"--no-herwig-env","--no-herwig-module"])
    proc=dest/"mg5_process"
    proc_card=dest/"process_card.dat"
    spec=mg.SampleSpec("dijet","QCD dijets","background","jj","LO","mg5","sm-no_b_mass",
                       "p p > j j QED=0 QCD=2","shower","native","dijet","dijet")
    mg.write_mg5_process_card(proc_card,spec,proc,cfg)
    mg.run_runtime_command(mg.mg5_command(cfg,proc_card),cfg,cwd=dest)
    mg.patch_mg5_toolchain(proc,cfg)
    mg.patch_mg5_hidden_pdf_defaults(proc,cfg)
    mg.patch_mg5_pdf_label_sync(proc,cfg)
    card=proc/"Cards/run_card.dat"
    mg.patch_run_card(card,sample=spec,nevents=a.events,seed=a.seed,cfg=cfg)
    cuts=dict(ptj=a.low,ptjmax=a.high,etaj=-1.,drjj=0.,maxjetflavor=5)
    set_cuts(card,cuts)
    mg.generate_and_patch_run_card_include(proc,cfg)
    include=proc/"Source/run_card.inc"
    check_include(include,cuts)
    mg.run_generate_events(proc,None,cfg)
    check_include(include,cuts)
    src=mg.find_lhe_file(proc,cfg.run_tag)
    audit_source(src,a.events,a.low,a.high)
    source=completion.source_weight_summary(src,a.events)
    if source["sum_weight"]!=a.events or source["negative_events"]:
        raise ValueError("dijet envelope requires unweighted positive hard events")
    record=dict(bin_GeV=[a.low,a.high],events=a.events,seed=a.seed,cross_section_pb=mg.parse_lhe_init_cross_section_pb(src),
                integration_error_pb=float(mg.parse_lhe_init_records(src)[1].split()[1]),lhe=str(src),
                lhe_sha256=mg.sha256_file(src),process_card_sha256=mg.sha256_file(proc_card),
                run_card_sha256=mg.sha256_file(card),effective_include_sha256=mg.sha256_file(include))
    (dest/"source.json").write_text(json.dumps(record,indent=2)+"\n")
    shower(dest,src,source,a.events,a.seed+1,spec.pdf_name)


def shower(dest,src,source,events,seed,pdf, *, disable_mpi=False):
    d=dest/"herwig"
    (d/"events").mkdir(parents=True)
    template=Path(__file__).resolve().parent/"HOAnalysis/HW-LHE.in"
    plan=completion.finite_shower_plan("dijet",events,events)
    (d/"shower_execution.json").write_text(json.dumps(plan))
    cardtext=Template(template.read_text()).substitute(lhe_file=src,pdf_name=pdf,
        nevents=events,seed=seed,max_errors=10000,matching_settings=("set /Herwig/Shower/ShowerHandler:MPIHandler NULL" if disable_mpi else ""),
        decay_settings="",sample_name="dijet")
    cardtext = cardtext.replace("saverun dijet", "get /Herwig/Shower/ShowerHandler:MPIHandler\nsaverun dijet")
    (d/"dijet.in").write_text(cardtext.replace("saverun dijet","set EventGenerator:DumpPeriod -1\nsaverun dijet"))
    for args,log,check in [(["Herwig","read","dijet.in"],"read.log",True),
                           (["Herwig","run","dijet.run",f"-N{plan['execution_events']}"],"run.log",False)]:
        with (d/log).open("w") as f:
            result=subprocess.run(args,cwd=d,stdout=f,stderr=subprocess.STDOUT,check=check)
    roots=sorted((d/"events").rglob("*.root"))
    metadata=json.loads(subprocess.check_output([sys.executable,"-c",ROOT_METADATA_SCRIPT,*roots],text=True).splitlines()[-1])
    proof=completion.build_completion(d,"dijet",events,events,metadata,
        "source_exhausted" if result.returncode else "requested_events",source_weights=source)
    (dest/"shower.json").write_text(json.dumps(proof,indent=2)+"\n")
    print("Finite shower complete",dest,flush=True)


if __name__=="__main__":
    main()
