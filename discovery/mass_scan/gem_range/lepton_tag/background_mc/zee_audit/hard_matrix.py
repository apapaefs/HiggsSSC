"""Evaluate the production Z-pole model on the original fully reconstructed eeγ collisions."""
import argparse
import gzip
import hashlib
import json
import math
import re
from pathlib import Path
import subprocess
import sys

import numpy as np


def matrix_sources(mg5, work):
    commands = ["import model sm-lepton_masses"]
    for model, final in (("full","e+ e- a QCD=0"), ("production_z","z a QCD=0, z > e+ e-")):
        commands += [f"generate u u~ > {final}",f"add process d d~ > {final}",
                     f"output standalone_cpp {work/model} -f"]
    script = work/"export.mg5"
    script.write_text("\n".join(commands)+"\n")
    with (work/"export.log").open("w") as stream:
        subprocess.run([sys.executable,mg5,script],stdout=stream,stderr=subprocess.STDOUT,check=True)


# In the MG5 3.5.15 standalone export the two members of each current
# are emission from the opposite ends of the same charged fermion line.
# Check the implementation itself before selecting a diagram pair.
CURRENT_DIAGRAMS = {"z_isr":(5,7), "gamma_isr":(4,6), "isr":(4,5,6,7), "fsr":(0,1,2,3)}
CURRENT_DIGESTS = {True:"930d1bc3d4b5ef23d4c8e099c0569b8c390fa028f02d1edb23882e6774930163",
                   False:"1a7ab5ac470cad05de719a84d64e7eeb8114d6b14369d0b4515499ffe2eafcf8"}


def evaluate(work, model, up, rows, parameters, current=None):
    # MG5 orders q, qbar, e+, e-, gamma in both independent exports.
    process, = (work/model/"SubProcesses").glob("P1*" if up else "P2*")
    lib, = (work/model/"lib").glob("libmodel_*.a")
    source = process/"CPPProcess.cc"
    if current is not None:
        original = source.read_text()
        definition = original.split("// Calculate all wavefunctions")[1].split("// Sum and square the color flows")[0]
        signature = hashlib.sha256(" ".join(definition.split()).encode()).hexdigest()
        if signature != CURRENT_DIGESTS[up] or model != "full":
            raise ValueError("Unrecognized generated helicity amplitudes")
        selected = " + ".join(f"amp[{i}]" for i in CURRENT_DIAGRAMS[current])
        modified, count = re.subn(r"jamp\[0\] = \+amp\[0\][^;]+;",f"jamp[0] = {selected};",original)
        if count != 1:
            raise ValueError("Unrecognized colour flow")
        source = process/(current+".cc")
        source.write_text(modified)
    exe = process/("evaluate_"+(current or "all"))
    subprocess.run(["c++","-std=c++11","-O2","-I",str(process),"-I",str(work/model/"src"),
                    Path(__file__).with_name("hard_evaluate.cc"),source,lib,"-o",exe],check=True)
    data = "".join(" ".join(format(x,".17g") for p in row["p4"] for x in p)+"\n" for row in rows)
    result = subprocess.run([exe,parameters],input=data,text=True,capture_output=True,check=True)
    values = [float(s.split()[1]) for s in result.stdout.splitlines() if s.startswith("M2 ")]
    if len(values) != len(rows) or any(not math.isfinite(x) or x < 0. for x in values):
        raise ValueError("Invalid MG5 matrix-element evaluation")
    return np.array(values)


def project(inputs, mg5, parameters, work, pdf_events=None):
    work.mkdir(parents=True,exist_ok=True)
    matrix_sources(mg5,work)
    result = {"regions":{}, "combined":{}}
    opener = gzip.open if inputs.suffix == ".gz" else open
    with opener(inputs, "rt") as stream:
        data = json.load(stream)
    weightsets, currents = {}, {}
    for up in (True, False):
        rows = [e for region in data.values() for e in region["events"] if (e["hard_flavour"] % 2 == 0) == up]
        full = evaluate(work,"full",up,rows,parameters)
        pole = evaluate(work,"production_z",up,rows,parameters)
        if min(full) <= 0.:
            raise ValueError("Zero full matrix element on an event with positive generation probability")
        components = {current:evaluate(work,"full",up,rows,parameters,current)/full for current in CURRENT_DIAGRAMS}
        # Independent, separately generated q qbar -> Z gamma, Z -> ee amplitude
        # must agree with both production-Z diagrams, including interference.
        np.testing.assert_allclose(components["z_isr"],pole/full,rtol=3.e-7,atol=1.e-10)
        for i,row in enumerate(rows):
            currents[row["integration"],row["event"]] = {k:col[i] for k,col in components.items()}
        for row,ratio in zip(rows,pole/full):
            weightsets[(row["integration"],row["event"])] = float(ratio)

    cases = ["full", "production_z", "production_z_pt10", "production_z_pt10_bw15", "full_pt10_bw15"]
    pdf = None
    if pdf_events is not None:
        pdf = json.loads(pdf_events.read_text())
        result["pdf_weights_sha256"] = hashlib.sha256(pdf_events.read_bytes()).hexdigest()
    elif all("pdf_weights_sha256" in region for region in data.values()):
        pdf = {name:region["events"] for name,region in data.items()}
        pdf_sha, = {region["pdf_weights_sha256"] for region in data.values()}
        result["pdf_weights_sha256"] = pdf_sha
    if pdf is not None:
        cases += ["full_ehlq", "production_z_pt10_ehlq", "production_z_pt10_bw15_ehlq",
                  "production_z_pt10_bw15_ehlq_plus20"]
    pieces = ("isr_z", "isr_virtual_photon", "isr_z_photon_interference", "all_isr",
              "all_fsr", "isr_fsr_interference")
    cases.extend(pieces)
    combined, covariance = np.zeros(len(cases)), np.zeros((len(cases),len(cases)))
    for name, region in data.items():
        if region["mapped"] != region["attempted"]:
            raise ValueError("Use full-event source mapping with hard_events.py --showers")
        events = region["events"]
        weights = np.array([row["contribution_fb"] for row in events])
        ratios = np.array([weightsets[row["integration"],row["event"]] for row in events])
        p4 = np.array([row["p4"] for row in events])
        pair = p4[:,2]+p4[:,3]
        mll = np.sqrt(np.maximum(0.,pair[:,0]**2-(pair[:,1:]**2).sum(axis=1)))
        ptcut = np.hypot(p4[:,4,1],p4[:,4,2]) > 10.
        bwcut = np.abs(mll-91.188) < 15.*2.441404
        factors = [np.ones(len(weights)), ratios, ratios*ptcut, ratios*ptcut*bwcut, ptcut*bwcut]
        if pdf is not None:
            keyed = {(r["integration"],r["event"]):r for r in pdf[name]}
            matches = [keyed.pop((r["integration"],r["event"])) for r in events]
            if keyed:
                raise ValueError("Unmatched PDF weights")
            np.testing.assert_allclose([r["contribution_fb"] for r in matches], weights,rtol=1.e-12,atol=0.)
            shifts = np.asarray([r["pdf_ratio"] for r in matches])
            factors.extend((shifts,shifts*factors[2],shifts*factors[3],1.2*shifts*factors[3]))
        amplitudes = [currents[r["integration"],r["event"]] for r in events]
        z, virtual, isr, fsr = (np.array([a[k] for a in amplitudes]) for k in ("z_isr","gamma_isr","isr","fsr"))
        factors.extend((z, virtual, isr-z-virtual, isr, fsr, 1.-isr-fsr))
        columns = weights[:,None]*np.column_stack(factors)
        groups = region["integrations"]
        ids = np.array([r["integration"] for r in events])
        cov = np.zeros_like(covariance)
        for g in groups:
            n = g["generated"]
            selected = columns[ids==g["integration_id"]]
            sums = selected.sum(axis=0)
            # Other generated collisions have zero selected rate. Hard integration
            # errors are fully correlated across the reweightings of one integration.
            cov += n/(n-1.)*(selected.T@selected-np.outer(sums,sums)/n)
            cov += np.outer(sums,sums)*(g["cross_section_error_fb"]/g["cross_section_fb"])**2
        total = columns.sum(axis=0)
        np.testing.assert_allclose(total[0],region["selected_fb"])
        combined += total
        covariance += cov
        coeff = np.asarray([g["generated"]/region["attempted"] for g in groups])
        sigma = sum(w*g["cross_section_fb"] for w,g in zip(coeff,groups))/1000
        se = math.sqrt(sum((w*g["cross_section_error_fb"])**2 for w,g in zip(coeff,groups)))/1000
        result["regions"][name] = dict(frozen_sha256=region["frozen_sha256"],
            attempted=region["attempted"], mapped=region["mapped"], hard_cross_section_pb=sigma,
            hard_integration_se_pb=se, frozen_window_fb=region["selected_fb"],
            selected_efficiency_including_fake=region["selected_fb"]/(1000*sigma),
            mapped_window_collisions=len(events), omitted_window_collisions=region["unmapped_window_collisions"],
            variations=measurements(cases,total,cov))
    result["combined"] = measurements(cases,combined,covariance)
    result["observable_order"] = cases
    result["covariance_fb2"] = covariance.tolist()
    result["mg5_version"]=(mg5.parent.parent/"VERSION").read_text().strip()
    result["parameters_sha256"]=hashlib.sha256(parameters.read_bytes()).hexdigest()
    result["event_records_sha256"]=hashlib.sha256(inputs.read_bytes()).hexdigest()
    return result, weightsets


def measurements(cases, total, covariance):
    result = {}
    for i,key in enumerate(cases):
        ratio = total[i]/total[0]
        gradient = np.eye(len(cases))[i]/total[0] - ratio*np.eye(len(cases))[0]/total[0]
        result[key] = dict(window_fb=float(total[i]),mc_se_fb=math.sqrt(covariance[i,i]),
                           relative_to_full=float(ratio),paired_ratio_mc_se=math.sqrt(max(0.,gradient@covariance@gradient)))
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("events","mg5","param_card","work","output"):
        parser.add_argument(key,type=Path)
    parser.add_argument("--event-weights",type=Path,help="Optional per-event matrix-element ratios for paired variations")
    parser.add_argument("--pdf-weights",type=Path,help="Exact per-event EHLQ/NNPDF ratios from pdf_reweight.py")
    a=parser.parse_args()
    result, ratios=project(a.events,a.mg5.resolve(),a.param_card.resolve(),a.work.resolve(), a.pdf_weights)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if a.event_weights:
        a.event_weights.write_text(json.dumps([dict(integration=k[0],event=k[1],production_z_over_full=r)
                                              for k,r in sorted(ratios.items())],indent=2)+"\n")
