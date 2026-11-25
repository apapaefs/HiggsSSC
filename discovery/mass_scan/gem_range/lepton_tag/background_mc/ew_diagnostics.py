"""Inspect flavour luminosities and generation-cut margins in accepted EW collisions.

The b and CKM variations change the incoming hard luminosities. The shower and
detector response remain those of the sampled light-quark collisions.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path


def mass(rows):
    px,py,pz,e=[sum(float(row[i]) for row in rows) for i in (6,7,8,9)]
    return math.sqrt(max(0.,e*e-px*px-py*py-pz*pz))


def geometry(p):
    px,py,pz=map(float,p[6:9])
    pt=math.hypot(px,py)
    return pt,math.asinh(pz/max(pt,1.e-100)), math.atan2(py,px)


def dr(a,b):
    return math.hypot(a[1]-b[1],math.remainder(a[2]-b[2],2.*math.pi))


def in_nominal(sample, rows):
    final=[r for r in rows if r[1]=="1"]
    photons=[r for r in final if r[0]=="22"]
    leptons=[r for r in final if abs(int(r[0])) in (11,13,15)]
    # The generated run-card eta is evaluated by MG5's lab rapidity function.
    def ygeo(row):
        pt,_,phi=geometry(row)
        e,pz=float(row[9]),float(row[8])
        return pt,math.atanh(pz/e) if abs(pz)<e else math.copysign(math.inf,pz),phi
    gp=list(map(ygeo,photons))
    lp=list(map(ygeo,leptons))
    if any(p[0]<14. or abs(p[1])>3. for p in gp):
        return False
    if sample=="z_ee_a":
        return all(p[0]>=14. and abs(p[1])<=3. and dr(p,gp[0])>=.2 for p in lp)
    if not 60.<=mass(photons)<=170. or max(p[0] for p in lp)<13.:
        return False
    return sample!="w_lnu_aa" or abs(lp[0][1])<=3.5


def verify_outside(lhe,sample):
    n=0
    with gzip.open(lhe,"rt") as stream:
        for line in stream:
            if line.strip()=="<event>":
                header=next(stream).split()
                rows=[next(stream).split() for _ in range(int(header[0]))]
                if in_nominal(sample,rows):
                    raise ValueError(f"complement event {n} passes nominal generator acceptance")
                n+=1
    return n


def incoming_shift(pdf, rows, scale, charged):
    first,second=rows
    def partons(row,flavours):
        x=float(row[9])/20000.
        sign=1 if int(row[0])>0 else -1
        return [pdf.xfxQ(sign*pid,x,scale) for pid in flavours]
    if charged:
        if abs(int(first[0])) not in (2,4):
            first,second=second,first
        u,c=partons(first,[2,4]);d,s,b=partons(second,[1,3,5])
        l=.2253
        A=.808
        rho,eta=.132,.341
        diag=(1.-l*l/2.)**2
        denominator=u*d+c*s
        light=(diag*(u*d+c*s)+l*l*(u*s+c*d))/denominator
        bottom=(A*A*l**6*(rho*rho+eta*eta)*u*b+A*A*l**4*c*b)/denominator
        return light-1.,bottom
    if abs(int(first[0])) not in (1,3):
        return 0.,0.
    a=partons(first,[1,3,5]);b=partons(second,[1,3,5])
    return 0.,a[2]*b[2]/(a[0]*b[0]+a[1]*b[1])


def measures(sample, rows, scale, pdf):
    inc=[r for r in rows if r[1]=="-1"]
    hard=[r for r in rows if r[1]=="1"]
    photons=[r for r in hard if r[0]=="22"]
    leptons=[r for r in hard if abs(int(r[0])) in (11,13,15)]
    gp=[geometry(row) for row in photons]
    lp=[geometry(row) for row in leptons]
    ckm,bottom=incoming_shift(pdf,inc,scale,sample=="w_lnu_aa")
    result={"all":1., "ckm_light_shift":ckm,"bottom_shift":bottom,
            "gen_photon_pt_below_16":float(min(g[0] for g in gp)<16.),
            "gen_photon_abseta_above_2p8":float(max(abs(g[1]) for g in gp)>2.8),
            "hard_tau":float(any(abs(int(r[0]))==15 for r in leptons))}
    if len(gp)==2:
        result["gen_pair_mass_outside_65_165"]=float(not 65.<mass(photons)<165.)
    if sample=="w_lnu_aa":
        result["gen_lepton_pt_below_15"]=float(lp[0][0]<15.)
        result["gen_lepton_abseta_above_3p3"]=float(abs(lp[0][1])>3.3)
    else:
        mll=mass(leptons)
        result.update(near_z_ll=float(abs(mll-91.188)<10.),
                      below_z_ll=float(mll<=81.188),above_z_ll=float(mll>=101.188),
                      near_z_all=float(abs(mass(hard)-91.188)<10.))
        if sample=="z_ll_aa":
            result["gen_leading_lepton_pt_below_15"]=float(max(p[0] for p in lp)<15.)
        else:
            result.update(gen_lepton_pt_below_16=float(min(p[0] for p in lp)<16.),
                          gen_lepton_abseta_above_2p8=float(max(abs(p[1]) for p in lp)>2.8),
                          gen_lepton_photon_dr_below_0p3=float(min(dr(p,gp[0]) for p in lp)<.3))
    result["any_generator_edge"]=float(any(v for k,v in result.items() if k.startswith("gen_")))
    return result


def calculate(ew):
    import lhapdf
    lhapdf.setVerbosity(0)
    pdf=lhapdf.mkPDF("NNPDF40_nlo_as_01180_qed",0)
    sums=defaultdict(Counter)
    counts=Counter()
    skipped=Counter()
    parents=defaultdict(list)
    for path in sorted((ew/"showers").glob("*/receipt.json")):
        rec=json.loads(path.read_text())
        parents[rec["integration_id"]].append((path,rec))
    for meta_path in sorted(ew.glob("*/hard.json")):
        meta=json.loads(meta_path.read_text())
        recs=parents.get(meta["integration_id"],[])
        sample=meta["sample"]
        window={}
        for path,rec in recs:
            if rec["saved_events"]!=rec["attempted_events"]:
                # Each saved index maps directly to its source only for lossless blocks.
                skipped[sample]+=rec["attempted_events"]
                continue
            counts[sample]+=rec["attempted_events"]
            with path.with_name("packet.csv").open() as stream:
                for row in csv.DictReader(stream):
                    if 123.3<=float(row["mass_gev"])<126.7:
                        idx=rec["source_first_event"]+int(row["source_event"])
                        window[idx]=window.get(idx,0.)+float(row["probability"])*float(row["evweight"])
        if not window:
            continue
        with gzip.open(meta_path.with_name("events.lhe.gz"),"rt") as stream:
            event=-1
            for line in stream:
                if line.strip()!="<event>":
                    continue
                event+=1
                header=next(stream).split()
                rows=[next(stream).split() for _ in range(int(header[0]))]
                if event in window:
                    factor=window.pop(event)*meta["cross_section_pb"]*1000.
                    for k,val in measures(sample,rows,float(header[3]),pdf).items():
                        sums[sample][k]+=factor*val
                    if not window:
                        break
        if window:
            raise ValueError("selected event absent from source")
    return {s:{"mapped_attempts":counts[s],"attempts_excluded_due_to_lost_showers":skipped[s],
               "selected_window_fb":dict((k,v/counts[s]) for k,v in vals.items()),
               "fraction_of_window":dict((k,v/vals["all"]) for k,v in vals.items())}
            for s,vals in sums.items()}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("ew",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    a.output.write_text(json.dumps(calculate(a.ew),indent=2)+"\n")
