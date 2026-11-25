"""Build physical 125 GeV mass templates and source-event uncertainties from the tagged showers."""
from pathlib import Path
from argparse import ArgumentParser
import hashlib
import json

import numpy as np
from scipy.ndimage import gaussian_filter1d

from transfer import PRODUCTION_PB, BR_H_GAMMAGAMMA

EDGES = np.round(np.linspace(70.,170.,2001), 8)
SMOOTH_GEV = .15


def positive(raw):
    # Signed MC can give negative individual 50 MeV bins. The signal resolution
    # exceeds the kernel width; restore each category's signed integral.
    smeared = gaussian_filter1d(raw, SMOOTH_GEV/.05, mode="reflect")
    smeared = np.clip(smeared, 0, None)
    if smeared.sum() == 0: return smeared
    return smeared * (raw.sum()/smeared.sum())


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("events",type=Path)
    args=parser.parse_args()
    with np.load(args.events,allow_pickle=False) as f:
        data=dict(zip(f["columns"],f["values"].T))
        meta=json.loads(str(f["metadata"]))
    events=data["evweight"]
    generated=sum(s["generated"] for s in meta["shards"])
    if len(meta["inputs"])!=9 or len(meta["shards"])!=45 or generated!=45000:
        raise ValueError("Expected all 45,000 independent hard collisions")
    raw={category: np.zeros(len(EDGES)-1) for category in ("tagged","untagged_associated","overlap_inclusive_mass")}
    mode_data={}
    for mode,pb in PRODUCTION_PB.items():
        shards=[i for i,s in enumerate(meta["shards"]) if s["parent"].startswith("signal_"+mode+"_")]
        keep=np.isin(data["shard"],shards)
        count=sum(meta["shards"][i]["generated"] for i in shards)
        norm=sum(meta["shards"][i]["signed_attempted"] for i in shards)
        weight=events[keep]
        if count!=15000 or norm<=0 or not np.all((weight==1)|(weight==-1)):
            raise ValueError(f"Unexpected source normalization for {mode}")
        factor=1000*pb*BR_H_GAMMAGAMMA/norm
        pi,pt,po=[data[k][keep] for k in ("p_inclusive","p_tag","p_overlap")]
        if np.any(po<0) or np.any(po>np.minimum(pi,pt)+1.e-9):
            raise ValueError("Invalid joint identification probability")
        inc,tag=[data[k][keep] for k in ("mass_inclusive","mass_tag")]
        moments={}
        def add(name,prob):
            mean=np.sum(weight*prob)/norm
            # Signed-weight ratio variance, including any unsaved shower
            # attempts at zero efficiency; generated weights have absolute one.
            var=np.sum(weight**2*(prob-mean)**2)+(count-len(weight))*mean**2
            moments[name+"_fb"]=float(np.sum(weight*prob)*factor)
            moments[name+"_mc_se_fb"]=float(np.sqrt(var)*factor)
        for name,prob in (("inclusive",pi),("tagged",pt),("overlap",po),("untagged_associated",pi-po),
                          ("tagged_pm2sigma",pt*(np.abs(tag-125)<=1.7))):
            add(name,prob)
        for origin in ("direct","tau","hf","other"):
            add("tag_origin_"+origin,data["p_"+origin][keep])
        for name,mass,prob in (("tagged",tag,pt),("untagged_associated",inc,pi-po),
                               ("overlap_inclusive_mass",inc,po)):
            raw[name]+=np.histogram(mass,bins=EDGES,weights=weight*prob*factor)[0]
        mode_data[mode]=dict(production_pb=pb,generated=count,saved=len(weight),negative_saved=int((weight<0).sum()),
                             signed_generated=norm,exactly_two_higgs_photons_above_2gev=int(np.sum(data["n_higgs_photons"][keep]==2)),**moments)
    totals={}
    for k in moments:
        values=np.array([s[k] for s in mode_data.values()])
        totals[k]=float(np.linalg.norm(values) if k.endswith("_mc_se_fb") else values.sum())
    physical={k: positive(v) for k,v in raw.items()}
    physical["inclusive"]=physical["untagged_associated"]+physical["overlap_inclusive_mass"]
    for k in physical:
        if k in raw and abs(physical[k].sum()-raw[k].sum())>1e-12:
            raise ValueError("Normalization shift in template "+k)
    audit={k:dict(negative_bins=int((v < -1.e-12).sum()),negative_weight_fb=float(np.minimum(v,0.).sum()),
                  signed_integral_in_edges_fb=float(v.sum())) for k,v in raw.items()}
    summary=dict(selected_fb=totals["tagged_fb"], selected_mc_se_fb=totals["tagged_mc_se_fb"],
                 before_nonnegative_smoothing=totals,by_mode=mode_data,
                 event_replica_sha256=hashlib.sha256(args.events.read_bytes()).hexdigest(),
                 generated=generated,saved=len(events),mass_edges_gev=EDGES.tolist(),
                 smoothing_kernel_sigma_gev=SMOOTH_GEV,template_audit=audit)
    for k,v in physical.items(): summary[k+"_spectrum_fb_per_bin"]=v.tolist()
    print(json.dumps(summary,allow_nan=False, separators=(",",":")))


if __name__=="__main__":
    main()
