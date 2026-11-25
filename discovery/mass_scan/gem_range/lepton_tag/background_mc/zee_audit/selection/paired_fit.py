"""Combine detector probabilities for mutually exclusive assignments on the same underlying eeγ collision."""
from collections import Counter
from pathlib import Path
import gzip
import hashlib
import json

import numpy as np


def read(path):
    path=Path(path)
    with (gzip.open if path.suffix==".gz" else open)(path,"rt") as f:
        return json.load(f)


def load(features_path, weights_path):
    """Load candidate rows in herwig_candidates.json.gz order, verified against the compact hard weights."""
    features=read(features_path)
    hard=read(weights_path)
    if hashlib.sha256(Path(features_path).read_bytes()).hexdigest()!=hard["features_sha256"]:
        raise ValueError("Candidate file differs from the paired hard input")
    events=hard["events"]
    lookup={(r["shard"],r["event"]):i for i,r in enumerate(events)}
    if len(lookup)!=len(events):
        raise ValueError("Two independent hard weights point to the same shower collision")
    if len({(r["integration"],r["hard_event"]) for r in events})!=len(events):
        raise ValueError("Underlying LHE event appears in multiple shower blocks")
    assignment_events=np.array([lookup[(r["shard"],r["event"])] for r in features])
    if len(set(assignment_events))!=len(events):
        raise ValueError("Underlying collision has no candidate assignments")
    return dict(assignments=features,events=events,groups=hard["groups"],assignment_events=assignment_events)


def error(variance):
    return float(np.sqrt(max(0.,variance)))


def combine(event_prob,packet,historical):
    events=packet["events"]
    # All generator integrations in one region contribute to its common attempted-event denominator.
    # The core and complementary regions describe disjoint hard phase space, so their cross sections add.
    ids=np.array([r["integration"] for r in events])
    hard=np.ones(len(events))
    if historical:
        hard=np.array([r["production_ratio"]*r["pdf_ratio"]*(r["hard_photon_pt_gev"]>10.)*1.20*.99 for r in events])
    nvar=event_prob.shape[1]
    total=np.zeros(nvar)
    cov=np.zeros((nvar,nvar))
    for key,g in packet["groups"].items():
        indices=np.flatnonzero(ids==key)
        n=g["generated"]
        norm=g["cross_section_fb"]/g["region_generated"]
        rows=event_prob[indices]*(hard[indices]*norm)[:,None]
        sums=rows.sum(axis=0)
        outer=np.outer(sums,sums)
        # The n-k unselected generated events have zero weights. The cross-section integration
        # uncertainty is common to all detector alternatives for a given integration.
        total+=sums
        cov+=n/(n-1.)*(rows.T@rows-outer/n) + outer*(g["cross_section_error_fb"]/g["cross_section_fb"])**2
    return total,cov


def evaluate(packet, variants, baseline="nominal", reference_fb=.0374):
    """Use a name -> array of probabilities in packet['assignments'] order. Sum assignments within each collision first.

    Output 'results' rows give rates, paired differences and ratios, with Monte Carlo errors. The fraction of the
    residual removed is (baseline - variation)/(baseline - historical table). The printed GEM rate is held fixed.
    'covariance_fb2' describes all the historical-input rates in 'observable_order'.
    """
    names=list(variants)
    original=names.index(baseline)
    values=np.column_stack([variants[s] for s in names])
    if values.shape!=(len(packet["assignments"]),len(names)) or not np.all(np.isfinite(values)) or np.min(values)<0.:
        raise ValueError("One finite positive detector probability required for each candidate and variation")
    assignments=packet["assignment_events"]
    event_prob=np.zeros((len(packet["events"]),len(names)))
    np.add.at(event_prob,assignments,values)
    if np.max(event_prob)>1.+1.e-10:
        raise ValueError("Exclusive detector outcomes for a collision exceed unit probability")
    modern,modern_cov=combine(event_prob,packet,False)
    total,cov=combine(event_prob,packet,True)
    base=total[original]
    excess=base-reference_fb
    if excess<=0:
        raise ValueError("Nominal model has no positive residual")
    results={}
    eye=np.eye(len(names))
    for i,name in enumerate(names):
        v=total[i]
        diff=eye[i]-eye[original]
        ratio=v/base
        ratio_gradient=(eye[i]-ratio*eye[original])/base
        residual_gradient=(v-reference_fb)/excess**2*eye[original]-eye[i]/excess
        results[name]=dict(rate_fb=float(v),mc_se_fb=error(cov[i,i]),
            change_from_baseline_fb=float(v-base),change_paired_se_fb=error(diff@cov@diff),
            relative_to_baseline=float(ratio),ratio_paired_se=error(ratio_gradient@cov@ratio_gradient),
            ratio_to_gem=float(v/reference_fb),ratio_to_gem_mc_se=error(cov[i,i])/reference_fb,
            fraction_of_residual_removed=float((base-v)/excess),residual_fraction_paired_se=error(residual_gradient@cov@residual_gradient),
            modern_hard_rate_fb=float(modern[i]),modern_hard_mc_se_fb=error(modern_cov[i,i]))
    return dict(mass_window_gev=[123.3,126.7],observable_order=names,baseline=baseline,gem_fb=reference_fb,
        generated=sum(g["generated"] for g in packet["groups"].values()),collision_count=len(event_prob),
        candidate_count=len(values),multiplicity=dict(Counter(Counter(assignments.tolist()).values())),
        results=results,covariance_fb2=cov.tolist(),modern_covariance_fb2=modern_cov.tolist())
