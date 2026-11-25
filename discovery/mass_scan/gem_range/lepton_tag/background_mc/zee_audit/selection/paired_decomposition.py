"""Attribute the paired isolation-rate change to each EM object by averaging all six orders."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from noise import pass_probability
from paired_fit import combine, error, load
from selection_rates import columns, model, RADII, ROLES


def decompose(packet, reference=.0374):
    rows=packet["assignments"]
    base_prob=np.array([r["probability"] for r in rows])
    before=columns(rows,"cone") < 5.
    after=pass_probability(columns(rows,"cell_iso_et"), RADII)
    variants,_=model(rows)
    np.testing.assert_allclose(before.prod(axis=1)*base_prob, variants["nominal"],atol=1.e-13)
    np.testing.assert_allclose(after.prod(axis=1)*base_prob, variants["noise"],atol=1.e-13)
    # Every path starts with nominal isolation for all three candidates and ends with
    # the updated marginal noise probabilities. Each order has the same 1/6 weight.
    changes=np.zeros_like(after)
    orders=list(itertools.permutations(range(len(ROLES))))
    for order in orders:
        roles=np.zeros(len(ROLES),dtype=bool)
        previous=variants["nominal"]
        for j in order:
            roles[j]=True
            current=base_prob*np.where(roles,after,before).prod(axis=1)
            changes[:,j]+= (current-previous)/len(orders)
            previous=current
    np.testing.assert_allclose(changes.sum(axis=1),variants["noise"]-variants["nominal"],atol=1.e-14)
    columns_out=np.column_stack((variants["nominal"], variants["noise"], changes))
    event=np.zeros((len(packet["events"]),columns_out.shape[1]))
    np.add.at(event,packet["assignment_events"],columns_out)
    rates,cov=combine(event,packet,True)
    np.testing.assert_allclose(rates[2:].sum(),rates[1]-rates[0],atol=1.e-14)
    diff=np.array([-1.,1.,0.,0.,0.])
    excess=rates[0]-reference
    results={}
    for i,role in enumerate(ROLES,2):
        gradient=np.zeros(len(rates))
        gradient[0]=rates[i]/excess**2
        gradient[i]=-1./excess
        results[role]=dict(change_fb=rates[i],paired_mc_se_fb=error(cov[i,i]),
                           fraction_of_original_residual_removed=-rates[i]/excess,
                           residual_fraction_paired_se=error(gradient@cov@gradient))
    return dict(nominal_fb=rates[0],with_cell_isolation_and_noise_fb=rates[1],
                combined_change_fb=diff@rates,combined_change_paired_se_fb=error(diff@cov@diff),
                roles=results,reference_fb=reference,covariance_order=["nominal","noise",*ROLES],
                covariance_fb2=cov.tolist(),assumption="Isolations of the three candidates have independent noise given their individual in-cone activities. The attribution averages all six orders of substituting the three isolation decisions.")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("features","hard","output"):
        parser.add_argument(key,type=Path)
    args=parser.parse_args()
    result=decompose(load(args.features,args.hard))
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
