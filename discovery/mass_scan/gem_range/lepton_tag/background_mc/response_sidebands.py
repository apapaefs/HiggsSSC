"""Blinded smooth-continuum cross-check of the independent MC window rate at 125 GeV."""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares
from scipy.stats import chi2

WINDOW = (123.3, 126.7)
NODES, WEIGHTS = leggauss(12)


def fractions(edges, lo, hi):
    return np.maximum(0., np.minimum(hi, edges[1:])-np.maximum(lo, edges[:-1]))/np.diff(edges)


def integrate(intervals, pars):
    lo, hi = np.asarray(intervals).T
    half = (hi-lo)/2
    m = (lo+hi)[:, None]/2 + half[:,None]*NODES
    basis = np.stack([((m-125.)/20.)**j for j in range(len(pars))], axis=-1)
    y = np.exp(np.einsum("ijk,k->ij",basis,pars))
    yield_fb = half*(y @ WEIGHTS)
    jac = np.einsum("i,ij,j,ijk->ik", half, y, WEIGHTS, basis)
    return yield_fb, jac


def fit(edges, rate, covariance, endpoints, order):
    intervals = [(a,b) for a,b in zip(endpoints,endpoints[1:]) if b<=WINDOW[0] or a>=WINDOW[1]]
    A = np.stack([fractions(edges,a,b) for a,b in intervals])
    y = A @ rate
    V = A @ covariance @ A.T
    if np.any(y<=0.) or np.any(np.diag(V)<=0.):
        return dict(usable=False, reason="One or more sidebands have no simulated rate or variance.")
    try:
        L = np.linalg.cholesky(V)
    except np.linalg.LinAlgError:
        return dict(usable=False, reason="The sideband covariance matrix is singular.")
    W = np.linalg.inv(L)
    def fun(p): return W @ (integrate(intervals,p)[0]-y)
    def jac(p): return W @ integrate(intervals,p)[1]
    opt = least_squares(fun,[np.log(y.sum()/sum(b-a for a,b in intervals))]+[0.]*order,jac=jac,
                        bounds=(-15.,15.))
    if not opt.success or np.any(opt.active_mask) or np.linalg.cond(opt.jac)>1.e7:
        return dict(usable=False, reason="No stable interior solution for this model.")
    dof = len(intervals)-len(opt.x)
    stat = float(opt.fun @ opt.fun)
    pvalue = float(chi2.sf(stat,dof))
    prediction, derivative = integrate([WINDOW],opt.x)
    # Linearized GLS propagation includes correlations between sideband MC and
    # any independent-event packets also contributing to the blinded window.
    gradient = (derivative @ np.linalg.pinv(opt.jac) @ W @ A)[0]
    direct = fractions(edges,*WINDOW)
    r = float(prediction[0])
    e = np.sqrt(max(0.,float(gradient@covariance@gradient)))
    diff_err = np.sqrt(max(0.,float((direct-gradient)@covariance@(direct-gradient))))
    observations = (y**2/np.diag(V)).tolist()
    return dict(usable=True, intervals_gev=intervals, rate_fb=r, error_fb=e, parameters=opt.x.tolist(),
                statistic=stat,dof=dof,asymptotic_chi_square_p=pvalue,
                sideband_effective_events=observations,min_effective_events=min(observations),
                direct_window_closure_pull=float((direct@rate-r)/diff_err) if diff_err else None)


def calculate(report):
    edges=np.array(report["edges_gev"])
    rate=np.array(report["rates_fb"])
    covariance=np.array(report["covariance_fb2"])
    fits={}
    for floor in (100.,110.):
        knots = [x for x in [100.,105.,110.,115.,120.,*WINDOW,130.,135.,140.,145.] if x>=floor]
        for order in (1,2,3):
            fits[f"{int(floor)}_145_order_{order}"] = fit(edges,rate,covariance,knots,order)
    sparse={name:source["window_125"] for name,source in report["sources"].items()
            if (source["window_125"]["effective_events"] or 0.)<5.}
    credible=[v for v in fits.values() if v.get("usable") and v["min_effective_events"]>=5.
              and v["asymptotic_chi_square_p"]>.05]
    values=[v["rate_fb"] for v in credible]
    # The envelope reports dependence on degree and range. Sparse large-weight
    # processes and unvalidated decay strata still require their own simulation.
    return dict(window_gev=WINDOW, nominal_direct_window=report["window_125"], fits=fits,
                fit_rate_envelope_fb=[min(values),max(values)] if values else None,
                envelope_criteria="Fitted sidebands each have at least five effective events and chi-square p > 0.05.",
                missing_sources=report["missing_sources"], zero_window_sources=report["zero_window_sources"],
                sparse_window_sources=sparse,top_decay_model=report["top_decay_model"],
                zero_window_components={name:s.get("zero_window_components",[]) for name,s in report["sources"].items()
                                        if s.get("zero_window_components")})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("spectrum",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    result=calculate(json.loads(a.spectrum.read_text()))
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("Sideband fit range and basis envelope [fb]",result["fit_rate_envelope_fb"])
