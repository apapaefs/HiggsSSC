"""Binned diphoton likelihood and independent sparse four-lepton counting experiment."""
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import brentq, minimize, root
from scipy.special import gammaincinv, ndtr
from scipy.stats import norm, poisson

EDGES = np.arange(105., 161.)
MASS = (EDGES[:-1] + EDGES[1:]) / 2
BLIND = (MASS >= 120) & (MASS < 130)
ALL = np.ones_like(BLIND)
NODES, WEIGHTS = np.polynomial.legendre.leggauss(6)
MODELS = (("exponential", 2), ("power", 2), ("exponential", 3))
ALPHA = norm.sf(5.)


class Shape:
    """Exponentiated polynomial, integrated over each 1 GeV bin by quadrature."""

    def __init__(self, family, degree):
        points = MASS[:, None] + NODES / 2
        x = (points - 125) / 25 if family == "exponential" else 5 * np.log(points / 125)
        self.basis = np.stack([x**i for i in range(degree + 1)], axis=-1)
        self.degree = degree

    def get(self, theta, hessian=False):
        terms = np.exp(self.basis @ theta) * WEIGHTS / 2
        means = terms.sum(axis=1)
        jac = np.einsum("nq,nqd->nd", terms, self.basis)
        if hessian:
            return means, jac, np.einsum("nq,nqi,nqj->nij", terms, self.basis, self.basis)
        return means, jac

    def fit(self, data, mask, signal=None):
        """Return parameters, continuum, continuum Jacobian and response to the input bins."""
        npars = self.degree + 1
        start = [np.log(np.median(data[mask])), -.65] + [0.] * (npars - 2)
        if signal is not None:
            start.append(0.)

        def likelihood(pars, curvature=False):
            continuum, jac, second = self.get(pars[:npars], True)
            mean = continuum.copy()
            if signal is not None:
                mean += pars[-1] * signal
                jac = np.column_stack((jac, signal))
            pred, obs, deriv = mean[mask], data[mask], jac[mask]
            if np.any(pred <= 0) or not np.isfinite(pred).all():
                if curvature:
                    raise ValueError("Nonpositive continuum fit")
                return 1.e32, np.full_like(pars, 1.e12)
            fraction = 1 - obs / pred
            gradient = deriv.T @ fraction
            if not curvature:
                return np.sum(pred - obs + obs * np.log(obs / pred)), gradient
            hess = deriv.T @ ((obs / pred**2)[:, None] * deriv)
            hess[:npars, :npars] += np.einsum("n,nij->ij", fraction, second[mask])
            return gradient, hess

        fit = minimize(likelihood, start, jac=True, method="BFGS", options={"gtol": 1.e-8})
        best = fit.x
        if np.linalg.norm(fit.jac, np.inf) > 1.e-4:
            best = root(likelihood, best, args=(True,), jac=True, tol=1.e-10).x
        grad, hess = likelihood(best, True)
        if np.linalg.norm(grad, np.inf) > 1.e-4 or np.linalg.eigvalsh(hess)[0] <= 0:
            raise ValueError("Continuum likelihood minimization failed")
        curve, jac = self.get(best[:npars])
        mean = curve if signal is None else curve + best[-1] * signal
        design = jac if signal is None else np.column_stack((jac, signal))
        response = np.zeros((len(best), len(mask)))
        response[:, mask] = np.linalg.solve(hess, (design[mask] / mean[mask, None]).T)
        return best, curve, jac, response


class Calibration:
    """Gaussian auxiliary measurement of a component with the Higgs mass shape."""

    def __init__(self, shape, theta, continuum, signal, center, width):
        self.shape, self.theta = shape, theta
        self.signal, self.center, self.width = signal, center, width
        self.truth = continuum + (1 + center) * signal
        if min(self.truth) <= 0 or min(continuum + center * signal) <= 0 or width <= 0:
            raise ValueError("Invalid background bias model")

    @lru_cache(maxsize=128)
    def q(self, luminosity):
        penalty = 1 / (self.width**2 * luminosity)
        def score(pars):
            mean, deriv, second = self.shape.get(pars[:-1], True)
            mean += (self.center + pars[-1]) * self.signal
            if min(mean) <= 0:
                raise ValueError("Nonpositive calibrated Poisson mean")
            deriv = np.column_stack((deriv, self.signal))
            fraction = (mean - self.truth) / mean
            grad = deriv.T @ fraction
            grad[-1] += pars[-1] * penalty
            hess = deriv.T @ ((self.truth / mean**2)[:, None] * deriv)
            hess[:-1, :-1] += np.einsum("n,nij->ij", fraction, second)
            hess[-1, -1] += penalty
            return grad, hess

        start = np.r_[self.theta, 1 / (1 + penalty / 8)]
        best = root(score, start, jac=True, tol=1.e-11).x
        grad, hess = score(best)
        if np.max(np.abs(grad)) > 1.e-6 or np.linalg.eigvalsh(hess)[0] <= 0:
            raise ValueError("Bias calibration likelihood minimization failed")
        null = self.shape.get(best[:-1])[0] + (self.center + best[-1]) * self.signal
        delta = (null - self.truth) / self.truth
        divergence = delta - np.log1p(delta)
        small = np.abs(delta) < 1.e-4
        t = delta[small]
        divergence[small] = t**2 * (.5 + t * (-1 / 3 + t * (.25 - .2 * t)))
        return float(2 * luminosity * (self.truth @ divergence) + (best[-1] / self.width)**2)


def crossing(difference):
    upper = 1.
    while difference(upper) < 0:
        upper *= 2
        if upper > 1.e7:
            raise ValueError("No discovery crossing in the supported range")
    return float(brentq(difference, 1.e-8, upper, xtol=2.e-9))


def fourlepton_median(signal, background):
    """First non-randomized 50% crossing, including the jumps in the critical count."""
    for count in range(1, 10000):
        luminosity = gammaincinv(count, .5) / (signal + background)
        tail = poisson.sf(count - 1, background * luminosity)
        if tail <= ALPHA:
            return float(luminosity), count
    raise ValueError("Four-lepton median exceeds count search range")


def joint_power(qgamma, luminosity, signal, background):
    """Threshold the joint likelihood ratio under H0, then calculate its H1 power."""
    q = qgamma(luminosity)
    sd = math.sqrt(q)
    counts = np.arange(2 + int(poisson.isf(1.e-14, (signal + background) * luminosity)))
    score = counts * math.log1p(signal / background)
    p0 = poisson.pmf(counts, background * luminosity)
    p1 = poisson.pmf(counts, (background + signal) * luminosity)
    threshold = brentq(lambda c: p0 @ ndtr((score - c) / sd) - ALPHA,
                       -10 * sd, 10 * sd + score[-1])
    return float(p1 @ ndtr((score + q - threshold) / sd))
