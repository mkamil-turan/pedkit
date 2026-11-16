"""Liability-threshold risk by Markov chain Monte Carlo (Gibbs sampling).

In the infinitesimal model the liability vector of the individuals is
multivariate normal: L ~ N(0, Σ), Σ = h²·2Φ + (1 - h²)·I (Φ the kinship
matrix; MZ twins have 2φ = 1). Observed phenotypes constrain L: L > T for
the affected, L ≤ T for the unaffected. The posterior P(L | phenotypes) is
a truncated multivariate normal and is sampled by Gibbs sampling (data
augmentation of Albert and Chib 1993): at each step one individual's
liability is drawn from its univariate normal full conditional, truncated
to the interval its phenotype allows. The conditional is read from the
precision matrix Q = Σ⁻¹:
    L_i | L_-i ~ N(-Σ_{j≠i} Q_ij L_j / Q_ii, 1 / Q_ii).

The target's risk is a Rao–Blackwell estimate: P(L_target > T | L_-target)
= 1 - Φ((T - μ)/s) is averaged over the steps. Its variance is much lower
than that of counting whether the target's sampled liability exceeds the
threshold.

Several independent chains run at once (as numpy vectors); convergence is
measured by the Gelman–Rubin R̂, uncertainty by the spread of the chain
means.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..inference.prepare import Prepared
from .model import MFParameters, kinship_matrix
from .normal import truncated_normal, upper_tail


@dataclass
class MCMCResult:
    risk: float                     # Rao–Blackwell risk estimate
    standard_error: float           # from the chain means
    risk_count: float               # simple estimate: share of L_target > T
    rhat_liability: float           # Gelman–Rubin R̂ for the target's liability
    rhat_risk: float                # R̂ for the Rao–Blackwell risk traces
    chains: int
    samples: int                    # steps per chain after burn-in
    burn_in: int
    risk_trace: np.ndarray          # [step, chain] Rao–Blackwell values
    liability_trace: np.ndarray     # [step, chain] target liability
    mean_liability: np.ndarray      # [individual] posterior mean
    liability_interval: np.ndarray  # [individual, 2] 5 % and 95 % quantiles
    target_samples: np.ndarray      # every sample of the target's liability


def mcmc_risk(h: Prepared, par: MFParameters, target: int, chains: int = 64,
              samples: int = 4000, burn_in: int = 1000, seed: int | None = None,
              thin: int = 4) -> MCMCResult:
    rng = np.random.default_rng(seed)
    n = h.n
    h2 = par.heritability
    sigma = h2 * 2 * kinship_matrix(h) + (1 - h2) * np.eye(n)
    Q = np.linalg.inv(sigma)
    Qd = np.diag(Q).copy()
    s = 1 / np.sqrt(Qd)
    thresholds = np.array([par.individual_threshold(c) for c in h.sex])

    # Start: each individual from N(0,1) truncated to its own constraint, independently
    L = np.empty((chains, n))
    for i in range(n):
        L[:, i] = _draw(rng, np.zeros(chains), 1.0, h.phenotype[i], thresholds[i])

    target_thresholds = _target_thresholds(h, par, target)
    total = burn_in + samples
    risk_trace = np.empty((samples, chains))
    liab_trace = np.empty((samples, chains))
    keep = max(1, thin)
    stored = []
    for step in range(total):
        for i in range(n):
            mu = -(L @ Q[:, i] - Qd[i] * L[:, i]) / Qd[i]
            if i == target and step >= burn_in:
                rb = np.mean([upper_tail((t - mu) / s[i]) for t in target_thresholds], axis=0)
                risk_trace[step - burn_in] = rb
            L[:, i] = _draw(rng, mu, s[i], h.phenotype[i], thresholds[i])
        if step >= burn_in:
            liab_trace[step - burn_in] = L[:, target]
            if (step - burn_in) % keep == 0:
                stored.append(L.copy())

    draws = np.concatenate(stored, axis=0)          # [samples*chains/thin, n]
    chain_means = risk_trace.mean(axis=0)
    count = np.mean([(liab_trace > t).mean() for t in target_thresholds])
    return MCMCResult(
        risk=float(chain_means.mean()),
        standard_error=float(chain_means.std(ddof=1) / math.sqrt(chains)) if chains > 1
        else math.nan,
        risk_count=float(count),
        rhat_liability=gelman_rubin(liab_trace),
        rhat_risk=gelman_rubin(risk_trace),
        chains=chains, samples=samples, burn_in=burn_in,
        risk_trace=risk_trace, liability_trace=liab_trace,
        mean_liability=draws.mean(axis=0),
        liability_interval=np.quantile(draws, [0.05, 0.95], axis=0).T,
        target_samples=liab_trace.ravel(),
    )


def _target_thresholds(h: Prepared, par: MFParameters, target: int) -> list[float]:
    """For a target of unknown sex the male and female thresholds count half each."""
    if h.sex[target] in (1, 2) or not par.sex_specific:
        return [par.individual_threshold(h.sex[target])]
    return [par.threshold_male, par.threshold_female]


def _draw(rng, mu, s, phenotype: int, threshold: float):
    if phenotype == 2:
        return truncated_normal(rng, mu, s, lower=threshold)
    if phenotype == 1:
        return truncated_normal(rng, mu, s, upper=threshold)
    return mu + s * rng.standard_normal(np.shape(mu))


def gelman_rubin(trace: np.ndarray) -> float:
    """Gelman–Rubin R̂; ``trace`` is [step, chain]. Values near 1 indicate convergence."""
    n_steps, m = trace.shape
    if m < 2 or n_steps < 2:
        return math.nan
    means = trace.mean(axis=0)
    B = n_steps * means.var(ddof=1)
    W = trace.var(axis=0, ddof=1).mean()
    if W <= 0:
        return 1.0
    var_plus = (n_steps - 1) / n_steps * W + B / n_steps
    return float(math.sqrt(var_plus / W))
