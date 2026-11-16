"""Elston–Stewart risk with the hypergeometric polygenic model (HPM).

The infinitesimal polygenic model cannot be peeled because its genotype
space is continuous. In the finite polygenic model (Cannings, Thompson and
Skolnick 1978; Fernando, Stricker and Elston 1994) the heritable part is
represented by the sum of 2n alleles of equal effect: an individual's
"polygenotype" is its number k ∈ {0..2n} of plus alleles, Binomial(2n, 1/2)
in founders. A parent transmits n of its k alleles without replacement
(hypergeometric); this preserves the parent–offspring covariance and the
Mendelian sampling variance exactly. As n grows the model converges to the
infinitesimal model.

At finite n the binomial tails are thinner than the normal ones, and for a
rare disease the risk is decided exactly in that tail; the approximation
error shrinks as 1/n. The risk is therefore computed with n and 2n loci
and carried to the infinitesimal limit by Richardson extrapolation
(R = 2·r(2n) - r(n)).

Breeding value A_k = h·(k - n)/√(n/2) (variance h²), environment
E ~ N(0, 1-h²):
    P(affected | k) = P(A_k + E > T) = 1 - Φ((T - A_k)/√(1-h²)).

MZ twins share one polygenotype (one variable) but have independent E.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from ..inference.prepare import Prepared
from ..peeling import Factor, eliminate
from .model import MFParameters
from .normal import phi, upper_tail

DEFAULT_LOCI = 20


@lru_cache(maxsize=8)
def _hpm(n: int) -> tuple[np.ndarray, np.ndarray]:
    """(founder prior [2n+1], transmission tensor T[k_mother, k_father, k_child])."""
    S = 2 * n + 1
    founder = np.array([math.comb(2 * n, k) for k in range(S)], dtype=float) / 2.0 ** (2 * n)
    # H[k, j]: probability of a gamete with j plus alleles from a parent with k
    H = np.zeros((S, n + 1))
    denom = math.comb(2 * n, n)
    for k in range(S):
        for j in range(n + 1):
            if j <= k and n - j <= 2 * n - k:
                H[k, j] = math.comb(k, j) * math.comb(2 * n - k, n - j) / denom
    T = np.zeros((S, S, S))
    for j1 in range(n + 1):
        for j2 in range(n + 1):
            T[:, :, j1 + j2] += np.outer(H[:, j1], H[:, j2])
    return founder, T


def breeding_values(n: int, h2: float) -> np.ndarray:
    k = np.arange(2 * n + 1)
    return math.sqrt(h2) * (k - n) / math.sqrt(n / 2)


def affection_probability(n: int, h2: float, threshold: float) -> np.ndarray:
    """P(affected | polygenotype k), k = 0..2n."""
    A = breeding_values(n, h2)
    return upper_tail((threshold - A) / math.sqrt(1 - h2))


def _observation(h: Prepared, par: MFParameters, n: int, i: int) -> np.ndarray:
    S = 2 * n + 1
    f = h.phenotype[i]
    if f == 0:
        return np.ones(S)
    p = affection_probability(n, par.heritability, par.individual_threshold(h.sex[i]))
    return p if f == 2 else 1 - p


def variable(h: Prepared, i: int) -> int:
    """Peeling variable of individual ``i`` (an MZ twin uses its co-twin's)."""
    return h.mz[i] if h.mz[i] >= 0 else i


def factors(h: Prepared, par: MFParameters, n: int,
            unobserved: int | None = None) -> list[Factor]:
    """Factors of the pedigree; the own phenotype of ``unobserved`` is left out."""
    founder, T = _hpm(n)
    out: list[Factor] = []
    for i in range(h.n):
        if h.mz[i] < 0:
            a, b = h.mother[i], h.father[i]
            a = variable(h, a) if a >= 0 else -1
            b = variable(h, b) if b >= 0 else -1
            if a >= 0 and b >= 0:
                out.append(Factor((a, b, i), T))
            elif a >= 0 or b >= 0:
                # One parent: the missing parent comes from the founder distribution
                half = np.einsum("y,xyz->xz", founder, T)
                out.append(Factor((a if a >= 0 else b, i), half))
            else:
                out.append(Factor((i,), founder))
        if i != unobserved and h.phenotype[i]:
            out.append(Factor((variable(h, i),), _observation(h, par, n, i)))
    return out


@dataclass
class ESResult:
    risk: float                       # risk carried to infinitely many loci (Richardson)
    risk_n: float                     # HPM risk with n loci
    risk_2n: float                    # HPM risk with 2n loci
    posterior: np.ndarray             # P(k_target | data), model with 2n loci
    log_likelihood: float             # model with 2n loci
    individual_risk: dict[int, float]  # i -> P(i affected | data of the others)
    loci: int                         # n


def target_risk(h: Prepared, par: MFParameters, target: int, n: int = DEFAULT_LOCI,
                own_phenotype: bool = False) -> tuple[float, np.ndarray]:
    """(risk, polygenotype posterior), HPM with n loci. If ``own_phenotype``
    is False the target's observed phenotype is ignored ("risk from the
    relatives")."""
    c = factors(h, par, n, unobserved=None if own_phenotype else target)
    table, _ = eliminate(c, keep=(variable(h, target),))
    total = table.sum()
    if total <= 0:
        raise ValueError("The observations are impossible under these parameters "
                         "(likelihood 0).")
    posterior = table / total
    return float(posterior @ _target_affection(h, par, n, target)), posterior


def richardson(r_n: float, r_2n: float) -> float:
    """Removes the error that shrinks as 1/n; clipped to the probability range."""
    return min(max(2 * r_2n - r_n, 0.0), 1.0)


def extrapolated_risk(h: Prepared, par: MFParameters, target: int,
                      n: int = DEFAULT_LOCI, own_phenotype: bool = False) -> float:
    r1 = target_risk(h, par, target, n, own_phenotype)[0]
    r2 = target_risk(h, par, target, 2 * n, own_phenotype)[0]
    return richardson(r1, r2)


def _target_affection(h: Prepared, par: MFParameters, n: int, i: int) -> np.ndarray:
    """For a target of unknown sex the male and female thresholds count half each."""
    if h.sex[i] in (1, 2) or not par.sex_specific:
        return affection_probability(n, par.heritability, par.individual_threshold(h.sex[i]))
    return 0.5 * (affection_probability(n, par.heritability, par.threshold_male)
                  + affection_probability(n, par.heritability, par.threshold_female))


def es_risk(h: Prepared, par: MFParameters, target: int,
            n: int = DEFAULT_LOCI, individual: bool = True) -> ESResult:
    r1, _ = target_risk(h, par, target, n, own_phenotype=True)
    r2, posterior = target_risk(h, par, target, 2 * n, own_phenotype=True)
    table, log_scale = eliminate(factors(h, par, 2 * n))
    log_l = math.log(float(table)) + log_scale if float(table) > 0 else -math.inf
    risks = {}
    if individual:
        for i in range(h.n):
            risks[i] = extrapolated_risk(h, par, i, n)
    return ESResult(richardson(r1, r2), r1, r2, posterior, log_l, risks, n)


def liability_density(posterior: np.ndarray, par: MFParameters, x: np.ndarray) -> np.ndarray:
    """Target's liability density from the HPM posterior: Σ_k P(k) · N(A_k, 1-h²)."""
    n = (len(posterior) - 1) // 2
    A = breeding_values(n, par.heritability)
    s = math.sqrt(1 - par.heritability)
    return (posterior[:, None] * phi((x[None, :] - A[:, None]) / s) / s).sum(axis=0)
