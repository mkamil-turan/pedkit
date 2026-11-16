"""Exact likelihoods for the Mendelian modes.

Two independent exact methods share the genetic model of ``modes``:

* ``peeling_likelihood``: generalised Elston–Stewart peeling (variable
  elimination, ``pedkit.peeling``). Its cost grows with the width of the
  pedigree, not with its size, so it has no individual limit. This is
  what ``infer_mode_exact`` uses.
* ``enumerate_genotypes`` / ``enumeration_likelihood``: sums over every
  joint genotype assignment. Its cost grows exponentially with the number
  of individuals (up to ``MAX_ENUMERATION``); the Bayes risk table needs
  it, and it serves as an independent check of peeling and Monte Carlo.

MZ twins share one genotype in both methods.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator

import numpy as np

from ..model import Disease, Mode, Pedigree
from ..peeling import Factor, eliminate
from . import modes as mm
from .prepare import Prepared, prepare
from .result import GenotypePool, ModeInference, ModeResult, compute_posteriors

MAX_ENUMERATION = 14


def infer_mode_exact(
    ped: Pedigree,
    disease: str | Disease | None = None,
    modes: Iterable[Mode] = mm.MENDELIAN_MODES,
    prior: dict[Mode, float] | None = None,
    default_allele_freq: float = 0.01,
) -> ModeInference:
    h = prepare(ped, disease, default_allele_freq)
    modes = list(modes)
    prior = prior or {m: 1 / len(modes) for m in modes}
    results = {}
    for m in modes:
        r = peeling_likelihood(h, m)
        r.prior = prior.get(m, 0.0)
        results[m] = r
    compute_posteriors(results)
    return ModeInference(h.disease, h.q, results, method="exact")


def states(h: Prepared, i: int, mode: Mode) -> list[tuple[bool, int, int, float]]:
    """(male, M, P, probability of the sex) options of individual ``i``."""
    s = h.sex[i]
    if s in (1, 2):
        sexes = [(s == 1, 1.0)]
    elif mm.sex_linked(mode):
        sexes = [(True, 0.5), (False, 0.5)]
    else:
        sexes = [(False, 1.0)]
    out = []
    for male, ps in sexes:
        ms = (0, 1) if mm.has_m_slot(mode) else (0,)
        pp = (0, 1) if mm.has_p_slot(mode, male) else (0,)
        out += [(male, m, p, ps) for m in ms for p in pp]
    return out


def _transmission(v: int, parent_state, slots: tuple[str, ...], q: float, mu: float) -> float:
    """P(child receives allele ``v`` | parent state); founder side if ``parent_state`` is None."""
    if parent_state is None:
        return q if v else 1 - q
    _, m, p = parent_state[:3]
    allele = {"M": m, "P": p}
    prob = 0.0
    for s in slots:
        t = allele[s]
        prob += (1.0 if v else 0.0) if t else (mu if v else 1 - mu)
    return prob / len(slots)


def _observation(h: Prepared, i: int, mode: Mode, male: bool, m: int, p: int) -> float:
    g = 1.0
    if h.phenotype[i]:
        g = float(mm.phenotype_probability(mode, m, p, male, h.phenotype[i],
                                           h.penetrance, h.phenocopy))
    if h.genotype[i] >= 0 and int(mm.code(mode, m, p, male)) != h.genotype[i]:
        g = 0.0
    return g


def _twin_male(h: Prepared, j: int, state_male: bool) -> bool:
    return (h.sex[j] == 1) if h.sex[j] in (1, 2) else state_male


# -- peeling -------------------------------------------------------------------------

def _factors(h: Prepared, mode: Mode, indicator: bool = False,
             skip_observations: frozenset[int] | set[int] = frozenset()
             ) -> tuple[list[Factor], dict[int, int], dict[int, list]]:
    """Factors of the pedigree. Returns (factors, variable of each individual,
    state list of each variable). MZ twins are aliased to their first co-twin.
    With ``indicator`` the observation factors are replaced by 1[weight > 0];
    individuals in ``skip_observations`` contribute no observation factor."""
    var = {i: (h.mz[i] if h.mz[i] >= 0 else i) for i in range(h.n)}
    st = {i: states(h, i, mode) for i in range(h.n) if h.mz[i] < 0}
    q, mu = h.q, h.mutation
    m_slots = mm.maternal_slots(mode)
    p_slots = mm.paternal_slots(mode)
    factors: list[Factor] = []

    for i in range(h.n):
        if h.mz[i] >= 0:
            continue
        mo = var[h.mother[i]] if h.mother[i] >= 0 and m_slots else None
        fa = var[h.father[i]] if h.father[i] >= 0 and p_slots else None
        parents = [v for v in (mo, fa) if v is not None]
        shape = [len(st[v]) for v in parents] + [len(st[i])]
        table = np.zeros(shape)
        for idx in np.ndindex(*shape):
            ps = dict(zip(parents, idx[:-1]))
            male, m, p, sp = st[i][idx[-1]]
            val = sp
            if m_slots:
                val *= _transmission(m, st[mo][ps[mo]] if mo is not None else None,
                                     m_slots, q, mu)
            if mm.has_p_slot(mode, male):
                val *= _transmission(p, st[fa][ps[fa]] if fa is not None else None,
                                     p_slots, q, mu)
            table[idx] = val
        factors.append(Factor(tuple(parents) + (i,), table))

    for i in range(h.n):
        if i in skip_observations:
            continue
        v = var[i]
        obs = np.ones(len(st[v]))
        for k, (male, m, p, _) in enumerate(st[v]):
            sm = _twin_male(h, i, male) if h.mz[i] >= 0 else male
            if h.mz[i] >= 0 and h.sex[i] in (1, 2) and h.sex[v] == 0 and mm.sex_linked(mode) \
                    and sm != male:
                obs[k] = 0.0
                continue
            w = _observation(h, i, mode, sm, m, p)
            obs[k] = (1.0 if w > 0 else 0.0) if indicator else w
        if not np.all(obs == 1.0):
            factors.append(Factor((v,), obs))
    return factors, var, st


def peeling_likelihood(h: Prepared, mode: Mode, marginals: bool = True) -> ModeResult:
    """Exact P(observations | mode), the prior probability of the supporting
    pool and the posterior genotype marginals, by peeling."""
    factors, var, st = _factors(h, mode)
    table, log_scale = eliminate(factors)
    likelihood = float(table) * math.exp(log_scale) if log_scale > -math.inf else 0.0

    ind_factors, _, _ = _factors(h, mode, indicator=True)
    t2, l2 = eliminate(ind_factors)
    support_prior = float(t2) * math.exp(l2) if l2 > -math.inf else 0.0

    marg = {}
    if marginals and likelihood > 0:
        cache: dict[int, np.ndarray] = {}
        for i, lab in enumerate(h.labels):
            v = var[i]
            if v not in cache:
                t, _ = eliminate(factors, keep=(v,))
                cache[v] = t / t.sum()
            dist: dict[str, float] = {}
            for k, (male, m, p, _) in enumerate(st[v]):
                sm = _twin_male(h, i, male) if h.mz[i] >= 0 else male
                c = mm.CODE_TEXT[int(mm.code(mode, m, p, sm))]
                dist[c] = dist.get(c, 0.0) + float(cache[v][k])
            marg[lab] = {c: pr for c, pr in sorted(dist.items(), key=lambda t: _code_order(t[0]))
                         if pr > 1e-15}
    return ModeResult(
        mode=mode, simulations=0,
        supporting=GenotypePool(probability=support_prior, marginal=marg),
        non_supporting=GenotypePool(probability=1 - support_prior),
        likelihood=likelihood, standard_error=0.0, effective_sample=float("inf"),
    )


def peeling_target_joints(h: Prepared, mode: Mode, target: int,
                          skip_observations: set[int]) -> np.ndarray:
    """Joint probability of each genotype code of ``target`` with the
    observations of everyone except ``skip_observations``, by peeling.
    Returns an array indexed by genotype code (see ``modes.CODE_TEXT``)."""
    factors, var, st = _factors(h, mode, skip_observations=skip_observations)
    v = var[target]
    table, log_scale = eliminate(factors, keep=(v,))
    out = np.zeros(mm.N_CODES)
    if log_scale == -math.inf:
        return out
    scale = math.exp(log_scale)
    for k, (male, m, p, _) in enumerate(st[v]):
        sm = _twin_male(h, target, male) if h.mz[target] >= 0 else male
        out[int(mm.code(mode, m, p, sm))] += float(table[k]) * scale
    return out


def _code_order(text: str) -> int:
    return {v: k for k, v in mm.CODE_TEXT.items()}[text]


# -- enumeration ----------------------------------------------------------------------

def enumerate_genotypes(h: Prepared, mode: Mode, prunable: set[int] | None = None
                        ) -> Iterator[tuple[dict[int, tuple[bool, int, int]], float, list[float]]]:
    """Yield every genotype assignment: (state, prior probability,
    per-individual observation weights). ``state[i] = (male, M, P)``.

    Branches with prior probability 0 are always pruned; branches with
    observation weight 0 only when the individual is in ``prunable``
    (default: everyone).
    """
    prunable = set(range(h.n)) if prunable is None else prunable
    q, mu = h.q, h.mutation
    m_slots = mm.maternal_slots(mode)
    p_slots = mm.paternal_slots(mode)
    state: dict[int, tuple[bool, int, int]] = {}
    weight = [1.0] * h.n

    def options(i: int):
        lead = h.mz[i]
        if lead >= 0:
            male_l, m, p = state[lead]
            male = _twin_male(h, i, male_l)
            if mm.sex_linked(mode) and h.sex[i] in (1, 2) and h.sex[lead] == 0 and male != male_l:
                return []
            return [(male, m, p, 1.0, True)]
        return [(male, m, p, ps, False) for male, m, p, ps in states(h, i, mode)]

    def rec(k: int, prior: float):
        if k == len(h.order):
            yield dict(state), prior, list(weight)
            return
        i = h.order[k]
        for male, m, p, ps, copied in options(i):
            o = ps
            if not copied:
                if m_slots:
                    o *= _transmission(m, state.get(h.mother[i]) if h.mother[i] >= 0 else None,
                                       m_slots, q, mu)
                if mm.has_p_slot(mode, male):
                    o *= _transmission(p, state.get(h.father[i]) if h.father[i] >= 0 else None,
                                       p_slots, q, mu)
            if o == 0:
                continue
            g = _observation(h, i, mode, male, m, p)
            if g == 0 and i in prunable:
                continue
            state[i] = (male, m, p)
            weight[i] = g
            yield from rec(k + 1, prior * o)
            del state[i]
            weight[i] = 1.0

    yield from rec(0, 1.0)


def enumeration_likelihood(h: Prepared, mode: Mode) -> ModeResult:
    if h.n > MAX_ENUMERATION:
        raise ValueError(f"Enumeration is limited to {MAX_ENUMERATION} individuals; "
                         "use peeling_likelihood.")
    total = 0.0
    support_prior = 0.0
    marginal = [[0.0] * mm.N_CODES for _ in range(h.n)]
    for state, prior, weight in enumerate_genotypes(h, mode):
        joint = prior * math.prod(weight)
        total += joint
        support_prior += prior
        for i, (male, m, p) in state.items():
            marginal[i][int(mm.code(mode, m, p, male))] += joint

    marg = {}
    if total > 0:
        for i, lab in enumerate(h.labels):
            marg[lab] = {mm.CODE_TEXT[c]: marginal[i][c] / total
                         for c in range(mm.N_CODES) if marginal[i][c] > 0}
    return ModeResult(
        mode=mode, simulations=0,
        supporting=GenotypePool(probability=support_prior, marginal=marg),
        non_supporting=GenotypePool(probability=1 - support_prior),
        likelihood=total, standard_error=0.0, effective_sample=float("inf"),
    )
