"""Monte Carlo mode-of-inheritance inference (gene dropping).

For each Mendelian mode, millions of times:
  1. Founders receive random alleles from the allele frequency (Hardy–Weinberg).
  2. Alleles are passed down the pedigree by Mendel's rules (autosomal, X,
     Y or mitochondrial, depending on the mode), with optional new
     mutation. MZ twins receive one shared genotype.
  3. The simulated genotypes are compared with the observations: tested
     genotypes must match exactly, phenotypes are weighted by the
     penetrance/phenocopy probabilities.

Simulations with weight above zero form the **supporting** genotype pool,
the rest the **non-supporting** pool. The mean weight is an unbiased
estimate of the likelihood P(observations | mode); Bayes' rule gives
P(mode | observations).

With rare alleles few simulations support a mode. If ``proposal_freq`` is
given, founder alleles are drawn with that higher frequency and every
simulation is reweighted by q/q' ratios (importance sampling); estimates
stay unbiased and the supporting pool grows.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np

from ..model import Disease, Mode, Pedigree
from . import modes as mm
from .prepare import Prepared, prepare
from .result import (
    GenotypeConfiguration, GenotypePool, ModeInference, ModeResult, compute_posteriors,
)

Progress = Callable[[Mode, int, int], None]  # (mode, done, total)


def infer_mode(
    ped: Pedigree,
    disease: str | Disease | None = None,
    modes: Iterable[Mode] = mm.MENDELIAN_MODES,
    n: int = 1_000_000,
    batch: int = 100_000,
    seed: int | None = None,
    prior: dict[Mode, float] | None = None,
    default_allele_freq: float = 0.01,
    proposal_freq: float | None = None,
    max_configurations: int | None = 20,
    progress: Progress | None = None,
) -> ModeInference:
    """Run ``n`` simulations per mode and return the posterior probabilities
    of the modes with their supporting/non-supporting genotype pools. With
    ``max_configurations`` None every configuration of each pool is kept."""
    h = prepare(ped, disease, default_allele_freq)
    modes = list(modes)
    prior = prior or {m: 1 / len(modes) for m in modes}
    rng = np.random.default_rng(seed)

    results: dict[Mode, ModeResult] = {}
    for mode in modes:
        acc = _Accumulator(h, mode, max_configurations)
        done = 0
        while done < n:
            b = min(batch, n - done)
            acc.add(*simulate_batch(h, mode, b, rng, proposal_freq))
            done += b
            if progress:
                progress(mode, done, n)
        res = acc.result()
        res.prior = prior.get(mode, 0.0)
        results[mode] = res
    compute_posteriors(results)
    return ModeInference(h.disease, h.q, results)


# -- one batch ----------------------------------------------------------------------

def simulate_batch(h: Prepared, mode: Mode, B: int, rng: np.random.Generator,
                   proposal: float | None):
    """``B`` simulations. Returns (codes[n,B], importance weight[B],
    observation weight[B], per-individual observation weight W[n,B])."""
    n = h.n
    M = np.zeros((n, B), dtype=np.uint8)
    P = np.zeros((n, B), dtype=np.uint8)
    MALE = np.zeros((n, B), dtype=bool)
    codes = np.zeros((n, B), dtype=np.int64)
    W = np.ones((n, B))
    log_imp = np.zeros(B)
    w_obs = np.ones(B)

    q = h.q
    qp = proposal if proposal is not None else q
    with np.errstate(divide="ignore"):
        lw1, lw0 = np.log(q / qp) if qp > 0 else 0.0, np.log((1 - q) / (1 - qp)) if qp < 1 else 0.0

    def founder_allele(mask=None) -> np.ndarray:
        nonlocal log_imp
        t = (rng.random(B) < qp).astype(np.uint8)
        if proposal is not None:
            lw = np.where(t == 1, lw1, lw0)
            log_imp = log_imp + (lw if mask is None else np.where(mask, lw, 0.0))
        return t

    def mutate(t: np.ndarray) -> np.ndarray:
        if h.mutation > 0:
            t = t | (rng.random(B) < h.mutation).astype(np.uint8)
        return t

    def transmit(parent: int, slots: tuple[str, ...]) -> np.ndarray:
        src = {"M": M[parent], "P": P[parent]}
        if len(slots) == 2:
            t = np.where(rng.random(B) < 0.5, src["M"], src["P"])
        else:
            t = src[slots[0]]
        return mutate(t)

    m_slots = mm.maternal_slots(mode)
    p_slots = mm.paternal_slots(mode)
    linked = mm.sex_linked(mode)

    for i in h.order:
        s = h.sex[i]
        twin = h.mz[i]
        if s == 1:
            male = np.ones(B, dtype=bool)
        elif s == 2:
            male = np.zeros(B, dtype=bool)
        elif twin >= 0:
            male = MALE[twin].copy()
        elif linked:
            male = rng.random(B) < 0.5  # sex unknown: one half each
        else:
            male = np.zeros(B, dtype=bool)
        MALE[i] = male

        if twin >= 0:
            M[i] = M[twin]
            P[i] = P[twin]
        else:
            if m_slots:
                M[i] = transmit(h.mother[i], m_slots) if h.mother[i] >= 0 else founder_allele()
            p_ok = mm.has_p_slot(mode, male)
            if p_ok.any():
                t = transmit(h.father[i], p_slots) if h.father[i] >= 0 else founder_allele(p_ok)
                P[i] = np.where(p_ok, t, 0)

        codes[i] = mm.code(mode, M[i], P[i], male)
        w = np.ones(B)
        if h.phenotype[i]:
            w = mm.phenotype_probability(mode, M[i], P[i], male, h.phenotype[i],
                                         h.penetrance, h.phenocopy)
        if h.genotype[i] >= 0:
            w = np.where(codes[i] == h.genotype[i], w, 0.0)
        W[i] = w
        w_obs *= w

    return codes, np.exp(log_imp), w_obs, W


# -- accumulation across batches -------------------------------------------------

class _Pool:
    def __init__(self, n: int):
        self.count = 0
        self.weight = 0.0
        self.configs: dict[int, list[float]] = {}  # key -> [count, weight]
        self.marginal = np.zeros((n, mm.N_CODES))

    def add(self, codes: np.ndarray, weight: np.ndarray, keys: np.ndarray | None):
        self.count += weight.size
        self.weight += float(weight.sum())
        for i in range(codes.shape[0]):
            self.marginal[i] += np.bincount(codes[i], weights=weight, minlength=mm.N_CODES)
        if keys is not None and keys.size:
            uniq, inv = np.unique(keys, return_inverse=True)
            counts = np.bincount(inv)
            weights = np.bincount(inv, weights=weight)
            for k, c, w in zip(uniq.tolist(), counts.tolist(), weights.tolist()):
                rec = self.configs.setdefault(k, [0, 0.0])
                rec[0] += c
                rec[1] += w


class _Accumulator:
    # base-6 key must fit in int64: 6**24 < 2**63
    MAX_KEYED = 24

    def __init__(self, h: Prepared, mode: Mode, limit: int | None):
        self.h = h
        self.mode = mode
        self.limit = limit
        self.total = 0
        self.total_imp = 0.0
        self.total_w = 0.0
        self.total_w2 = 0.0
        self.support = _Pool(h.n)
        self.against = _Pool(h.n)
        self.conflict = np.zeros(h.n)
        self.base = (mm.N_CODES ** np.arange(h.n)).astype(np.int64) \
            if h.n <= self.MAX_KEYED else None

    def add(self, codes, imp, w_obs, W) -> None:
        zero = W == 0
        w = imp * w_obs
        self.total += w.size
        self.total_imp += float(imp.sum())
        self.total_w += float(w.sum())
        self.total_w2 += float((w * w).sum())
        keys = self.base @ codes if self.base is not None else None
        d = w_obs > 0
        self.support.add(codes[:, d], w[d], keys[d] if keys is not None else None)
        k = ~d
        self.against.add(codes[:, k], imp[k], keys[k] if keys is not None else None)
        self.conflict += (zero[:, k] * imp[k]).sum(axis=1)

    def _pool(self, pool: _Pool) -> GenotypePool:
        labels = self.h.labels
        ranked = sorted(pool.configs.items(), key=lambda t: -t[1][1])
        configs = []
        for k, (count, weight) in ranked[: self.limit]:
            cs = [(k // mm.N_CODES ** i) % mm.N_CODES for i in range(self.h.n)]
            configs.append(GenotypeConfiguration(
                {labels[i]: mm.CODE_TEXT[c] for i, c in enumerate(cs)},
                int(count), weight / pool.weight if pool.weight > 0 else 0.0))
        marginal = {}
        if pool.weight > 0:
            for i, lab in enumerate(labels):
                marginal[lab] = {mm.CODE_TEXT[c]: float(pool.marginal[i, c] / pool.weight)
                                 for c in range(mm.N_CODES) if pool.marginal[i, c] > 0}
        return GenotypePool(
            count=pool.count,
            frequency=pool.count / self.total if self.total else 0.0,
            probability=0.0,
            n_configurations=len(pool.configs),
            configurations=configs,
            marginal=marginal,
        )

    def result(self) -> ModeResult:
        n = self.total
        mean = self.total_w / n
        var = max(self.total_w2 / n - mean ** 2, 0.0)
        support = self._pool(self.support)
        against = self._pool(self.against)
        # prior probability masses: importance weights of the supporting pool
        support_imp = self.total_imp - self.against.weight
        support.probability = support_imp / self.total_imp
        against.probability = self.against.weight / self.total_imp
        conflict = {lab: float(self.conflict[i] / self.against.weight)
                    if self.against.weight > 0 else 0.0
                    for i, lab in enumerate(self.h.labels)}
        ess = self.total_w ** 2 / self.total_w2 if self.total_w2 > 0 else 0.0
        return ModeResult(
            mode=self.mode, simulations=n, supporting=support, non_supporting=against,
            likelihood=mean, standard_error=(var / n) ** 0.5, effective_sample=ess,
            conflict=conflict,
        )
