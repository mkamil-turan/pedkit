"""Bayes risk for a given individual.

The classical genetic-counselling table is built; the hypotheses are the
possible genotypes of the target individual:

    prior(g)        = P(g | upper data)
    conditional(g)  = P(lower data | g, upper data)
    joint(g)        = prior × conditional
    posterior(g)    = joint / Σ joint = P(g | all data)

The "lower data" are the observations of the target's descendants and of
their married-in (founder) parents; everything else, including the target's
own phenotype and test result, is "upper data". This is the general form of
the textbook split "affected brother, prior 2/3; two healthy children,
conditional (1/2)²".

Without a given mode a separate table is built for each mode and the risk
is averaged with the posterior probabilities of the modes.

Methods: "exact" uses peeling (no size limit), "enumeration" sums over all
genotype assignments (at most ``MAX_ENUMERATION`` individuals; an
independent check), "monte-carlo" uses gene dropping. "auto" is "exact".
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np

from ..inference import modes as mm
from ..inference.diagnosis import MODE_NAME
from ..inference.exact import MAX_ENUMERATION, enumerate_genotypes, peeling_target_joints
from ..inference.montecarlo import simulate_batch
from ..inference.prepare import Prepared, prepare
from ..model import Disease, Mode, Pedigree
from .events import event_name, event_probability, valid_event

METHODS = ("auto", "exact", "enumeration", "monte-carlo")


# -- result structures ---------------------------------------------------------------

@dataclass
class BayesRow:
    genotype: str
    prior: float
    conditional: float | None    # undefined when the prior is 0
    joint: float
    posterior: float
    event: float                 # P(event | this genotype)


@dataclass
class ModeRisk:
    mode: Mode
    weight: float | None             # P(mode | data)
    likelihood: float                # P(data | mode)
    rows: list[BayesRow] = field(default_factory=list)
    risk: float | None = None        # P(event | data, mode)
    supporting: int | None = None    # Monte Carlo: compatible simulations


@dataclass
class RiskResult:
    individual: str
    disease: Disease
    event: str
    method: str
    lower_individuals: list[str]      # individuals giving the conditional information
    modes: list[ModeRisk]
    risk: float | None                # mode-averaged P(event | data)
    simulations: int = 0

    def report(self) -> str:
        return report(self)


# -- main function ---------------------------------------------------------------------

def bayes_risk(
    ped: Pedigree,
    individual: str,
    event: str = "carrier",
    disease: str | Disease | None = None,
    modes: Mode | Iterable[Mode] | None = None,
    method: str = "auto",
    n: int = 1_000_000,
    batch: int = 100_000,
    seed: int | None = None,
    proposal_freq: float | None = None,
    prior: dict[Mode, float] | None = None,
    default_allele_freq: float = 0.01,
) -> RiskResult:
    """Bayes risk of ``event`` for ``individual``.

    Without ``modes``: only the disease's ``mode`` if it is given in the
    file, every Mendelian mode (averaged) if it is ``unknown``.
    """
    if not valid_event(event):
        raise ValueError(f"Unknown event {event!r}; valid: carrier, mutant, affected, "
                         "intermediate, genotype:0/1 ...")
    h = prepare(ped, disease, default_allele_freq)
    if h.disease.mode == Mode.MF:
        raise ValueError(f"{h.disease.id} is multifactorial (mode=MF); a Bayes risk over "
                         "Mendelian genotypes is not defined. Use: pedkit mf FILE.pedx "
                         f"--individual {individual}")
    target = target_index(h, ped, individual)
    modes = _mode_list(h, modes)
    prior = prior or {m: 1 / len(modes) for m in modes}
    if method not in METHODS:
        raise ValueError(f"Unknown method {method!r}; one of {', '.join(METHODS)}")
    if method == "auto":
        method = "exact"
    if method == "enumeration" and h.n > MAX_ENUMERATION:
        raise ValueError(f"Enumeration is limited to {MAX_ENUMERATION} individuals; "
                         "use --method exact.")

    lower = lower_set(h, target)
    rng = np.random.default_rng(seed)
    results: list[ModeRisk] = []
    for m in modes:
        support = None
        if method == "exact":
            upper = peeling_target_joints(h, m, target, lower)
            full = peeling_target_joints(h, m, target, set())
        elif method == "enumeration":
            upper, full = _enumeration_joints(h, m, target, lower)
        else:
            upper, full, support = _mc_joints(h, m, target, lower, n, batch, rng, proposal_freq)
        results.append(_mode_table(h, m, target, event, upper, full, support))

    total = sum(prior.get(r.mode, 0) * r.likelihood for r in results)
    for r in results:
        r.weight = prior.get(r.mode, 0) * r.likelihood / total if total > 0 else None
    risk = (sum(r.weight * r.risk for r in results if r.weight) if total > 0 else None)
    results.sort(key=lambda r: -(r.weight or 0))
    return RiskResult(
        individual=h.labels[target], disease=h.disease, event=event, method=method,
        lower_individuals=[h.labels[i] for i in sorted(lower)], modes=results, risk=risk,
        simulations=n if method == "monte-carlo" else 0,
    )


def lower_set(h: Prepared, target: int) -> set[int]:
    """Individuals giving the conditional information: the target's
    descendants and the married-in (parentless) other parents of them."""
    children: dict[int, list[int]] = {i: [] for i in range(h.n)}
    for i in range(h.n):
        for p in (h.mother[i], h.father[i]):
            if p >= 0:
                children[p].append(i)
    desc: set[int] = set()
    stack = list(children[target])
    while stack:
        c = stack.pop()
        if c not in desc:
            desc.add(c)
            stack += children[c]
    spouses = {p for c in desc for p in (h.mother[c], h.father[c])
               if p >= 0 and p != target and p not in desc and h.mother[p] < 0 and h.father[p] < 0}
    return desc | spouses


# -- engines -------------------------------------------------------------------------------

def _enumeration_joints(h: Prepared, mode: Mode, target: int, lower: set[int]):
    upper = np.zeros(mm.N_CODES)
    full = np.zeros(mm.N_CODES)
    upper_set = [i for i in range(h.n) if i not in lower]
    for state, prior, w in enumerate_genotypes(h, mode, prunable=set(upper_set)):
        male, m, p = state[target]
        c = int(mm.code(mode, m, p, male))
        w_up = prior * math.prod(w[i] for i in upper_set)
        upper[c] += w_up
        full[c] += w_up * math.prod(w[i] for i in lower)
    return upper, full


def _mc_joints(h: Prepared, mode: Mode, target: int, lower: set[int], n: int, batch: int,
               rng: np.random.Generator, proposal: float | None):
    upper = np.zeros(mm.N_CODES)
    full = np.zeros(mm.N_CODES)
    upper_set = [i for i in range(h.n) if i not in lower]
    support = 0
    done = 0
    while done < n:
        b = min(batch, n - done)
        codes, imp, w_obs, W = simulate_batch(h, mode, b, rng, proposal)
        w_up = imp * W[upper_set].prod(axis=0)
        w_all = imp * w_obs
        upper += np.bincount(codes[target], weights=w_up, minlength=mm.N_CODES)
        full += np.bincount(codes[target], weights=w_all, minlength=mm.N_CODES)
        support += int((w_obs > 0).sum())
        done += b
    return upper / n, full / n, support


def _mode_table(h: Prepared, mode: Mode, target: int, event: str,
                upper: np.ndarray, full: np.ndarray, support: int | None) -> ModeRisk:
    likelihood = float(full.sum())
    upper_total = float(upper.sum())
    r = ModeRisk(mode, None, likelihood, supporting=support)
    if upper_total <= 0:
        return r
    # The joint column is prior × conditional = P(g, lower data | upper data)
    for c in range(mm.N_CODES):
        if upper[c] <= 0 and full[c] <= 0:
            continue
        prior_g = upper[c] / upper_total
        cond = full[c] / upper[c] if upper[c] > 0 else None
        r.rows.append(BayesRow(
            genotype=mm.CODE_TEXT[c], prior=float(prior_g),
            conditional=None if cond is None else float(cond),
            joint=float(full[c] / upper_total),
            posterior=float(full[c] / likelihood) if likelihood > 0 else 0.0,
            event=event_probability(event, mode, c, h.penetrance, h.phenocopy,
                                    h.phenotype[target]),
        ))
    if likelihood > 0:
        r.risk = sum(row.posterior * row.event for row in r.rows)
    return r


def target_index(h: Prepared, ped: Pedigree, individual: str) -> int:
    if individual in h.labels:
        return h.labels.index(individual)
    candidates = [i for i, ind in enumerate(ped.individuals) if ind.id == individual]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise KeyError(f"'{individual}' is not in the pedigree.")
    raise KeyError(f"'{individual}' occurs in more than one family; write it as 'FAMILY/ID'.")


def _mode_list(h: Prepared, modes) -> list[Mode]:
    if isinstance(modes, Mode):
        return [modes]
    if modes is not None:
        return list(modes)
    if h.disease.mode is not Mode.UNKNOWN:
        return [h.disease.mode]
    return list(mm.MENDELIAN_MODES)


# -- report -------------------------------------------------------------------------------------

def fraction(x: float | None, max_denominator: int = 200) -> str:
    """Write the value as a decimal; add the fraction when it is very close to a simple one."""
    if x is None:
        return "—"
    if x == 0:
        return "0"
    if x == 1:
        return "1"
    f = Fraction(x).limit_denominator(max_denominator)
    dec = f"{x:.3e}" if x < 1e-3 else f"{x:.4f}".rstrip("0").rstrip(".")
    if f.denominator > 1 and abs(float(f) - x) < 1e-9 * max(1.0, abs(x)):
        return f"{f.numerator}/{f.denominator} ({dec})"
    return dec


def report(r: RiskResult) -> str:
    # A Monte Carlo estimate is not shown as a fraction (69/100 is not an exact ratio)
    def _f(x):
        return fraction(x) if r.method in ("exact", "enumeration") else _decimal(x)

    d = r.disease
    lines = [
        f"Bayes risk · {r.individual} · {d.id}" + (f" ({d.name})" if d.name else ""),
        f"Event: {event_name(r.event)} · method: {r.method}"
        + (f" ({r.simulations:,} simulations/mode)" if r.simulations else ""),
        "Prior information: the individual's own observations, ancestors and collateral relatives",
        "Conditional information: " + (", ".join(r.lower_individuals) if r.lower_individuals
                                       else "none (no descendants or married-in spouses)"),
        "",
    ]
    if r.risk is None:
        lines.append("No mode is compatible with the observations; the risk cannot be computed.")
        return "\n".join(lines)

    several = len(r.modes) > 1
    if several:
        lines.append("Mode weights P(mode | data):")
        for m in r.modes:
            lines.append(f"  {m.mode.value:<4} {_f(m.weight):>22}   "
                         f"P(event | data, {m.mode.value}) = {_f(m.risk)}")
        lines.append("")

    for m in r.modes:
        if not m.weight:
            continue
        title = f"── {m.mode.value} · {MODE_NAME[m.mode]}"
        if several:
            title += f" · weight {_f(m.weight)}"
        if m.supporting is not None:
            title += f" · {m.supporting:,} compatible simulations"
        lines += [title + " ──",
                  f"  {'Genotype':<9}{'Prior':>20}{'Conditional':>20}{'Joint':>20}"
                  f"{'Posterior':>20}{'P(event|g)':>12}"]
        for row in m.rows:
            lines.append(f"  {row.genotype:<9}{_f(row.prior):>20}{_f(row.conditional):>20}"
                         f"{_f(row.joint):>20}{_f(row.posterior):>20}{_f(row.event):>12}")
        lines += [f"  P({event_name(r.event)} | data, {m.mode.value}) = {_f(m.risk)}", ""]

    lines.append(f"Bayes risk: P({event_name(r.event)}) = {_f(r.risk)}"
                 + (" (averaged over modes)" if several else ""))
    return "\n".join(lines)


def _decimal(x: float | None) -> str:
    if x is None:
        return "—"
    if x in (0, 1):
        return str(int(x))
    return f"{x:.3e}" if x < 1e-3 else f"{x:.4f}".rstrip("0").rstrip(".")
