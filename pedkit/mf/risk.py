"""Individual risk for a multifactorial disease: Elston–Stewart and MCMC together."""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass

import numpy as np

from ..inference.prepare import Prepared
from ..model import Disease, Mode, Pedigree
from .elston import DEFAULT_LOCI, ESResult, es_risk
from .mcmc import MCMCResult, mcmc_risk
from .model import MFParameters, kinship_matrix, parameters, prepare_mf, target_index

_SEX = {0: "unknown", 1: "male", 2: "female"}
_SEX_SHORT = {0: "?", 1: "male", 2: "female"}
_PHENOTYPE = {0: "unknown", 1: "unaffected", 2: "affected"}

METHODS = ("es", "mcmc", "both")


@dataclass
class MFRiskResult:
    individual: str
    target: int
    disease: Disease
    par: MFParameters
    h: Prepared                       # preparation including the target's own phenotype
    es: ESResult | None
    mcmc: MCMCResult | None
    kinship: np.ndarray               # kinship coefficients φ with the target

    @property
    def risk(self) -> float:
        return self.es.risk if self.es else self.mcmc.risk

    @property
    def population_risk(self) -> float:
        """Prevalence in the target's sex (their average if the sex is unknown)."""
        c = self.h.sex[self.target]
        if c == 1:
            return self.par.prevalence_male
        if c == 2:
            return self.par.prevalence_female
        return (self.par.prevalence_male + self.par.prevalence_female) / 2

    @property
    def relative_risk(self) -> float:
        return self.risk / self.population_risk

    def report(self) -> str:
        return report(self)

    def plots(self, folder: str, fmt: str = "png") -> list[str]:
        from .plots import save_plots
        return save_plots(self, folder, fmt)


def mf_risk(
    ped: Pedigree,
    individual: str,
    disease: str | Disease | None = None,
    method: str = "both",
    prevalence: float | None = None,
    heritability: float | None = None,
    prevalence_male: float | None = None,
    prevalence_female: float | None = None,
    loci: int = DEFAULT_LOCI,
    chains: int = 64,
    samples: int = 4000,
    burn_in: int = 1000,
    seed: int | None = None,
) -> MFRiskResult:
    """Risk of ``individual`` being affected, given the phenotypes of the relatives.

    ``method``: "es" (Elston–Stewart, hypergeometric polygenic model),
    "mcmc" (Gibbs sampling) or "both". An observed phenotype of the
    individual is ignored: the result is the risk that would be expected
    if the individual's status were unknown.
    """
    if method not in METHODS:
        raise ValueError(f"Unknown method {method!r}: es, mcmc or both")
    h = prepare_mf(ped, disease)
    if h.disease.mode not in (Mode.MF, Mode.UNKNOWN):
        raise ValueError(f"{h.disease.id} has mode {h.disease.mode.value} in the file; "
                         "a multifactorial risk needs mode=MF.")
    par = parameters(h.disease, prevalence, heritability, prevalence_male, prevalence_female)
    target = target_index(h, ped, individual)
    h0 = dataclasses.replace(h, phenotype=[0 if i == target else f
                                           for i, f in enumerate(h.phenotype)])

    es = es_risk(h0, par, target, n=loci, individual=False) if method in ("es", "both") else None
    if es is not None:
        es.individual_risk = es_risk(h, par, target, n=loci, individual=True).individual_risk
    mc = (mcmc_risk(h0, par, target, chains=chains, samples=samples, burn_in=burn_in, seed=seed)
          if method in ("mcmc", "both") else None)
    return MFRiskResult(h.labels[target], target, h.disease, par, h, es, mc,
                        kinship_matrix(h)[target])


# -- report -----------------------------------------------------------------------------

def _f(x: float | None, digits: int = 4) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    if x == 0:
        return "0"
    if abs(x) < 1e-3:
        return f"{x:.3e}"
    return f"{x:.{digits}f}".rstrip("0").rstrip(".")


def report(r: MFRiskResult) -> str:
    """Lines stay within 76 characters (to fit a narrow panel and a terminal)."""
    p, h, i = r.par, r.h, r.target
    d = r.disease
    lines = [f"Multifactorial risk · {r.individual} · {d.id}" + (f" ({d.name})" if d.name else ""),
             f"Model: liability-threshold (Falconer) · h² = {_f(p.heritability)}"]
    if p.sex_specific:
        lines.append(f"  male: K = {_f(p.prevalence_male)}, threshold T = "
                     f"{p.threshold_male:.3f} · female: K = {_f(p.prevalence_female)}, "
                     f"T = {p.threshold_female:.3f}")
    else:
        lines.append(f"  K = {_f(p.prevalence)}, threshold T = {p.threshold:.3f}")
    lines.append(f"Individual: {r.individual} · sex {_SEX[h.sex[i]]} · phenotype "
                 f"{_PHENOTYPE[h.phenotype[i]]}")
    if h.phenotype[i]:
        lines += ["  The individual's own phenotype was ignored: the risk below is",
                  "  the one expected from the relatives if the status were unknown."]
    if h.sex[i] == 0 and p.sex_specific:
        lines.append("  Sex unknown: the male and female risks were averaged.")
    lines += ["", "Risk of being affected given the relatives' phenotypes:"]
    if r.es:
        lines += [
            f"  Elston–Stewart (HPM, Richardson n→∞): {_f(r.es.risk)}",
            f"      n={r.es.loci} loci: {_f(r.es.risk_n)} · "
            f"n={2 * r.es.loci} loci: {_f(r.es.risk_2n)}",
        ]
    if r.mcmc:
        m = r.mcmc
        lines += [
            f"  MCMC (Gibbs, {m.chains} chains × {m.samples} steps, "
            f"burn-in {m.burn_in}): {_f(m.risk)}",
            f"      standard error {_f(m.standard_error)} · count estimate "
            f"(L > T) {_f(m.risk_count)}",
            f"      R̂ liability {m.rhat_liability:.4f} · R̂ risk {m.rhat_risk:.4f}",
        ]
    if r.es and r.mcmc and r.mcmc.standard_error > 0:
        diff = (r.mcmc.risk - r.es.risk) / r.mcmc.standard_error
        lines.append(f"  Difference of the two methods: {diff:+.1f} standard errors")
    lines += [
        f"  Population risk (prevalence): {_f(r.population_risk)}",
        f"  Relative risk (λ = risk / prevalence): {r.relative_risk:.2f}",
        "",
        "Individuals:",
        f"  {'Indiv.':<8}{'Sex':<7}{'Phenotype':<12}{'φ':>7}{'ES risk*':>11}"
        f"{'L mean':>8}{'L 90% interval':>17}",
    ]
    for j, lab in enumerate(h.labels):
        es_r = _f(r.es.individual_risk.get(j)) if r.es and r.es.individual_risk else "—"
        if r.mcmc:
            ml = f"{r.mcmc.mean_liability[j]:+.2f}"
            a = r.mcmc.liability_interval[j]
            iv = f"[{a[0]:+.2f}, {a[1]:+.2f}]"
        else:
            ml = iv = "—"
        lines.append(f"  {lab:<8}{_SEX_SHORT[h.sex[j]]:<7}{_PHENOTYPE[h.phenotype[j]]:<12}"
                     f"{_f(r.kinship[j], 4):>7}{es_r:>11}{ml:>8}{iv:>17}")
    lines += [
        "  φ: kinship coefficient with the target.",
        "  * ES risk: the individual's risk from the relatives if its own",
        "    phenotype were unknown.",
        "  L mean, L 90% interval: MCMC posterior mean and 5-95% quantiles of",
        "  the liability given all data (population mean 0, SD 1); the",
        f"  target's ({r.individual}) own phenotype is ignored here as well.",
    ]
    return "\n".join(lines)
