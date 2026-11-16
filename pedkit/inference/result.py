"""Result structures and report of mode-of-inheritance inference."""

from __future__ import annotations

from dataclasses import dataclass, field

from ..model import Disease, Mode


@dataclass
class GenotypeConfiguration:
    """A joint genotype assignment of every individual in the pedigree."""

    genotypes: dict[str, str]   # individual -> "0/1", "1", "-" ...
    count: int                  # number of simulations with this assignment
    probability: float          # (weighted) relative probability within the pool


@dataclass
class GenotypePool:
    """Aggregate of the supporting or the non-supporting simulations."""

    count: int = 0              # simulations in the pool
    frequency: float = 0.0      # count / total simulations
    probability: float = 0.0    # prior probability of the pool under the mode
    n_configurations: int = 0
    configurations: list[GenotypeConfiguration] = field(default_factory=list)  # most frequent
    marginal: dict[str, dict[str, float]] = field(default_factory=dict)  # individual -> genotype -> P


@dataclass
class ModeResult:
    mode: Mode
    simulations: int
    supporting: GenotypePool
    non_supporting: GenotypePool
    likelihood: float            # P(observations | mode)
    standard_error: float
    effective_sample: float      # effective size of the weighted supporting pool
    conflict: dict[str, float] = field(default_factory=dict)  # individual -> share of conflicts
    prior: float = 0.0
    posterior: float | None = None   # P(mode | observations)
    bayes_factor: float | None = None  # likelihood ratio to the most likely mode


@dataclass
class ModeInference:
    disease: Disease
    allele_freq: float
    results: dict[Mode, ModeResult]
    method: str = "monte-carlo"      # "monte-carlo" or "exact"

    @property
    def ranked(self) -> list[ModeResult]:
        return sorted(self.results.values(), key=lambda r: -r.likelihood * r.prior)

    @property
    def most_likely(self) -> ModeResult | None:
        r = self.ranked
        return r[0] if r and r[0].likelihood > 0 else None

    @property
    def posteriors(self) -> dict[Mode, float | None]:
        return {m: r.posterior for m, r in self.results.items()}

    @property
    def diagnosis(self):
        """Diagnosis at the autosomal/gonosomal/mitochondrial and mode levels."""
        from .diagnosis import hierarchical_diagnosis
        return hierarchical_diagnosis(self)

    def report(self, details: bool = True, limit: int | None = None) -> str:
        return report(self, details, limit)


def compute_posteriors(results: dict[Mode, ModeResult]) -> None:
    total = sum(r.prior * r.likelihood for r in results.values())
    best = max((r.likelihood for r in results.values()), default=0.0)
    for r in results.values():
        r.posterior = r.prior * r.likelihood / total if total > 0 else None
        r.bayes_factor = r.likelihood / best if best > 0 else None


# -- report ----------------------------------------------------------------------

def fmt(x: float | None, digits: int = 4) -> str:
    if x is None:
        return "—"
    if x == 0:
        return "0"
    if abs(x) < 10 ** -(digits - 1) or abs(x) >= 1e5:
        return f"{x:.{digits - 2}e}"
    return f"{x:.{digits}f}".rstrip("0").rstrip(".")


def report(c: ModeInference, details: bool = True, limit: int | None = None) -> str:
    """With ``limit`` None every collected configuration is listed."""
    from .decision import summary_line

    d = c.disease
    first = next(iter(c.results.values()), None)
    n = first.simulations if first else 0
    lines = [
        c.diagnosis.report(),
        "",
        "═" * 72,
        f"Mode-of-inheritance inference · {d.id}" + (f" ({d.name})" if d.name else ""),
        f"q={fmt(c.allele_freq)} · penetrance={fmt(d.penetrance)} · phenocopy="
        f"{fmt(d.phenocopy)} · mutation rate={fmt(d.mutation_rate)} · "
        + (f"{n:,} simulations/mode" if c.method == "monte-carlo"
           else "exact computation (peeling)"),
        "",
        f"{'Mode':<6}{'Supporting':>12}{'Frequency':>11}{'P(data|mode)':>15}{'±SE':>11}"
        f"{'Prior':>8}{'Posterior':>11}",
    ]
    for r in c.ranked:
        lines.append(
            f"{r.mode.value:<6}{r.supporting.count:>12,}"
            f"{fmt(r.supporting.frequency):>11}{fmt(r.likelihood):>15}{fmt(r.standard_error):>11}"
            f"{fmt(r.prior, 3):>8}{fmt(r.posterior, 4):>11}")
    lines += ["", f"Posterior ranking: {summary_line(c)}"]

    if not details:
        return "\n".join(lines)
    exact = c.method == "exact"
    for r in c.ranked:
        lines += ["", f"── {r.mode.value} ──"]
        if exact and r.likelihood > 0:
            lines.append(f"Supporting pool: prior probability {fmt(r.supporting.probability)} "
                         "(computed exactly)")
        elif r.supporting.count:
            lines.append(f"Supporting pool: {r.supporting.count:,} simulations, "
                         f"{r.supporting.n_configurations} distinct configurations, "
                         f"effective sample {r.effective_sample:.0f}")
        if r.supporting.count or (exact and r.likelihood > 0):
            lines += _configurations(r.supporting, limit)
            lines.append("  Individual genotype probabilities (given the observations):")
            for ind, dist in r.supporting.marginal.items():
                part = ", ".join(f"{g}: {fmt(p, 3)}" for g, p in dist.items() if p > 0)
                lines.append(f"    {ind:<8} {part}")
        else:
            lines.append("Supporting pool is empty: this mode is incompatible with the observations.")
        if exact:
            lines.append(f"Non-supporting pool: prior probability "
                         f"{fmt(r.non_supporting.probability)}")
        elif r.non_supporting.count:
            lines.append(f"Non-supporting pool: {r.non_supporting.count:,} simulations "
                         f"(prior probability {fmt(r.non_supporting.probability)})")
            if r.conflict:
                part = ", ".join(f"{i}: {fmt(o, 3)}" for i, o in
                                 sorted(r.conflict.items(), key=lambda t: -t[1]) if o > 0)
                lines.append(f"  Individuals causing the conflict (share): {part}")
            lines += _configurations(r.non_supporting, limit)
            lines.append("  Individual genotype probabilities (in the non-supporting pool):")
            for ind, dist in r.non_supporting.marginal.items():
                part = ", ".join(f"{g}: {fmt(p, 3)}" for g, p in dist.items() if p > 0)
                lines.append(f"    {ind:<8} {part}")
    return "\n".join(lines)


def _configurations(pool: GenotypePool, limit: int | None) -> list[str]:
    o = []
    for cfg in pool.configurations[:limit]:
        g = " ".join(f"{i}={v}" for i, v in cfg.genotypes.items())
        o.append(f"  {fmt(cfg.probability, 3):>7}  ({cfg.count:>9,})  {g}")
    return o
