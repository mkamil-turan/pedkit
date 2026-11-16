"""Multifactorial risk plots (matplotlib).

Colours: Elston–Stewart blue, MCMC orange, population grey. The categorical
order is fixed: a method has the same colour in every plot.
"""

from __future__ import annotations

import dataclasses
import math
import os

import numpy as np

from .elston import extrapolated_risk, liability_density
from .normal import phi

ES = "#2a78d6"
MC = "#eb6834"
THIRD = "#1baf7a"
FOURTH = "#eda100"
POPULATION = "#8c8b86"
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
GRID = "#e6e5e0"
SURFACE = "#fcfcfb"


STYLE = {
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT, "axes.titlecolor": TEXT,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "xtick.color": TEXT2, "ytick.color": TEXT2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "legend.frameon": False, "lines.linewidth": 2.0,
}

PLOT_TITLES = {
    "liability": "Liability distribution",
    "heritability": "Heritability sensitivity",
    "prevalence": "Prevalence sensitivity",
    "relatives": "Risk of each individual",
    "mcmc": "MCMC convergence",
}


def figures(r) -> dict:
    """{name: matplotlib Figure}. pyplot is not used: the figures can be
    saved to files or embedded in a Qt canvas, and built in a worker thread."""
    try:
        import matplotlib
        from matplotlib.figure import Figure
    except ImportError as e:
        raise RuntimeError("The plots need matplotlib: pip install 'pedkit[mf]'") from e
    makers = {"liability": _liability, "heritability": _heritability,
              "prevalence": _prevalence, "relatives": _relatives, "mcmc": _mcmc}
    out = {}
    with matplotlib.rc_context(STYLE):
        for name, fn in makers.items():
            if name == "mcmc" and r.mcmc is None:
                continue
            if name in ("heritability", "prevalence", "relatives") and r.es is None:
                continue
            fig = fn(Figure, r)
            if fig.get_layout_engine() is None:
                fig.set_layout_engine("tight")    # keep labels inside a Qt canvas
            out[name] = fig
    return out


def save_plots(r, folder: str, fmt: str = "png") -> list[str]:
    os.makedirs(folder, exist_ok=True)
    paths = []
    for name, fig in figures(r).items():
        path = os.path.join(folder, f"{name}.{fmt}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        paths.append(path)
    return paths


def _target_thresholds(r) -> list[tuple[str, float]]:
    p, c = r.par, r.h.sex[r.target]
    if c in (1, 2) or not p.sex_specific:
        return [("threshold", p.individual_threshold(c))]
    return [("threshold (male)", p.threshold_male), ("threshold (female)", p.threshold_female)]


def _liability(Figure, r):
    fig = Figure(figsize=(7.5, 4.2))
    ax = fig.subplots()
    x = np.linspace(-4, 5.5, 800)
    thresholds = _target_thresholds(r)
    T = thresholds[0][1]
    ax.plot(x, phi(x), color=POPULATION, label="Population: N(0, 1)")
    ax.fill_between(x, 0, phi(x), where=x >= T, color=POPULATION, alpha=0.35, linewidth=0)
    if r.mcmc is not None:
        ax.hist(r.mcmc.target_samples, bins=120, range=(-4, 5.5), density=True,
                histtype="step", color=MC, linewidth=1.6, label=f"{r.individual}, MCMC samples")
    if r.es is not None:
        y = liability_density(r.es.posterior, r.par, x)
        ax.plot(x, y, color=ES, label=f"{r.individual}, Elston–Stewart")
        ax.fill_between(x, 0, y, where=x >= T, color=ES, alpha=0.25, linewidth=0)
    for name, t in thresholds:
        ax.axvline(t, color=TEXT2, linewidth=1.2, linestyle="--")
        ax.text(t, ax.get_ylim()[1] * 0.97, f" {name}\n T = {t:.2f}", va="top", ha="left",
                fontsize=9, color=TEXT2)
    ax.set_xlabel("Liability (in population standard deviations)")
    ax.set_ylabel("Density")
    ax.set_title(f"{r.individual} liability · area above threshold = risk {r.risk:.3g} "
                 f"(population {r.population_risk:.3g})")
    ax.legend(loc="upper left")
    ax.set_ylim(bottom=0)
    return fig


def _heritability(Figure, r):
    fig = Figure(figsize=(7.5, 4.2))
    ax = fig.subplots()
    h0 = without_target(r)
    grid = np.linspace(0.0, 0.95, 20)
    risks = [extrapolated_risk(h0, r.par.replace(heritability=float(g)), r.target, r.es.loci)
             for g in grid]
    ax.plot(grid, risks, color=ES, label="Elston–Stewart")
    ax.axhline(r.population_risk, color=POPULATION, linestyle="--", linewidth=1.5)
    ax.text(grid[-1], r.population_risk, "population risk ", va="bottom", ha="right",
            color=TEXT2, fontsize=9)
    ax.plot([r.par.heritability], [r.es.risk], "o", color=ES, markersize=8,
            markeredgecolor=SURFACE, markeredgewidth=2)
    if r.mcmc is not None:
        ax.errorbar([r.par.heritability], [r.mcmc.risk], yerr=[[2 * r.mcmc.standard_error]] * 2,
                    fmt="s", color=MC, markersize=7, capsize=4, label="MCMC (±2 SE)",
                    markeredgecolor=SURFACE, markeredgewidth=1.5)
    ax.annotate(f"h² = {r.par.heritability:.2f}\nrisk = {r.es.risk:.3g}",
                (r.par.heritability, r.es.risk), textcoords="offset points", xytext=(-70, 12),
                fontsize=9, color=TEXT)
    ax.set_xlabel("Heritability h²")
    ax.set_ylabel(f"Risk for {r.individual}")
    ax.set_title("Sensitivity of the risk to heritability (prevalence fixed)")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left")
    return fig


def _prevalence(Figure, r):
    fig = Figure(figsize=(7.5, 4.2))
    ax = fig.subplots()
    h0 = without_target(r)
    p = r.par
    factors = np.logspace(-2, math.log10(min(0.2 / max(p.prevalence_male, p.prevalence_female),
                                              30)), 18)
    K, risks = [], []
    for c in factors:
        new = p.replace(prevalence=p.prevalence * c, prevalence_male=p.prevalence_male * c,
                        prevalence_female=p.prevalence_female * c)
        K.append(r.population_risk * c)
        risks.append(extrapolated_risk(h0, new, r.target, r.es.loci))
    ax.loglog(K, risks, color=ES, label="Elston–Stewart")
    ax.loglog(K, K, color=POPULATION, linestyle="--", linewidth=1.5)
    ax.text(K[0], K[0], " population risk (risk = prevalence)", va="bottom", ha="left",
            color=TEXT2, fontsize=9)
    ax.plot([r.population_risk], [r.es.risk], "o", color=ES, markersize=8,
            markeredgecolor=SURFACE, markeredgewidth=2)
    if r.mcmc is not None:
        ax.errorbar([r.population_risk], [r.mcmc.risk], yerr=[[2 * r.mcmc.standard_error]] * 2,
                    fmt="s", color=MC, markersize=7, capsize=4, label="MCMC (±2 SE)",
                    markeredgecolor=SURFACE, markeredgewidth=1.5)
    ax.annotate(f"λ = {r.relative_risk:.1f}", (r.population_risk, r.es.risk),
                textcoords="offset points", xytext=(10, -14), fontsize=9, color=TEXT)
    ax.set_xlabel("Population prevalence (log)")
    ax.set_ylabel(f"Risk for {r.individual} (log)")
    ax.set_title("Sensitivity of the risk to prevalence (h² fixed)")
    ax.legend(loc="upper left")
    return fig


def _relatives(Figure, r):
    h = r.h
    label = {0: "unknown", 1: "unaffected", 2: "affected"}
    order = list(range(h.n))
    risks = [r.es.individual_risk[i] for i in order]
    fig = Figure(figsize=(7.5, 0.42 * h.n + 1.4))
    ax = fig.subplots()
    y = np.arange(len(order))[::-1]
    ax.barh(y, risks, color=ES, height=0.62)
    for yy, i, rk in zip(y, order, risks):
        ax.text(rk, yy, f"  {rk:.3g}", va="center", ha="left", fontsize=9, color=TEXT)
    ax.set_yticks(y, [f"{h.labels[i]} ({label[h.phenotype[i]]})" for i in order])
    # The prevalence in each individual's own sex: a grey vertical tick
    population = [_sex_prevalence(r.par, h.sex[i]) for i in order]
    for yy, t in zip(y, population):
        ax.plot([t, t], [yy - 0.42, yy + 0.42], color=TEXT, linewidth=2)
    ax.plot([], [], color=TEXT, linewidth=2, label="prevalence in own sex")
    ax.legend(loc="lower right")
    ax.set_xlabel("Risk of being affected from the relatives, own phenotype unknown")
    ax.set_title("Risk from the relatives (Elston–Stewart)")
    ax.set_xlim(0, max(risks + population) * 1.25)
    ax.grid(axis="y", visible=False)
    return fig


def _mcmc(Figure, r):
    m = r.mcmc
    fig = Figure(figsize=(7.5, 6.2))
    a1, a2 = fig.subplots(2, 1)
    step = np.arange(1, m.samples + 1)
    cumulative = np.cumsum(m.risk_trace, axis=0) / step[:, None]   # [step, chain]
    overall = cumulative.mean(axis=1)
    se = cumulative.std(axis=1, ddof=1) / math.sqrt(m.chains)
    a1.plot(step, overall, color=MC, label="MCMC cumulative mean")
    a1.fill_between(step, overall - 2 * se, overall + 2 * se, color=MC, alpha=0.2, linewidth=0,
                    label="±2 standard errors")
    if r.es is not None:
        a1.axhline(r.es.risk, color=ES, linestyle="--", linewidth=1.5, label="Elston–Stewart")
    a1.set_ylabel("Risk estimate")
    a1.set_xlabel("Step after burn-in")
    a1.set_title(f"MCMC convergence · R̂ {m.rhat_liability:.3f} (liability), "
                 f"{m.rhat_risk:.3f} (risk)")
    a1.legend(loc="upper right")
    colours = [MC, ES, THIRD]
    window = min(m.samples, 400)
    for z in range(min(3, m.chains)):
        a2.plot(step[:window], m.liability_trace[:window, z], color=colours[z], linewidth=1.0,
                label=f"chain {z + 1}")
    for k, (_, t) in enumerate(_target_thresholds(r)):
        a2.axhline(t, color=TEXT2, linestyle="--", linewidth=1.2,
                   label="threshold" if k == 0 else None)
    a2.set_ylabel(f"{r.individual} liability")
    a2.set_xlabel("Step after burn-in")
    a2.set_title(f"Trace plot · first {window} steps")
    a2.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=4, fontsize=8)
    fig.set_layout_engine("tight")
    return fig


def _sex_prevalence(p, sex: int) -> float:
    if sex == 1:
        return p.prevalence_male
    if sex == 2:
        return p.prevalence_female
    return (p.prevalence_male + p.prevalence_female) / 2


def without_target(r):
    """The preparation with the target's own phenotype removed."""
    return dataclasses.replace(r.h, phenotype=[0 if i == r.target else f
                                               for i, f in enumerate(r.h.phenotype)])
