"""Validation against independent computations.

Each test compares pedkit with a result obtained by a different route:
textbook first-principles answers, a second exact algorithm, unbiased
simulation, or direct numerical integration of the liability model.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from pedkit import Level, Mode, read, read_text, validate
from pedkit.inference import (
    MENDELIAN_MODES, enumeration_likelihood, infer_mode, peeling_likelihood, prepare,
)
from pedkit.mf import kinship_matrix, mf_risk
from pedkit.mf.model import prepare_mf
from pedkit.mf.normal import threshold, upper_tail
from pedkit.risk import bayes_risk

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
HEAD = "#family\tindividual\tfather\tmother\tsex\tD1"


def pedx(disease: str, rows: str, extra_cols: str = "") -> str:
    return f"##format=PEDX/0.1\n##disease=<ID=D1,{disease}>\n{HEAD}{extra_cols}\n{rows}"


# -- 1. textbook answers derived by hand -------------------------------------------------

def test_cf_healthy_sib_carrier_two_thirds():
    """Healthy sib of an AR-affected child of carrier parents: 2/3."""
    ped = read(EXAMPLES / "cf_family.pedx")
    r = bayes_risk(ped, "II-2", "carrier")
    assert Fraction(r.risk).limit_denominator(100) == Fraction(2, 3)


def test_cf_child_risk_exact_fraction():
    """III-1: mother 2/3 carrier × father 2q/(1+q) carrier × 1/4, q = 0.02.
    2/3 · (0.04/1.02) · 1/4 = 1/153."""
    q = 0.02
    expected = Fraction(2, 3) * Fraction(2 * q).limit_denominator() / \
        Fraction(1 + q).limit_denominator() * Fraction(1, 4)
    ped = read(EXAMPLES / "cf_family.pedx")
    assert bayes_risk(ped, "III-1", "affected").risk == pytest.approx(float(expected), rel=1e-12)
    assert expected == Fraction(1, 153)


def test_dmd_one_fifth_and_one_twentieth():
    """Bayes table by hand: prior 1/2 vs 1/2, conditional (1/2)² vs 1 → 1/5;
    affected son of the pregnancy: 1/5 · 1/2 · 1/2 = 1/20 (new mutation ignored)."""
    ped = read(EXAMPLES / "dmd_bayes.pedx")
    assert bayes_risk(ped, "II-2", "carrier").risk == pytest.approx(1 / 5, abs=1e-4)
    assert bayes_risk(ped, "III-3", "affected").risk == pytest.approx(1 / 20, abs=1e-4)


def test_ad_child_of_affected_heterozygote_one_half():
    ped = read_text(pedx("mode=AD,allele_freq=0.0001",
                         "F\tP\t0\t0\t1\t2\nF\tM\t0\t0\t2\t1\nF\tC\tP\tM\t0\t0\n"))
    assert bayes_risk(ped, "C", "affected").risk == pytest.approx(0.5, abs=1e-3)


def test_reduced_penetrance_bayes():
    """AD, penetrance f: an unaffected child of an affected heterozygote is a
    carrier with probability (1/2)(1-f) / ((1/2)(1-f) + 1/2) = (1-f)/(2-f)."""
    f = 0.7
    ped = read_text(pedx(f"mode=AD,allele_freq=1e-6,penetrance={f}",
                         "F\tP\t0\t0\t1\t2\nF\tM\t0\t0\t2\t1\nF\tC\tP\tM\t2\t1\n"))
    assert bayes_risk(ped, "C", "carrier").risk == pytest.approx((1 - f) / (2 - f), abs=1e-4)


# -- 2. two exact algorithms and Monte Carlo agree ---------------------------------------------

@pytest.mark.parametrize("file", ["cf_family.pedx", "dmd_bayes.pedx", "codominant.pedx",
                                  "mitochondrial.pedx"])
def test_peeling_equals_enumeration(file):
    h = prepare(read(EXAMPLES / file))
    for mode in MENDELIAN_MODES:
        a = peeling_likelihood(h, mode)
        b = enumeration_likelihood(h, mode)
        assert a.likelihood == pytest.approx(b.likelihood, rel=1e-10, abs=1e-300)
        assert a.supporting.probability == pytest.approx(b.supporting.probability, rel=1e-10,
                                                         abs=1e-15)
        for ind, dist in b.supporting.marginal.items():
            for g, p in dist.items():
                assert a.supporting.marginal[ind].get(g, 0) == pytest.approx(p, abs=1e-10)


@pytest.mark.parametrize("individual, event", [("II-2", "carrier"), ("III-1", "affected"),
                                                ("II-3", "mutant")])
def test_bayes_exact_equals_enumeration(individual, event):
    text = (EXAMPLES / "cf_family.pedx").read_text(encoding="utf-8").replace("\t1/1\t", "\t.\t")
    ped = read_text(text)
    a = bayes_risk(ped, individual, event, method="exact")
    b = bayes_risk(ped, individual, event, method="enumeration")
    assert a.risk == pytest.approx(b.risk, abs=1e-12)
    for ma, mb in zip(a.modes, b.modes):
        assert ma.mode is mb.mode
        assert ma.weight == pytest.approx(mb.weight, abs=1e-12)


def test_monte_carlo_likelihood_within_four_standard_errors():
    ped = read(EXAMPLES / "mitochondrial.pedx")
    h = prepare(ped)
    mc = infer_mode(ped, n=400_000, seed=13, proposal_freq=0.2)
    for mode in MENDELIAN_MODES:
        exact = peeling_likelihood(h, mode).likelihood
        r = mc.results[mode]
        assert abs(r.likelihood - exact) <= 4 * r.standard_error + 1e-15, mode


def _descent_line(generations: int) -> str:
    """AD, an affected founder and a line of descendants of unknown phenotype,
    each with an unrelated spouse of unknown phenotype."""
    rows = ["F\tG0\t0\t0\t1\t2", "F\tS0\t0\t0\t2\t0"]
    father, mother = "G0", "S0"
    for g in range(1, generations + 1):
        rows.append(f"F\tG{g}\t{father}\t{mother}\t1\t0")
        if g < generations:
            rows.append(f"F\tS{g}\t0\t0\t2\t0")
        father, mother = f"G{g}", f"S{g}"
    return pedx("mode=AD,allele_freq=1e-9", "\n".join(rows) + "\n")


def test_peeling_has_no_size_limit_descent_line():
    """The carrier probability halves each generation: (1/2)^g as q -> 0."""
    g = 12
    ped = read_text(_descent_line(g))
    assert len(ped.individuals) > 2 * 12          # far above MAX_ENUMERATION
    r = bayes_risk(ped, f"G{g}", "carrier", method="exact")
    assert r.risk == pytest.approx(0.5 ** g, rel=1e-5)


def test_unrelated_family_does_not_change_the_risk():
    """Adding an unrelated 12-person family leaves III-1 at exactly 1/153."""
    text = (EXAMPLES / "cf_family.pedx").read_text(encoding="utf-8")
    cells = "\t.\t.\t.\t0\t.\tbirth\t."
    extra = ["FAM09\tA0\t0\t0\t1\t1\t1" + cells, "FAM09\tB0\t0\t0\t2\t2\t1" + cells]
    for k in range(1, 11):
        extra.append(f"FAM09\tK{k}\tA0\tB0\t{1 + k % 2}\t{1 + (k % 3 == 0)}\t1" + cells)
    ped = read_text(text + "\n".join(extra) + "\n")
    assert not [f for f in validate(ped) if f.level is Level.ERROR]
    assert len(ped.individuals) == 18
    r = bayes_risk(ped, "FAM01/III-1", "affected", modes=Mode.AR)
    assert r.risk == pytest.approx(1 / 153, rel=1e-10)


# -- 3. monozygotic twins share one genotype ---------------------------------------------------

def test_mz_twin_of_affected_is_affected_under_full_penetrance():
    ped = read_text(pedx("mode=AD,allele_freq=0.001",
                         "F\tP\t0\t0\t1\t1\t.\nF\tM\t0\t0\t2\t2\t.\n"
                         "F\tA\tP\tM\t1\t2\tMZ:t\nF\tB\tP\tM\t1\t0\tMZ:t\n", "\ttwin"))
    assert bayes_risk(ped, "B", "affected").risk == pytest.approx(1.0)
    dz = read_text(pedx("mode=AD,allele_freq=0.001",
                        "F\tP\t0\t0\t1\t1\t.\nF\tM\t0\t0\t2\t2\t.\n"
                        "F\tA\tP\tM\t1\t2\tDZ:t\nF\tB\tP\tM\t1\t0\tDZ:t\n", "\ttwin"))
    assert bayes_risk(dz, "B", "affected").risk == pytest.approx(0.5, abs=1e-3)


def test_mz_twins_discordant_genotypes_are_impossible():
    ped = read_text(pedx("mode=AR,allele_freq=0.01",
                         "F\tP\t0\t0\t1\t1\t.\t.\nF\tM\t0\t0\t2\t1\t.\t.\n"
                         "F\tA\tP\tM\t1\t1\tMZ:t\t0/0\nF\tB\tP\tM\t1\t2\tMZ:t\t.\n",
                         "\ttwin\tG:D1"))
    h = prepare(ped)
    assert peeling_likelihood(h, Mode.AR).likelihood == 0
    assert enumeration_likelihood(h, Mode.AR).likelihood == 0


def test_mz_twin_kinship_and_mf_risk():
    """MZ twins share A, so their liabilities correlate by h²."""
    K, h2 = 0.01, 0.6
    ped = read_text(pedx(f"mode=MF,prevalence={K},heritability={h2}",
                         "F\tP\t0\t0\t1\t0\t.\nF\tM\t0\t0\t2\t0\t.\n"
                         "F\tA\tP\tM\t2\t2\tMZ:t\nF\tB\tP\tM\t2\t0\tMZ:t\n", "\ttwin"))
    h = prepare_mf(ped)
    phi = kinship_matrix(h)
    assert phi[2, 3] == pytest.approx(0.5)
    r = mf_risk(ped, "B", seed=5, chains=32, samples=2000, burn_in=300)
    exact = _bivariate_conditional(K, h2)
    assert r.es.risk == pytest.approx(exact, rel=0.01)
    assert abs(r.mcmc.risk - exact) < 5 * r.mcmc.standard_error + 1e-3


# -- 4. validation rule W2 -------------------------------------------------------------------

def test_w2_genotype_contradicts_mode():
    ped = read_text(pedx("mode=AR,allele_freq=0.01",
                         "F\tP\t0\t0\t1\t1\t.\nF\tM\t0\t0\t2\t1\t.\nF\tC\tP\tM\t1\t1\t1/1\n",
                         "\tG:D1"))
    w2 = [f for f in validate(ped) if f.rule == "W2"]
    assert len(w2) == 1 and w2[0].level is Level.WARNING and w2[0].individual == "C"
    reduced = read_text(pedx("mode=AR,allele_freq=0.01,penetrance=0.9",
                             "F\tP\t0\t0\t1\t1\t.\nF\tM\t0\t0\t2\t1\t.\nF\tC\tP\tM\t1\t1\t1/1\n",
                             "\tG:D1"))
    assert not [f for f in validate(reduced) if f.rule == "W2"]


# -- 5. the liability-threshold model against direct computation --------------------------------

def _bivariate_conditional(K: float, r: float) -> float:
    """P(L2 > T | L1 > T) for a standard bivariate normal with correlation r."""
    T = threshold(K)
    x = np.linspace(T, T + 12, 200_001)
    dens = np.exp(-x ** 2 / 2) / math.sqrt(2 * math.pi)
    cond = upper_tail((T - r * x) / math.sqrt(1 - r * r))
    return float(np.trapezoid(dens * cond, x)) / K


@pytest.mark.parametrize("relation, r_factor", [("parent", 0.5), ("sib", 0.5)])
def test_first_degree_relative_risk_bivariate_integral(relation, r_factor):
    K, h2 = 0.02, 0.7
    if relation == "parent":
        rows = "F\tP\t0\t0\t1\t2\nF\tM\t0\t0\t2\t0\nF\tC\tP\tM\t1\t0\n"
    else:
        rows = "F\tP\t0\t0\t1\t0\nF\tM\t0\t0\t2\t0\nF\tS\tP\tM\t2\t2\nF\tC\tP\tM\t1\t0\n"  # parents unknown
    ped = read_text(pedx(f"mode=MF,prevalence={K},heritability={h2}", rows))
    r = mf_risk(ped, "C", method="es")
    assert r.es.risk == pytest.approx(_bivariate_conditional(K, r_factor * h2), rel=0.01)


def test_mf_against_direct_multivariate_normal_simulation():
    """Draw liabilities from N(0, Σ) directly, keep the draws that match the
    phenotypes, and count the target above the threshold."""
    K, h2 = 0.05, 0.6
    ped = read(EXAMPLES / "multifactorial.pedx")
    text = (EXAMPLES / "multifactorial.pedx").read_text(encoding="utf-8")
    text = text.replace("prevalence_male=0.005,prevalence_female=0.001,heritability=0.8",
                        f"prevalence={K},heritability={h2}")
    ped = read_text(text)
    h = prepare_mf(ped)
    target = h.labels.index("III-2")
    sigma = h2 * 2 * kinship_matrix(h) + (1 - h2) * np.eye(h.n)
    rng = np.random.default_rng(11)
    T = threshold(K)
    chol = np.linalg.cholesky(sigma)
    hits = kept = 0
    for _ in range(20):
        L = rng.standard_normal((500_000, h.n)) @ chol.T
        ok = np.ones(len(L), dtype=bool)
        for i, f in enumerate(h.phenotype):
            if i == target or f == 0:
                continue
            ok &= (L[:, i] > T) if f == 2 else (L[:, i] <= T)
        kept += int(ok.sum())
        hits += int((L[ok, target] > T).sum())
    direct = hits / kept
    se = math.sqrt(direct * (1 - direct) / kept)
    r = mf_risk(ped, "III-2", seed=9, chains=32, samples=2000, burn_in=300)
    assert abs(r.es.risk - direct) < 4 * se + 0.002
    assert abs(r.mcmc.risk - direct) < 4 * math.hypot(se, r.mcmc.standard_error) + 0.002
