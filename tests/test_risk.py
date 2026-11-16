"""Bayes risk tests: textbook examples, mode averaging, Monte Carlo."""

from __future__ import annotations

from pathlib import Path

import pytest

from pedkit import Mode, read, read_text
from pedkit.risk import answer_questions, bayes_risk

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


@pytest.fixture(scope="module")
def dmd():
    return read(EXAMPLES / "dmd_bayes.pedx")


@pytest.fixture(scope="module")
def cf():
    return read(EXAMPLES / "cf_family.pedx")


def _row(r, mode, genotype):
    k = next(k for k in r.modes if k.mode is mode)
    return next(s for s in k.rows if s.genotype == genotype)


def test_dmd_classic_bayes_table(dmd):
    """Affected brother: prior 1/2; two healthy sons: conditional 1/4; posterior 1/5."""
    r = bayes_risk(dmd, "II-2", "carrier")
    assert r.method == "exact"
    assert [k.mode for k in r.modes] == [Mode.XR]  # mode=XR in the file
    carrier = _row(r, Mode.XR, "0/1")
    normal = _row(r, Mode.XR, "0/0")
    assert carrier.prior == pytest.approx(1 / 2)
    assert carrier.conditional / normal.conditional == pytest.approx(1 / 4)
    assert r.risk == pytest.approx(1 / 5)
    assert sorted(r.lower_individuals) == ["II-3", "III-1", "III-2", "III-3"]


def test_dmd_pregnancy_affection_risk(dmd):
    # 1/5 carrier × 1/2 male × 1/2 mutant X
    assert bayes_risk(dmd, "III-3", "affected").risk == pytest.approx(1 / 20)


def test_cf_family_risk_questions(cf):
    assert bayes_risk(cf, "II-2", "carrier").risk == pytest.approx(2 / 3)
    assert bayes_risk(cf, "III-1", "affected").risk == pytest.approx(1 / 153)


def test_observed_phenotype_is_taken_as_given(cf):
    assert bayes_risk(cf, "II-1", "affected").risk == pytest.approx(1.0)
    assert bayes_risk(cf, "II-2", "affected").risk == pytest.approx(0.0)


def test_mode_averaging():
    text = (EXAMPLES / "cf_family.pedx").read_text(encoding="utf-8").replace("\t1/1\t", "\t.\t")
    r = bayes_risk(read_text(text), "II-2", "carrier")
    weights = {k.mode: k.weight for k in r.modes}
    assert sum(weights.values()) == pytest.approx(1.0)
    assert weights[Mode.XR] > 0.9
    expected = sum(k.weight * k.risk for k in r.modes if k.weight)
    assert r.risk == pytest.approx(expected)
    # Under XR only: the mother is an obligate carrier, II-2 is 1/2
    xr = bayes_risk(read_text(text), "II-2", "carrier", modes=Mode.XR)
    assert xr.risk == pytest.approx(0.5, abs=1e-3)


def test_monte_carlo_agrees_with_exact(dmd):
    mc = bayes_risk(dmd, "II-2", "carrier", method="monte-carlo", n=300_000, seed=4,
                    proposal_freq=0.2)
    assert mc.method == "monte-carlo"
    assert mc.modes[0].supporting > 1000
    assert mc.risk == pytest.approx(1 / 5, abs=0.02)


def test_questions_are_answered(cf, dmd):
    answers = answer_questions(cf)
    assert [a.question.id for a in answers] == ["Q1", "Q2", "Q3"]
    assert all(a.error is None for a in answers)
    assert answers[0].result.diagnosis.diagnosis == "autosomal recessive"
    assert answers[1].result.risk == pytest.approx(2 / 3)
    assert "1/5" in answer_questions(dmd)[0].report()


def test_invalid_inputs(cf):
    with pytest.raises(ValueError, match="Unknown event"):
        bayes_risk(cf, "II-2", "sick")
    with pytest.raises(KeyError):
        bayes_risk(cf, "IV-9")
