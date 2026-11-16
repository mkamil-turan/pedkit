"""Monte Carlo mode inference tests; comparison with the exact computation."""

from __future__ import annotations

from pathlib import Path

import pytest

from pedkit import Mode, read, read_text
from pedkit.inference import infer_mode, infer_mode_exact

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "cf_family.pedx"


@pytest.fixture(scope="module")
def ped():
    return read(EXAMPLE)


@pytest.fixture(scope="module")
def untested():
    """The example with II-1's genotype test removed (G:D1 = .)."""
    text = EXAMPLE.read_text(encoding="utf-8").replace("\t1/1\t", "\t.\t")
    return read_text(text)


def _close(mc, exact, mode, sigma=4.0):
    m, k = mc.results[mode], exact.results[mode]
    assert abs(m.likelihood - k.likelihood) <= sigma * m.standard_error + 1e-12, \
        (mode, m.likelihood, m.standard_error, k.likelihood)


def test_exact_classic_carrier_probabilities(ped):
    k = infer_mode_exact(ped, "D1")
    assert k.most_likely.mode is Mode.AR
    assert k.results[Mode.AR].posterior == pytest.approx(1.0)
    marg = k.results[Mode.AR].supporting.marginal
    assert marg["II-2"]["0/1"] == pytest.approx(2 / 3)       # healthy sib of an affected
    q = 0.02
    assert marg["II-3"]["0/1"] == pytest.approx(2 * q / (1 + q))  # healthy, married in


def test_biallelic_test_in_a_male_excludes_x_modes(ped):
    k = infer_mode_exact(ped, "D1")
    for mode in (Mode.XD, Mode.XR, Mode.Y, Mode.MT, Mode.AD):
        assert k.results[mode].likelihood == 0


def test_monte_carlo_agrees_with_exact(ped):
    mc = infer_mode(ped, "D1", n=400_000, seed=7)
    k = infer_mode_exact(ped, "D1")
    _close(mc, k, Mode.AR)
    s = mc.results[Mode.AR]
    assert s.supporting.count + s.non_supporting.count == 400_000
    assert s.supporting.probability + s.non_supporting.probability == pytest.approx(1.0)
    assert s.supporting.frequency == pytest.approx(s.supporting.count / 400_000)
    assert mc.most_likely.mode is Mode.AR
    # Almost every conflict comes from II-1 (drawing 1/1 is rare); the rest
    # arise when II-1 fits but another individual (e.g. II-2 also 1/1) conflicts.
    assert s.conflict["II-1"] > 0.999
    assert max(s.conflict, key=s.conflict.get) == "II-1"


def test_without_test_xr_beats_ar_and_agrees_with_exact(untested):
    k = infer_mode_exact(untested, "D1")
    assert k.most_likely.mode is Mode.XR  # an X from the mother is likelier than two carrier parents
    mc = infer_mode(untested, "D1", n=500_000, seed=3)
    for mode in (Mode.AR, Mode.XR, Mode.MT):
        _close(mc, k, mode)
    assert mc.most_likely.mode is Mode.XR
    assert mc.results[Mode.XR].posterior == pytest.approx(k.results[Mode.XR].posterior, abs=0.03)


def test_importance_sampling_is_unbiased_and_grows_the_pool(untested):
    k = infer_mode_exact(untested, "D1", modes=[Mode.AR])
    plain = infer_mode(untested, "D1", modes=[Mode.AR], n=200_000, seed=5)
    imp = infer_mode(untested, "D1", modes=[Mode.AR], n=200_000, seed=5, proposal_freq=0.25)
    _close(imp, k, Mode.AR)
    assert imp.results[Mode.AR].supporting.count > 20 * plain.results[Mode.AR].supporting.count
    assert imp.results[Mode.AR].standard_error < plain.results[Mode.AR].standard_error
    assert imp.results[Mode.AR].supporting.marginal["II-2"]["0/1"] == \
        pytest.approx(k.results[Mode.AR].supporting.marginal["II-2"]["0/1"], abs=0.02)


def test_new_mutation_and_reduced_penetrance_open_ad():
    text = EXAMPLE.read_text(encoding="utf-8").replace("\t1/1\t", "\t.\t").replace(
        "penetrance=1.0,mutation_rate=0", "penetrance=0.8,mutation_rate=0.001")
    ped2 = read_text(text)
    k = infer_mode_exact(ped2, "D1")
    assert k.results[Mode.AD].likelihood > 0
    mc = infer_mode(ped2, "D1", n=400_000, seed=11, proposal_freq=0.2)
    for mode in (Mode.AD, Mode.AR, Mode.XR):
        _close(mc, k, mode)


def test_invalid_pedigree_is_rejected():
    text = EXAMPLE.read_text(encoding="utf-8").replace("FAM01\tII-1\tI-1", "FAM01\tII-1\tI-9")
    with pytest.raises(ValueError, match="E4"):
        infer_mode(read_text(text), n=10)


def test_seed_is_reproducible(ped):
    a = infer_mode(ped, "D1", n=50_000, seed=42, modes=[Mode.AR])
    b = infer_mode(ped, "D1", n=50_000, seed=42, modes=[Mode.AR])
    assert a.results[Mode.AR].supporting.count == b.results[Mode.AR].supporting.count


def test_report_is_produced(ped):
    text = infer_mode(ped, "D1", n=20_000, seed=1).report()
    assert "Supporting" in text and "Posterior" in text and "Posterior ranking:" in text
