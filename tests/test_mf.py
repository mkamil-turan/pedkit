"""Multifactorial (liability-threshold) risk tests."""

import math
from pathlib import Path

import numpy as np
import pytest

from pedkit.cli import main
from pedkit.io import read, read_text
from pedkit.mf import kinship_matrix, mf_risk
from pedkit.mf.model import prepare_mf
from pedkit.mf.normal import Phi_inv, threshold, upper_tail
from pedkit.risk.bayes import bayes_risk

EXAMPLE = Path(__file__).parent.parent / "examples" / "multifactorial.pedx"


def _sibs(K=0.01, h2=0.8, km=None, kf=None, target_sex=1):
    extra = f",prevalence_male={km},prevalence_female={kf}" if km else f",prevalence={K}"
    return read_text(f"""##format=PEDX/0.1
##disease=<ID=D1,mode=MF,heritability={h2}{extra}>
#family\tindividual\tfather\tmother\tsex\tD1
A\tB\t0\t0\t1\t0
A\tA\t0\t0\t2\t0
A\tC1\tB\tA\t2\t2
A\tC2\tB\tA\t{target_sex}\t0
""")


def sib_risk_exact(K, h2):
    """P(L2 > T | L1 > T) with correlation h²/2: a bivariate normal integral."""
    T = threshold(K)
    r = h2 / 2
    x = np.linspace(T, T + 12, 200_001)
    dens = np.exp(-x ** 2 / 2) / math.sqrt(2 * math.pi)
    cond = upper_tail((T - r * x) / math.sqrt(1 - r * r))
    return float(np.trapezoid(dens * cond, x)) / K


def test_threshold_and_inverse():
    assert threshold(0.01) == pytest.approx(2.326348, abs=1e-5)
    assert Phi_inv(0.975) == pytest.approx(1.959964, abs=1e-6)


def test_kinship():
    h = prepare_mf(_sibs())
    phi = kinship_matrix(h)
    i = {e.split("/")[-1]: k for k, e in enumerate(h.labels)}
    assert phi[i["C1"], i["C2"]] == pytest.approx(0.25)
    assert phi[i["B"], i["C1"]] == pytest.approx(0.25)
    assert phi[i["A"], i["B"]] == 0
    assert phi[i["C1"], i["C1"]] == pytest.approx(0.5)


@pytest.mark.parametrize("K,h2", [(0.01, 0.8), (0.05, 0.5)])
def test_sib_risk_matches_exact_integral(K, h2):
    exact = sib_risk_exact(K, h2)
    r = mf_risk(_sibs(K, h2), "C2", seed=3, chains=32, samples=2000, burn_in=300)
    assert r.es.risk == pytest.approx(exact, rel=0.01)
    assert abs(r.mcmc.risk - exact) < 5 * r.mcmc.standard_error + 1e-3
    assert r.mcmc.rhat_liability < 1.05


def test_zero_heritability_gives_population_risk():
    r = mf_risk(_sibs(0.02, 0.0), "C2", method="es")
    assert r.es.risk == pytest.approx(0.02, rel=1e-6)
    assert r.relative_risk == pytest.approx(1.0, rel=1e-6)


def test_own_phenotype_is_ignored():
    r1 = mf_risk(_sibs(), "C1", method="es")       # C1 affected, relatives unknown
    assert r1.es.risk == pytest.approx(0.01, rel=1e-3)   # HPM + Richardson approximation


def test_carter_effect():
    """An affected individual of the less often affected sex (female) passes more risk."""
    ped = read(EXAMPLE)
    son = mf_risk(ped, "III-2", method="es")
    daughter = mf_risk(ped, "III-3", method="es")
    assert son.es.risk > daughter.es.risk
    assert daughter.relative_risk > son.relative_risk
    # Compared with an affected father: disease in the mother indicates more liability
    head = ("##format=PEDX/0.1\n##disease=<ID=D1,mode=MF,heritability=0.8,prevalence_male=0.005,"
            "prevalence_female=0.001>\n#family\tindividual\tfather\tmother\tsex\tD1\n")
    father = mf_risk(read_text(head + "A\tB\t0\t0\t1\t2\nA\tA\t0\t0\t2\t0\nA\tC\tB\tA\t1\t0\n"),
                     "C", method="es")
    mother = mf_risk(read_text(head + "A\tB\t0\t0\t1\t0\nA\tA\t0\t0\t2\t2\nA\tC\tB\tA\t1\t0\n"),
                     "C", method="es")
    assert mother.es.risk > father.es.risk


def test_example_methods_agree():
    r = mf_risk(read(EXAMPLE), "III-2", seed=7, chains=32, samples=2000, burn_in=300)
    assert abs(r.mcmc.risk - r.es.risk) < 5 * r.mcmc.standard_error + 2e-3
    assert "Relative risk" in r.report()
    assert max(len(s) for s in r.report().splitlines()) <= 76


def test_parameter_errors():
    ped = read_text("##format=PEDX/0.1\n##disease=<ID=D1,mode=MF,prevalence=0.01>\n"
                    "#family\tindividual\tfather\tmother\tsex\tD1\nA\tX\t0\t0\t1\t2\n")
    with pytest.raises(ValueError, match="heritability"):
        mf_risk(ped, "X")
    assert mf_risk(ped, "X", heritability=0.5, method="es").es.risk == pytest.approx(0.01, rel=1e-3)


def test_bayes_points_to_mf():
    with pytest.raises(ValueError, match="pedkit mf"):
        bayes_risk(read(EXAMPLE), "III-2")


def test_plots(tmp_path):
    pytest.importorskip("matplotlib")
    r = mf_risk(read(EXAMPLE), "III-2", seed=1, chains=8, samples=300, burn_in=100)
    paths = r.plots(str(tmp_path), "svg")
    assert {Path(p).stem for p in paths} == {"liability", "heritability", "prevalence",
                                             "relatives", "mcmc"}
    assert all(Path(p).stat().st_size > 1000 for p in paths)


def test_cli(capsys):
    code = main(["mf", str(EXAMPLE), "--questions", "--method", "es"])
    out = capsys.readouterr().out
    assert code == 0
    assert "Question Q1" in out and "Question Q2" in out and "Elston–Stewart" in out


def test_viewer_mf_tab(monkeypatch):
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("matplotlib")
    QtWidgets = pytest.importorskip("PyQt5.QtWidgets")
    from PyQt5 import QtCore

    from pedkit.vis.qt import create_window
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    assert create_window(read(EXAMPLE.parent / "cf_family.pedx")).mf_panel is None
    w = create_window(read(EXAMPLE))
    m = w.mf_panel
    assert m.individual_box.currentText() == "III-2"
    m.method_box.setCurrentIndex(1)            # Elston–Stewart only: fast
    m.compute()
    m.worker.wait(60_000)
    for _ in range(10):
        app.processEvents()
        QtCore.QThread.msleep(10)
    assert m.plot_tabs.count() == 4            # no MCMC plot
    assert "Elston–Stewart" in m.report.toPlainText()
