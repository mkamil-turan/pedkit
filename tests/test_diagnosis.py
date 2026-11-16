"""Hierarchical diagnosis, codominance (phenotype 3) and mitochondrial inheritance."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from pedkit import Mode, Phenotype, fam_text, read, read_text, write_text
from pedkit.inference import infer_mode, infer_mode_exact
from pedkit.vis import draw

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def _path(c) -> list[str]:
    return [d[0].name for d in c.diagnosis.path]


@pytest.mark.parametrize("file, path, diagnosis", [
    ("cf_family.pedx", ["Autosomal", "autosomal recessive"], "autosomal recessive"),
    ("codominant.pedx", ["Autosomal", "autosomal codominant"], "autosomal codominant"),
    ("mitochondrial.pedx", ["Mitochondrial (maternal)"], "mitochondrial (maternal)"),
])
def test_diagnosis_of_the_example_families(file, path, diagnosis):
    c = infer_mode_exact(read(EXAMPLES / file))
    assert _path(c) == path
    assert c.diagnosis.diagnosis == diagnosis


def _untested():
    text = (EXAMPLES / "cf_family.pedx").read_text(encoding="utf-8").replace("\t1/1\t", "\t.\t")
    return read_text(text)


def test_x_linked_recessive_without_genotype_test():
    c = infer_mode_exact(_untested())
    assert _path(c) == ["Gonosomal", "X-linked", "X-linked recessive"]


def test_node_probabilities_add_up():
    c = infer_mode_exact(_untested())
    root = c.diagnosis.root
    assert root.probability == pytest.approx(1.0)
    assert sum(a.probability for a in root.children) == pytest.approx(1.0)
    for d in root.children:
        if d.children:
            assert sum(a.probability for a in d.children) == pytest.approx(d.probability)
            assert sum(a.conditional for a in d.children) == pytest.approx(
                1.0 if d.probability else 0.0)
    gonosomal = next(d for d in root.children if d.name == "Gonosomal")
    assert gonosomal.probability == pytest.approx(c.results[Mode.XR].posterior
                                                  + c.results[Mode.XD].posterior
                                                  + c.results[Mode.Y].posterior)


def test_no_diagnosis_when_no_mode_fits():
    # The intermediate phenotype is possible only under codominance; without ACD nothing fits
    c = infer_mode_exact(read(EXAMPLES / "codominant.pedx"), modes=[Mode.AD, Mode.AR, Mode.XR])
    assert c.diagnosis.path == []
    assert "no diagnosis" in c.diagnosis.diagnosis


def test_codominant_monte_carlo_agrees_with_exact():
    ped = read(EXAMPLES / "codominant.pedx")
    k = infer_mode_exact(ped, modes=[Mode.ACD])
    mc = infer_mode(ped, modes=[Mode.ACD], n=300_000, seed=2, proposal_freq=0.3)
    s = mc.results[Mode.ACD]
    assert abs(s.likelihood - k.results[Mode.ACD].likelihood) <= 4 * s.standard_error
    # Individuals with the intermediate phenotype are heterozygous, II-1 (full phenotype) homozygous
    assert s.supporting.marginal["I-1"] == {"0/1": pytest.approx(1.0)}
    assert s.supporting.marginal["II-1"] == {"1/1": pytest.approx(1.0)}
    assert s.supporting.marginal["II-3"] == {"0/0": pytest.approx(1.0)}


def test_mitochondrial_father_does_not_transmit():
    ped = read(EXAMPLES / "mitochondrial.pedx")
    k = infer_mode_exact(ped, modes=[Mode.MT])
    marg = k.results[Mode.MT].supporting.marginal
    assert marg["III-1"] == {"0": pytest.approx(1.0)}  # not from the affected father (II-1)
    assert marg["III-3"] == {"1": pytest.approx(1.0)}  # from the affected mother (II-3)


def test_intermediate_phenotype_is_read_written_and_drawn():
    text = (EXAMPLES / "codominant.pedx").read_text(encoding="utf-8")
    ped = read_text(text)
    assert ped.individual("I-1").phenotype("D1") is Phenotype.INTERMEDIATE
    assert write_text(ped) == text
    svg = draw(ped)
    ET.fromstring(svg.encode("utf-8"))
    assert "Intermediate phenotype" in svg
    fam = [s.split("\t") for s in fam_text(ped).splitlines()]
    assert [f[5] for f in fam] == ["2", "2", "2", "2", "1", "2"]


def test_report_starts_with_the_diagnosis():
    text = infer_mode(read(EXAMPLES / "cf_family.pedx"), n=10_000, seed=1).report()
    assert text.startswith("Hierarchical diagnosis")
    assert text.index("Diagnosis:") < text.index("Mode-of-inheritance inference")
