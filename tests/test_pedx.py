"""PEDX reading, validation, writing and drawing tests."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from pedkit import (
    Disease, Individual, Level, Mode, Pedigree, Phenotype, Sex, read, read_text, validate,
    write, write_fam, write_text,
)
from pedkit.vis import draw

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "cf_family.pedx"
SVG = "{http://www.w3.org/2000/svg}"

# line numbers in cf_family.pedx
FORMAT, DISEASE_D2, HEADER = 1, 4, 9
I1, I2, II1, II2, II3, III1 = 10, 11, 12, 13, 14, 15


@pytest.fixture
def text() -> str:
    return EXAMPLE.read_text(encoding="utf-8")


def corrupt(text: str, line: int, old: str, new: str) -> str:
    """Replace ``old`` by ``new`` on line ``line``."""
    lines = text.split("\n")
    assert old in lines[line - 1], (old, lines[line - 1])
    lines[line - 1] = lines[line - 1].replace(old, new, 1)
    return "\n".join(lines)


def set_cell(text: str, line: int, column: int, new: str) -> str:
    lines = text.split("\n")
    cells = lines[line - 1].split("\t")
    cells[column] = new
    lines[line - 1] = "\t".join(cells)
    return "\n".join(lines)


# -- 1. round trip ------------------------------------------------------------------

def test_read_then_write_is_byte_identical(tmp_path):
    out = tmp_path / "copy.pedx"
    write(read(EXAMPLE), out)
    assert out.read_bytes() == EXAMPLE.read_bytes()


def test_unknown_columns_and_metadata_are_kept():
    text = ("##format=PEDX/0.1\n"
            "##source=lecture notes, p. 12\n"
            "##disease=<ID=D1,name=\"A, B\",custom_field=xyz>\n"
            "#family\tindividual\tfather\tmother\tsex\tD1\tnote\n"
            "\n"
            "F\tA\t0\t0\t1\t2\ta note\n")
    ped = read_text(text)
    assert ped.diseases[0].name == "A, B"
    assert ped.diseases[0].extra == {"custom_field": "xyz"}
    assert ped.individuals[0].extra == {"note": "a note"}
    assert write_text(ped) == text


def test_a_change_touches_only_its_own_line(text):
    ped = read_text(text)
    ped.individual("II-2").age = 28
    ped.disease("D1").mode = Mode.AR
    new = write_text(ped).split("\n")
    old = text.split("\n")
    changed = [i + 1 for i, (a, b) in enumerate(zip(old, new)) if a != b]
    assert changed == [3, II2]
    assert new[II2 - 1].split("\t")[11] == "28"
    assert "mode=AR" in new[2] and "mutation_rate=0," in new[2]


def test_invalid_value_is_kept_on_write(text):
    bad = set_cell(text, II2, 4, "3")
    assert write_text(read_text(bad)) == bad


def test_pedigree_built_in_memory_round_trips():
    ped = Pedigree(diseases=[Disease("D1", name="Haemophilia A", mode=Mode.XR)])
    ped.individuals = [
        Individual("F", "father", sex=Sex.MALE, phenotypes={"D1": Phenotype.UNAFFECTED}),
        Individual("F", "mother", sex=Sex.FEMALE, phenotypes={"D1": Phenotype.UNAFFECTED}),
        Individual("F", "son", "father", "mother", Sex.MALE, {"D1": Phenotype.AFFECTED},
                   deceased=True),
    ]
    back = read_text(write_text(ped))
    assert validate(back) == []
    assert back.columns == ["family", "individual", "father", "mother", "sex", "D1", "deceased"]
    assert back.disease("D1").mode is Mode.XR
    assert back.individual("son").deceased and back.individual("son").affected("D1")


# -- 2. example file ----------------------------------------------------------------

def test_example_file_validates_without_findings():
    ped = read(EXAMPLE)
    assert validate(ped) == []
    assert len(ped.individuals) == 6
    assert [d.id for d in ped.diseases] == ["D1", "D2"]
    assert ped.individual("II-1").genotype("D1").mutant_count == 2
    assert ped.individual("II-3").genotype("D2").biallelic is False
    assert len(ped.questions) == 3 and ped.questions[1].event == "carrier"


# -- 3. every ERROR rule ------------------------------------------------------------

def _has_error(text: str, rule: str, line: int) -> None:
    findings = validate(read_text(text))
    match = [f for f in findings if f.rule == rule and f.line == line]
    assert match, f"expected {rule} on line {line}; found: {[str(f) for f in findings]}"
    assert match[0].level is Level.ERROR
    assert match[0].message


def test_e1_format_line(text):
    _has_error(corrupt(text, FORMAT, "##format=", "##fmt="), "E1", FORMAT)


def test_e1_no_header(text):
    _has_error(corrupt(text, HEADER, "#family", "family"), "E1", HEADER)


def test_e2_order_of_first_five_columns(text):
    _has_error(corrupt(text, HEADER, "father\tmother", "mother\tfather"), "E2", HEADER)


def test_e2_required_column_missing(text):
    _has_error(corrupt(text, HEADER, "\tsex", "\tgender_"), "E2", HEADER)


def test_e3_duplicate_individual(text):
    _has_error(set_cell(text, III1, 1, "II-3"), "E3", III1)


def test_e4_missing_parent(text):
    _has_error(set_cell(text, II1, 2, "I-9"), "E4", II1)


def test_e5_father_is_female(text):
    _has_error(set_cell(text, I1, 4, "2"), "E5", I1)


def test_e5_mother_is_male(text):
    _has_error(set_cell(text, I2, 4, "1"), "E5", I2)


def test_e6_loop(text):
    bad = set_cell(set_cell(text, I1, 2, "II-1"), I1, 3, "II-2")
    _has_error(bad, "E6", I1)


def test_e7_disease_without_column(text):
    _has_error(corrupt(text, DISEASE_D2, "ID=D2", "ID=D3"), "E7", DISEASE_D2)


def test_e7_column_without_disease(text):
    lines = text.split("\n")
    lines = [s + "\t1" if s.startswith("FAM01") else s for s in lines]
    lines[HEADER - 1] += "\tD9"
    _has_error("\n".join(lines), "E7", HEADER)


def test_e8_value_outside_allowed_set(text):
    _has_error(set_cell(text, II2, 4, "3"), "E8", II2)


def test_e8_metadata_value(text):
    _has_error(corrupt(text, DISEASE_D2, "mode=XR", "mode=XY"), "E8", DISEASE_D2)


def test_e8_autosomal_hemizygote(text):
    # With an unknown mode a hemizygous genotype in a male is allowed (it may
    # be X-linked); under an autosomal mode it is invalid.
    autosomal = corrupt(text, 3, "mode=unknown", "mode=AR")
    _has_error(set_cell(autosomal, II1, 7, "1"), "E8", II1)


def test_e9_x_linked_male_with_two_alleles(text):
    _has_error(set_cell(text, II3, 8, "0/1"), "E9", II3)


# -- warnings -------------------------------------------------------------------------

def _has_warning(text: str, rule: str, line: int) -> None:
    findings = validate(read_text(text))
    assert any(f.rule == rule and f.line == line and f.level is Level.WARNING
               for f in findings), [str(f) for f in findings]


def test_w3_question_about_missing_individual(text):
    _has_warning(corrupt(text, 7, "individual=II-2", "individual=II-9"), "W3", 7)


def test_w4_single_parent(text):
    _has_warning(set_cell(text, II1, 3, "0"), "W4", II1)


# -- 4. PLINK .fam --------------------------------------------------------------------

def test_fam_six_columns_and_phenotype(tmp_path):
    ped = read(EXAMPLE)
    for disease, expected in ((None, ["1", "1", "2", "1", "1", "0"]),
                              ("D2", ["1", "1", "1", "1", "2", "0"])):
        path = tmp_path / "family.fam"
        write_fam(ped, path, disease)
        lines = path.read_text().splitlines()
        assert len(lines) == 6
        assert all(len(s.split()) == 6 for s in lines)
        assert [s.split()[5] for s in lines] == expected
        assert "-9" not in path.read_text()
    assert lines[2].split() == ["FAM01", "II-1", "I-1", "I-2", "1", "1"]
    assert lines[5].split()[2:5] == ["II-3", "II-2", "0"]


# -- 5. SVG ---------------------------------------------------------------------------

def test_draw_valid_svg_with_six_symbols():
    svg = draw(read(EXAMPLE))
    root = ET.fromstring(svg.encode("utf-8"))
    assert root.tag == SVG + "svg"
    symbols = [g for g in root.iter(SVG + "g") if g.get("class") == "individual"]
    assert len(symbols) == 6
    assert "individual-FAM01-III-1" in {g.get("id") for g in symbols}


def test_svg_is_qt_tiny_compatible():
    """Elements QtSvg (SVG 1.2 Tiny) does not support must not be used."""
    svg = draw(read(EXAMPLE))
    for tag in ("<clipPath", "<pattern", "<marker", "url(#"):
        assert tag not in svg


def test_draw_twins_and_special_kinds():
    text = ("##format=PEDX/0.1\n"
            "##disease=<ID=D1,pattern=dotted>\n"
            "##disease=<ID=D2,pattern=quadrant>\n"
            "#family\tindividual\tfather\tmother\tsex\tD1\tD2\ttwin\tkind\tadopted\tdeceased\n"
            "F\tA\t0\t0\t1\t1\t1\t.\tbirth\t.\t0\n"
            "F\tB\t0\t0\t2\t1\t1\t.\tbirth\t.\t0\n"
            "F\tC\tA\tB\t1\t2\t2\tMZ:a\tbirth\t.\t0\n"
            "F\tD\tA\tB\t1\t2\t0\tMZ:a\tbirth\t.\t0\n"
            "F\tE\tA\tB\t0\t0\t0\t.\tmiscarriage\t.\t0\n"
            "F\tG\tA\tB\t2\t1\t1\t.\tstillbirth\t.\t1\n"
            "F\tK\tA\tB\t2\t1\t1\t.\tbirth\tin\t0\n")
    ped = read_text(text)
    assert [f for f in validate(ped) if f.level is Level.ERROR] == []
    svg = draw(ped)
    ET.fromstring(svg.encode("utf-8"))
    assert svg.count('class="individual"') == 7
    assert ">SB<" in svg and ">?<" in svg


# -- command line ---------------------------------------------------------------------

def test_cli_validate_draw_fam(tmp_path, text, capsys):
    from pedkit.cli import main
    assert main(["validate", str(EXAMPLE)]) == 0
    assert "0 errors, 0 warnings" in capsys.readouterr().out
    bad = tmp_path / "bad.pedx"
    bad.write_text(set_cell(text, II1, 2, "I-9"), encoding="utf-8")
    assert main(["validate", str(bad)]) == 1
    assert "ERROR E4 (line 12)" in capsys.readouterr().out
    assert main(["draw", str(EXAMPLE), "-o", str(tmp_path / "a.svg")]) == 0
    ET.parse(tmp_path / "a.svg")
    assert main(["fam", str(EXAMPLE), "-o", str(tmp_path / "a.fam"), "--disease", "D2"]) == 0
    assert (tmp_path / "a.fam").read_text().splitlines()[4].endswith("\t2")
