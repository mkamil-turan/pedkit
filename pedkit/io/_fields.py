"""PEDX field definitions: how each field is parsed from text and how it is
written. The reader and the writer share this table, so the format rules
live in one place."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Callable

from ..model import (
    Adoption, Genotype, Individual, Kind, Mode, Pattern, Phenotype, Role, Sex, Twin,
    TwinType, UnionStatus,
)

MISSING = "."
FIXED_COLUMNS = ["family", "individual", "father", "mother", "sex"]
GENOTYPE_PREFIX = "G:"


@dataclass(frozen=True)
class Field:
    parse: Callable[[str], Any]     # raises ValueError on an invalid value
    format: Callable[[Any], str]
    allowed: str                    # value set shown in error messages


# -- helper parsers --------------------------------------------------------

def _choice(enum_cls, missing: str | None = None, default=None) -> Field:
    values = [e.value for e in enum_cls]

    def parse(t: str):
        if missing is not None and t == missing:
            return default
        return enum_cls(t)

    def fmt(v) -> str:
        return missing if v is None else str(v.value)

    allowed = ", ".join(str(v) for v in values) + (f" or {missing}" if missing else "")
    return Field(parse, fmt, allowed)


def _int_choice(enum_cls) -> Field:
    def parse(t: str):
        if t not in {str(e.value) for e in enum_cls}:
            raise ValueError(t)
        return enum_cls(int(t))
    return Field(parse, lambda v: str(int(v)), ", ".join(str(e.value) for e in enum_cls))


def _parent() -> Field:
    def parse(t: str):
        if not t:
            raise ValueError(t)
        return None if t == "0" else t
    return Field(parse, lambda v: "0" if v is None else v, "an individual ID or 0")


def _text(missing: str | None = MISSING) -> Field:
    def parse(t: str):
        if missing is not None and t == missing:
            return None
        if not t:
            raise ValueError(t)
        return t
    return Field(parse, lambda v: missing if v is None else v, "non-empty text")


def _natural() -> Field:
    def parse(t: str):
        if t == MISSING:
            return None
        if not re.fullmatch(r"\d+", t):
            raise ValueError(t)
        return int(t)
    return Field(parse, lambda v: MISSING if v is None else str(v),
                 "a non-negative integer or .")


def _deceased() -> Field:
    def parse(t: str):
        if t not in ("0", "1"):
            raise ValueError(t)
        return t == "1"
    return Field(parse, lambda v: "1" if v else "0", "0, 1")


def _twin() -> Field:
    def parse(t: str):
        if t == MISSING:
            return None
        m = re.fullmatch(r"(MZ|DZ|\?):(\w+)", t)
        if not m:
            raise ValueError(t)
        return Twin(TwinType(m.group(1)), m.group(2))
    return Field(parse, lambda v: MISSING if v is None else str(v), "MZ:a, DZ:a, ?:a or .")


def _genotype() -> Field:
    def parse(t: str):
        return None if t == MISSING else Genotype.parse(t)
    return Field(parse, lambda v: MISSING if v is None else str(v),
                 "0/0, 0/1, 1/1, (hemizygous) 0, 1 or .")


def _probability() -> Field:
    def parse(t: str):
        v = float(t)
        if not 0.0 <= v <= 1.0:
            raise ValueError(t)
        return v
    return Field(parse, _format_number, "a number between 0 and 1")


def _format_number(v: float | None) -> str:
    if v is None:
        return ""
    return str(int(v)) if float(v).is_integer() else repr(float(v))


# -- individual columns ----------------------------------------------------

PHENOTYPE_FIELD = _int_choice(Phenotype)
GENOTYPE_FIELD = _genotype()

INDIVIDUAL_FIELDS: dict[str, tuple[str, Field]] = {
    # column name: (Individual attribute, field)
    "family": ("family", _text(missing=None)),
    "individual": ("id", _text(missing=None)),
    "father": ("father", _parent()),
    "mother": ("mother", _parent()),
    "sex": ("sex", _int_choice(Sex)),
    "role": ("role", _choice(Role, MISSING)),
    "deceased": ("deceased", _deceased()),
    "age": ("age", _natural()),
    "kind": ("kind", _choice(Kind)),
    "twin": ("twin", _twin()),
    "adopted": ("adoption", _choice(Adoption, MISSING)),
    "order": ("order", _natural()),
    "gender": ("gender", _text()),
}

DEFAULT_CELL = {"father": "0", "mother": "0", "sex": "0", "deceased": "0", "kind": "birth"}


def column_kind(column: str, disease_ids: set[str]) -> str:
    """'fixed', 'phenotype', 'genotype' or 'extra'."""
    if column in INDIVIDUAL_FIELDS:
        return "fixed"
    if column in disease_ids or re.fullmatch(r"D\d+", column):
        return "phenotype"
    if column.startswith(GENOTYPE_PREFIX):
        return "genotype"
    return "extra"


def get_value(ind: Individual, column: str, kind: str) -> Any:
    if kind == "fixed":
        return getattr(ind, INDIVIDUAL_FIELDS[column][0])
    if kind == "phenotype":
        return ind.phenotypes.get(column, Phenotype.UNKNOWN)
    if kind == "genotype":
        return ind.genotypes.get(column[len(GENOTYPE_PREFIX):])
    return ind.extra.get(column, MISSING)


def set_value(ind: Individual, column: str, kind: str, value: Any) -> None:
    if kind == "fixed":
        setattr(ind, INDIVIDUAL_FIELDS[column][0], value)
    elif kind == "phenotype":
        ind.phenotypes[column] = value
    elif kind == "genotype":
        ind.genotypes[column[len(GENOTYPE_PREFIX):]] = value
    else:
        ind.extra[column] = value


def field_for(column: str, kind: str) -> Field | None:
    if kind == "fixed":
        return INDIVIDUAL_FIELDS[column][1]
    if kind == "phenotype":
        return PHENOTYPE_FIELD
    if kind == "genotype":
        return GENOTYPE_FIELD
    return None


# -- metadata fields -------------------------------------------------------

def _meta_text() -> Field:
    def parse(t: str):
        if t == "":
            raise ValueError(t)
        return t
    return Field(parse, lambda v: v, "text")


META_FIELDS: dict[str, dict[str, tuple[str, Field, Any]]] = {
    # key: {field name: (attribute, field, default)}
    "family": {
        "ID": ("id", _meta_text(), None),
        "description": ("description", _meta_text(), None),
    },
    "disease": {
        "ID": ("id", _meta_text(), None),
        "name": ("name", _meta_text(), None),
        "mode": ("mode", _choice(Mode), Mode.UNKNOWN),
        "allele_freq": ("allele_freq", _probability(), None),
        "penetrance": ("penetrance", _probability(), 1.0),
        "mutation_rate": ("mutation_rate", _probability(), 0.0),
        "phenocopy": ("phenocopy", _probability(), 0.0),
        "pattern": ("pattern", _choice(Pattern), Pattern.SOLID),
        "prevalence": ("prevalence", _probability(), None),
        "heritability": ("heritability", _probability(), None),
        "prevalence_male": ("prevalence_male", _probability(), None),
        "prevalence_female": ("prevalence_female", _probability(), None),
    },
    "union": {
        "p1": ("p1", _meta_text(), None),
        "p2": ("p2", _meta_text(), None),
        "consanguinity": ("consanguinity", _meta_text(), None),
        "status": ("status", _choice(UnionStatus), UnionStatus.MARRIED),
    },
    "question": {
        "ID": ("id", _meta_text(), None),
        "type": ("type", _meta_text(), None),
        "individual": ("individual", _meta_text(), None),
        "disease": ("disease", _meta_text(), None),
        "event": ("event", _meta_text(), None),
    },
}


# -- ##key=<k=v,...> syntax ------------------------------------------------

def parse_struct(content: str) -> list[tuple[str, str]]:
    """``k1=v1,k2="v,2"`` -> [(k1, 'v1'), (k2, '"v,2"')]. Values are returned
    raw (quoted); use ``unquote`` to strip quotes."""
    items: list[str] = []
    buf: list[str] = []
    quoted = False
    for ch in content:
        if ch == '"':
            quoted = not quoted
        if ch == "," and not quoted:
            items.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if quoted:
        raise ValueError("unterminated quote")
    items.append("".join(buf))
    result = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"expected 'key=value': {item!r}")
        k, v = item.split("=", 1)
        result.append((k.strip(), v.strip()))
    return result


def unquote(raw: str) -> str:
    if len(raw) >= 2 and raw[0] == raw[-1] == '"':
        return raw[1:-1]
    return raw


def quote(value: str) -> str:
    if re.fullmatch(r"[\w.\-/:+]*", value) and value:
        return value
    return '"' + value.replace('"', "'") + '"'
