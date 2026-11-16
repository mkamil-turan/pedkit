"""PEDX writer and PLINK ``.fam`` export.

Lossless round trip: if a field has not changed since it was read, its raw
text from the file is written; if no field of a line has changed, the line
itself is written verbatim. A file saved without changes is therefore
byte-for-byte identical.
"""

from __future__ import annotations

import os
from typing import Any

from ..model import (
    Disease, FamilyInfo, Individual, MetaLine, Pedigree, Phenotype, Question, Union,
)
from ._fields import (
    FIXED_COLUMNS, GENOTYPE_PREFIX, INDIVIDUAL_FIELDS, META_FIELDS, column_kind, field_for,
    get_value, quote,
)

_META_KEY = {FamilyInfo: "family", Disease: "disease", Union: "union", Question: "question"}
_OPTIONAL_COLUMNS = ["role", "deceased", "age", "kind", "twin", "adopted", "order", "gender"]


def write(ped: Pedigree, path: str | os.PathLike) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(write_text(ped))


def write_text(ped: Pedigree) -> str:
    lines: list[str] = []
    layout = ped.source.layout if ped.source else []
    columns = active_columns(ped)
    raw_line_ok = bool(ped.source) and tuple(columns) == ped.source.columns

    metas = _meta_objects(ped)
    written: set[int] = set()

    def format_line(raw: str | None = None) -> str:
        new = f"##format={ped.format}"
        return raw if raw == new else new

    def write_pending_meta() -> None:
        for m in metas:
            if id(m) not in written:
                lines.append(_meta_line(m))
                written.add(id(m))

    if not any(item[0] == "format" for item in layout):
        lines.append(format_line())

    header_written = False
    live_meta = {id(m) for m in metas}
    live_inds = {id(i) for i in ped.individuals}
    for item in layout:
        kind = item[0]
        if kind == "format":
            lines.append(format_line(item[1]))
        elif kind == "meta":
            if id(item[1]) in live_meta:
                lines.append(_meta_line(item[1]))
                written.add(id(item[1]))
        elif kind == "header":
            write_pending_meta()
            lines.append("#" + "\t".join(columns))
            header_written = True
        elif kind == "individual":
            ind = item[1]
            if ind is not None and id(ind) in live_inds:
                lines.append(_individual_line(ped, ind, columns, raw_line_ok))
                written.add(id(ind))
        else:  # raw
            lines.append(item[1])

    if not header_written:
        write_pending_meta()
        lines.append("#" + "\t".join(columns))
    for ind in ped.individuals:
        if id(ind) not in written:
            lines.append(_individual_line(ped, ind, columns, raw_line_ok))

    end = "\n" if (ped.source is None or ped.source.final_newline) else ""
    return "\n".join(lines) + end


def active_columns(ped: Pedigree) -> list[str]:
    """Columns to write: the file's header plus any missing disease columns.
    Without a header (pedigree built in memory) the needed columns are chosen."""
    if ped.columns:
        cols = list(ped.columns)
        for d in ped.diseases:
            if d.id not in cols:
                cols.append(d.id)
        return cols

    blank = Individual(family="", id="")
    cols = FIXED_COLUMNS + [d.id for d in ped.diseases]
    for d in ped.diseases:
        if any(i.genotype(d.id) is not None for i in ped.individuals):
            cols.append(GENOTYPE_PREFIX + d.id)
    for c in _OPTIONAL_COLUMNS:
        attr = INDIVIDUAL_FIELDS[c][0]
        if any(getattr(i, attr) != getattr(blank, attr) for i in ped.individuals):
            cols.append(c)
    for i in ped.individuals:
        for c in i.extra:
            if c not in cols:
                cols.append(c)
    return cols


# -- line builders -----------------------------------------------------------

def _individual_line(ped: Pedigree, ind: Individual, columns: list[str], raw_ok: bool) -> str:
    disease_ids = {d.id for d in ped.diseases}
    blank = Individual(family="", id="")
    src_fields = ind.source.fields if ind.source else {}
    changed = not (raw_ok and ind.source and ind.source.raw_line is not None)
    cells = []
    for c in columns:
        kind = column_kind(c, disease_ids)
        current = get_value(ind, c, kind)
        before = src_fields.get(c)
        if before is not None and current == before[1]:
            cells.append(before[0])
            continue
        if before is not None or current != get_value(blank, c, kind):
            changed = True
        fld = field_for(c, kind)
        cells.append(fld.format(current) if fld else str(current))
    if not changed:
        return ind.source.raw_line
    return "\t".join(cells)


def _meta_objects(ped: Pedigree) -> list:
    objs: list = []
    if ped.family is not None:
        objs.append(ped.family)
    return objs + ped.diseases + ped.unions + ped.questions + ped.other_meta


def _meta_line(obj) -> str:
    if isinstance(obj, MetaLine):
        return f"##{obj.key}={obj.value}"

    key = _META_KEY[type(obj)]
    defs = META_FIELDS[key]
    src_fields = obj.source.fields if obj.source else {}
    items: list[str] = []
    changed = obj.source is None or obj.source.raw_line is None

    def current(name: str) -> Any:
        if name in defs:
            return getattr(obj, defs[name][0])
        return obj.extra.get(name)

    def fmt(name: str, value: Any) -> str:
        text = defs[name][1].format(value) if name in defs else str(value)
        return quote(text)

    for name, (raw, before) in src_fields.items():
        v = current(name)
        if v == before:
            items.append(f"{name}={raw}")
        else:
            changed = True
            if v is not None:
                items.append(f"{name}={fmt(name, v)}")
    to_add = [n for n in defs if n not in src_fields] + \
             [n for n in obj.extra if n not in src_fields]
    for name in to_add:
        v = current(name)
        default = defs[name][2] if name in defs else None
        if v is not None and v != default:
            changed = True
            items.append(f"{name}={fmt(name, v)}")

    if not changed:
        return obj.source.raw_line
    return f"##{key}=<{','.join(items)}>"


# -- PLINK .fam ---------------------------------------------------------------

def write_fam(ped: Pedigree, path: str | os.PathLike, disease: str | Disease | None = None) -> None:
    """Write the first six columns as a PLINK ``.fam`` file (missing phenotype ``0``)."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(fam_text(ped, disease))


def fam_text(ped: Pedigree, disease: str | Disease | None = None) -> str:
    if isinstance(disease, Disease):
        did: str | None = disease.id
    elif disease is not None:
        if ped.disease(disease) is None:
            raise KeyError(f"Disease '{disease}' is not in the pedigree.")
        did = disease
    else:
        did = ped.primary_disease.id if ped.primary_disease else None

    lines = []
    for ind in ped.individuals:
        fields = [ind.family, ind.id, ind.father or "0", ind.mother or "0", str(int(ind.sex)),
                  str(_fam_phenotype(ind.phenotype(did) if did else Phenotype.UNKNOWN))]
        if any(not f or any(ch.isspace() for ch in f) for f in fields[:4]):
            raise ValueError(f"PLINK IDs cannot contain whitespace: {ind.family}/{ind.id}")
        lines.append("\t".join(fields))
    return "".join(line + "\n" for line in lines)


def _fam_phenotype(p: Phenotype) -> int:
    """PLINK binary phenotype: the intermediate phenotype (codominant
    heterozygote) is expressed, so it is counted as a case (2)."""
    return 2 if p is Phenotype.INTERMEDIATE else int(p)
