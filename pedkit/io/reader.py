"""PEDX reader.

``read`` parses the file structurally and never stops on a content error:
syntax and value errors are collected in ``pedigree.source.read_findings``;
semantic checks are done by ``pedkit.validate``.
"""

from __future__ import annotations

import os

from ..model import (
    Disease, FamilyInfo, FileSource, Finding, Individual, Level, MetaLine, Pedigree,
    Question, Source, Union,
)
from ._fields import (
    FIXED_COLUMNS, META_FIELDS, column_kind, field_for, get_value, parse_struct, set_value,
    unquote,
)

_META_CLASSES = {
    "family": FamilyInfo,
    "disease": Disease,
    "union": Union,
    "question": Question,
}


def read(path: str | os.PathLike) -> Pedigree:
    """Read a PEDX file."""
    with open(path, encoding="utf-8") as f:
        return read_text(f.read(), path=os.fspath(path))


def read_text(text: str, path: str | None = None) -> Pedigree:
    """Parse PEDX text (for tests and in-memory use)."""
    return _Parser(text, path).run()


class _Parser:
    def __init__(self, text: str, path: str | None):
        self.src = FileSource(path=path, final_newline=text.endswith("\n"))
        self.ped = Pedigree(format="", source=self.src)
        lines = text.split("\n")
        if self.src.final_newline:
            lines.pop()
        self.lines = lines

    def finding(self, rule: str, message: str, line: int | None, individual: str | None = None,
                level: Level = Level.ERROR) -> None:
        self.src.read_findings.append(Finding(level, rule, message, line, individual))

    # -- main loop ------------------------------------------------------------

    def run(self) -> Pedigree:
        data_lines: list[tuple[int, str, int]] = []  # (line no, text, layout index)
        for no, line in enumerate(self.lines, start=1):
            if line.startswith("##"):
                self._meta(no, line)
            elif line.startswith("#"):
                if self.src.header_line is None:
                    self.ped.columns = line[1:].split("\t")
                    self.src.columns = tuple(self.ped.columns)
                    self.src.header_line = no
                    self.src.layout.append(("header",))
                else:
                    self.src.layout.append(("raw", line))  # comment
            elif not line.strip():
                self.src.layout.append(("raw", line))
            elif self.src.header_line is None:
                self.src.layout.append(("raw", line))
                data_lines.append((no, line, -1))
            else:
                data_lines.append((no, line, len(self.src.layout)))
                self.src.layout.append(("individual", None))

        self._check_format()
        if self.src.header_line is None:
            self.finding("E1", "Column header (#family\tindividual\t...) not found.",
                         data_lines[0][0] if data_lines else None)
            return self.ped
        self._check_header()
        for no, line, idx in data_lines:
            ind = self._individual(no, line)
            self.ped.individuals.append(ind)
            self.src.layout[idx] = ("individual", ind)
        return self.ped

    # -- metadata -------------------------------------------------------------

    def _meta(self, no: int, line: str) -> None:
        body = line[2:]
        if "=" not in body:
            self.finding("E8", f"Metadata line is not of the form '##key=value': {line!r}", no)
            self.src.layout.append(("raw", line))
            return
        key, value = body.split("=", 1)

        if key == "format":
            if self.src.format_line is None:
                self.ped.format = value
                self.src.format_line = no
                self.src.layout.append(("format", line))
            else:
                self.finding("E8", "More than one ##format line.", no)
                self.src.layout.append(("raw", line))
            return

        if key in _META_CLASSES and value.startswith("<") and value.endswith(">"):
            obj = self._structured_meta(no, line, key, value[1:-1])
            if obj is not None:
                self.src.layout.append(("meta", obj))
                return

        meta = MetaLine(key, value, Source(no, line))
        self.ped.other_meta.append(meta)
        self.src.layout.append(("meta", meta))

    def _structured_meta(self, no: int, line: str, key: str, content: str):
        try:
            items = parse_struct(content)
        except ValueError as e:
            self.finding("E8", f"Could not parse ##{key} line: {e}", no)
            return None

        cls = _META_CLASSES[key]
        defs = META_FIELDS[key]
        required = {"disease": ["ID"], "union": ["p1", "p2"]}.get(key, [])
        values: dict = {}
        source = Source(no, line)
        extra: dict[str, str] = {}

        for name, raw in items:
            if name in defs:
                attr, fld, default = defs[name]
                try:
                    parsed = fld.parse(unquote(raw))
                except ValueError:
                    self.finding("E8", f"Invalid value {raw!r} for '{name}' in ##{key} line; "
                                       f"allowed: {fld.allowed}.", no)
                    parsed = default
                values[attr] = parsed
                source.fields[name] = (raw, parsed)
            else:
                extra[name] = unquote(raw)
                source.fields[name] = (raw, extra[name])

        for name in required:
            if name not in source.fields:
                self.finding("E8", f"Required field '{name}' missing in ##{key} line.", no)
                values.setdefault(defs[name][0], "")

        obj = cls(**values, extra=extra, source=source)
        if key == "family":
            if self.ped.family is not None:
                self.finding("E8", "More than one ##family line.", no, level=Level.WARNING)
            self.ped.family = obj
        elif key == "disease":
            if self.ped.disease(obj.id) is not None:
                self.finding("E8", f"Disease '{obj.id}' is defined more than once.", no)
            self.ped.diseases.append(obj)
        elif key == "union":
            self.ped.unions.append(obj)
        else:
            self.ped.questions.append(obj)
        return obj

    # -- header ---------------------------------------------------------------

    def _check_format(self) -> None:
        first = self.lines[0] if self.lines else ""
        if not first.startswith("##format="):
            self.finding("E1", "The first line must be '##format=PEDX/...'.", 1)
        elif not self.ped.format.startswith("PEDX/"):
            self.finding("E1", f"Unrecognised format '{self.ped.format}'; expected 'PEDX/...'.", 1)

    def _check_header(self) -> None:
        no = self.src.header_line
        cols = self.ped.columns
        missing = [c for c in FIXED_COLUMNS if c not in cols]
        if missing:
            self.finding("E2", f"Required column missing: {', '.join(missing)}.", no)
        elif cols[:5] != FIXED_COLUMNS:
            self.finding("E2", "The first five columns must be family, individual, father, "
                               f"mother, sex in this order; found: {', '.join(cols[:5])}.", no)
        primary = self.ped.primary_disease
        if primary is not None and primary.id in cols and (len(cols) < 6 or cols[5] != primary.id):
            self.finding("E2", f"The sixth column must be the phenotype of the primary disease "
                               f"({primary.id}).", no)
        dup = {c for c in cols if cols.count(c) > 1}
        if dup:
            self.finding("E2", f"Duplicate column name: {', '.join(sorted(dup))}.", no)

    # -- individual line --------------------------------------------------------

    def _individual(self, no: int, line: str) -> Individual:
        cols = self.ped.columns
        cells = line.split("\t")
        disease_ids = {d.id for d in self.ped.diseases}
        row = dict(zip(cols, cells))
        ind = Individual(family=row.get("family", ""), id=row.get("individual", f"?line{no}"))
        ind.source = Source(no, line)

        if len(cells) != len(cols):
            self.finding("E8", f"Line has {len(cells)} cells, the header has {len(cols)} columns.",
                         no, ind.id)

        for j, col in enumerate(cols):
            kind = column_kind(col, disease_ids)
            cell = cells[j] if j < len(cells) else None
            fld = field_for(col, kind)
            if cell is None:
                continue  # missing cell: the default value stays
            if fld is None:  # unknown column
                set_value(ind, col, kind, cell)
                ind.source.fields[col] = (cell, cell)
                continue
            try:
                value = fld.parse(cell)
            except ValueError:
                self.finding("E8", f"Invalid value {cell!r} in column '{col}'; "
                                   f"allowed: {fld.allowed}.", no, ind.id)
                # Invalid cell: the model keeps its default, the raw text is
                # kept so it can be written back.
                ind.source.fields[col] = (cell, get_value(ind, col, kind))
                continue
            set_value(ind, col, kind, value)
            ind.source.fields[col] = (cell, value)
        return ind
