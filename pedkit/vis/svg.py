"""Draws a pedigree as SVG with NSGC standard symbols (no external dependencies).

The output is SVG 1.2 Tiny compatible (supported by QtSvg): fill patterns,
clipping and arrowheads are drawn as plain geometry.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from xml.sax.saxutils import escape, quoteattr

from ..model import (
    Individual, Kind, Mode, Pattern, Pedigree, Phenotype, Role, Sex, TwinType, UnionStatus,
)
from . import geometry as geo
from .layout import Layout, compute_layout, key

SIZE = 40.0             # symbol side / diameter
SMALL_RATIO = 0.6       # planned, pregnancy, miscarriage
GEN_SPACING = 140.0
LEFT_MARGIN = 70.0      # for generation numbers
TOP_MARGIN = 40.0
MARGIN = 30.0
STROKE = 1.5
FONT = "font-family=\"DejaVu Sans, Arial, sans-serif\""

MODE_NAME = {
    Mode.AD: "autosomal dominant", Mode.AR: "autosomal recessive",
    Mode.ACD: "autosomal codominant", Mode.MF: "multifactorial",
    Mode.XD: "X-linked dominant", Mode.XR: "X-linked recessive",
    Mode.Y: "Y-linked", Mode.MT: "mitochondrial", Mode.UNKNOWN: "mode unknown",
}


def draw(ped: Pedigree, layout: Layout | None = None) -> str:
    """Return the pedigree as SVG text."""
    return _Renderer(ped, layout or compute_layout(ped)).run()


def save_svg(ped: Pedigree, path: str | os.PathLike) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(draw(ped))


def roman(n: int) -> str:
    values = [(1000, "M"), (900, "CM"), (500, "D"), (400, "CD"), (100, "C"), (90, "XC"),
              (50, "L"), (40, "XL"), (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")]
    result = ""
    for v, s in values:
        while n >= v:
            result += s
            n -= v
    return result


@dataclass
class _Shape:
    """Geometry of a symbol: a polygon for clipping, an SVG element for the outline."""

    polygon: geo.Polygon
    outline: str      # start of the unfilled outline element (e.g. '<rect ...')
    top: float        # y of the top edge (the sibship line comes down to it)
    bottom: float
    radius: float     # horizontal half-width


def _n(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


def _line(x1, y1, x2, y2, width=STROKE) -> str:
    return (f'<line x1="{_n(x1)}" y1="{_n(y1)}" x2="{_n(x2)}" y2="{_n(y2)}" '
            f'stroke="#000" stroke-width="{width}"/>')


def _polygon(points: geo.Polygon, fill="#000", stroke: str | None = None) -> str:
    p = " ".join(f"{_n(x)},{_n(y)}" for x, y in points)
    s = f' stroke="{stroke}" stroke-width="{STROKE}"' if stroke else ' stroke="none"'
    return f'<polygon points="{p}" fill="{fill}"{s}/>'


def _text(x, y, text, size=11, anchor="middle", bold=False) -> str:
    b = ' font-weight="bold"' if bold else ""
    return (f'<text x="{_n(x)}" y="{_n(y)}" font-size="{size}" text-anchor="{anchor}" '
            f'{FONT}{b} fill="#000">{escape(text)}</text>')


class _Renderer:
    def __init__(self, ped: Pedigree, layout: Layout):
        self.ped = ped
        self.lay = layout
        self.lines: list[str] = []
        self.symbols: list[str] = []
        self.right = 0.0
        self.bottom = 0.0
        self.has_carrier = False
        self.has_intermediate = False

    # -- coordinates -------------------------------------------------------------

    def pos(self, ind: Individual) -> tuple[float, float]:
        k = key(ind)
        return (self.lay.x[k] + LEFT_MARGIN + SIZE,
                self.lay.generation[k] * GEN_SPACING + TOP_MARGIN + SIZE / 2)

    def shape(self, ind: Individual) -> _Shape:
        x, y = self.pos(ind)
        small = ind.kind in (Kind.PLANNED, Kind.PREGNANCY, Kind.MISCARRIAGE)
        s = SIZE * (SMALL_RATIO if small else 1.0)
        r = s / 2
        if ind.kind is Kind.MISCARRIAGE:
            h = r * 1.1
            pts = [(x, y - h), (x + h, y + h * 0.75), (x - h, y + h * 0.75)]
            return _Shape(pts, _polygon_outline(pts), y - h, y + h * 0.75, h)
        if ind.sex is Sex.MALE:
            pts = [(x - r, y - r), (x + r, y - r), (x + r, y + r), (x - r, y + r)]
            outline = (f'<rect x="{_n(x - r)}" y="{_n(y - r)}" width="{_n(s)}" '
                       f'height="{_n(s)}"')
            return _Shape(pts, outline, y - r, y + r, r)
        if ind.sex is Sex.FEMALE:
            outline = f'<circle cx="{_n(x)}" cy="{_n(y)}" r="{_n(r)}"'
            return _Shape(geo.circle(x, y, r), outline, y - r, y + r, r)
        h = r * 1.2
        pts = [(x, y - h), (x + h, y), (x, y + h), (x - h, y)]
        return _Shape(pts, _polygon_outline(pts), y - h, y + h, h)

    # -- main ---------------------------------------------------------------------

    def run(self) -> str:
        for ind in self.lay.individuals:
            x, yy = self.pos(ind)
            self.right = max(self.right, x + SIZE)
            self.bottom = max(self.bottom, yy + SIZE / 2 + 50)
        self._couples()
        self._sibships()
        for ind in self.lay.individuals:
            self._symbol(ind)
        gen_labels = [
            _text(MARGIN, g * GEN_SPACING + TOP_MARGIN + SIZE / 2 + 5, roman(g + 1), 14,
                  "start", True)
            for g in range(self.lay.n_generations)]
        legend, lw, lh = self._legend()
        width = max(self.right + MARGIN, lw + 2 * MARGIN + LEFT_MARGIN)
        height = self.bottom + lh + MARGIN
        legend_x = width - MARGIN - lw
        legend_y = self.bottom
        parts = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            f'<svg xmlns="http://www.w3.org/2000/svg" version="1.2" baseProfile="tiny" '
            f'width="{_n(width)}" height="{_n(height)}" '
            f'viewBox="0 0 {_n(width)} {_n(height)}">',
            f'<title>{escape(self._title())}</title>',
            f'<rect x="0" y="0" width="{_n(width)}" height="{_n(height)}" fill="#fff"/>',
            '<g class="generations">', *gen_labels, '</g>',
            '<g class="lines">', *self.lines, '</g>',
            '<g class="symbols">', *self.symbols, '</g>',
            f'<g class="legend" transform="translate({_n(legend_x)},{_n(legend_y)})">',
            *legend, '</g>',
            '</svg>',
        ]
        return "\n".join(parts) + "\n"

    def _title(self) -> str:
        if self.ped.family and (self.ped.family.description or self.ped.family.id):
            return self.ped.family.description or self.ped.family.id
        return "Pedigree"

    # -- lines ----------------------------------------------------------------------

    def _couples(self) -> None:
        for a, b in self.lay.couples:
            (xa, ya), (xb, yb) = self.pos(a), self.pos(b)
            if xa > xb:
                a, b, xa, ya, xb, yb = b, a, xb, yb, xa, ya
            x1 = xa + self.shape(a).radius
            x2 = xb - self.shape(b).radius
            u = self.ped.union(a.id, b.id)
            consanguineous = (u is not None and u.consanguineous) or bool(
                self.ped.ancestors(a) & self.ped.ancestors(b))
            if consanguineous:
                self.lines += [_line(x1, ya - 3, x2, yb - 3), _line(x1, ya + 3, x2, yb + 3)]
            else:
                self.lines.append(_line(x1, ya, x2, yb))
            mx, my = (x1 + x2) / 2, (ya + yb) / 2
            if u is not None and u.status is UnionStatus.SEPARATED:
                for dx in (-4, 4):
                    self.lines.append(_line(mx + dx - 5, my + 8, mx + dx + 5, my - 8))
            elif u is not None and u.status is UnionStatus.CHILDLESS:
                self.lines += [_line(mx, my, mx, my + 25),
                               _line(mx - 10, my + 25, mx + 10, my + 25)]

    def _sibships(self) -> None:
        for sib in self.lay.sibships:
            if not sib.children:
                continue
            first_child_y = min(self.pos(c)[1] for c in sib.children)
            bar_y = first_child_y - SIZE / 2 - 22
            if len(sib.parents) == 2:
                (xa, ya), (xb, yb) = (self.pos(p) for p in sib.parents)
                ra, rb = self.shape(sib.parents[0]).radius, self.shape(sib.parents[1]).radius
                if xa <= xb:
                    drop_x = (xa + ra + xb - rb) / 2
                else:
                    drop_x = (xb + rb + xa - ra) / 2
                drop_y = (ya + yb) / 2
            else:
                p = sib.parents[0]
                drop_x, _ = self.pos(p)
                drop_y = self.shape(p).bottom
            self.lines.append(_line(drop_x, drop_y, drop_x, bar_y))

            anchors = [drop_x]
            for group in self._twin_groups(sib.children):
                xs = [self.pos(c)[0] for c in group]
                if len(group) == 1:
                    c = group[0]
                    self.lines.append(_line(xs[0], bar_y, xs[0], self.shape(c).top))
                    anchors.append(xs[0])
                    continue
                apex = sum(xs) / len(xs)
                anchors.append(apex)
                ends = []
                for c, cx in zip(group, xs):
                    top = self.shape(c).top
                    self.lines.append(_line(apex, bar_y, cx, top))
                    ends.append((cx, top))
                if group[0].twin.type is TwinType.MZ:
                    t = 0.5
                    (x1, y1), (x2, y2) = ends[0], ends[-1]
                    self.lines.append(_line(apex + t * (x1 - apex), bar_y + t * (y1 - bar_y),
                                            apex + t * (x2 - apex), bar_y + t * (y2 - bar_y)))
                elif group[0].twin.type is TwinType.UNKNOWN:
                    self.lines.append(_text(apex, bar_y + 18, "?", 11))
            self.lines.append(_line(min(anchors), bar_y, max(anchors), bar_y))

    @staticmethod
    def _twin_groups(children: list[Individual]) -> list[list[Individual]]:
        groups: list[list[Individual]] = []
        for c in children:
            if c.twin and groups and groups[-1][0].twin and \
                    groups[-1][0].twin.group == c.twin.group:
                groups[-1].append(c)
            else:
                groups.append([c])
        return groups

    # -- symbol ---------------------------------------------------------------------

    def _symbol(self, ind: Individual) -> None:
        x, y = self.pos(ind)
        sh = self.shape(ind)
        o: list[str] = []
        o.append(f"{sh.outline} fill=\"#fff\" stroke=\"none\"/>")
        o += self._fill(ind, sh)
        o.append(f"{sh.outline} fill=\"none\" stroke=\"#000\" stroke-width=\"{STROKE}\"/>")

        if ind.kind in (Kind.PLANNED, Kind.PREGNANCY):
            o.append(_text(x, y + 4, "P", 11, bold=True))
        if ind.deceased or ind.kind is Kind.STILLBIRTH:
            e = sh.radius + 7
            o.append(_line(x - e, sh.bottom + 7, x + e, sh.top - 7))
        if ind.adoption is not None:
            e, t, b = sh.radius + 5, sh.top - 4, sh.bottom + 4
            o.append(f'<polyline points="{_n(x - e + 5)},{_n(t)} {_n(x - e)},{_n(t)} '
                     f'{_n(x - e)},{_n(b)} {_n(x - e + 5)},{_n(b)}" fill="none" stroke="#000" '
                     f'stroke-width="{STROKE}"/>')
            o.append(f'<polyline points="{_n(x + e - 5)},{_n(t)} {_n(x + e)},{_n(t)} '
                     f'{_n(x + e)},{_n(b)} {_n(x + e - 5)},{_n(b)}" fill="none" stroke="#000" '
                     f'stroke-width="{STROKE}"/>')
        if ind.role is not None:
            o += self._arrow(ind, x, sh)

        label_y = sh.bottom + 14
        labels = [ind.id]
        if ind.age is not None:
            labels.append(f"d. {ind.age}" if ind.deceased else str(ind.age))
        if ind.kind is Kind.STILLBIRTH:
            labels.append("SB")
        if ind.gender:
            labels.append(ind.gender)
        for i, lab in enumerate(labels):
            o.append(_text(x, label_y + 13 * i, lab, 11))
        self.bottom = max(self.bottom, label_y + 13 * len(labels) + 10)

        title = escape(self._summary(ind))
        self.symbols.append(
            f'<g class="individual" id={quoteattr(_xml_id(ind))}><title>{title}</title>'
            + "".join(o) + "</g>")

    def _fill(self, ind: Individual, sh: _Shape) -> list[str]:
        diseases = self.ped.diseases
        if not diseases:
            return []
        x0, y0, x1, y1 = geo.bounds(sh.polygon)
        w = (x1 - x0) / len(diseases)
        out: list[str] = []
        born = ind.kind in (Kind.BIRTH, Kind.STILLBIRTH)
        for i, d in enumerate(diseases):
            a, c = x0 + i * w, x0 + (i + 1) * w
            part = geo.clip_box(sh.polygon, xmin=a, xmax=c)
            if len(part) < 3:
                continue
            mid_x = (a + c) / 2
            mid_y = (sh.top + sh.bottom) / 2 + (0 if ind.kind is not Kind.MISCARRIAGE
                                                else sh.radius * 0.2)
            p = ind.phenotype(d.id)
            if p is Phenotype.AFFECTED:
                out += _pattern(part, d.pattern, a, c, (sh.top + sh.bottom) / 2)
            elif p is Phenotype.INTERMEDIATE:
                # intermediate phenotype (codominant heterozygote): lower half filled
                lower = geo.clip_box(part, ymin=(sh.top + sh.bottom) / 2)
                if len(lower) >= 3:
                    self.has_intermediate = True
                    out += _pattern(lower, d.pattern, a, c, (sh.top + sh.bottom) / 2)
            elif p is Phenotype.UNKNOWN and born:
                out.append(_text(mid_x, mid_y + 5, "?", 14 if len(diseases) == 1 else 11))
            g = ind.genotype(d.id)
            if p is Phenotype.UNAFFECTED and g is not None and g.heterozygous:
                self.has_carrier = True
                out.append(f'<circle cx="{_n(mid_x)}" cy="{_n(mid_y)}" r="4" fill="#000"/>')
        if len(diseases) > 1 and born:
            for i in range(1, len(diseases)):
                kx = x0 + i * w
                seg = geo.clip_segment((kx, sh.top - 1), (kx, sh.bottom + 1), sh.polygon)
                if seg:
                    out.append(_line(*seg[0], *seg[1], width=0.8))
        return out

    def _arrow(self, ind: Individual, x: float, sh: _Shape) -> list[str]:
        # arrow pointing at the symbol from the lower left
        ux, uy = x - sh.radius * 0.75, sh.bottom - sh.radius * 0.25 + 2
        bx, by = ux - 22, uy + 18
        ang = math.atan2(uy - by, ux - bx)
        length, half = 8.0, 4.0
        tx, ty = ux - length * math.cos(ang), uy - length * math.sin(ang)
        head = [(ux, uy), (tx - half * math.sin(ang), ty + half * math.cos(ang)),
                (tx + half * math.sin(ang), ty - half * math.cos(ang))]
        o = [_line(bx, by, tx, ty), _polygon(head)]
        if ind.role is Role.PROBAND:
            o.append(_text(bx - 3, by + 5, "P", 12, "end", bold=True))
        return o

    def _summary(self, ind: Individual) -> str:
        parts = [ind.id]
        for d in self.ped.diseases:
            p = ind.phenotype(d.id)
            g = ind.genotype(d.id)
            state = {Phenotype.AFFECTED: "affected", Phenotype.UNAFFECTED: "unaffected",
                     Phenotype.UNKNOWN: "unknown", Phenotype.INTERMEDIATE: "intermediate"}[p]
            parts.append(f"{d.id}: {state}" + (f" ({g})" if g else ""))
        return " · ".join(parts)

    # -- legend -----------------------------------------------------------------------

    def _legend(self) -> tuple[list[str], float, float]:
        rows: list[tuple[str, list[str]]] = []
        k = 16.0
        for d in self.ped.diseases:
            square = [(0, 0), (k, 0), (k, k), (0, k)]
            art = [_polygon(square, "#fff")] + _pattern(square, d.pattern, 0, k, k / 2) + \
                  [_polygon(square, "none", "#000")]
            name = d.name or d.id
            text = f"{d.id}: {name}" + (f" ({MODE_NAME[d.mode]})" if d.mode is not Mode.UNKNOWN
                                        else "")
            rows.append((text, art))
        if self.has_intermediate:
            square = [(0, 0), (k, 0), (k, k), (0, k)]
            rows.append(("Intermediate phenotype (codominant heterozygote)",
                         [_polygon(square, "#fff"),
                          _polygon([(0, k / 2), (k, k / 2), (k, k), (0, k)]),
                          _polygon(square, "none", "#000")]))
        if self.has_carrier:
            rows.append(("Tested heterozygous carrier",
                         [_polygon([(0, 0), (k, 0), (k, k), (0, k)], "#fff", "#000"),
                          f'<circle cx="{_n(k / 2)}" cy="{_n(k / 2)}" r="3.5" fill="#000"/>']))
        if not rows:
            return [], 0.0, 0.0
        o: list[str] = []
        width = 0.0
        for i, (text, art) in enumerate(rows):
            yy = i * 24
            o.append(f'<g transform="translate(0,{_n(yy)})">' + "".join(art) + "</g>")
            o.append(_text(k + 8, yy + 12.5, text, 12, "start"))
            width = max(width, k + 8 + 7.0 * len(text))
        return o, width, 24.0 * len(rows)


def _pattern(part: geo.Polygon, pattern: Pattern, a: float, c: float, mid_y: float) -> list[str]:
    """Fill of an affected part; ``a``-``c`` is the part's horizontal range."""
    if pattern is Pattern.SOLID:
        return [_polygon(part)]
    if pattern is Pattern.QUADRANT:
        quad = geo.clip_box(part, xmax=(a + c) / 2, ymax=mid_y)
        return [_polygon(quad)] if len(quad) >= 3 else []
    x0, y0, x1, y1 = geo.bounds(part)
    if pattern is Pattern.HATCHED:
        o = []
        step = 5.0 * math.sqrt(2)
        const = x0 + y0
        while const < x1 + y1:
            const += step
            seg = geo.clip_segment((const - y0 + 1, y0 - 1), (const - y1 - 1, y1 + 1), part)
            if seg:
                o.append(_line(*seg[0], *seg[1], width=1.2))
        return o
    # dotted
    o = []
    step = 5.0
    yy = y0 + step / 2
    row = 0
    while yy < y1:
        xx = x0 + step / 2 + (step / 2 if row % 2 else 0)
        while xx < x1:
            if geo.inside((xx, yy), part):
                o.append(f'<circle cx="{_n(xx)}" cy="{_n(yy)}" r="1.1" fill="#000"/>')
            xx += step
        yy += step
        row += 1
    return o


def _polygon_outline(points: geo.Polygon) -> str:
    p = " ".join(f"{_n(x)},{_n(y)}" for x, y in points)
    return f'<polygon points="{p}"'


def _xml_id(ind: Individual) -> str:
    safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in f"{ind.family}-{ind.id}")
    return "individual-" + safe
