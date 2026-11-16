"""Convex polygon clipping helpers.

QtSvg supports only SVG 1.2 Tiny (no ``clipPath``, ``pattern`` or
``marker``), so fill patterns and split symbols are drawn as plain clipped
geometry.
"""

from __future__ import annotations

import math

Point = tuple[float, float]
Polygon = list[Point]


def circle(cx: float, cy: float, r: float, n: int = 72) -> Polygon:
    return [(cx + r * math.cos(2 * math.pi * i / n), cy + r * math.sin(2 * math.pi * i / n))
            for i in range(n)]


def clip_box(poly: Polygon, xmin=-math.inf, xmax=math.inf,
             ymin=-math.inf, ymax=math.inf) -> Polygon:
    """Sutherland–Hodgman: clip a polygon to an axis-aligned rectangle."""
    edges = [
        (lambda p: p[0] >= xmin, lambda a, b: _x_cross(a, b, xmin)),
        (lambda p: p[0] <= xmax, lambda a, b: _x_cross(a, b, xmax)),
        (lambda p: p[1] >= ymin, lambda a, b: _y_cross(a, b, ymin)),
        (lambda p: p[1] <= ymax, lambda a, b: _y_cross(a, b, ymax)),
    ]
    result = list(poly)
    for inside_fn, cross in edges:
        if not result:
            break
        src, result = result, []
        prev = src[-1]
        for p in src:
            if inside_fn(p):
                if not inside_fn(prev):
                    result.append(cross(prev, p))
                result.append(p)
            elif inside_fn(prev):
                result.append(cross(prev, p))
            prev = p
    return result


def _x_cross(a: Point, b: Point, x: float) -> Point:
    t = (x - a[0]) / (b[0] - a[0])
    return (x, a[1] + t * (b[1] - a[1]))


def _y_cross(a: Point, b: Point, y: float) -> Point:
    t = (y - a[1]) / (b[1] - a[1])
    return (a[0] + t * (b[0] - a[0]), y)


def _orientation(poly: Polygon) -> float:
    """Positive (signed area) for counter-clockwise polygons."""
    return sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(poly, poly[1:] + poly[:1]))


def clip_segment(p1: Point, p2: Point, poly: Polygon) -> tuple[Point, Point] | None:
    """Cyrus–Beck: clip a line segment to the inside of a convex polygon."""
    if len(poly) < 3:
        return None
    sign = 1.0 if _orientation(poly) > 0 else -1.0
    t0, t1 = 0.0, 1.0
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    for a, b in zip(poly, poly[1:] + poly[:1]):
        # inward normal
        nx, ny = -(b[1] - a[1]) * sign, (b[0] - a[0]) * sign
        num = nx * (p1[0] - a[0]) + ny * (p1[1] - a[1])
        den = nx * dx + ny * dy
        if abs(den) < 1e-12:
            if num < 0:
                return None
            continue
        t = -num / den
        if den > 0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
        if t0 > t1:
            return None
    return ((p1[0] + t0 * dx, p1[1] + t0 * dy), (p1[0] + t1 * dx, p1[1] + t1 * dy))


def inside(p: Point, poly: Polygon) -> bool:
    sign = 1.0 if _orientation(poly) > 0 else -1.0
    for a, b in zip(poly, poly[1:] + poly[:1]):
        if ((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])) * sign < 0:
            return False
    return True


def bounds(poly: Polygon) -> tuple[float, float, float, float]:
    xs = [p[0] for p in poly]
    ys = [p[1] for p in poly]
    return min(xs), min(ys), max(xs), max(ys)
