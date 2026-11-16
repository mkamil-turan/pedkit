"""Generalised Elston–Stewart peeling (variable elimination).

A pedigree likelihood is the sum, over all genotypes, of a product of
factors: founder priors P(g_i), transmissions P(g_c | g_mother, g_father) and
penetrances P(y_i | g_i). The Elston–Stewart algorithm evaluates this sum by
"peeling" nuclear families; Cannings, Thompson and Skolnick (1978)
generalised it to any elimination order. Here, at each step the variable
whose elimination produces the smallest intermediate table is summed out
(greedy minimum-size order). On loop-free pedigrees this is equivalent to
classic nuclear-family peeling; with consanguinity loops it is still exact,
only the intermediate tables grow.

The module is shared by the Mendelian inference (``pedkit.inference``) and
the multifactorial model (``pedkit.mf``). To avoid underflow each
intermediate table is divided by its maximum and the log scale is
accumulated separately.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

_LETTERS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"


@dataclass
class Factor:
    variables: tuple[int, ...]
    table: np.ndarray


def eliminate(factors: list[Factor], keep: tuple[int, ...] = ()) -> tuple[np.ndarray, float]:
    """Sum out every variable not in ``keep``.

    Returns (table over the kept variables, log scale). The true value is
    ``table * exp(log_scale)``. With an empty ``keep`` the table is a scalar.
    """
    factors = list(factors)
    log_scale = 0.0
    every = {v for f in factors for v in f.variables}
    remaining = every - set(keep)
    size = {v: f.table.shape[k] for f in factors for k, v in enumerate(f.variables)}

    while remaining:
        def cost(v: int) -> int:
            union = {u for f in factors if v in f.variables for u in f.variables}
            return math.prod(size[u] for u in union)

        v = min(remaining, key=cost)
        involved = [f for f in factors if v in f.variables]
        others = [f for f in factors if v not in f.variables]
        new = _multiply_sum(involved, drop=v)
        largest = float(new.table.max()) if new.table.size else 0.0
        if largest > 0:
            new.table = new.table / largest
            log_scale += math.log(largest)
        else:
            return (np.zeros([size[t] for t in keep]) if keep else np.zeros(())), -math.inf
        factors = others + [new]
        remaining.discard(v)

    result = _multiply_sum(factors, drop=None, order=keep)
    return result.table, log_scale


def _multiply_sum(factors: list[Factor], drop: int | None,
                  order: tuple[int, ...] | None = None) -> Factor:
    variables: list[int] = []
    for f in factors:
        for v in f.variables:
            if v not in variables:
                variables.append(v)
    if len(variables) > len(_LETTERS):
        raise ValueError("Elimination produced a table with more than 52 variables; "
                         "the pedigree is too densely looped for exact peeling.")
    letter = {v: _LETTERS[i] for i, v in enumerate(variables)}
    out = list(order) if order is not None else [v for v in variables if v != drop]
    inputs = ",".join("".join(letter[v] for v in f.variables) for f in factors)
    expr = f"{inputs}->{''.join(letter[v] for v in out)}"
    table = np.einsum(expr, *[f.table for f in factors], optimize=True) if factors \
        else np.ones(())
    return Factor(tuple(out), np.asarray(table))
