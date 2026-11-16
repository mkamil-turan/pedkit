"""One-line summary of the results (no threshold or interpretation is applied)."""

from __future__ import annotations

from .result import ModeInference


def summary_line(c: ModeInference) -> str:
    """List every mode in decreasing order of posterior probability."""
    if c.most_likely is None:
        return "No Mendelian mode has a likelihood above zero."
    return " > ".join(f"{r.mode.value} ({r.posterior:.4g})" for r in c.ranked)
