"""Standard normal distribution functions (numpy only, no scipy).

Relative accuracy matters in the tails: for a disease with 0.1 % prevalence
the liability of an affected individual is drawn above the threshold
(≈ 3.09), where 1 - Φ(x) is very small. Φ is therefore computed through
the complementary error function (erfc), not directly.
"""

from __future__ import annotations

import numpy as np

_SQRT2 = np.sqrt(2.0)


def erfc(x):
    """Complementary error function; relative error < 1.2e-7 everywhere
    (Numerical Recipes 'erfcc' Chebyshev approximation)."""
    x = np.asarray(x, dtype=float)
    z = np.abs(x)
    t = 1.0 / (1.0 + 0.5 * z)
    y = t * np.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (
        0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (
            1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    return np.where(x >= 0, y, 2.0 - y)


def phi(x):
    """Density φ(x)."""
    x = np.asarray(x, dtype=float)
    return np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)


def Phi(x):
    """Cumulative distribution Φ(x) = P(Z ≤ x)."""
    return 0.5 * erfc(-np.asarray(x, dtype=float) / _SQRT2)


def upper_tail(x):
    """1 - Φ(x), accurate in relative terms also for large x."""
    return 0.5 * erfc(np.asarray(x, dtype=float) / _SQRT2)


# Acklam's inverse normal approximation (relative error < 1.15e-9)
_A = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
      1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
_B = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
      6.680131188771972e+01, -1.328068155288572e+01)
_C = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
      -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
_D = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
      3.754408661907416e+00)


def Phi_inv(p):
    """Φ⁻¹(p), 0 < p < 1."""
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1 - 1e-16)
    out = np.empty_like(p)
    low = p < 0.02425
    high = p > 1 - 0.02425
    mid = ~(low | high)

    q = p[mid] - 0.5
    r = q * q
    out[mid] = ((((((_A[0] * r + _A[1]) * r + _A[2]) * r + _A[3]) * r + _A[4]) * r + _A[5]) * q
                / (((((_B[0] * r + _B[1]) * r + _B[2]) * r + _B[3]) * r + _B[4]) * r + 1))

    def tail(pp):
        q = np.sqrt(-2 * np.log(pp))
        return ((((((_C[0] * q + _C[1]) * q + _C[2]) * q + _C[3]) * q + _C[4]) * q + _C[5])
                / ((((_D[0] * q + _D[1]) * q + _D[2]) * q + _D[3]) * q + 1))

    out[low] = tail(p[low])
    out[high] = -tail(1 - p[high])
    return out


def threshold(prevalence: float) -> float:
    """Threshold on the standard liability scale of a disease with prevalence K."""
    return float(-Phi_inv(np.array([prevalence]))[0])


def truncated_normal(rng: np.random.Generator, mean, sd, lower=None, upper=None):
    """Draw from N(mean, sd²) truncated to [lower, upper]
    (inverse-CDF method, computed on the tail side)."""
    mean = np.asarray(mean, dtype=float)
    u = rng.random(mean.shape)
    if lower is not None:
        za = (lower - mean) / sd
        rest = upper_tail(za)                  # P(Z > za)
        z = -Phi_inv(np.maximum(u * rest, 1e-300))
        z = np.maximum(z, za)                  # rounding guard
    elif upper is not None:
        zu = (upper - mean) / sd
        z = Phi_inv(np.maximum(u * Phi(zu), 1e-300))
        z = np.minimum(z, zu)
    else:
        z = Phi_inv(u)
    return mean + sd * z
