"""Liability-threshold (Falconer) model: parameters and pedigree preparation.

Each individual's liability is L = A + E, on the standard scale L ~ N(0, 1):
  A  heritable (additive polygenic) part, variance h²,
  E  environmental part, variance 1 - h², independent between individuals.
An individual whose L exceeds the threshold T = Φ⁻¹(1 - K) is affected.
The covariance of the A values of two relatives is 2φ·h² (φ the kinship
coefficient); MZ twins share A (2φ = 1). With sex-specific prevalences a
separate threshold is used for males and females (Carter effect: affected
individuals of the less often affected sex carry a higher liability and
pass more risk to their relatives).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..inference.prepare import Prepared, prepare
from ..model import Disease, Pedigree
from .normal import threshold


@dataclass(frozen=True)
class MFParameters:
    prevalence: float          # K (sex-independent or the average)
    heritability: float        # h²
    prevalence_male: float
    prevalence_female: float

    @property
    def threshold(self) -> float:
        return threshold(self.prevalence)

    @property
    def threshold_male(self) -> float:
        return threshold(self.prevalence_male)

    @property
    def threshold_female(self) -> float:
        return threshold(self.prevalence_female)

    @property
    def sex_specific(self) -> bool:
        return self.prevalence_male != self.prevalence_female

    def individual_threshold(self, sex: int) -> float:
        """An individual of unknown sex gets the general threshold."""
        if sex == 1:
            return self.threshold_male
        if sex == 2:
            return self.threshold_female
        return self.threshold

    def replace(self, **fields) -> MFParameters:
        d = dict(prevalence=self.prevalence, heritability=self.heritability,
                 prevalence_male=self.prevalence_male, prevalence_female=self.prevalence_female)
        d.update(fields)
        if "prevalence" in fields and not self.sex_specific:
            d["prevalence_male"] = d["prevalence_female"] = d["prevalence"]
        return MFParameters(**d)


def parameters(disease: Disease, prevalence: float | None = None,
               heritability: float | None = None, prevalence_male: float | None = None,
               prevalence_female: float | None = None) -> MFParameters:
    """Build the parameters from the ``##disease`` fields of the file; given
    values override those of the file."""
    h2 = heritability if heritability is not None else disease.heritability
    km = prevalence_male if prevalence_male is not None else disease.prevalence_male
    kf = prevalence_female if prevalence_female is not None else disease.prevalence_female
    k = prevalence if prevalence is not None else disease.prevalence
    if km is not None and kf is not None:
        k = k if k is not None else (km + kf) / 2
    elif k is not None:
        km = km if km is not None else k
        kf = kf if kf is not None else k
    else:
        raise ValueError(f"No prevalence for {disease.id}: add 'prevalence=' (or "
                         "prevalence_male/prevalence_female) to the ##disease line or "
                         "pass --prevalence.")
    if h2 is None:
        raise ValueError(f"No heritability for {disease.id}: add 'heritability=' to the "
                         "##disease line or pass --heritability.")
    if not 0 <= h2 < 1:
        raise ValueError(f"Heritability must be in [0, 1): {h2}")
    for name, value in (("prevalence", k), ("prevalence_male", km), ("prevalence_female", kf)):
        if not 0 < value < 1:
            raise ValueError(f"{name} must be between 0 and 1: {value}")
    return MFParameters(k, h2, km, kf)


def prepare_mf(ped: Pedigree, disease: str | Disease | None = None) -> Prepared:
    h = prepare(ped, disease)
    if any(f == 3 for f in h.phenotype):
        raise ValueError("The intermediate phenotype (3) is not defined in the "
                         "liability-threshold model.")
    return h


def kinship_matrix(h: Prepared) -> np.ndarray:
    """Kinship coefficients φ; 2φ is the additive relationship matrix.
    Parents are processed before children (``Prepared.order``). An MZ twin
    has with its co-twin the kinship of an individual with itself."""
    n = h.n
    phi = np.zeros((n, n))
    done: list[int] = []
    for i in h.order:
        a, b = h.mother[i], h.father[i]
        phi[i, i] = 0.5 * (1 + (phi[a, b] if a >= 0 and b >= 0 else 0.0))
        for j in done:
            value = 0.5 * ((phi[j, a] if a >= 0 else 0.0) + (phi[j, b] if b >= 0 else 0.0))
            phi[i, j] = phi[j, i] = value
        if h.mz[i] >= 0:
            t = h.mz[i]
            phi[i, t] = phi[t, i] = phi[t, t]
        done.append(i)
    return phi


def target_index(h: Prepared, ped: Pedigree, individual: str) -> int:
    if individual in h.labels:
        return h.labels.index(individual)
    candidates = [i for i, ind in enumerate(ped.individuals) if ind.id == individual]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise KeyError(f"'{individual}' is not in the pedigree.")
    raise KeyError(f"'{individual}' occurs in more than one family; write it as 'FAMILY/ID'.")
