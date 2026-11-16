"""Transmission rules of the Mendelian modes.

Every individual carries at most two allele slots:
  M  allele from the mother  (autosomal, X and MT modes)
  P  allele from the father  (autosomal; for X only in daughters; for Y only in sons)
Allele value 1 is the disease allele, 0 the wild type.

The functions in this module work both on scalars (exact computation) and
on numpy arrays (Monte Carlo), so both methods share one genetic model.
"""

from __future__ import annotations

import numpy as np

from ..model import Genotype, Mode

MENDELIAN_MODES = (Mode.AD, Mode.AR, Mode.ACD, Mode.XD, Mode.XR, Mode.Y, Mode.MT)
AUTOSOMAL = (Mode.AD, Mode.AR, Mode.ACD)

# Genotype codes: biallelic dosage 0-2, haploid 3+allele, no allele 5
CODE_TEXT = {0: "0/0", 1: "0/1", 2: "1/1", 3: "0", 4: "1", 5: "-"}
N_CODES = 6


def sex_linked(mode: Mode) -> bool:
    return mode in (Mode.XD, Mode.XR, Mode.Y)


def maternal_slots(mode: Mode) -> tuple[str, ...]:
    """Slots the mother can transmit (one chosen with equal probability)."""
    if mode in AUTOSOMAL or mode in (Mode.XD, Mode.XR):
        return ("M", "P")
    if mode is Mode.MT:
        return ("M",)
    return ()


def paternal_slots(mode: Mode) -> tuple[str, ...]:
    """Slots the father can transmit: one of two autosomal copies, his single
    X (the maternal M slot) for X-linked modes, the Y chromosome (P) for Y."""
    if mode in AUTOSOMAL:
        return ("M", "P")
    if mode in (Mode.XD, Mode.XR):
        return ("M",)
    if mode is Mode.Y:
        return ("P",)
    return ()


def has_m_slot(mode: Mode) -> bool:
    return bool(maternal_slots(mode))


def has_p_slot(mode: Mode, male):
    """Is an allele received from the father? ``male`` is a bool or a bool array."""
    if mode in AUTOSOMAL:
        return np.ones_like(male, dtype=bool) if isinstance(male, np.ndarray) else True
    if mode in (Mode.XD, Mode.XR):
        return ~male if isinstance(male, np.ndarray) else not male
    if mode is Mode.Y:
        return male
    return np.zeros_like(male, dtype=bool) if isinstance(male, np.ndarray) else False


def at_risk(mode: Mode, m, p, male):
    """Is the genotype predisposed to the disease (penetrance applies to it)?
    Under codominance any expression (intermediate or full) counts."""
    if mode in (Mode.AD, Mode.ACD):
        return m | p
    if mode is Mode.AR:
        return m & p
    if mode is Mode.XD:
        return np.where(male, m, m | p)
    if mode is Mode.XR:
        return np.where(male, m, m & p)
    if mode is Mode.Y:
        return np.where(male, p, 0)
    if mode is Mode.MT:
        return m
    raise ValueError(f"Not a Mendelian mode: {mode}")


def phenotype_probability(mode: Mode, m, p, male, observed: int,
                          penetrance: float, phenocopy: float):
    """P(observed phenotype | genotype). ``observed``: 1 unaffected,
    2 affected (full), 3 intermediate; 0 (unknown) returns 1.

    Codominant (ACD): 0/0 -> 1, 0/1 -> 3 (intermediate), 1/1 -> 2; the
    penetrance is the probability of the expected expression, the phenocopy
    rate a full phenotype in the wild type. In the other modes the
    phenotype is binary and an intermediate observation is impossible.
    """
    if observed == 0:
        return np.ones_like(np.asarray(m), dtype=float) if isinstance(m, np.ndarray) else 1.0
    if mode is Mode.ACD:
        dose = np.asarray(m, dtype=np.int64) + np.asarray(p, dtype=np.int64)
        if observed == 1:
            return np.where(dose == 0, 1 - phenocopy, 1 - penetrance)
        if observed == 2:
            return np.where(dose == 2, penetrance, np.where(dose == 0, phenocopy, 0.0))
        return np.where(dose == 1, penetrance, 0.0)
    r = np.asarray(at_risk(mode, m, p, male)).astype(bool)
    if observed == 1:
        return np.where(r, 1 - penetrance, 1 - phenocopy)
    if observed == 2:
        return np.where(r, penetrance, phenocopy)
    return np.zeros(r.shape)


def code(mode: Mode, m, p, male):
    """Code of a simulated genotype (see ``CODE_TEXT``)."""
    m_ok = has_m_slot(mode)
    p_ok = has_p_slot(mode, male)
    m = np.asarray(m, dtype=np.int64)
    p = np.asarray(p, dtype=np.int64)
    if m_ok:
        return np.where(p_ok, m + p, 3 + m)
    return np.where(p_ok, 3 + p, 5)


def observed_code(g: Genotype | None) -> int:
    """Code of a test result; -1 if untested."""
    if g is None:
        return -1
    return g.mutant_count if g.biallelic else 3 + g.alleles[0]
