"""Events of a risk question: the probability of an event given a genotype.

Genotype codes are those of ``pedkit.inference.modes.CODE_TEXT``:
0 "0/0", 1 "0/1", 2 "1/1", 3 "0" (haploid), 4 "1" (haploid), 5 "-" (no allele).
"""

from __future__ import annotations

from ..model import Mode

EVENT_NAME = {
    "carrier": "heterozygous carrier (0/1)",
    "mutant": "carrying at least one disease allele",
    "affected": "being affected",
    "intermediate": "intermediate phenotype (codominant heterozygote)",
}


def event_name(event: str) -> str:
    if event.startswith("genotype:"):
        return f"genotype {event.split(':', 1)[1]}"
    return EVENT_NAME.get(event, event)


def valid_event(event: str) -> bool:
    return event in EVENT_NAME or (event.startswith("genotype:") and
                                   event.split(":", 1)[1] in ("0/0", "0/1", "1/1", "0", "1"))


def expression(mode: Mode, code: int) -> int:
    """Expected expression of a genotype: 0 none, 1 intermediate (ACD only), 2 full."""
    if code == 5:
        return 0
    if code in (3, 4):  # haploid: X in males, Y, MT
        return 2 if code == 4 else 0
    if mode in (Mode.AD, Mode.XD):
        return 2 if code >= 1 else 0
    if mode in (Mode.AR, Mode.XR):
        return 2 if code == 2 else 0
    if mode is Mode.ACD:
        return code  # 0 -> none, 1 -> intermediate, 2 -> full
    return 0  # no biallelic code arises under Y or MT


def event_probability(event: str, mode: Mode, code: int, penetrance: float, phenocopy: float,
                      observed_phenotype: int = 0) -> float:
    """P(event | genotype). For phenotype events an observed phenotype of the
    individual is taken as given (probability 0 or 1)."""
    if event == "carrier":
        return 1.0 if code == 1 else 0.0
    if event == "mutant":
        return 1.0 if code in (1, 2, 4) else 0.0
    if event.startswith("genotype:"):
        target = event.split(":", 1)[1]
        text = {0: "0/0", 1: "0/1", 2: "1/1", 3: "0", 4: "1", 5: "-"}[code]
        return 1.0 if text == target else 0.0
    if event in ("affected", "intermediate"):
        wanted = 2 if event == "affected" else 3
        if observed_phenotype:
            return 1.0 if observed_phenotype == wanted else 0.0
        e = expression(mode, code)
        if event == "intermediate":
            return penetrance if e == 1 else 0.0
        if e == 2:
            return penetrance
        if e == 1:
            return 0.0  # a codominant heterozygote does not show the full phenotype
        return phenocopy
    raise ValueError(f"Unknown event: {event!r}")
