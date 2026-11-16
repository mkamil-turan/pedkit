"""Pedigree validator.

``validate`` collects every finding and does not stop at the first error.
Syntax errors found while reading (``source.read_findings``) are included.

Rules (docs/format.md):
  E1 format line/header missing        E6 loop in the pedigree (individual is own ancestor)
  E2 required column missing/misplaced E7 disease <-> column mismatch
  E3 duplicate individual              E8 value outside the allowed set
  E4 reference to a missing individual E9 biallelic genotype of a male for an X-linked disease
  E5 contradictory parent sex
  W1 MZ twins differ                   W3 question refers to a missing individual/disease
  W2 genotype contradicts the mode     W4 only one parent given
"""

from __future__ import annotations

from .io._fields import GENOTYPE_PREFIX, column_kind
from .model import (
    Finding, Genotype, Individual, Level, Mode, Pedigree, Phenotype, Sex, TwinType,
)


def validate(ped: Pedigree) -> list[Finding]:
    findings: list[Finding] = list(ped.source.read_findings) if ped.source else []
    for check in (_duplicates, _missing_individuals, _parent_sex, _loops, _disease_columns,
                  _genotype_ploidy, _mz_twins, _mode_contradiction, _questions,
                  _single_parent):
        findings.extend(check(ped))
    return sorted(findings, key=lambda f: (f.line is None, f.line or 0))


def errors(findings: list[Finding]) -> list[Finding]:
    return [f for f in findings if f.level is Level.ERROR]


def _error(rule, message, line=None, individual=None) -> Finding:
    return Finding(Level.ERROR, rule, message, line, individual)


def _warning(rule, message, line=None, individual=None) -> Finding:
    return Finding(Level.WARNING, rule, message, line, individual)


# -- ERROR rules ----------------------------------------------------------------

def _duplicates(ped: Pedigree):
    seen: dict[tuple[str, str], Individual] = {}
    for ind in ped.individuals:
        key = (ind.family, ind.id)
        if key in seen:
            first = seen[key].line
            yield _error("E3", f"ID '{ind.id}' is repeated in family '{ind.family}' "
                               f"(first defined on line {first}).", ind.line, ind.id)
        else:
            seen[key] = ind


def _missing_individuals(ped: Pedigree):
    for ind in ped.individuals:
        for role, pid in (("father", ind.father), ("mother", ind.mother)):
            if pid is not None and ped.individual(pid, ind.family) is None:
                yield _error("E4", f"'{pid}', given as {role}, is not in the file.",
                             ind.line, ind.id)
    for u in ped.unions:
        for who in (u.p1, u.p2):
            if who and ped.individual(who) is None:
                yield _error("E4", f"'{who}' in a ##union line is not in the file.", u.line, who)


def _parent_sex(ped: Pedigree):
    for ind in ped.individuals:
        if ind.father is not None and ind.father == ind.mother:
            yield _error("E5", f"'{ind.father}' is given as both father and mother.",
                         ind.line, ind.id)
    for p in ped.individuals:
        children = ped.children(p)
        as_father = [c.id for c in children if c.father == p.id]
        as_mother = [c.id for c in children if c.mother == p.id]
        if as_father and p.sex is Sex.FEMALE:
            yield _error("E5", f"'{p.id}' is female (sex=2) but is the father of: "
                               f"{', '.join(as_father)}.", p.line, p.id)
        if as_mother and p.sex is Sex.MALE:
            yield _error("E5", f"'{p.id}' is male (sex=1) but is the mother of: "
                               f"{', '.join(as_mother)}.", p.line, p.id)


def _loops(ped: Pedigree):
    loops: list[list[Individual]] = []
    ancestors = {id(i): ped.ancestors(i) for i in ped.individuals}
    for ind in ped.individuals:
        if ind.id not in ancestors[id(ind)]:
            continue
        if any(ind in loop for loop in loops):
            continue
        # individuals in the same loop as ind: ancestors of ind whose ancestor is ind
        members = [x for x in ped.individuals if x.family == ind.family and
                   (x is ind or (x.id in ancestors[id(ind)] and ind.id in ancestors[id(x)]))]
        loops.append(members)
        yield _error("E6", "The pedigree contains a loop; these individuals are their own "
                           f"ancestors: {', '.join(m.id for m in members)}.", ind.line, ind.id)


def _disease_columns(ped: Pedigree):
    if not ped.columns:
        return
    header = ped.source.header_line if ped.source else None
    ids = {d.id for d in ped.diseases}
    for d in ped.diseases:
        if d.id and d.id not in ped.columns:
            yield _error("E7", f"Disease '{d.id}' has no phenotype column.", d.line)
    for c in ped.columns:
        kind = column_kind(c, ids)
        if kind == "phenotype" and c not in ids:
            yield _error("E7", f"Column '{c}' has no ##disease line.", header)
        elif kind == "genotype" and c[len(GENOTYPE_PREFIX):] not in ids:
            yield _error("E7", f"Column '{c}' has no ##disease line.", header)


def _genotype_ploidy(ped: Pedigree):
    """E8/E9: is the number of alleles consistent with the mode and the sex?"""
    for ind in ped.individuals:
        for d in ped.diseases:
            g = ind.genotype(d.id)
            if g is None:
                continue
            col = GENOTYPE_PREFIX + d.id
            if d.mode.x_linked and ind.sex is Sex.MALE and g.biallelic:
                yield _error("E9", f"{d.id} is X-linked ({d.mode.value}); '{col}' of a male "
                                   f"must be haploid (0 or 1), found '{g}'.", ind.line, ind.id)
            elif not g.biallelic and not _haploid_allowed(d.mode, ind.sex):
                yield _error("E8", f"Haploid value '{g}' in column '{col}' is not valid for "
                                   f"this individual; use 0/0, 0/1 or 1/1.", ind.line, ind.id)
            elif g.biallelic and d.mode in (Mode.Y, Mode.MT):
                yield _error("E8", f"{d.id} has mode {d.mode.value}; '{col}' must be haploid, "
                                   f"found '{g}'.", ind.line, ind.id)


def _haploid_allowed(mode: Mode, sex: Sex) -> bool:
    if mode is Mode.MT:
        return True
    if mode in (Mode.XD, Mode.XR, Mode.Y, Mode.UNKNOWN):
        return sex is Sex.MALE
    return False


# -- WARNING rules --------------------------------------------------------------

def _mz_twins(ped: Pedigree):
    for members in ped.mz_groups():
        first = members[0]
        for other in members[1:]:
            if other.sex != first.sex:
                yield _warning("W1", f"MZ twins '{first.id}' and '{other.id}' differ in sex.",
                               other.line, other.id)
            for d in ped.diseases:
                g1, g2 = first.genotype(d.id), other.genotype(d.id)
                if g1 is not None and g2 is not None and g1 != g2:
                    yield _warning("W1", f"MZ twins '{first.id}' ({g1}) and '{other.id}' ({g2}) "
                                         f"have different {d.id} genotypes.", other.line, other.id)


def _mode_contradiction(ped: Pedigree):
    for d in ped.diseases:
        if d.mode in (Mode.UNKNOWN, Mode.MF):
            continue
        for ind in ped.individuals:
            g = ind.genotype(d.id)
            p = ind.phenotype(d.id)
            if g is None or p is Phenotype.UNKNOWN:
                continue
            if contradicts_mode(d.mode, ind.sex, g, p, d.penetrance, d.phenocopy):
                state = {Phenotype.AFFECTED: "affected", Phenotype.UNAFFECTED: "unaffected",
                         Phenotype.INTERMEDIATE: "intermediate"}[p]
                yield _warning("W2", f"Contradicts mode {d.mode.value} of {d.id}: {state} "
                                     f"individual with genotype {g}.", ind.line, ind.id)


def expected_phenotype(mode: Mode, sex: Sex, genotype: Genotype) -> Phenotype | None:
    """Phenotype a tested genotype produces under ``mode`` with full penetrance
    and no phenocopies; ``None`` if the combination is undefined (e.g. a Y
    genotype in a female or a haploid autosomal genotype)."""
    m = genotype.mutant_count
    if mode.autosomal:
        if not genotype.biallelic:
            return None
        if mode is Mode.AD:
            return Phenotype.AFFECTED if m >= 1 else Phenotype.UNAFFECTED
        if mode is Mode.AR:
            return Phenotype.AFFECTED if m == 2 else Phenotype.UNAFFECTED
        return (Phenotype.UNAFFECTED, Phenotype.INTERMEDIATE, Phenotype.AFFECTED)[m]
    if mode.x_linked:
        if not genotype.biallelic:          # hemizygous male
            return Phenotype.AFFECTED if m == 1 else Phenotype.UNAFFECTED
        if mode is Mode.XD:
            return Phenotype.AFFECTED if m >= 1 else Phenotype.UNAFFECTED
        return Phenotype.AFFECTED if m == 2 else Phenotype.UNAFFECTED
    if mode is Mode.Y:
        if sex is Sex.FEMALE or genotype.biallelic:
            return None
        return Phenotype.AFFECTED if m == 1 else Phenotype.UNAFFECTED
    if mode is Mode.MT:
        if genotype.biallelic:
            return None
        return Phenotype.AFFECTED if m == 1 else Phenotype.UNAFFECTED
    return None


def contradicts_mode(mode: Mode, sex: Sex, genotype: Genotype, phenotype: Phenotype,
                     penetrance: float, phenocopy: float) -> bool:
    """Is a tested genotype incompatible with the observed phenotype under ``mode``?

    Incomplete penetrance (``penetrance < 1``) can explain an unaffected
    individual with a disease genotype, and phenocopies (``phenocopy > 0``)
    an affected individual with a wild-type genotype; these are not counted
    as contradictions. An intermediate phenotype is defined only for the
    codominant mode, and under it only the heterozygote is intermediate.
    """
    if phenotype is Phenotype.UNKNOWN:
        return False
    expected = expected_phenotype(mode, sex, genotype)
    if expected is None:
        return False  # ploidy problems are reported by E8/E9
    if phenotype == expected:
        return False
    if phenotype is Phenotype.INTERMEDIATE and mode is not Mode.ACD:
        return True
    if phenotype is Phenotype.UNAFFECTED:
        # a disease genotype that did not penetrate
        return penetrance >= 1.0
    if expected is Phenotype.UNAFFECTED:
        # an affected (or intermediate) wild-type: only a phenocopy explains it
        return phenocopy <= 0.0
    # observed affected vs expected intermediate, or the reverse (codominance)
    return True


def _questions(ped: Pedigree):
    for q in ped.questions:
        if q.individual is not None and ped.individual(q.individual) is None:
            yield _warning("W3", f"Question {q.id or ''} refers to a missing individual: "
                                 f"'{q.individual}'.", q.line, q.individual)
        if q.disease is not None and ped.disease(q.disease) is None:
            yield _warning("W3", f"Question {q.id or ''} refers to a missing disease: "
                                 f"'{q.disease}'.", q.line)


def _single_parent(ped: Pedigree):
    for ind in ped.individuals:
        if (ind.father is None) != (ind.mother is None):
            missing = "father" if ind.father is None else "mother"
            yield _warning("W4", f"Only one parent is given ({missing} is 0).", ind.line, ind.id)
