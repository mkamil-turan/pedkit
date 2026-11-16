"""Turns a pedigree into flat arrays for inference."""

from __future__ import annotations

from dataclasses import dataclass

from ..model import Disease, Pedigree
from ..validate import errors, validate
from .modes import observed_code


@dataclass
class Prepared:
    disease: Disease
    q: float                 # allele frequency
    penetrance: float
    phenocopy: float
    mutation: float
    labels: list[str]        # individual names ("FAMILY/ID" when there are several families)
    mother: list[int]        # index or -1
    father: list[int]
    sex: list[int]           # 0, 1, 2
    phenotype: list[int]     # 0 unknown, 1 unaffected, 2 affected, 3 intermediate
    genotype: list[int]      # observed genotype code or -1
    order: list[int]         # parents before children
    mz: list[int]            # index of the MZ co-twin processed first, or -1

    @property
    def n(self) -> int:
        return len(self.labels)


def resolve_disease(ped: Pedigree, disease: str | Disease | None) -> Disease:
    if isinstance(disease, Disease):
        return disease
    if disease is not None:
        d = ped.disease(disease)
        if d is None:
            raise KeyError(f"Disease '{disease}' is not in the pedigree.")
        return d
    d = ped.primary_disease
    if d is None:
        raise ValueError("No disease is defined in the pedigree.")
    return d


def prepare(ped: Pedigree, disease: str | Disease | None = None,
            default_allele_freq: float = 0.01) -> Prepared:
    errs = errors(validate(ped))
    if errs:
        raise ValueError("The pedigree has errors; inference is not possible:\n"
                         + "\n".join(str(e) for e in errs))
    d = resolve_disease(ped, disease)

    inds = ped.individuals
    multi = len(ped.families) > 1
    index = {(i.family, i.id): k for k, i in enumerate(inds)}

    def find(family: str, pid: str | None) -> int:
        return index.get((family, pid), -1) if pid else -1

    mother = [find(i.family, i.mother) for i in inds]
    father = [find(i.family, i.father) for i in inds]

    # Topological order (validation guarantees there is no loop)
    order: list[int] = []
    placed: set[int] = set()
    while len(order) < len(inds):
        for k in range(len(inds)):
            if k not in placed and all(p < 0 or p in placed for p in (mother[k], father[k])):
                order.append(k)
                placed.add(k)

    # MZ twins share one genotype: each later twin points to the first one in ``order``
    mz = [-1] * len(inds)
    position = {k: r for r, k in enumerate(order)}
    for group in ped.mz_groups():
        ks = sorted((index[(i.family, i.id)] for i in group), key=position.get)
        for k in ks[1:]:
            if mother[k] == mother[ks[0]] and father[k] == father[ks[0]]:
                mz[k] = ks[0]

    q = d.allele_freq if d.allele_freq is not None else default_allele_freq
    return Prepared(
        disease=d, q=q, penetrance=d.penetrance, phenocopy=d.phenocopy, mutation=d.mutation_rate,
        labels=[f"{i.family}/{i.id}" if multi else i.id for i in inds],
        mother=mother, father=father,
        sex=[int(i.sex) for i in inds],
        phenotype=[int(i.phenotype(d.id)) for i in inds],
        genotype=[observed_code(i.genotype(d.id)) for i in inds],
        order=order, mz=mz,
    )
