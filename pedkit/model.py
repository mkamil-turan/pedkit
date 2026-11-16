"""pedkit data model.

The classes in this module are independent of the file format: the
inference engine, the risk calculations and the renderer all work on them
directly. Parser-specific information (line numbers, raw text) lives only in
the optional ``source`` fields, which are ``None`` for a pedigree built in
memory.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, IntEnum
from typing import Any, Iterator


# --------------------------------------------------------------------------
# Enumerations
# --------------------------------------------------------------------------

class Sex(IntEnum):
    UNKNOWN = 0
    MALE = 1
    FEMALE = 2


class Phenotype(IntEnum):
    UNKNOWN = 0
    UNAFFECTED = 1
    AFFECTED = 2
    INTERMEDIATE = 3  # distinct heterozygote expression under codominance


class Mode(str, Enum):
    """Mode of inheritance."""

    AD = "AD"    # autosomal dominant
    AR = "AR"    # autosomal recessive
    ACD = "ACD"  # autosomal codominant
    XD = "XD"    # X-linked dominant
    XR = "XR"    # X-linked recessive
    Y = "Y"      # Y-linked
    MT = "MT"    # mitochondrial
    MF = "MF"    # multifactorial (liability-threshold model)
    UNKNOWN = "unknown"

    @property
    def x_linked(self) -> bool:
        return self in (Mode.XD, Mode.XR)

    @property
    def autosomal(self) -> bool:
        return self in (Mode.AD, Mode.AR, Mode.ACD)


class Pattern(str, Enum):
    """Fill pattern of an affected symbol."""

    SOLID = "solid"
    HATCHED = "hatched"
    DOTTED = "dotted"
    QUADRANT = "quadrant"


class Kind(str, Enum):
    """Kind of individual (born individual, pregnancy, miscarriage, ...)."""

    BIRTH = "birth"
    PLANNED = "planned"
    PREGNANCY = "pregnancy"
    MISCARRIAGE = "miscarriage"
    STILLBIRTH = "stillbirth"


class Role(str, Enum):
    PROBAND = "P"
    CONSULTAND = "C"


class TwinType(str, Enum):
    MZ = "MZ"
    DZ = "DZ"
    UNKNOWN = "?"


class Adoption(str, Enum):
    """``in``: adopted into the family; ``out``: placed for adoption."""

    IN = "in"
    OUT = "out"


class Level(str, Enum):
    ERROR = "ERROR"
    WARNING = "WARNING"


class UnionStatus(str, Enum):
    MARRIED = "married"
    SEPARATED = "separated"
    CHILDLESS = "childless"


# --------------------------------------------------------------------------
# Small value objects
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Genotype:
    """Tested genotype. ``0`` is the wild-type allele, ``1`` the disease allele.

    Either biallelic (``0/1``) or hemizygous/haploid (``1``).
    """

    alleles: tuple[int, ...]

    @classmethod
    def parse(cls, text: str) -> Genotype:
        parts = text.split("/")
        if len(parts) not in (1, 2) or any(p not in ("0", "1") for p in parts):
            raise ValueError(text)
        alleles = tuple(int(p) for p in parts)
        if len(alleles) == 2 and alleles[0] > alleles[1]:
            raise ValueError(text)  # only 0/0, 0/1, 1/1 are valid spellings
        return cls(alleles)

    @property
    def biallelic(self) -> bool:
        return len(self.alleles) == 2

    @property
    def mutant_count(self) -> int:
        return sum(self.alleles)

    @property
    def heterozygous(self) -> bool:
        return self.biallelic and self.alleles[0] != self.alleles[1]

    def __str__(self) -> str:
        return "/".join(str(a) for a in self.alleles)


@dataclass(frozen=True)
class Twin:
    type: TwinType
    group: str

    def __str__(self) -> str:
        return f"{self.type.value}:{self.group}"


@dataclass(eq=False)
class Source:
    """Where a record came from in the file; used only for lossless writing
    and error messages. ``fields``: field name -> (raw text, parsed value)."""

    line: int | None = None
    raw_line: str | None = None
    fields: dict[str, tuple[str, Any]] = field(default_factory=dict)


@dataclass
class Finding:
    """A validation finding."""

    level: Level
    rule: str  # e.g. "E4", "W2"
    message: str
    line: int | None = None
    individual: str | None = None

    def __str__(self) -> str:
        where = f"line {self.line}" if self.line is not None else "line ?"
        who = f" [{self.individual}]" if self.individual else ""
        return f"{self.level.value} {self.rule} ({where}){who}: {self.message}"


# --------------------------------------------------------------------------
# Records
# --------------------------------------------------------------------------

@dataclass(eq=False)
class Individual:
    family: str
    id: str
    father: str | None = None
    mother: str | None = None
    sex: Sex = Sex.UNKNOWN
    phenotypes: dict[str, Phenotype] = field(default_factory=dict)
    genotypes: dict[str, Genotype | None] = field(default_factory=dict)
    role: Role | None = None
    deceased: bool = False
    age: int | None = None
    kind: Kind = Kind.BIRTH
    twin: Twin | None = None
    adoption: Adoption | None = None
    order: int | None = None
    gender: str | None = None
    extra: dict[str, str] = field(default_factory=dict)  # unknown columns
    source: Source | None = field(default=None, repr=False)

    @property
    def line(self) -> int | None:
        return self.source.line if self.source else None

    @property
    def founder(self) -> bool:
        """An individual whose parents are not given in the file."""
        return self.father is None and self.mother is None

    def phenotype(self, disease: str) -> Phenotype:
        return self.phenotypes.get(disease, Phenotype.UNKNOWN)

    def affected(self, disease: str) -> bool:
        return self.phenotype(disease) is Phenotype.AFFECTED

    def genotype(self, disease: str) -> Genotype | None:
        return self.genotypes.get(disease)


@dataclass(eq=False)
class Disease:
    id: str
    name: str | None = None
    mode: Mode = Mode.UNKNOWN
    allele_freq: float | None = None
    penetrance: float = 1.0
    mutation_rate: float = 0.0
    phenocopy: float = 0.0
    pattern: Pattern = Pattern.SOLID
    # Multifactorial (liability-threshold) parameters
    prevalence: float | None = None          # population prevalence K
    heritability: float | None = None        # heritability of liability h²
    prevalence_male: float | None = None     # sex-specific K (Carter effect)
    prevalence_female: float | None = None
    extra: dict[str, str] = field(default_factory=dict)
    source: Source | None = field(default=None, repr=False)

    @property
    def line(self) -> int | None:
        return self.source.line if self.source else None


@dataclass(eq=False)
class Union:
    p1: str
    p2: str
    consanguinity: str | None = None  # "none", "cousin1", "cousin2", "uncle_niece", ...
    status: UnionStatus = UnionStatus.MARRIED
    extra: dict[str, str] = field(default_factory=dict)
    source: Source | None = field(default=None, repr=False)

    @property
    def line(self) -> int | None:
        return self.source.line if self.source else None

    @property
    def consanguineous(self) -> bool:
        return self.consanguinity not in (None, "", "none")


@dataclass(eq=False)
class Question:
    """A ``##question`` record (answered by ``pedkit infer/risk/mf --questions``)."""

    id: str | None = None
    type: str | None = None
    individual: str | None = None
    disease: str | None = None
    event: str | None = None
    extra: dict[str, str] = field(default_factory=dict)
    source: Source | None = field(default=None, repr=False)

    @property
    def line(self) -> int | None:
        return self.source.line if self.source else None


@dataclass(eq=False)
class FamilyInfo:
    id: str | None = None
    description: str | None = None
    extra: dict[str, str] = field(default_factory=dict)
    source: Source | None = field(default=None, repr=False)


@dataclass(eq=False)
class MetaLine:
    """An unrecognised ``##key=value`` line; written back verbatim."""

    key: str
    value: str
    source: Source | None = field(default=None, repr=False)


@dataclass(eq=False)
class FileSource:
    """Layout of the file that was read: line order, header, read findings."""

    path: str | None = None
    format_line: int | None = None
    header_line: int | None = None
    columns: tuple[str, ...] = ()  # header as read
    # Order of items in the file: ("format",), ("meta", obj), ("header",),
    # ("individual", ind), ("raw", text). Comments/blank lines are "raw".
    layout: list[tuple] = field(default_factory=list)
    final_newline: bool = True
    read_findings: list[Finding] = field(default_factory=list)


# --------------------------------------------------------------------------
# Pedigree
# --------------------------------------------------------------------------

@dataclass(eq=False)
class Pedigree:
    format: str = "PEDX/0.1"
    family: FamilyInfo | None = None
    diseases: list[Disease] = field(default_factory=list)
    unions: list[Union] = field(default_factory=list)
    questions: list[Question] = field(default_factory=list)
    other_meta: list[MetaLine] = field(default_factory=list)
    columns: list[str] = field(default_factory=list)  # column order of the header
    individuals: list[Individual] = field(default_factory=list)
    source: FileSource | None = field(default=None, repr=False)

    # -- lookup -------------------------------------------------------------

    def individual(self, id: str, family: str | None = None) -> Individual | None:
        for ind in self.individuals:
            if ind.id == id and (family is None or ind.family == family):
                return ind
        return None

    def __iter__(self) -> Iterator[Individual]:
        return iter(self.individuals)

    def __len__(self) -> int:
        return len(self.individuals)

    def disease(self, id: str) -> Disease | None:
        return next((d for d in self.diseases if d.id == id), None)

    @property
    def primary_disease(self) -> Disease | None:
        return self.diseases[0] if self.diseases else None

    @property
    def families(self) -> list[str]:
        return list(dict.fromkeys(ind.family for ind in self.individuals))

    # -- relationships ------------------------------------------------------

    def parents(self, ind: Individual) -> list[Individual]:
        return [p for p in (self.individual(ind.father, ind.family) if ind.father else None,
                            self.individual(ind.mother, ind.family) if ind.mother else None)
                if p is not None]

    def children(self, ind: Individual) -> list[Individual]:
        return [c for c in self.individuals
                if c.family == ind.family and ind.id in (c.father, c.mother)]

    def spouses(self, ind: Individual) -> list[Individual]:
        """Individuals with a common child or linked by a ``##union`` line."""
        ids: dict[str, None] = {}
        for c in self.children(ind):
            other = c.mother if c.father == ind.id else c.father
            if other:
                ids[other] = None
        for u in self.unions:
            if u.p1 == ind.id:
                ids[u.p2] = None
            elif u.p2 == ind.id:
                ids[u.p1] = None
        return [s for i in ids if (s := self.individual(i, ind.family)) is not None]

    def sibships(self, family: str | None = None
                 ) -> dict[tuple[str | None, str | None], list[Individual]]:
        """(father, mother) -> children, sorted by ``order``, else file order."""
        groups: dict[tuple[str | None, str | None], list[Individual]] = {}
        for ind in self.individuals:
            if ind.founder or (family is not None and ind.family != family):
                continue
            groups.setdefault((ind.father, ind.mother), []).append(ind)
        for key, members in groups.items():
            groups[key] = _sort(members)
        return groups

    def ancestors(self, ind: Individual) -> set[str]:
        """IDs of all ancestors of an individual (robust to loops)."""
        seen: set[str] = set()
        stack = list(self.parents(ind))
        while stack:
            p = stack.pop()
            if p.id in seen:
                continue
            seen.add(p.id)
            stack.extend(self.parents(p))
        return seen

    def union(self, a: str, b: str) -> Union | None:
        for u in self.unions:
            if {u.p1, u.p2} == {a, b}:
                return u
        return None

    def mz_groups(self, family: str | None = None) -> list[list[Individual]]:
        """Groups of monozygotic twins (two or more members)."""
        groups: dict[tuple[str, str], list[Individual]] = {}
        for ind in self.individuals:
            if ind.twin is not None and ind.twin.type is TwinType.MZ and (
                    family is None or ind.family == family):
                groups.setdefault((ind.family, ind.twin.group), []).append(ind)
        return [g for g in groups.values() if len(g) > 1]


def _sort(members: list[Individual]) -> list[Individual]:
    inf = float("inf")
    return [ind for _, _, ind in sorted(
        ((ind.order if ind.order is not None else inf, i, ind) for i, ind in enumerate(members)),
        key=lambda t: (t[0], t[1]))]
