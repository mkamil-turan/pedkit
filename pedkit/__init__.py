"""pedkit: pedigrees in the PEDX format.

Reading/writing, validation and NSGC-standard drawing need only the standard
library. Inference (``pedkit.inference``), Bayes risk (``pedkit.risk``) and
multifactorial risk (``pedkit.mf``) need numpy; the plots need matplotlib
and the viewer PyQt5.

    from pedkit import read, validate
    ped = read("family.pedx")
    for finding in validate(ped):
        print(finding)
"""

__version__ = "0.1.0"

from .io import fam_text, read, read_text, write, write_fam, write_text
from .model import (
    Adoption, Disease, FamilyInfo, Finding, Genotype, Individual, Kind, Level, Mode, Pattern,
    Pedigree, Phenotype, Question, Role, Sex, Twin, TwinType, Union, UnionStatus,
)
from .validate import errors, validate

__all__ = [
    "__version__", "read", "read_text", "write", "write_text", "write_fam", "fam_text",
    "validate", "errors", "Pedigree", "Individual", "Disease", "Union", "Question",
    "FamilyInfo", "Finding", "Genotype", "Twin", "Sex", "Phenotype", "Mode", "Pattern", "Kind",
    "Role", "TwinType", "Adoption", "Level", "UnionStatus",
]
