"""Mode-of-inheritance inference (requires numpy).

    from pedkit import read
    from pedkit.inference import infer_mode
    result = infer_mode(read("family.pedx"), "D1", n=1_000_000, seed=1)
    print(result.report())
"""

from .decision import summary_line
from .exact import (
    MAX_ENUMERATION, enumeration_likelihood, infer_mode_exact, peeling_likelihood,
)
from .modes import MENDELIAN_MODES
from .montecarlo import infer_mode
from .prepare import Prepared, prepare
from .result import (
    GenotypeConfiguration, GenotypePool, ModeInference, ModeResult,
)

__all__ = ["infer_mode", "infer_mode_exact", "peeling_likelihood", "enumeration_likelihood",
           "summary_line", "MENDELIAN_MODES", "MAX_ENUMERATION", "ModeInference", "ModeResult",
           "GenotypePool", "GenotypeConfiguration", "prepare", "Prepared"]
