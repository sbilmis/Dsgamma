"""Isolated dimension-five completion of the existing two-point matrix.

Adds the heavy-line background-gluon topology to the light-propagator x^2
contribution already in production. Tree level; ms=0 Wilson coefficients,
with the specified d=mQ+ms current normalization retained. No alpha_s,
ms*C5 dimension-six terms, or gluon-only condensates are added here.
"""
from pathlib import Path
import math
import sys

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from twopoint_ds1_matrix_sumrule import Inputs, ope_matrix


def heavy_line_mixed_component(m2: float, inp: Inputs) -> np.ndarray:
    m, d = inp.mc, inp.mc + inp.ms
    factor = inp.mixed_ss * math.exp(-m*m/m2) * (1-m*m/m2)
    return factor * np.array([[0.0, 1/(12*d)], [1/(12*d), m/(6*d*d)]])


def completed_ope_matrix(m2: float, s0: float, inp: Inputs):
    original, pieces = ope_matrix(m2, s0, inp)
    extra = heavy_line_mixed_component(m2, inp)
    return original + extra, pieces | {"d5_heavy_line_gluon": extra}
