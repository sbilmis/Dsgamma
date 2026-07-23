"""Authoritative mixing-angle inputs for the Ds/Bs radiative paper.

The nominal values are external inputs from the dedicated QCD-sum-rule
analysis arXiv:2408.08014.  The heavy-quark-limit value and the broad interval
are diagnostics only; they are not statistical alternatives to the nominal
Gaussian inputs.
"""

from __future__ import annotations

import numpy as np


THETA_DS_DEG = 26.6
THETA_DS_SIGMA_DEG = 0.6
THETA_BS_DEG = 38.5
THETA_BS_SIGMA_DEG = 0.1

THETA_HQET_DEG = 35.3
THETA_SENSITIVITY_MIN_DEG = 25.0
THETA_SENSITIVITY_MAX_DEG = 45.0

ANGLE_SOURCE = "Aliev-Bilmis-Savci, arXiv:2408.08014v2"


def sample_gaussian_angle(
    rng: np.random.Generator,
    central_deg: float,
    sigma_deg: float,
    nsigma_clip: float = 4.0,
) -> float:
    """Sample a quoted angle uncertainty with a harmless far-tail clip."""

    low = central_deg - nsigma_clip * sigma_deg
    high = central_deg + nsigma_clip * sigma_deg
    return float(np.clip(rng.normal(central_deg, sigma_deg), low, high))
