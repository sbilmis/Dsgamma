"""Split the published axial F1 three-particle DA combination.

The Colangelo axial Stage-2 term uses

    F1 = S + St - T1 - T2 + T3 + T4 + 2 v (-S - T3 + T2).

For the tensor-current Stage-2 derivation, this script separates the known
axial result into

    F1_sigma_like  = S + St - T1 - T2 + T3 + T4,
    F1_xgamma_like = 2 v (-S - T3 + T2),

and integrates both pieces over the standard two domains.  This does not yet
produce the final tensor x_alpha gamma_beta contribution; instead it provides
the axial calibration target that the tensor derivative kernels must reproduce
when current B is replaced by current A.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import photon_da as pda
from stage1_axial_g1_baseline import gauss_legendre_integral
from stage2_axial_g1_three_particle import F1_integral


def f1_sigma_like(aq, aqb):
    return (
        pda.S_3p(aq, aqb)
        + pda.St_3p(aq, aqb)
        - pda.T1_3p(aq, aqb)
        - pda.T2_3p(aq, aqb)
        + pda.T3_3p(aq, aqb)
        + pda.T4_3p(aq, aqb)
    )


def f1_xgamma_like(aq, aqb, v):
    return 2.0 * v * (
        -pda.S_3p(aq, aqb)
        - pda.T3_3p(aq, aqb)
        + pda.T2_3p(aq, aqb)
    )


def integrate_piece(piece, u0=0.5, n=100):
    def integrand(v, ag):
        aq = u0 - (1.0 - v) * ag
        aqb = 1.0 - u0 - v * ag
        return piece(aq, aqb, v)

    def integrate_v_scalar(func, lo, hi):
        if hi <= lo:
            return 0.0
        x, w = np.polynomial.legendre.leggauss(n)
        points = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
        values = np.array([func(float(v)) for v in points], dtype=float)
        return float(0.5 * (hi - lo) * np.sum(w * values))

    def domain1_v(v):
        upper = u0 / (1.0 - v)
        return gauss_legendre_integral(lambda ag: integrand(v, ag), 0.0, upper, n=n)

    def domain2_v(v):
        upper = (1.0 - u0) / v
        return gauss_legendre_integral(lambda ag: integrand(v, ag), 0.0, upper, n=n)

    part1 = integrate_v_scalar(domain1_v, 0.0, 1.0 - u0)
    part2 = integrate_v_scalar(domain2_v, 1.0 - u0, 1.0)
    return part1 + part2, part1, part2


def main():
    total, total_d1, total_d2 = F1_integral(u0=0.5)
    sigma_total, sigma_d1, sigma_d2 = integrate_piece(
        lambda aq, aqb, v: f1_sigma_like(aq, aqb), u0=0.5
    )
    xg_total, xg_d1, xg_d2 = integrate_piece(f1_xgamma_like, u0=0.5)
    residual = total - sigma_total - xg_total

    lines = [
        "Axial F1 split into sigma-like and xgamma-like pieces",
        "====================================================",
        f"F1 total              = {total:+.10e}",
        f"  domain 1            = {total_d1:+.10e}",
        f"  domain 2            = {total_d2:+.10e}",
        f"F1 sigma-like         = {sigma_total:+.10e}",
        f"  domain 1            = {sigma_d1:+.10e}",
        f"  domain 2            = {sigma_d2:+.10e}",
        f"F1 xgamma-like        = {xg_total:+.10e}",
        f"  domain 1            = {xg_d1:+.10e}",
        f"  domain 2            = {xg_d2:+.10e}",
        f"closure residual      = {residual:+.3e}",
        "",
        "Interpretation:",
        "  The xgamma-like piece is the axial calibration target for the",
        "  derivative-reduced x_alpha G^{alpha beta} gamma_beta kernels.",
        "  The tensor-current version must be matched to this normalization",
        "  before it can be added numerically.",
    ]
    path = OUT / "step10_stage2_f1_split.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
