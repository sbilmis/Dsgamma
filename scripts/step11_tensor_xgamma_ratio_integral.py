"""Partial tensor-current matching of the x_alpha G gamma_beta term.

The derivative-reduced xgamma kernels show that the tensor-current result is
not obtained by multiplying the full axial F1 integral by one constant.  For
the structures that do have axial partners, S, T2 and T3, this script applies
the double-pole residue ratios

    B_i / [(m_c/(m_c+m_s)) A_i]

to the axial xgamma combination

    F1_xgamma = 2 v (-S - T3 + T2).

The T4 tensor kernel has no axial partner in this E1 projection
(A_T4 = 0, B_T4 != 0), so it is reported as unresolved and is not added to the
matched integral here.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import photon_da as pda
from stage1_axial_g1_baseline import central_inputs, gauss_legendre_integral
from step10_stage2_f1_split import f1_xgamma_like, integrate_piece


def residue_ratios(a: float, p2: float, pq: float, mc: float) -> dict[str, float]:
    """Double-pole residue ratios from step9 after k=p+a q.

    These are the on-residue ratios for S, T2, and T3.  T4 is absent because
    its axial kernel is zero.
    """
    mc2 = mc * mc
    return {
        "S": 1.0 + (1.0 - a) * pq / mc2,
        "T2": -2.0 + (3.0 * p2 + 2.0 * pq + 3.0 * a * pq) / mc2,
        "T3": (p2 + 2.0 * a * pq) / mc2,
    }


def tensor_xgamma_matched_piece(aq: float, aqb: float, v: float, p2: float, pq: float, mc: float):
    ag = 1.0 - aq - aqb
    a = aqb + v * ag
    ratios = residue_ratios(a, p2, pq, mc)
    return 2.0 * v * (
        -ratios["S"] * pda.S_3p(aq, aqb)
        -ratios["T3"] * pda.T3_3p(aq, aqb)
        +ratios["T2"] * pda.T2_3p(aq, aqb)
    )


def t4_residue_indicator(aq: float, aqb: float, v: float, p2: float, pq: float, mc: float):
    """Dimensionless shape of the unresolved T4 tensor-only residue.

    The derivative reduction gives B_T4 proportional to -(1-a) p.q after
    dropping the single-pole D/D^2 part.  This is not normalized to an axial
    partner, so it is diagnostic only.
    """
    ag = 1.0 - aq - aqb
    a = aqb + v * ag
    return 2.0 * v * (-(1.0 - a) * pq / (mc * mc)) * pda.T4_3p(aq, aqb)


def integrate_xgamma(piece, u0=0.5, n=100):
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


def kinematics_options(inputs: dict[str, float]):
    mc = inputs["mc"]
    physical_p2 = inputs["m_ds"] ** 2
    physical_pq = (inputs["m_ds1"] ** 2 - inputs["m_ds"] ** 2) / 2.0
    a0 = 0.5
    return {
        "hadron_p2_physical_pq": (physical_p2, physical_pq),
        "heavy_residue_p2_physical_pq": (mc * mc - 2.0 * a0 * physical_pq, physical_pq),
    }


def main():
    inputs = central_inputs()
    mc = inputs["mc"]
    axial_total, axial_d1, axial_d2 = integrate_piece(f1_xgamma_like, u0=0.5)

    lines = [
        "Partial tensor xgamma matched integral",
        "=======================================",
        "Matched structures: S, T2, T3.",
        "Unresolved structure: T4 has A_T4 = 0 but B_T4 != 0.",
        f"Axial xgamma calibration integral = {axial_total:+.10e}",
        f"  domain 1 = {axial_d1:+.10e}",
        f"  domain 2 = {axial_d2:+.10e}",
    ]

    for label, (p2, pq) in kinematics_options(inputs).items():
        matched, d1, d2 = integrate_xgamma(
            lambda aq, aqb, v, p2=p2, pq=pq: tensor_xgamma_matched_piece(
                aq, aqb, v, p2, pq, mc
            ),
            u0=0.5,
        )
        t4_diag, t4_d1, t4_d2 = integrate_xgamma(
            lambda aq, aqb, v, p2=p2, pq=pq: t4_residue_indicator(
                aq, aqb, v, p2, pq, mc
            ),
            u0=0.5,
        )
        ratios_mid = residue_ratios(0.5, p2, pq, mc)
        lines.extend(
            [
                "",
                f"Kinematics option: {label}",
                f"  p2 = {p2:+.8e} GeV^2",
                f"  pq = {pq:+.8e} GeV^2",
                "  midpoint residue ratios:",
                f"    R_S  = {ratios_mid['S']:+.8e}",
                f"    R_T2 = {ratios_mid['T2']:+.8e}",
                f"    R_T3 = {ratios_mid['T3']:+.8e}",
                f"  matched S,T2,T3 integral = {matched:+.10e}",
                f"    domain 1 = {d1:+.10e}",
                f"    domain 2 = {d2:+.10e}",
                f"  diagnostic T4 shape integral = {t4_diag:+.10e}",
                f"    domain 1 = {t4_d1:+.10e}",
                f"    domain 2 = {t4_d2:+.10e}",
            ]
        )

    path = OUT / "step11_tensor_xgamma_ratio_integral.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
