"""Stage-2 axial-current three-particle photon DA correction.

This implements the Colangelo-De Fazio-Ozpineci three-particle DA term F1 in
the axial-current g1 sum rule.  It is a calibration baseline for Stage 2 before
deriving the corresponding tensor-current three-particle projection.

The full twist-4 bracket in the published g1 formula contains

  -1/4 (A - 8 Hbar) (1 + m_c^2/M^2)
  - int int F1.

Our Stage-1 axial script already contains the two-particle part.  This script
adds only the three-particle correction

  Delta_QCD = + e_s <sbar s> exp(-m_c^2/M^2) int int F1,

with the project sign convention for chi left untouched.
"""

from __future__ import annotations

import csv
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
from stage1_axial_g1_baseline import central_inputs, gauss_legendre_integral, g1_stage1


def F1_colangelo(aq, aqb, ag, v):
    """Colangelo F1 combination evaluated from BBK/Rohrwild photon DAs."""
    # photon_da functions take (alpha_q, alpha_qbar); alpha_g is implicit.
    S = pda.S_3p(aq, aqb)
    St = pda.St_3p(aq, aqb)
    T1 = pda.T1_3p(aq, aqb)
    T2 = pda.T2_3p(aq, aqb)
    T3 = pda.T3_3p(aq, aqb)
    T4 = pda.T4_3p(aq, aqb)
    return S + St - T1 - T2 + T3 + T4 + 2.0 * v * (-S - T3 + T2)


def integrate_ag(func, lo, hi, n=120):
    return gauss_legendre_integral(func, lo, hi, n=n)


def integrate_v_scalar(func, lo, hi, n=120):
    """Gauss-Legendre integral for scalar functions with variable inner limits."""
    if hi <= lo:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(n)
    points = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
    values = np.array([func(float(v)) for v in points], dtype=float)
    return float(0.5 * (hi - lo) * np.sum(w * values))


def F1_integral(u0=0.5):
    """Two-domain F1 integral appearing in the axial g1 sum rule."""

    def integrand(v, ag):
        aq = u0 - (1.0 - v) * ag
        aqb = 1.0 - u0 - v * ag
        return F1_colangelo(aq, aqb, ag, v)

    def domain1_v(v):
        upper = u0 / (1.0 - v)
        return integrate_ag(lambda ag: integrand(v, ag), 0.0, upper, n=100)

    def domain2_v(v):
        upper = (1.0 - u0) / v
        return integrate_ag(lambda ag: integrand(v, ag), 0.0, upper, n=100)

    part1 = integrate_v_scalar(domain1_v, 0.0, 1.0 - u0, n=100)
    part2 = integrate_v_scalar(domain2_v, 1.0 - u0, 1.0, n=100)
    return part1 + part2, part1, part2


def prefactor_a(M2, inputs):
    mc = inputs["mc"]
    ms = inputs["ms"]
    return (
        math.exp((inputs["m_ds1"]**2 + inputs["m_ds"]**2) / (2.0 * M2))
        * (mc + ms)
        / (inputs["m_ds1"] * inputs["f_ds1"] * inputs["m_ds"]**2 * inputs["f_ds"])
    )


def width_keV(m_initial, m_final, G):
    alpha = 1.0 / 137.036
    qgamma = (m_initial**2 - m_final**2) / (2.0 * m_initial)
    return alpha / 3.0 * G**2 * qgamma**3 * 1.0e6


def g1_stage2(M2, s0, inputs, f1_total):
    stage1 = g1_stage1(M2, s0, **inputs)
    delta_qcd = (
        inputs["es"]
        * inputs["ss"]
        * math.exp(-(inputs["mc"] ** 2) / M2)
        * f1_total
    )
    delta_g1 = prefactor_a(M2, inputs) * delta_qcd
    g1_full = stage1["g1_GeV_inv"] + delta_g1
    return {
        **stage1,
        "F1_integral": f1_total,
        "F1_delta_qcd": delta_qcd,
        "F1_delta_g1": delta_g1,
        "g1_stage2_GeV_inv": g1_full,
        "width_stage2_keV": width_keV(inputs["m_ds1"], inputs["m_ds"], g1_full),
    }


def main():
    inputs = central_inputs()
    f1_total, f1_part1, f1_part2 = F1_integral(u0=0.5)

    rows = []
    for s0_root in (2.50, 2.55, 2.60):
        s0 = s0_root**2
        for M2 in np.linspace(3.0, 6.0, 13):
            row = g1_stage2(float(M2), s0, inputs, f1_total)
            row["s0_root"] = s0_root
            rows.append(row)

    csv_path = OUT / "stage2_axial_g1_three_particle.csv"
    keys = [
        "M2",
        "s0_root",
        "s0",
        "g1_GeV_inv",
        "g1_stage2_GeV_inv",
        "F1_integral",
        "F1_delta_g1",
        "width_keV",
        "width_stage2_keV",
        "pert",
        "heavy_local",
        "twist2",
        "twist4_2p",
        "twist3_2p",
        "F1_delta_qcd",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in keys})

    g1_s1 = np.array([r["g1_GeV_inv"] for r in rows])
    g1_s2 = np.array([r["g1_stage2_GeV_inv"] for r in rows])
    d_g1 = np.array([r["F1_delta_g1"] for r in rows])
    w_s2 = np.array([r["width_stage2_keV"] for r in rows])
    lines = [
        "Stage-2 axial-current three-particle DA correction",
        "==================================================",
        f"F1 integral total = {f1_total:+.6e}",
        f"  domain 1 = {f1_part1:+.6e}",
        f"  domain 2 = {f1_part2:+.6e}",
        f"Stage-1 g1 range = {g1_s1.min():+.4f} to {g1_s1.max():+.4f} GeV^-1",
        f"F1 delta g1 range = {d_g1.min():+.4f} to {d_g1.max():+.4f} GeV^-1",
        f"Stage-2 axial g1 range = {g1_s2.min():+.4f} to {g1_s2.max():+.4f} GeV^-1",
        f"Stage-2 axial width range = {w_s2.min():.2f} to {w_s2.max():.2f} keV",
        f"Wrote grid to {csv_path}",
    ]
    summary_path = OUT / "stage2_axial_g1_three_particle_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
