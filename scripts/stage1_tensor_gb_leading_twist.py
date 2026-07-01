"""Partial Stage-1 tensor-current estimate for Ds1 -> Ds gamma.

This script computes only the leading-twist two-particle photon-DA contribution
to the tensor-current form factor G_B.  It is not the full tensor-current sum
rule yet.

The trace-level result from scripts/step3_soft_2p_e1_projection.wl is

  E1[A current, tensor DA] = 8 (k.q)/(p.q)
  E1[B current, tensor DA] = 8 m_c/(m_c+m_s)

At the two-particle DA saddle after the double Borel transform,
k = p + \bar u q, so k.q = p.q.  Therefore the tensor-current OPE contribution
from chi phi_gamma is m_c/(m_c+m_s) times the axial-current OPE contribution,
before dividing by the corresponding axial/tensor decay constants.

We quote two f_B conventions:
  1. direct:      f_B = f_T
  2. converted:  f_B = f_T * m_Ds1/(m_c+m_s)

The second is motivated by contracting a conventional tensor decay constant
with P^alpha/(m_c+m_s), but it remains a normalization choice to verify.
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

from stage1_axial_g1_baseline import central_inputs, g1_stage1


def gb_twist2_from_axial_twist2(row, inputs, fB):
    """Convert the axial chi-phi OPE term into the tensor-current G_B term."""
    mc = inputs["mc"]
    ms = inputs["ms"]
    M2 = row["M2"]
    m_ds1 = inputs["m_ds1"]
    m_ds = inputs["m_ds"]
    f_ds = inputs["f_ds"]

    qcd_b_twist2 = (mc / (mc + ms)) * row["twist2"]
    prefactor_b = (
        math.exp((m_ds1**2 + m_ds**2) / (2.0 * M2))
        * (mc + ms)
        / (m_ds1 * fB * m_ds**2 * f_ds)
    )
    return prefactor_b * qcd_b_twist2


def width_keV(m_initial, m_final, G):
    alpha = 1.0 / 137.036
    qgamma = (m_initial**2 - m_final**2) / (2.0 * m_initial)
    return alpha / 3.0 * G**2 * qgamma**3 * 1.0e6


def main():
    inputs = central_inputs()
    fT = 0.256
    fB_options = {
        "fB_equal_fT": fT,
        "fB_converted_fT_M_over_msum": fT * inputs["m_ds1"] / (inputs["mc"] + inputs["ms"]),
    }

    theta = math.radians(35.3)

    rows = []
    for s0_root in (2.50, 2.55, 2.60):
        s0 = s0_root**2
        for M2 in np.linspace(3.0, 6.0, 13):
            axial = g1_stage1(float(M2), s0, **inputs)
            GA = axial["g1_GeV_inv"]
            for label, fB in fB_options.items():
                GB = gb_twist2_from_axial_twist2(axial, inputs, fB)
                G2460 = math.sin(theta) * GA + math.cos(theta) * GB
                G2536 = math.cos(theta) * GA - math.sin(theta) * GB
                rows.append({
                    "fB_scheme": label,
                    "M2": float(M2),
                    "s0_root": s0_root,
                    "s0": s0,
                    "GA_stage1_full": GA,
                    "GB_twist2_only": GB,
                    "G2460_using_partial_GB": G2460,
                    "G2536_using_partial_GB": G2536,
                    "Gamma2460_keV_partial": width_keV(inputs["m_ds1"], inputs["m_ds"], G2460),
                    "Gamma2536_keV_partial": width_keV(2.53511, inputs["m_ds"], G2536),
                    "fB": fB,
                })

    csv_path = OUT / "stage1_tensor_gb_leading_twist_partial.csv"
    keys = list(rows[0].keys())
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Partial Stage-1 tensor-current leading-twist estimate",
        "====================================================",
        "Only the chi phi_gamma contribution to G_B is included.",
        "These are not final physical widths.",
    ]
    for scheme in fB_options:
        subset = [r for r in rows if r["fB_scheme"] == scheme]
        gb = np.array([r["GB_twist2_only"] for r in subset])
        g2460 = np.array([r["G2460_using_partial_GB"] for r in subset])
        w2460 = np.array([r["Gamma2460_keV_partial"] for r in subset])
        lines.append(
            f"{scheme}: GB_tw2 {gb.min():+.4f} to {gb.max():+.4f} GeV^-1; "
            f"G2460(partial) {g2460.min():+.4f} to {g2460.max():+.4f}; "
            f"Gamma2460(partial) {w2460.min():.2f} to {w2460.max():.2f} keV"
        )
    lines.append(f"Wrote grid to {csv_path}")
    summary_path = OUT / "stage1_tensor_gb_leading_twist_partial_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
