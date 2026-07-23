"""Soft two-particle tensor-current estimate for G_B.

This script extends the previous leading-twist-only estimate by using the
trace-level relation for the full two-particle tensor DA family:

  E1[A, tensor DA] = 8 (k.q)/(p.q)
  E1[B, tensor DA] = 8 m_c/(m_c+m_s)

At the two-particle DA Borel saddle k = p + ubar q, so k.q = p.q.
Therefore the tensor-current OPE contribution from the tensor photon DA
structures (chi phi_gamma and the A/h_gamma twist-4 family) is

  Pi_B^tensor-DA = [m_c/(m_c+m_s)] Pi_A^tensor-DA

before hadronic decay-constant normalization.

The vector DA contribution to J_B contains extra p^2/(p.q) kinematics.  With
the explicitly adopted physical-residue numerator prescription
``p^2=m_P^2`` and ``p.q=(m_A^2-m_P^2)/2`` its trace ratio is fixed, so it is
included in the complete two-particle result.  The hard-photon tensor-current
contribution is not included here.
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


def prefactor_b(M2, inputs, fB):
    mc = inputs["mc"]
    ms = inputs["ms"]
    return (
        math.exp((inputs["m_ds1"]**2 + inputs["m_ds"]**2) / (2.0 * M2))
        * (mc + ms)
        / (inputs["m_ds1"] * fB * inputs["m_ds"]**2 * inputs["f_ds"])
    )


def width_keV(m_initial, m_final, G):
    alpha = 1.0 / 137.036
    qgamma = (m_initial**2 - m_final**2) / (2.0 * m_initial)
    return alpha / 3.0 * G**2 * qgamma**3 * 1.0e6


def gb_soft_tensor_da(axial_row, inputs, fB):
    """Tensor-current G_B from tensor-DA two-particle terms only."""
    mc = inputs["mc"]
    ms = inputs["ms"]
    qcd_b = (mc / (mc + ms)) * (axial_row["twist2"] + axial_row["twist4_2p"])
    return prefactor_b(axial_row["M2"], inputs, fB) * qcd_b


def gb_vector_da(axial_row, inputs, fB, ubar=0.5):
    """Vector-DA B-current term from the physical-residue trace ratio.

    The relation is estimated by comparing the B-current vector-DA E1 trace to
    the A-current vector-DA trace and then multiplying the known axial
    f3gamma*Psi^v OPE term.  The B trace contains p^2/(p.q), evaluated with
    the stated physical-residue numerator prescription.
    """
    mc = inputs["mc"]
    ms = inputs["ms"]
    p2 = inputs["m_ds"]**2
    pq = (inputs["m_ds1"]**2 - inputs["m_ds"]**2) / 2.0
    ratio_trace = (3.0 * p2 + (2.0 * ubar + 3.0) * pq) / (3.0 * mc * (mc + ms))
    qcd_b = ratio_trace * axial_row["twist3_2p"]
    return prefactor_b(axial_row["M2"], inputs, fB) * qcd_b


def gb_vector_da_diagnostic(axial_row, inputs, fB, ubar=0.5):
    """Backward-compatible alias for older audit scripts."""
    return gb_vector_da(axial_row, inputs, fB, ubar=ubar)


def gb_soft_two_particle(axial_row, inputs, fB):
    """Complete tensor-current two-particle photon-DA contribution."""
    return gb_soft_tensor_da(axial_row, inputs, fB) + gb_vector_da(
        axial_row, inputs, fB
    )


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
                GB_tensor = gb_soft_tensor_da(axial, inputs, fB)
                GB_vec_diag = gb_vector_da(axial, inputs, fB)
                for include_vec in (False, True):
                    GB = GB_tensor + (GB_vec_diag if include_vec else 0.0)
                    G2460 = math.sin(theta) * GA + math.cos(theta) * GB
                    G2536 = math.cos(theta) * GA - math.sin(theta) * GB
                    rows.append({
                        "fB_scheme": label,
                        "include_vector_da_diagnostic": include_vec,
                        "M2": float(M2),
                        "s0_root": s0_root,
                        "s0": s0,
                        "GA_stage1_full": GA,
                        "GB_tensor_DA_soft2p": GB_tensor,
                        "GB_vector_DA_diagnostic": GB_vec_diag,
                        "GB_soft2p_estimate": GB,
                        "G2460_partial": G2460,
                        "G2536_partial": G2536,
                        "Gamma2460_keV_partial": width_keV(inputs["m_ds1"], inputs["m_ds"], G2460),
                        "Gamma2536_keV_partial": width_keV(2.53511, inputs["m_ds"], G2536),
                        "fB": fB,
                    })

    csv_path = OUT / "stage1_tensor_gb_soft2p_estimate.csv"
    keys = list(rows[0].keys())
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Stage-1 tensor-current soft two-particle estimate",
        "=================================================",
        "Complete result includes tensor-DA family and vector f3gamma psi^v DA.",
        "Vector DA uses the stated physical-residue p^2/(p.q) prescription.",
        "Hard tensor-current contribution is not included.",
    ]
    for scheme in fB_options:
        for include_vec in (False, True):
            subset = [
                r for r in rows
                if r["fB_scheme"] == scheme
                and r["include_vector_da_diagnostic"] == include_vec
            ]
            gb = np.array([r["GB_soft2p_estimate"] for r in subset])
            w2460 = np.array([r["Gamma2460_keV_partial"] for r in subset])
            suffix = "with vector diagnostic" if include_vec else "tensor DA only"
            lines.append(
                f"{scheme}, {suffix}: GB {gb.min():+.4f} to {gb.max():+.4f}; "
                f"Gamma2460(partial) {w2460.min():.2f} to {w2460.max():.2f} keV"
            )
    lines.append(f"Wrote grid to {csv_path}")
    summary_path = OUT / "stage1_tensor_gb_soft2p_estimate_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
