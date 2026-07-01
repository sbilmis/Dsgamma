"""Stage-2 tensor-current three-particle DA estimate.

The first FeynCalc trace pass for the sigma_{ab} part of the background-gluon
heavy-quark propagator finds the same trace ratio as the two-particle tensor
DA sector:

    kernel_B / kernel_A = m_c/(m_c+m_s)

for the local S, T1 and T2 structures after the Borel saddle k.q = p.q.  This
script uses that relation to estimate the tensor-current three-particle
correction from the calibrated axial F1 contribution.

Important limitation:
  The sigma_{ab} part of S_Q^G has now been checked also for the explicit
  x-dependent T3/T4 structures.  The separate x_alpha G^{alpha beta} gamma_beta
  term has been derivative-reduced and contains double-pole pieces, but it has
  not yet been matched to the axial F1 normalization and integrated.  Therefore
  this is still an estimate, not the final tensor-current Stage-2 result.
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

from stage1_axial_g1_baseline import central_inputs
from stage1_tensor_gb_hard_candidate import hard_tensor_gb
from stage1_tensor_gb_soft2p_estimate import gb_soft_tensor_da, width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2


def prefactor_b(M2, inputs, fB):
    mc = inputs["mc"]
    ms = inputs["ms"]
    return (
        math.exp((inputs["m_ds1"]**2 + inputs["m_ds"]**2) / (2.0 * M2))
        * (mc + ms)
        / (inputs["m_ds1"] * fB * inputs["m_ds"]**2 * inputs["f_ds"])
    )


def gb_three_particle_estimate(axial_stage2_row, inputs, fB):
    """Estimate J_B three-particle term from axial F1 QCD-side correction."""
    mc = inputs["mc"]
    ms = inputs["ms"]
    qcd_b = (mc / (mc + ms)) * axial_stage2_row["F1_delta_qcd"]
    return prefactor_b(axial_stage2_row["M2"], inputs, fB) * qcd_b


def main():
    inputs = central_inputs()
    f1_total, _, _ = F1_integral(u0=0.5)
    fT = 0.256
    fB_options = {
        "fB_equal_fT_diagnostic": fT,
        "fB_converted_fT_M_over_msum_central": fT
        * inputs["m_ds1"]
        / (inputs["mc"] + inputs["ms"]),
    }
    theta = math.radians(35.3)

    rows = []
    for s0_root in (2.50, 2.55, 2.60):
        s0 = s0_root**2
        for M2 in np.linspace(3.0, 6.0, 13):
            axial2 = g1_stage2(float(M2), s0, inputs, f1_total)
            GA2 = axial2["g1_stage2_GeV_inv"]
            for label, fB in fB_options.items():
                GB_soft2p = gb_soft_tensor_da(axial2, inputs, fB)
                GB_hard, hard_qcd = hard_tensor_gb(float(M2), s0, inputs, fB)
                GB_3p = gb_three_particle_estimate(axial2, inputs, fB)
                GB_stage2 = GB_hard + GB_soft2p + GB_3p
                G2460 = math.sin(theta) * GA2 + math.cos(theta) * GB_stage2
                G2536 = math.cos(theta) * GA2 - math.sin(theta) * GB_stage2
                rows.append(
                    {
                        "fB_scheme": label,
                        "M2": float(M2),
                        "s0_root": s0_root,
                        "s0": s0,
                        "GA_stage2": GA2,
                        "GB_hard": GB_hard,
                        "GB_soft2p": GB_soft2p,
                        "GB_3p_estimate": GB_3p,
                        "GB_stage2_estimate": GB_stage2,
                        "hard_qcd_side": hard_qcd,
                        "G2460_stage2_estimate": G2460,
                        "G2536_stage2_estimate": G2536,
                        "Gamma2460_keV_stage2_estimate": width_keV(
                            inputs["m_ds1"], inputs["m_ds"], G2460
                        ),
                        "Gamma2536_keV_stage2_estimate": width_keV(
                            2.53511, inputs["m_ds"], G2536
                        ),
                        "fB": fB,
                    }
                )

    csv_path = OUT / "stage2_tensor_gb_three_particle_estimate.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Stage-2 tensor-current three-particle DA estimate",
        "=================================================",
        "Uses B/A trace ratio m_c/(m_c+m_s) for the sigma_ab part of S_Q^G.",
        "This includes local S,T1,T2 and explicit-x T3,T4 kernels.",
        "The separate x_alpha G^{alpha beta} gamma_beta term has double-pole derivative kernels.",
        "It remains to match those kernels to the axial F1 normalization and integrate them.",
    ]
    for scheme in fB_options:
        subset = [r for r in rows if r["fB_scheme"] == scheme]
        gb3 = np.array([r["GB_3p_estimate"] for r in subset])
        gb = np.array([r["GB_stage2_estimate"] for r in subset])
        w2460 = np.array([r["Gamma2460_keV_stage2_estimate"] for r in subset])
        w2536 = np.array([r["Gamma2536_keV_stage2_estimate"] for r in subset])
        lines.append(
            f"{scheme}: GB_3p {gb3.min():+.5f} to {gb3.max():+.5f}; "
            f"GB_stage2 {gb.min():+.4f} to {gb.max():+.4f}; "
            f"Gamma2460 {w2460.min():.2f} to {w2460.max():.2f} keV; "
            f"Gamma2536 {w2536.min():.2f} to {w2536.max():.2f} keV"
        )
    lines.append(f"Wrote grid to {csv_path}")

    summary_path = OUT / "stage2_tensor_gb_three_particle_estimate_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
