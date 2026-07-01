"""Stage-2 tensor-current three-particle DA estimate.

The first FeynCalc trace pass for the sigma_{ab} part of the background-gluon
heavy-quark propagator finds the same trace ratio as the two-particle tensor
DA sector:

    kernel_B / kernel_A = m_c/(m_c+m_s)

for the local S, T1 and T2 structures after the Borel saddle k.q = p.q.

This script reports two versions:
  - sigma-part matched: applies the trace ratio only to the sigma_ab part of
    S_Q^G, whose local and explicit-x kernels have been checked;
  - full-ratio diagnostic: the older estimate that multiplies the full axial F1
    correction by m_c/(m_c+m_s).

Important limitation:
  The sigma_{ab} part of S_Q^G has now been checked also for the explicit
  x-dependent T3/T4 structures.  The separate x_alpha G^{alpha beta} gamma_beta
  term has been derivative-reduced and contains double-pole pieces, but it has
  not yet been matched to the axial F1 normalization and integrated.  Therefore
  the full Stage-2 tensor result is not final.
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
from step10_stage2_f1_split import f1_sigma_like, f1_xgamma_like, integrate_piece


def prefactor_b(M2, inputs, fB):
    mc = inputs["mc"]
    ms = inputs["ms"]
    return (
        math.exp((inputs["m_ds1"]**2 + inputs["m_ds"]**2) / (2.0 * M2))
        * (mc + ms)
        / (inputs["m_ds1"] * fB * inputs["m_ds"]**2 * inputs["f_ds"])
    )


def gb_three_particle_from_integral(axial_stage2_row, inputs, fB, f1_integral):
    """Convert one axial-normalized F1 integral into the J_B trace-ratio term."""
    mc = inputs["mc"]
    ms = inputs["ms"]
    delta_qcd = (
        inputs["es"]
        * inputs["ss"]
        * math.exp(-(inputs["mc"] ** 2) / axial_stage2_row["M2"])
        * f1_integral
    )
    qcd_b = (mc / (mc + ms)) * delta_qcd
    return prefactor_b(axial_stage2_row["M2"], inputs, fB) * qcd_b


def main():
    inputs = central_inputs()
    f1_total, _, _ = F1_integral(u0=0.5)
    f1_sigma, _, _ = integrate_piece(lambda aq, aqb, v: f1_sigma_like(aq, aqb), u0=0.5)
    f1_xgamma, _, _ = integrate_piece(f1_xgamma_like, u0=0.5)
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
                GB_3p_sigma = gb_three_particle_from_integral(axial2, inputs, fB, f1_sigma)
                GB_3p_xgamma_unmatched = gb_three_particle_from_integral(
                    axial2, inputs, fB, f1_xgamma
                )
                GB_3p_full_ratio = gb_three_particle_from_integral(axial2, inputs, fB, f1_total)

                GB_sigma_matched = GB_hard + GB_soft2p + GB_3p_sigma
                GB_full_ratio = GB_hard + GB_soft2p + GB_3p_full_ratio

                G2460_sigma = math.sin(theta) * GA2 + math.cos(theta) * GB_sigma_matched
                G2536_sigma = math.cos(theta) * GA2 - math.sin(theta) * GB_sigma_matched
                G2460_full = math.sin(theta) * GA2 + math.cos(theta) * GB_full_ratio
                G2536_full = math.cos(theta) * GA2 - math.sin(theta) * GB_full_ratio
                rows.append(
                    {
                        "fB_scheme": label,
                        "M2": float(M2),
                        "s0_root": s0_root,
                        "s0": s0,
                        "GA_stage2": GA2,
                        "GB_hard": GB_hard,
                        "GB_soft2p": GB_soft2p,
                        "GB_3p_sigma_matched": GB_3p_sigma,
                        "GB_3p_xgamma_unmatched_axial_ratio": GB_3p_xgamma_unmatched,
                        "GB_3p_full_ratio_diagnostic": GB_3p_full_ratio,
                        "GB_stage2_sigma_matched": GB_sigma_matched,
                        "GB_stage2_full_ratio_diagnostic": GB_full_ratio,
                        "hard_qcd_side": hard_qcd,
                        "G2460_sigma_matched": G2460_sigma,
                        "G2536_sigma_matched": G2536_sigma,
                        "G2460_full_ratio_diagnostic": G2460_full,
                        "G2536_full_ratio_diagnostic": G2536_full,
                        "Gamma2460_keV_sigma_matched": width_keV(
                            inputs["m_ds1"], inputs["m_ds"], G2460_sigma
                        ),
                        "Gamma2536_keV_sigma_matched": width_keV(
                            2.53511, inputs["m_ds"], G2536_sigma
                        ),
                        "Gamma2460_keV_full_ratio_diagnostic": width_keV(
                            inputs["m_ds1"], inputs["m_ds"], G2460_full
                        ),
                        "Gamma2536_keV_full_ratio_diagnostic": width_keV(
                            2.53511, inputs["m_ds"], G2536_full
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
        f"F1 sigma-like integral = {f1_sigma:+.6e}",
        f"F1 xgamma-like axial calibration integral = {f1_xgamma:+.6e}",
        f"F1 total integral = {f1_total:+.6e}",
        "The separate x_alpha G^{alpha beta} gamma_beta term has double-pole derivative kernels.",
        "It remains to match those kernels to the axial F1 normalization and integrate them.",
    ]
    for scheme in fB_options:
        subset = [r for r in rows if r["fB_scheme"] == scheme]
        gb3_sigma = np.array([r["GB_3p_sigma_matched"] for r in subset])
        gb3_full = np.array([r["GB_3p_full_ratio_diagnostic"] for r in subset])
        gb_sigma = np.array([r["GB_stage2_sigma_matched"] for r in subset])
        gb_full = np.array([r["GB_stage2_full_ratio_diagnostic"] for r in subset])
        w2460_sigma = np.array([r["Gamma2460_keV_sigma_matched"] for r in subset])
        w2536_sigma = np.array([r["Gamma2536_keV_sigma_matched"] for r in subset])
        w2460_full = np.array([r["Gamma2460_keV_full_ratio_diagnostic"] for r in subset])
        w2536_full = np.array([r["Gamma2536_keV_full_ratio_diagnostic"] for r in subset])
        lines.append(
            f"{scheme}, sigma-matched only: GB_3p {gb3_sigma.min():+.5f} to {gb3_sigma.max():+.5f}; "
            f"GB {gb_sigma.min():+.4f} to {gb_sigma.max():+.4f}; "
            f"Gamma2460 {w2460_sigma.min():.2f} to {w2460_sigma.max():.2f} keV; "
            f"Gamma2536 {w2536_sigma.min():.2f} to {w2536_sigma.max():.2f} keV"
        )
        lines.append(
            f"{scheme}, full-ratio diagnostic: GB_3p {gb3_full.min():+.5f} to {gb3_full.max():+.5f}; "
            f"GB {gb_full.min():+.4f} to {gb_full.max():+.4f}; "
            f"Gamma2460 {w2460_full.min():.2f} to {w2460_full.max():.2f} keV; "
            f"Gamma2536 {w2536_full.min():.2f} to {w2536_full.max():.2f} keV"
        )
    lines.append(f"Wrote grid to {csv_path}")

    summary_path = OUT / "stage2_tensor_gb_three_particle_estimate_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
