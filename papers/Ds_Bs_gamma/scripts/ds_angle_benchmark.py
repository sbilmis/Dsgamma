"""Central Ds interference benchmarks at the three documented angles."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from final_stage2_uncertainty_scan import precompute_f1_basis
from lattice_photon_normalization_comparison import (
    central_inputs_for_scenario,
    evaluate,
)
from mixing_angle_inputs import THETA_DS_DEG, THETA_HQET_DEG
from stage1_tensor_gb_soft2p_estimate import width_keV
from stage2_axial_g1_three_particle import F1_integral
from twopoint_ds1_matrix_sumrule import (
    Inputs,
    fitted_threshold,
    ope_matrix,
    projected_ope,
    rotation_from_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"

# Midpoint of the selected final Ds transition window used in the publication
# Monte Carlo: M^2=[3.0,4.5] GeV^2 and s0=[7.5,8.5] GeV^2.
TRANSITION_M2 = 3.75
TRANSITION_S0_ROOT = math.sqrt(8.0)
TWOPOINT_M2 = 2.3
TWOPOINT_S0_MIX = 10.0


def projected_residues(theta_deg: float) -> tuple[float, float, float, float]:
    inp = Inputs()
    s01, _ = fitted_threshold(TWOPOINT_M2, theta_deg, 0, inp.mass_low, inp)
    s02, _ = fitted_threshold(TWOPOINT_M2, theta_deg, 1, inp.mass_high, inp)
    pi1, _ = projected_ope(TWOPOINT_M2, s01, theta_deg, 0, inp)
    pi2, _ = projected_ope(TWOPOINT_M2, s02, theta_deg, 1, inp)
    f1 = math.sqrt(pi1 * math.exp(inp.mass_low**2 / TWOPOINT_M2) / inp.mass_low**2)
    f2 = math.sqrt(pi2 * math.exp(inp.mass_high**2 / TWOPOINT_M2) / inp.mass_high**2)
    return f1, f2, s01, s02


def main() -> None:
    inp = Inputs()
    matrix_angle, _ = rotation_from_matrix(
        ope_matrix(TWOPOINT_M2, TWOPOINT_S0_MIX, inp)[0]
    )
    angle_cases = [
        ("previous_study_nominal", THETA_DS_DEG),
        ("truncated_matrix_diagnostic", matrix_angle),
        ("hqet_benchmark", THETA_HQET_DEG),
    ]

    vals = central_inputs_for_scenario("lattice_fperp_s")
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows = []
    for label, theta_deg in angle_cases:
        basis = evaluate(
            vals,
            TRANSITION_M2,
            TRANSITION_S0_ROOT,
            theta_deg,
            f1_axial,
            f1_basis,
        )
        f1, f2, s01, s02 = projected_residues(theta_deg)
        theta = math.radians(theta_deg)
        f_a = vals["f_ds1"]
        f_b = vals["fT"] * vals["m_ds1"] / (vals["mc"] + vals["ms"])
        low_a = math.sin(theta) * f_a * basis["GA"] / f1
        low_b = math.cos(theta) * f_b * basis["GB"] / f1
        high_a = math.cos(theta) * f_a * basis["GA"] / f2
        high_b = -math.sin(theta) * f_b * basis["GB"] / f2
        g_low = low_a + low_b
        g_high = high_a + high_b
        rows.append(
            {
                "case": label,
                "theta_deg": theta_deg,
                "f1_GeV": f1,
                "f2_GeV": f2,
                "s01_GeV2": s01,
                "s02_GeV2": s02,
                "GA_GeV_inv": basis["GA"],
                "GB_GeV_inv": basis["GB"],
                "G_low_A_component_GeV_inv": low_a,
                "G_low_B_component_GeV_inv": low_b,
                "G_low_GeV_inv": g_low,
                "Gamma_low_keV": width_keV(vals["m_ds1"], vals["m_ds"], g_low),
                "G_high_A_component_GeV_inv": high_a,
                "G_high_B_component_GeV_inv": high_b,
                "G_high_GeV_inv": g_high,
                "Gamma_high_keV": width_keV(vals["m_ds1_2536"], vals["m_ds"], g_high),
                "low_interference": "constructive" if low_a * low_b > 0.0 else "destructive",
                "high_interference": "constructive" if high_a * high_b > 0.0 else "destructive",
            }
        )

    cancellation_angle = math.degrees(
        math.atan2(vals["f_ds1"] * rows[0]["GA_GeV_inv"], f_b * rows[0]["GB_GeV_inv"])
    )
    if cancellation_angle < 0.0:
        cancellation_angle += 180.0

    path = OUT / "ds_mixing_angle_benchmark_table.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (OUT / "ds1_interference_diagnostic_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Ds mixing-angle and interference benchmark",
        "=============================================",
        f"Transition point: M^2={TRANSITION_M2} GeV^2, sqrt(s0)={TRANSITION_S0_ROOT} GeV.",
        f"Two-point projection point: M^2={TWOPOINT_M2} GeV^2.",
        f"High-state central cancellation angle: {cancellation_angle:.3f} deg.",
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"{row['case']}: theta={row['theta_deg']:.4f} deg",
                f"  f1={row['f1_GeV']:.6f} GeV; f2={row['f2_GeV']:.6f} GeV",
                "  low:  A={:+.6f}, B={:+.6f}, G={:+.6f} GeV^-1; Gamma={:.6f} keV ({})".format(
                    row["G_low_A_component_GeV_inv"], row["G_low_B_component_GeV_inv"],
                    row["G_low_GeV_inv"], row["Gamma_low_keV"], row["low_interference"]
                ),
                "  high: A={:+.6f}, B={:+.6f}, G={:+.6f} GeV^-1; Gamma={:.6f} keV ({})".format(
                    row["G_high_A_component_GeV_inv"], row["G_high_B_component_GeV_inv"],
                    row["G_high_GeV_inv"], row["Gamma_high_keV"], row["high_interference"]
                ),
                "",
            ]
        )
    text = "\n".join(lines) + "\n"
    (OUT / "ds_mixing_angle_benchmark_summary.txt").write_text(text)
    (OUT / "ds1_interference_diagnostic_summary.txt").write_text(text)
    print(text)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
