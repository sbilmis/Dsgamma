"""Central Ds interference benchmarks at the three documented angles."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from lattice_photon_normalization_comparison import (
    central_inputs_for_scenario,
)
from mixing_angle_inputs import THETA_DS_DEG, THETA_HQET_DEG
from rohrwild_transition_exact import (
    physical_couplings,
    precompute_convolutions,
    transition_invariants,
)
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
    invariants = transition_invariants(
        TRANSITION_M2,
        TRANSITION_S0_ROOT**2,
        vals,
        convolutions=precompute_convolutions(),
    )
    rows = []
    for label, theta_deg in angle_cases:
        f1, f2, s01, s02 = projected_residues(theta_deg)
        physical = physical_couplings(
            invariants,
            theta_deg=theta_deg,
            m_state_1=vals["m_ds1"],
            m_state_2=vals["m_ds1_2536"],
            f_1=f1,
            f_2=f2,
            m_p=vals["m_ds"],
            f_p=vals["f_ds"],
            m_q=vals["mc"],
            m_s=vals["ms"],
        )
        low_a = physical["g_1_A"]
        low_b = physical["g_1_B"]
        high_a = physical["g_2_A"]
        high_b = physical["g_2_B"]
        g_low = physical["g_1"]
        g_high = physical["g_2"]
        rows.append(
            {
                "case": label,
                "theta_deg": theta_deg,
                "f1_GeV": f1,
                "f2_GeV": f2,
                "s01_GeV2": s01,
                "s02_GeV2": s02,
                "T_A_GeV3": invariants["T_A"],
                "T_B_GeV3": invariants["T_B"],
                "G_low_A_component_GeV_inv": low_a,
                "G_low_B_component_GeV_inv": low_b,
                "G_low_GeV_inv": g_low,
                "Gamma_low_keV": physical["Gamma_1_keV"],
                "G_high_A_component_GeV_inv": high_a,
                "G_high_B_component_GeV_inv": high_b,
                "G_high_GeV_inv": g_high,
                "Gamma_high_keV": physical["Gamma_2_keV"],
                "low_interference": "constructive" if low_a * low_b > 0.0 else "destructive",
                "high_interference": "constructive" if high_a * high_b > 0.0 else "destructive",
            }
        )

    cancellation_angle = math.degrees(
        math.atan2(float(invariants["T_A"]), float(invariants["T_B"]))
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
