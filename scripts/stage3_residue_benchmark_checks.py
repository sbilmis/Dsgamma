"""Diagnostic checks for the Ds1 mixed-current residue model.

This script does not define the headline normalization.  It checks how the
Stage-3 two-point residue model behaves in limiting current choices and when it
is compared with an external pure-axial local-QCDSR benchmark.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

THETA_PHYS_DEG = 35.3

# Basis-current values used in the present Stage-3 bookkeeping.
F_A_WORKING = 0.225
F_A_WORKING_ERR = 0.025
F_T_WORKING = 0.256
F_T_WORKING_ERR = 0.017
M_DS1 = 2.4595
M_C = 1.27
M_S = 0.093
F_B_WORKING = F_T_WORKING * M_DS1 / (M_C + M_S)

# External pure-axial benchmark read from Wang's Table 1.  We use it only as a
# comparison target for the pure axial current, not as an authority for the
# physical mixed-current residues.
F_A_WANG_TABLE1 = 0.245
F_A_WANG_TABLE1_ERR = 0.017

# Current headline Stage-3 anchor in the project.  It is kept here as a
# scenario to compare against the pure-axial benchmark.
F1_STAGE3_ANCHOR = 0.345
F1_STAGE3_ANCHOR_ERR = 0.017


def residue_pair(theta_deg: float, f_a: float, f_b: float, rho: float) -> tuple[float, float]:
    theta = math.radians(theta_deg)
    s = math.sin(theta)
    c = math.cos(theta)
    f1_sq = s * s * f_a * f_a + c * c * f_b * f_b + 2.0 * s * c * rho * f_a * f_b
    f2_sq = c * c * f_a * f_a + s * s * f_b * f_b - 2.0 * s * c * rho * f_a * f_b
    return math.sqrt(max(f1_sq, 0.0)), math.sqrt(max(f2_sq, 0.0))


def rho_for_f1(theta_deg: float, f_a: float, f_b: float, f1_target: float) -> tuple[float, bool]:
    theta = math.radians(theta_deg)
    s = math.sin(theta)
    c = math.cos(theta)
    denom = 2.0 * s * c * f_a * f_b
    if abs(denom) < 1.0e-14:
        return float("nan"), False
    rho = (f1_target * f1_target - s * s * f_a * f_a - c * c * f_b * f_b) / denom
    clipped = False
    if rho < -1.0:
        rho = -1.0
        clipped = True
    if rho > 1.0:
        rho = 1.0
        clipped = True
    return rho, clipped


def fmt(x: float) -> str:
    if math.isnan(x):
        return "nan"
    return f"{x:.6g}"


def main() -> None:
    rows: list[dict[str, object]] = []

    for theta_deg, rho, description in [
        (90.0, 0.0, "J1 is pure axial; f1 must equal fA"),
        (0.0, 0.0, "J2 is pure axial; f2 must equal fA"),
        (THETA_PHYS_DEG, -1.0, "maximally anticorrelated AA/BB overlap"),
        (THETA_PHYS_DEG, 0.0, "orthogonal AA/BB basis residues"),
        (THETA_PHYS_DEG, 1.0, "maximally correlated AA/BB overlap"),
    ]:
        f1, f2 = residue_pair(theta_deg, F_A_WORKING, F_B_WORKING, rho)
        rows.append(
            {
                "check": "fixed_rho_working_basis",
                "theta_deg": theta_deg,
                "rho_AB": rho,
                "fA_input_GeV": F_A_WORKING,
                "fB_input_GeV": F_B_WORKING,
                "target_f1_GeV": "",
                "target_description": description,
                "rho_clipped": "",
                "f1_GeV": f1,
                "f2_GeV": f2,
            }
        )

    scenarios = [
        (
            "physical_angle_fit_wang_pure_axial_with_working_fA",
            F_A_WORKING,
            F_B_WORKING,
            F_A_WANG_TABLE1,
            "physical angle; choose rho_AB so f1 matches Wang Table-1 pure axial value",
        ),
        (
            "physical_angle_fit_wang_pure_axial_with_wang_fA",
            F_A_WANG_TABLE1,
            F_B_WORKING,
            F_A_WANG_TABLE1,
            "same check but use Wang value also as fA basis input",
        ),
        (
            "physical_angle_fit_stage3_anchor_with_working_fA",
            F_A_WORKING,
            F_B_WORKING,
            F1_STAGE3_ANCHOR,
            "current headline Stage-3 f1-anchor scenario",
        ),
    ]

    for check, f_a, f_b, f1_target, description in scenarios:
        rho, clipped = rho_for_f1(THETA_PHYS_DEG, f_a, f_b, f1_target)
        f1, f2 = residue_pair(THETA_PHYS_DEG, f_a, f_b, rho)
        rows.append(
            {
                "check": check,
                "theta_deg": THETA_PHYS_DEG,
                "rho_AB": rho,
                "fA_input_GeV": f_a,
                "fB_input_GeV": f_b,
                "target_f1_GeV": f1_target,
                "target_description": description,
                "rho_clipped": int(clipped),
                "f1_GeV": f1,
                "f2_GeV": f2,
            }
        )

    out_csv = OUT / "stage3_residue_benchmark_checks.csv"
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Stage-3 mixed-current residue benchmark checks",
        "================================================",
        f"theta_phys = {THETA_PHYS_DEG:.1f} deg",
        f"working fA = {F_A_WORKING:.3f} +- {F_A_WORKING_ERR:.3f} GeV",
        f"working fB = fT*M/(mc+ms) = {F_B_WORKING:.3f} GeV",
        f"Wang Table-1 pure axial benchmark = {F_A_WANG_TABLE1:.3f} +- {F_A_WANG_TABLE1_ERR:.3f} GeV",
        f"current Stage-3 f1 anchor = {F1_STAGE3_ANCHOR:.3f} +- {F1_STAGE3_ANCHOR_ERR:.3f} GeV",
        "",
    ]
    for row in rows:
        formatted = {k: fmt(v) if isinstance(v, float) else v for k, v in row.items()}
        lines.append(
            "{check}: theta={theta_deg} deg, rho={rho_AB}, fA={fA_input_GeV}, "
            "fB={fB_input_GeV}, target={target_f1_GeV}, f1={f1_GeV}, f2={f2_GeV}, "
            "clipped={rho_clipped}".format(
                **formatted
            )
        )
        lines.append(f"  {row['target_description']}")
    lines.append("")
    lines.append(
        "Interpretation: the pure-current limits reproduce fA by construction. "
        "At the physical angle, a Wang-like pure-axial target lies near the lower "
        "edge of the mixed-current residue range and requires a strongly negative "
        "AA/BB overlap.  The current Stage-3 anchor is therefore a separate "
        "normalization scenario, not a direct consequence of Wang's pure axial "
        "current."
    )

    out_txt = OUT / "stage3_residue_benchmark_checks.txt"
    out_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
