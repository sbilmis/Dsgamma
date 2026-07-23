"""Rohrwild-nonlocal Stage-2 baseline for Bs1 -> Bs gamma.

This script reuses the completed Ds1 LCSR machinery with the heavy-quark
replacement c -> b.  The imported functions still use the internal key ``mc``
for the heavy-quark mass; in this file that key should be read as m_b.

The numerical inputs are intentionally labeled as a working baseline.  The
observed narrow state is B_s1(5830).  A lower 1+ partner has not been firmly
established experimentally, so the lower-state row uses a model mass only as a
diagnostic.
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

from final_stage2_uncertainty_scan import matched_f1_integral, precompute_f1_basis
from stage1_tensor_gb_hard_candidate import hard_tensor_gb
from stage1_tensor_gb_soft2p_estimate import gb_soft_two_particle, width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from stage2_tensor_gb_three_particle_estimate import gb_em_da, gb_three_particle_from_integral


FPERP_S_2GEV = -0.0510
TENSOR_RUNNING_2_TO_1 = 1.08
FPERP_S_1GEV = FPERP_S_2GEV * TENSOR_RUNNING_2_TO_1


def hq_scaled_decay_constant(f_charm: float, m_charm_state: float, m_bottom_state: float) -> float:
    """Naive HQ scaling f_H sqrt(M_H) = constant for a working baseline."""
    return f_charm * math.sqrt(m_charm_state / m_bottom_state)


def bs1_inputs(m_initial: float, scenario: str) -> dict[str, float]:
    """Working central inputs for Bs1 -> Bs gamma.

    ``scenario`` chooses the leading-twist photon normalization:
      - legacy_chi_condensate: chi_s(1 GeV) = 3.15 GeV^-2;
      - lattice_fperp_s: replace chi_s <sbar s> by f_gamma,s^perp(1 GeV).
    """
    qq = -(0.240) ** 3
    ss = 0.8 * qq
    chi = 3.15
    if scenario == "lattice_fperp_s":
        chi = FPERP_S_1GEV / ss
    elif scenario != "legacy_chi_condensate":
        raise ValueError(f"unknown scenario: {scenario}")

    # f_Bs1 and fT are working estimates obtained from the Ds1 central values
    # by f sqrt(M) scaling.  They should be replaced by dedicated Bs1 sum-rule
    # or lattice inputs before a final paper table.
    f_bs1_a = hq_scaled_decay_constant(0.225, 2.4595, m_initial)
    fT_bs1 = hq_scaled_decay_constant(0.256, 2.4595, m_initial)

    return {
        "mc": 4.18,  # read as m_b(m_b)
        "ms": 0.093,
        "m_ds1": m_initial,  # read as m_Bs1
        "m_ds": 5.36692,  # read as m_Bs
        "f_ds1": f_bs1_a,  # read as f_Bs1,A
        "f_ds": 0.2303,  # read as f_Bs
        "fT": fT_bs1,
        "ss": ss,
        "chi": chi,
        "f3g": -0.0039,
        "ec": -1.0 / 3.0,  # read as e_b
        "es": -1.0 / 3.0,
        "fperp_s_used": chi * ss,
    }


def evaluate_point(inputs: dict[str, float], M2: float, s0_root: float, theta_deg: float, f1_axial: float, f1_basis: dict[str, float]):
    calc_inputs = dict(inputs)
    fT = calc_inputs.pop("fT")
    fperp = calc_inputs.pop("fperp_s_used")
    s0 = s0_root * s0_root
    theta = math.radians(theta_deg)

    axial = g1_stage2(M2, s0, calc_inputs, f1_axial)
    GA = axial["g1_stage2_GeV_inv"]
    fB = fT * calc_inputs["m_ds1"] / (calc_inputs["mc"] + calc_inputs["ms"])
    GB_hard, hard_qcd = hard_tensor_gb(M2, s0, calc_inputs, fB)
    GB_soft = gb_soft_two_particle(axial, calc_inputs, fB)
    f1_tensor = matched_f1_integral(calc_inputs, f1_basis)
    GB_3p = gb_three_particle_from_integral(axial, calc_inputs, fB, f1_tensor)
    GB_em, _ = gb_em_da(axial, calc_inputs, fB)
    GB = GB_hard + GB_soft + GB_3p + GB_em

    G_low = math.sin(theta) * GA + math.cos(theta) * GB
    G_high = math.cos(theta) * GA - math.sin(theta) * GB

    return {
        "M2": M2,
        "s0_root": s0_root,
        "s0": s0,
        "theta_deg": theta_deg,
        "GA": GA,
        "GB": GB,
        "GB_hard": GB_hard,
        "GB_soft2p": GB_soft,
        "GB_3p": GB_3p,
        "GB_em": GB_em,
        "G_low_combo": G_low,
        "G_high_combo": G_high,
        "Gamma_low_combo_keV": width_keV(calc_inputs["m_ds1"], calc_inputs["m_ds"], G_low),
        "Gamma_high_combo_keV": width_keV(calc_inputs["m_ds1"], calc_inputs["m_ds"], G_high),
        "fB_effective": fB,
        "hard_qcd_side": hard_qcd,
        "f1_tensor_matched": f1_tensor,
        "chi_effective": calc_inputs["chi"],
        "fperp_s_used_GeV": fperp,
    }


def summarize(values):
    arr = np.array(values, dtype=float)
    return {
        "min": float(np.min(arr)),
        "median": float(np.percentile(arr, 50.0)),
        "max": float(np.max(arr)),
    }


def fmt_range(stats, unit=""):
    return f"{stats['median']:.4g}{unit} [{stats['min']:.4g}, {stats['max']:.4g}]"


def main():
    theta_deg = 35.3
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()

    targets = [
        {
            "state": "Bs1_5830_observed_high_combo",
            "m_initial": 5.82870,
            "quoted_combo": "high",
            "mass_comment": "observed narrow B_s1(5830)",
            "M2_min": 8.0,
            "M2_max": 14.0,
            "s0_roots": (6.05, 6.15, 6.25),
        },
        {
            "state": "Bs1_lower_model_5720_low_combo",
            "m_initial": 5.720,
            "quoted_combo": "low",
            "mass_comment": "model lower 1+ partner diagnostic",
            "M2_min": 8.0,
            "M2_max": 14.0,
            "s0_roots": (5.95, 6.05, 6.15),
        },
    ]
    scenarios = ("legacy_chi_condensate", "lattice_fperp_s")

    rows = []
    for target in targets:
        for scenario in scenarios:
            inputs = bs1_inputs(target["m_initial"], scenario)
            for s0_root in target["s0_roots"]:
                for M2 in np.linspace(target["M2_min"], target["M2_max"], 13):
                    row = evaluate_point(inputs, float(M2), float(s0_root), theta_deg, f1_axial, f1_basis)
                    row.update(
                        {
                            "state": target["state"],
                            "m_initial": target["m_initial"],
                            "m_final": inputs["m_ds"],
                            "quoted_combo": target["quoted_combo"],
                            "mass_comment": target["mass_comment"],
                            "scenario": scenario,
                            "mb_MSbar": inputs["mc"],
                            "ms_MSbar": inputs["ms"],
                            "f_Bs": inputs["f_ds"],
                            "f_Bs1_A_working": inputs["f_ds1"],
                            "fT_Bs1_working": inputs["fT"],
                        }
                    )
                    row["G_quoted"] = row[f"G_{target['quoted_combo']}_combo"]
                    row["Gamma_quoted_keV"] = row[f"Gamma_{target['quoted_combo']}_combo_keV"]
                    rows.append(row)

    csv_path = OUT / "bs1_stage2_baseline.csv"
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Bs1 -> Bs gamma Rohrwild-nonlocal Stage-2 baseline",
        "===================================================",
        "Ordinary transition local condensate: excluded.",
        "S_gamma/T4^gamma electromagnetic sector: included.",
        "Internal key mc is used as m_b in this script.",
        "Hard photon and complete two-/three-particle terms are included within",
        "the stated tensor physical-residue prescription.",
        "Borel window: M^2 in [8,14] GeV^2.",
        "The Bs1 decay constants are working HQ-scaled estimates, not final inputs.",
        "",
    ]
    for target in targets:
        lines.append(f"State target: {target['state']} ({target['mass_comment']})")
        for scenario in scenarios:
            subset = [r for r in rows if r["state"] == target["state"] and r["scenario"] == scenario]
            gamma_stats = summarize([r["Gamma_quoted_keV"] for r in subset])
            g_stats = summarize([r["G_quoted"] for r in subset])
            ga_stats = summarize([r["GA"] for r in subset])
            gb_stats = summarize([r["GB"] for r in subset])
            chi_stats = summarize([r["chi_effective"] for r in subset])
            lines.append(
                f"  {scenario}: G_quoted {fmt_range(g_stats, ' GeV^-1')}; "
                f"Gamma {fmt_range(gamma_stats, ' keV')}; "
                f"GA {fmt_range(ga_stats, ' GeV^-1')}; "
                f"GB {fmt_range(gb_stats, ' GeV^-1')}; "
                f"chi_eff {fmt_range(chi_stats, ' GeV^-2')}"
            )
        lines.append("")
    lines.append(f"Wrote grid to {csv_path}")

    summary_path = OUT / "bs1_stage2_baseline_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
