"""One-at-a-time uncertainty budget for the matched Stage-2 result.

The global Monte Carlo scan is the preferred quoted uncertainty.  This script
is a diagnostic companion: it varies one input group at a time around the
central point so that the source of the width spread is visible.
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
sys.path.insert(0, str(ROOT.parents[1] / "shared"))
sys.path.insert(0, str(ROOT / "scripts"))

import photon_da as pda
from final_stage2_uncertainty_scan import (
    clipped_normal,
    matched_f1_integral,
    precompute_f1_basis,
    set_photon_shape_params,
)
from stage1_axial_g1_baseline import central_inputs
from stage1_tensor_gb_hard_candidate import hard_tensor_gb
from stage1_tensor_gb_soft2p_estimate import gb_soft_two_particle, width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from stage2_tensor_gb_three_particle_estimate import gb_em_da, gb_three_particle_from_integral


N_POINTS = 500
SEED = 20260702


def central_sample():
    inputs = central_inputs()
    return {
        **inputs,
        "m_ds1_2536": 2.53511,
        "fT": 0.256,
        "omegaA": -2.1,
        "omegaV": 3.8,
    }


def vary_group(rng, group):
    vals = central_sample()
    qq_scale = clipped_normal(rng, 0.240, 0.010, 0.200, 0.280)
    ss_ratio = clipped_normal(rng, 0.8, 0.1, 0.4, 1.2)
    qq = -(qq_scale**3)

    if group == "borel_threshold":
        pass
    elif group == "theta":
        pass
    elif group == "quark_and_hadron_masses":
        vals.update(
            {
                "mc": clipped_normal(rng, 1.27, 0.02, 1.15, 1.40),
                "ms": clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
                "m_ds1": clipped_normal(rng, 2.4595, 0.0006),
                "m_ds1_2536": clipped_normal(rng, 2.53511, 0.00006),
                "m_ds": clipped_normal(rng, 1.96835, 0.00007),
            }
        )
    elif group == "decay_constants":
        vals.update(
            {
                "f_ds1": clipped_normal(rng, 0.225, 0.025, 0.150, 0.320),
                "f_ds": clipped_normal(rng, 0.2499, 0.0005, 0.245, 0.255),
                "fT": clipped_normal(rng, 0.256, 0.017, 0.190, 0.330),
            }
        )
    elif group == "condensates_and_chi":
        vals.update(
            {
                "ss": ss_ratio * qq,
                "chi": clipped_normal(rng, 3.15, 0.30, 2.0, 4.2),
                "f3g": clipped_normal(rng, -0.0039, 0.0020, -0.010, 0.004),
            }
        )
    elif group == "photon_DA_shape":
        vals.update(
            {
                "omegaA": clipped_normal(rng, -2.1, 1.0, -6.0, 2.0),
                "omegaV": clipped_normal(rng, 3.8, 1.8, -2.0, 10.0),
            }
        )
    elif group == "all_inputs_fixed_theta":
        vals.update(
            {
                "mc": clipped_normal(rng, 1.27, 0.02, 1.15, 1.40),
                "ms": clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
                "m_ds1": clipped_normal(rng, 2.4595, 0.0006),
                "m_ds1_2536": clipped_normal(rng, 2.53511, 0.00006),
                "m_ds": clipped_normal(rng, 1.96835, 0.00007),
                "f_ds1": clipped_normal(rng, 0.225, 0.025, 0.150, 0.320),
                "f_ds": clipped_normal(rng, 0.2499, 0.0005, 0.245, 0.255),
                "fT": clipped_normal(rng, 0.256, 0.017, 0.190, 0.330),
                "ss": ss_ratio * qq,
                "chi": clipped_normal(rng, 3.15, 0.30, 2.0, 4.2),
                "f3g": clipped_normal(rng, -0.0039, 0.0020, -0.010, 0.004),
                "omegaA": clipped_normal(rng, -2.1, 1.0, -6.0, 2.0),
                "omegaV": clipped_normal(rng, 3.8, 1.8, -2.0, 10.0),
            }
        )
    elif group == "all_inputs_and_theta":
        vals = vary_group(rng, "all_inputs_fixed_theta")
    else:
        raise ValueError(f"Unknown group: {group}")
    return vals


def evaluate(vals, M2, s0_root, theta_deg, f1_axial_total, f1_basis):
    set_photon_shape_params(vals["omegaA"], vals["omegaV"])
    inputs = dict(vals)
    m2536 = inputs.pop("m_ds1_2536")
    fT = inputs.pop("fT")
    inputs.pop("omegaA")
    inputs.pop("omegaV")

    s0 = s0_root * s0_root
    theta = math.radians(theta_deg)
    axial = g1_stage2(M2, s0, inputs, f1_axial_total)
    GA = axial["g1_stage2_GeV_inv"]
    fB = fT * inputs["m_ds1"] / (inputs["mc"] + inputs["ms"])
    GB_hard, _ = hard_tensor_gb(M2, s0, inputs, fB)
    GB_soft = gb_soft_two_particle(axial, inputs, fB)
    f1_tensor = matched_f1_integral(inputs, f1_basis)
    GB_3p = gb_three_particle_from_integral(axial, inputs, fB, f1_tensor)
    GB_em, _ = gb_em_da(axial, inputs, fB)
    GB = GB_hard + GB_soft + GB_3p + GB_em
    G2460 = math.sin(theta) * GA + math.cos(theta) * GB
    G2536 = math.cos(theta) * GA - math.sin(theta) * GB
    return {
        "M2": M2,
        "s0_root": s0_root,
        "theta_deg": theta_deg,
        "GA": GA,
        "GB": GB,
        "G2460": G2460,
        "G2536": G2536,
        "Gamma2460_keV": width_keV(inputs["m_ds1"], inputs["m_ds"], G2460),
        "Gamma2536_keV": width_keV(m2536, inputs["m_ds"], G2536),
    }


def summarize(values):
    arr = np.array(values, dtype=float)
    return {
        "median": float(np.percentile(arr, 50.0)),
        "p16": float(np.percentile(arr, 16.0)),
        "p84": float(np.percentile(arr, 84.0)),
        "half_width": float(0.5 * (np.percentile(arr, 84.0) - np.percentile(arr, 16.0))),
    }


def main():
    rng = np.random.default_rng(SEED)
    f1_axial_total, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    groups = [
        "borel_threshold",
        "theta",
        "quark_and_hadron_masses",
        "decay_constants",
        "condensates_and_chi",
        "photon_DA_shape",
        "all_inputs_fixed_theta",
        "all_inputs_and_theta",
    ]

    rows = []
    for group in groups:
        for _ in range(N_POINTS):
            vals = vary_group(rng, group)
            M2 = 4.5
            s0_root = 2.55
            theta_deg = 35.3
            if group in {"borel_threshold", "all_inputs_fixed_theta", "all_inputs_and_theta"}:
                M2 = float(rng.uniform(3.0, 6.0))
                s0_root = float(rng.uniform(2.50, 2.60))
            if group in {"theta", "all_inputs_and_theta"}:
                theta_deg = float(rng.uniform(25.0, 45.0))
            rows.append({"group": group, **evaluate(vals, M2, s0_root, theta_deg, f1_axial_total, f1_basis)})

    csv_path = OUT / "uncertainty_budget_one_at_time.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "One-at-a-time Stage-2 uncertainty budget",
        "========================================",
        f"{N_POINTS} random points per group.",
        "Central point unless the group name is varied.",
        "",
    ]
    summary_rows = []
    for group in groups:
        subset = [r for r in rows if r["group"] == group]
        stats2460 = summarize([r["Gamma2460_keV"] for r in subset])
        stats2536 = summarize([r["Gamma2536_keV"] for r in subset])
        summary = {
            "group": group,
            "Gamma2460_median": stats2460["median"],
            "Gamma2460_p16": stats2460["p16"],
            "Gamma2460_p84": stats2460["p84"],
            "Gamma2460_half_width": stats2460["half_width"],
            "Gamma2536_median": stats2536["median"],
            "Gamma2536_p16": stats2536["p16"],
            "Gamma2536_p84": stats2536["p84"],
            "Gamma2536_half_width": stats2536["half_width"],
        }
        summary_rows.append(summary)
        lines.append(
            f"{group}: Gamma2460 {stats2460['median']:.3g} "
            f"[{stats2460['p16']:.3g}, {stats2460['p84']:.3g}] keV; "
            f"Gamma2536 {stats2536['median']:.3g} "
            f"[{stats2536['p16']:.3g}, {stats2536['p84']:.3g}] keV"
        )
    lines.append("")
    lines.append(f"Wrote scan rows to {csv_path}")

    summary_csv = OUT / "uncertainty_budget_one_at_time_summary.csv"
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    summary_txt = OUT / "uncertainty_budget_one_at_time_summary.txt"
    summary_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
