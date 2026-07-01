"""Final matched Stage-2 uncertainty scan for Ds1 -> Ds gamma.

The scan uses the completed tensor-current Stage-2 baseline:

  G_B = hard photon + two-particle tensor DA
        + matched three-particle sigma and xgamma terms.

It reports two ensembles:
  - theta_fixed: theta = 35.3 deg;
  - theta_scan_25_45: theta sampled uniformly in [25,45] deg.

The purpose is not to replace a full correlated error analysis.  It is a
transparent first uncertainty budget over the input table, Borel window, and
continuum threshold.
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
from stage1_axial_g1_baseline import central_inputs
from stage1_tensor_gb_hard_candidate import hard_tensor_gb
from stage1_tensor_gb_soft2p_estimate import gb_soft_tensor_da, width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from stage2_tensor_gb_three_particle_estimate import gb_three_particle_from_integral
from step10_stage2_f1_split import f1_sigma_like, integrate_piece


def clipped_normal(rng, mean, sigma, lo=None, hi=None):
    value = rng.normal(mean, sigma)
    if lo is not None:
        value = max(lo, value)
    if hi is not None:
        value = min(hi, value)
    return float(value)


def set_photon_shape_params(omega_a, omega_v):
    pda.PARAMS["omegaA"] = (omega_a, pda.PARAMS["omegaA"][1])
    pda.PARAMS["omegaV"] = (omega_v, pda.PARAMS["omegaV"][1])


def sample_inputs(rng):
    qq_scale = clipped_normal(rng, 0.240, 0.010, 0.200, 0.280)
    ss_ratio = clipped_normal(rng, 0.8, 0.1, 0.4, 1.2)
    qq = -(qq_scale**3)
    return {
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
        "ec": 2.0 / 3.0,
        "es": -1.0 / 3.0,
        "omegaA": clipped_normal(rng, -2.1, 1.0, -6.0, 2.0),
        "omegaV": clipped_normal(rng, 3.8, 1.8, -2.0, 10.0),
    }


def precompute_f1_basis():
    """Precompute DA integrals entering the matched tensor three-particle term."""
    f1_sigma, _, _ = integrate_piece(lambda aq, aqb, v: f1_sigma_like(aq, aqb), u0=0.5)
    i_s, _, _ = integrate_piece(lambda aq, aqb, v: 2.0 * v * (-pda.S_3p(aq, aqb)), u0=0.5)
    i_t2, _, _ = integrate_piece(lambda aq, aqb, v: 2.0 * v * pda.T2_3p(aq, aqb), u0=0.5)
    i_t3, _, _ = integrate_piece(lambda aq, aqb, v: 2.0 * v * (-pda.T3_3p(aq, aqb)), u0=0.5)
    i_t4, _, _ = integrate_piece(lambda aq, aqb, v: 2.0 * v * pda.T4_3p(aq, aqb), u0=0.5)
    return {
        "sigma": f1_sigma,
        "S": i_s,
        "T2": i_t2,
        "T3": i_t3,
        "T4": i_t4,
    }


def matched_f1_integral(inputs, basis):
    """Matched tensor integral at u0=1/2 using residue ratios.

    At the double-pole residue with a=1/2,
      R_S = R_T2 = 1 + (p.q)/(2 m_c^2),
      R_T3 = 1,
      R_T4 = -(p.q)/(2 m_c^2),
    where T4 is normalized to the axial T2 unit.
    """
    pq = (inputs["m_ds1"] ** 2 - inputs["m_ds"] ** 2) / 2.0
    half_pq_over_mc2 = 0.5 * pq / (inputs["mc"] ** 2)
    r_s_t2 = 1.0 + half_pq_over_mc2
    r_t4 = -half_pq_over_mc2
    return (
        basis["sigma"]
        + r_s_t2 * basis["S"]
        + r_s_t2 * basis["T2"]
        + basis["T3"]
        + r_t4 * basis["T4"]
    )


def evaluate_point(rng, ensemble, f1_axial_total, f1_basis):
    sampled = sample_inputs(rng)
    set_photon_shape_params(sampled["omegaA"], sampled["omegaV"])

    inputs = dict(sampled)
    m2536 = inputs.pop("m_ds1_2536")
    fT = inputs.pop("fT")
    inputs.pop("omegaA")
    inputs.pop("omegaV")

    M2 = float(rng.uniform(3.0, 6.0))
    s0_root = float(rng.uniform(2.50, 2.60))
    s0 = s0_root * s0_root
    theta_deg = 35.3 if ensemble == "theta_fixed" else float(rng.uniform(25.0, 45.0))
    theta = math.radians(theta_deg)

    axial = g1_stage2(M2, s0, inputs, f1_axial_total)
    GA = axial["g1_stage2_GeV_inv"]
    fB = fT * inputs["m_ds1"] / (inputs["mc"] + inputs["ms"])
    GB_hard, hard_qcd = hard_tensor_gb(M2, s0, inputs, fB)
    GB_soft = gb_soft_tensor_da(axial, inputs, fB)
    f1_tensor = matched_f1_integral(inputs, f1_basis)
    GB_3p = gb_three_particle_from_integral(axial, inputs, fB, f1_tensor)
    GB = GB_hard + GB_soft + GB_3p

    G2460 = math.sin(theta) * GA + math.cos(theta) * GB
    G2536 = math.cos(theta) * GA - math.sin(theta) * GB
    return {
        "ensemble": ensemble,
        "M2": M2,
        "s0_root": s0_root,
        "theta_deg": theta_deg,
        "GA": GA,
        "GB": GB,
        "GB_hard": GB_hard,
        "GB_soft2p": GB_soft,
        "GB_3p": GB_3p,
        "G2460": G2460,
        "G2536": G2536,
        "Gamma2460_keV": width_keV(inputs["m_ds1"], inputs["m_ds"], G2460),
        "Gamma2536_keV": width_keV(m2536, inputs["m_ds"], G2536),
        "fB": fB,
        "hard_qcd_side": hard_qcd,
        "f1_tensor_matched": f1_tensor,
        **{f"input_{k}": v for k, v in sampled.items()},
    }


def summarize(values):
    arr = np.array(values, dtype=float)
    return {
        "min": float(np.min(arr)),
        "p16": float(np.percentile(arr, 16.0)),
        "median": float(np.percentile(arr, 50.0)),
        "p84": float(np.percentile(arr, 84.0)),
        "max": float(np.max(arr)),
    }


def main():
    rng = np.random.default_rng(20260701)
    set_photon_shape_params(-2.1, 3.8)
    f1_axial_total, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()

    rows = []
    for ensemble in ("theta_fixed", "theta_scan_25_45"):
        for _ in range(400):
            rows.append(evaluate_point(rng, ensemble, f1_axial_total, f1_basis))

    csv_path = OUT / "final_stage2_uncertainty_scan.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Final matched Stage-2 uncertainty scan",
        "======================================",
        "400 random points per ensemble.",
        "Borel window: M^2 in [3,6] GeV^2.",
        "Continuum window: sqrt(s0) in [2.50,2.60] GeV.",
        "Input errors sampled independently from inputs_table.csv.",
        "",
    ]
    for ensemble in ("theta_fixed", "theta_scan_25_45"):
        subset = [r for r in rows if r["ensemble"] == ensemble]
        lines.append(f"Ensemble: {ensemble}")
        for key in ("G2460", "G2536", "Gamma2460_keV", "Gamma2536_keV", "GA", "GB"):
            stats = summarize([r[key] for r in subset])
            lines.append(
                f"  {key}: median {stats['median']:+.4g}; "
                f"68% [{stats['p16']:+.4g}, {stats['p84']:+.4g}]; "
                f"min/max [{stats['min']:+.4g}, {stats['max']:+.4g}]"
            )
        lines.append("")
    lines.append(f"Wrote scan grid to {csv_path}")

    summary_path = OUT / "final_stage2_uncertainty_scan_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
