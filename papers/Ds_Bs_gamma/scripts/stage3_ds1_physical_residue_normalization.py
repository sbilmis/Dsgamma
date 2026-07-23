"""Legacy Stage-3 overlap normalization for Ds1 -> Ds gamma.

The Stage-1/2 scripts compute basis-current form factors GA and GB.  This
script converts those results to a physical-current normalization by inserting
the residues f1 and f2 of the mixed currents

    J1 = sin(theta) JA + cos(theta) JB,
    J2 = cos(theta) JA - sin(theta) JB.

This file preserves the former comparison scheme.  The preferred calculation
is now ``twopoint_ds1_matrix_sumrule.py``, which evaluates the full normalized
AA/AB/BB matrix.  Here, for legacy reproducibility only, we calibrate the
off-diagonal two-point overlap rho_AB from a working Stage-3 anchor for the low
state, f1 ~= f_Ds1(2460) = 0.345(17) GeV.  The same rho_AB then predicts f2 and
the physical-current G2536.  The separate script
stage3_residue_benchmark_checks.py compares this anchor with pure-axial
benchmarks such as Wang's local-QCDSR f_A value.

Do not use the outputs of this script as the current paper normalization.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from stage1_tensor_gb_soft2p_estimate import width_keV


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

F1_STAGE3_ANCHOR = 0.345
F1_STAGE3_ANCHOR_ERR = 0.017


def clipped_normal(rng, mean, sigma, lo=None, hi=None):
    x = float(rng.normal(mean, sigma))
    if lo is not None:
        x = max(lo, x)
    if hi is not None:
        x = min(hi, x)
    return x


def summarize(values):
    arr = np.asarray(values, dtype=float)
    return {
        "p16": float(np.percentile(arr, 16.0)),
        "median": float(np.percentile(arr, 50.0)),
        "p84": float(np.percentile(arr, 84.0)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
    }


def calibrate_from_f1(theta, fA, fB, f1_target):
    """Return rho_AB and f2 implied by the target f1.

    We model the two-point residue matrix as
      f1^2 = s^2 fA^2 + c^2 fB^2 + 2 s c rho fA fB,
      f2^2 = c^2 fA^2 + s^2 fB^2 - 2 s c rho fA fB.
    """
    s = math.sin(theta)
    c = math.cos(theta)
    denom = 2.0 * s * c * fA * fB
    if abs(denom) < 1.0e-12:
        return None

    lower = abs(s * fA - c * fB)
    upper = s * fA + c * fB
    clipped = False
    f1_used = f1_target
    if f1_used < lower:
        f1_used = lower
        clipped = True
    if f1_used > upper:
        f1_used = upper
        clipped = True

    rho = (f1_used * f1_used - s * s * fA * fA - c * c * fB * fB) / denom
    rho = max(-1.0, min(1.0, rho))
    f2_sq = c * c * fA * fA + s * s * fB * fB - 2.0 * s * c * rho * fA * fB
    if f2_sq <= 0.0:
        return None
    return rho, f1_used, math.sqrt(f2_sq), clipped


def physical_couplings(row, rng, scenario):
    theta = math.radians(float(row["theta_deg"]))
    s = math.sin(theta)
    c = math.cos(theta)
    GA = float(row["GA"])
    GB = float(row["GB"])

    if "input_f_ds1" in row and row.get("input_f_ds1"):
        fA = float(row["input_f_ds1"])
    else:
        fA = clipped_normal(rng, 0.225, 0.025, 0.150, 0.320)

    if "fB" in row and row.get("fB"):
        fB = float(row["fB"])
    else:
        fT = clipped_normal(rng, 0.256, 0.017, 0.190, 0.330)
        mc = clipped_normal(rng, 1.27, 0.02, 1.15, 1.40)
        ms = clipped_normal(rng, 0.093, 0.011, 0.050, 0.140)
        m_initial = float(row.get("input_m_ds1", 2.4595))
        fB = fT * m_initial / (mc + ms)

    f1_target = clipped_normal(rng, F1_STAGE3_ANCHOR, F1_STAGE3_ANCHOR_ERR, 0.250, 0.450)
    calibrated = calibrate_from_f1(theta, fA, fB, f1_target)
    if calibrated is None:
        return None
    rho, f1, f2, clipped = calibrated

    g2460 = (s * fA * GA + c * fB * GB) / f1
    g2536 = (c * fA * GA - s * fB * GB) / f2

    m2460 = float(row.get("input_m_ds1", 2.4595))
    m2536 = float(row.get("input_m_ds1_2536", 2.53511))
    mds = float(row.get("input_m_ds", 1.96835))

    return {
        "scenario": scenario,
        "ensemble": row["ensemble"],
        "theta_deg": float(row["theta_deg"]),
        "GA_basis": GA,
        "GB_basis": GB,
        "G2460_basis": float(row["G2460"]),
        "G2536_basis": float(row["G2536"]),
        "Gamma2460_basis_keV": float(row["Gamma2460_keV"]),
        "Gamma2536_basis_keV": float(row["Gamma2536_keV"]),
        "fA_GeV": fA,
        "fB_eff_GeV": fB,
        "f1_GeV": f1,
        "f2_GeV": f2,
        "rho_AB": rho,
        "f1_target_GeV": f1_target,
        "f1_target_clipped": int(clipped),
        "G2460_phys": g2460,
        "G2536_phys": g2536,
        "Gamma2460_phys_keV": width_keV(m2460, mds, g2460),
        "Gamma2536_phys_keV": width_keV(m2536, mds, g2536),
    }


def read_rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    rng = np.random.default_rng(20260705)
    rows = []

    # Legacy chi scan with full input columns.
    for row in read_rows(OUT / "final_stage2_uncertainty_scan.csv"):
        out = physical_couplings(row, rng, "legacy_chi_condensate")
        if out is not None:
            rows.append(out)

    # Lattice-normalized scan; rows do not store fA/fB, so we sample them from
    # the same input distributions used in the original scan.
    for row in read_rows(OUT / "lattice_fperp_comparison_scan.csv"):
        if row["scenario"] != "lattice_fperp_s" or row["ensemble"] == "central":
            continue
        out = physical_couplings(row, rng, "lattice_fperp_s")
        if out is not None:
            rows.append(out)

    write_csv(OUT / "stage3_ds1_physical_residue_scan.csv", rows)

    summary_rows = []
    lines = [
        "Stage-3 Ds1 physical-current residue normalization",
        "===================================================",
        "Transition amplitudes use the completed Rohrwild-nonlocal Stage-2 scan;",
        "ordinary local condensates are excluded and S_gamma/T4^gamma is included.",
        "The off-diagonal two-point overlap rho_AB is calibrated from",
        f"the working Stage-3 anchor f1 = {F1_STAGE3_ANCHOR:.3f} +- {F1_STAGE3_ANCHOR_ERR:.3f} GeV.",
        "The calibrated rho_AB then predicts f2 and the physical G2536.",
        "",
    ]
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for ensemble in ("theta_fixed", "theta_scan_25_45"):
            subset = [r for r in rows if r["scenario"] == scenario and r["ensemble"] == ensemble]
            if not subset:
                continue
            stats = {key: summarize([r[key] for r in subset]) for key in (
                "f1_GeV",
                "f2_GeV",
                "rho_AB",
                "G2460_phys",
                "G2536_phys",
                "Gamma2460_phys_keV",
                "Gamma2536_phys_keV",
            )}
            clipped_fraction = sum(r["f1_target_clipped"] for r in subset) / len(subset)
            summary = {
                "scenario": scenario,
                "ensemble": ensemble,
                "n": len(subset),
                "f1_median": stats["f1_GeV"]["median"],
                "f1_p16": stats["f1_GeV"]["p16"],
                "f1_p84": stats["f1_GeV"]["p84"],
                "f2_median": stats["f2_GeV"]["median"],
                "f2_p16": stats["f2_GeV"]["p16"],
                "f2_p84": stats["f2_GeV"]["p84"],
                "rho_median": stats["rho_AB"]["median"],
                "rho_p16": stats["rho_AB"]["p16"],
                "rho_p84": stats["rho_AB"]["p84"],
                "G2460_median": stats["G2460_phys"]["median"],
                "G2460_p16": stats["G2460_phys"]["p16"],
                "G2460_p84": stats["G2460_phys"]["p84"],
                "G2536_median": stats["G2536_phys"]["median"],
                "G2536_p16": stats["G2536_phys"]["p16"],
                "G2536_p84": stats["G2536_phys"]["p84"],
                "Gamma2460_median": stats["Gamma2460_phys_keV"]["median"],
                "Gamma2460_p16": stats["Gamma2460_phys_keV"]["p16"],
                "Gamma2460_p84": stats["Gamma2460_phys_keV"]["p84"],
                "Gamma2536_median": stats["Gamma2536_phys_keV"]["median"],
                "Gamma2536_p16": stats["Gamma2536_phys_keV"]["p16"],
                "Gamma2536_p84": stats["Gamma2536_phys_keV"]["p84"],
                "f1_target_clipped_fraction": clipped_fraction,
            }
            summary_rows.append(summary)
            lines.append(f"{scenario}, {ensemble}, n={len(subset)}")
            lines.append(
                "  f1={:.3f} [{:.3f},{:.3f}] GeV; f2={:.3f} [{:.3f},{:.3f}] GeV; "
                "rho={:+.3f} [{:+.3f},{:+.3f}]".format(
                    summary["f1_median"], summary["f1_p16"], summary["f1_p84"],
                    summary["f2_median"], summary["f2_p16"], summary["f2_p84"],
                    summary["rho_median"], summary["rho_p16"], summary["rho_p84"],
                )
            )
            lines.append(
                "  G2460={:+.4g} [{:+.4g},{:+.4g}] GeV^-1; "
                "Gamma2460={:.4g} [{:.4g},{:.4g}] keV".format(
                    summary["G2460_median"], summary["G2460_p16"], summary["G2460_p84"],
                    summary["Gamma2460_median"], summary["Gamma2460_p16"], summary["Gamma2460_p84"],
                )
            )
            lines.append(
                "  G2536={:+.4g} [{:+.4g},{:+.4g}] GeV^-1; "
                "Gamma2536={:.4g} [{:.4g},{:.4g}] keV".format(
                    summary["G2536_median"], summary["G2536_p16"], summary["G2536_p84"],
                    summary["Gamma2536_median"], summary["Gamma2536_p16"], summary["Gamma2536_p84"],
                )
            )
            lines.append(f"  clipped f1-target fraction: {clipped_fraction:.3f}")
            lines.append("")

    write_csv(OUT / "stage3_ds1_physical_residue_summary.csv", summary_rows)
    (OUT / "stage3_ds1_physical_residue_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
