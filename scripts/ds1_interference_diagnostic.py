"""Summarize constructive/destructive Ds1 amplitude interference.

The Stage-3 scan stores the basis-current amplitudes GA and GB together with
the sampled residues.  This diagnostic decomposes the physical couplings into
their normalized axial and tensor pieces,

    G2460 = (sin(theta) fA GA)/f1 + (cos(theta) fB GB)/f1,
    G2536 = (cos(theta) fA GA)/f2 - (sin(theta) fB GB)/f2.

The output makes the Bondar-Milstein cancellation argument visible in the
LCSR variables.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def summarize(values):
    arr = np.asarray(values, dtype=float)
    return {
        "median": float(np.percentile(arr, 50.0)),
        "p16": float(np.percentile(arr, 16.0)),
        "p84": float(np.percentile(arr, 84.0)),
    }


def read_scan():
    with (OUT / "stage3_ds1_physical_residue_scan.csv").open() as f:
        return list(csv.DictReader(f))


def decompose(row):
    theta = math.radians(float(row["theta_deg"]))
    s = math.sin(theta)
    c = math.cos(theta)
    ga = float(row["GA_basis"])
    gb = float(row["GB_basis"])
    f_a = float(row["fA_GeV"])
    f_b = float(row["fB_eff_GeV"])
    f1 = float(row["f1_GeV"])
    f2 = float(row["f2_GeV"])
    a2460 = s * f_a * ga / f1
    b2460 = c * f_b * gb / f1
    a2536 = c * f_a * ga / f2
    b2536 = -s * f_b * gb / f2
    return {
        "scenario": row["scenario"],
        "ensemble": row["ensemble"],
        "theta_deg": float(row["theta_deg"]),
        "A2460": a2460,
        "B2460": b2460,
        "G2460": a2460 + b2460,
        "A2536": a2536,
        "B2536": b2536,
        "G2536": a2536 + b2536,
        "abs_ratio_G2460_over_G2536": abs(a2460 + b2460) / max(abs(a2536 + b2536), 1.0e-12),
    }


def main():
    rows = [decompose(row) for row in read_scan()]
    with (OUT / "ds1_interference_diagnostic_scan.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []
    lines = [
        "Ds1 amplitude interference diagnostic",
        "======================================",
        "G2460 = A2460 + B2460; G2536 = A2536 + B2536.",
        "The 2536 pieces have opposite signs in the physical combination.",
        "",
    ]
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for ensemble in ("theta_fixed", "theta_scan_25_45"):
            subset = [r for r in rows if r["scenario"] == scenario and r["ensemble"] == ensemble]
            if not subset:
                continue
            out = {"scenario": scenario, "ensemble": ensemble, "n": len(subset)}
            for key in ("A2460", "B2460", "G2460", "A2536", "B2536", "G2536", "abs_ratio_G2460_over_G2536"):
                stats = summarize([r[key] for r in subset])
                out[f"{key}_median"] = stats["median"]
                out[f"{key}_p16"] = stats["p16"]
                out[f"{key}_p84"] = stats["p84"]
            summary_rows.append(out)
            lines.append(f"{scenario}, {ensemble}, n={len(subset)}")
            lines.append(
                "  2460: A={:+.4f} [{:+.4f},{:+.4f}], "
                "B={:+.4f} [{:+.4f},{:+.4f}], "
                "sum={:+.4f} [{:+.4f},{:+.4f}]".format(
                    out["A2460_median"], out["A2460_p16"], out["A2460_p84"],
                    out["B2460_median"], out["B2460_p16"], out["B2460_p84"],
                    out["G2460_median"], out["G2460_p16"], out["G2460_p84"],
                )
            )
            lines.append(
                "  2536: A={:+.4f} [{:+.4f},{:+.4f}], "
                "B={:+.4f} [{:+.4f},{:+.4f}], "
                "sum={:+.4f} [{:+.4f},{:+.4f}]".format(
                    out["A2536_median"], out["A2536_p16"], out["A2536_p84"],
                    out["B2536_median"], out["B2536_p16"], out["B2536_p84"],
                    out["G2536_median"], out["G2536_p16"], out["G2536_p84"],
                )
            )
            lines.append(
                "  |G2460/G2536| median={:.2f} [{:.2f},{:.2f}]".format(
                    out["abs_ratio_G2460_over_G2536_median"],
                    out["abs_ratio_G2460_over_G2536_p16"],
                    out["abs_ratio_G2460_over_G2536_p84"],
                )
            )
            lines.append("")

    with (OUT / "ds1_interference_diagnostic_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    (OUT / "ds1_interference_diagnostic_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
