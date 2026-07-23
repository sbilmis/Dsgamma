"""Audit how the numerical mixed-state decay constants are obtained.

This is deliberately separate from the transition-form-factor calculation.
It records whether f1 and f2 are external inputs, closure-derived quantities,
or genuine predictions of a completed two-point QCD sum rule.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

from mixing_angle_inputs import THETA_BS_DEG, THETA_DS_DEG
from stage3_ds1_physical_residue_normalization import (
    F1_STAGE3_ANCHOR,
    calibrate_from_f1,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


def main() -> None:
    theta_legacy = math.radians(35.3)

    # Ds: f1 is an external anchor; f2 follows from the overlap model.
    f_a_ds = 0.225
    f_t_ds = 0.256
    f_b_ds = f_t_ds * 2.4595 / (1.27 + 0.093)
    ds = calibrate_from_f1(theta_legacy, f_a_ds, f_b_ds, F1_STAGE3_ANCHOR)
    if ds is None:
        raise RuntimeError("Ds central overlap calibration failed")
    rho_ds, f1_ds, f2_ds, clipped = ds

    with (OUT / "twopoint_ds1_matrix_mc.csv").open() as f:
        ds_matrix_rows = [
            row for row in csv.DictReader(f) if int(row["accepted"]) == 1
        ]
    theta_matrix = sorted(float(row["theta_deg"]) for row in ds_matrix_rows)
    theta_diagnostic = sorted(float(row["theta_matrix_deg"]) for row in ds_matrix_rows)
    pi12_residual = sorted(float(row["Pi12_normalized"]) for row in ds_matrix_rows)
    f1_matrix = sorted(float(row["f1_GeV"]) for row in ds_matrix_rows)
    f2_matrix = sorted(float(row["f2_GeV"]) for row in ds_matrix_rows)

    def percentile(values: list[float], fraction: float) -> float:
        position = fraction * (len(values) - 1)
        lo = int(position)
        hi = min(lo + 1, len(values) - 1)
        weight = position - lo
        return values[lo] * (1.0 - weight) + values[hi] * weight

    # Bs: the two basis constants are inputs; f1/f2 follow from Pi12=0.
    f_a_bs = 0.305
    f_t_bs = 0.285
    f_b_bs = f_t_bs * 5.82870 / (4.18 + 0.093)
    with (OUT / "final_window_mc_summary.csv").open() as f:
        bs_rows = [
            row for row in csv.DictReader(f)
            if row["sector"] == "Bs"
            and row["state"] == "B_{s1}(5830)"
            and row["scenario"] == "lattice_fperp_s"
            and row["ensemble"] == "theta_prior_gaussian"
            and row["window_id"] == "central_10_14"
        ]
    if len(bs_rows) != 1:
        raise RuntimeError("preferred Bs summary row not found")
    bs_row = bs_rows[0]

    lines = [
        "Decay-constant provenance audit",
        "================================",
        "",
        "Ds1 sector",
        "  preferred method: normalized-current AA/AB/BB two-point QCDSR",
        "  OPE: exact-mass LO perturbative + local d=3 through ms^2 + local d=5",
        "  complete matrix at the stated truncation used: yes",
        "  external f1/f2 or overlap input used: no",
        "  external theta input = {:.6f} [{:.6f},{:.6f}] deg".format(
            percentile(theta_matrix, 0.50), percentile(theta_matrix, 0.16),
            percentile(theta_matrix, 0.84)),
        "  truncated-matrix diagnostic theta = {:.6f} [{:.6f},{:.6f}] deg".format(
            percentile(theta_diagnostic, 0.50), percentile(theta_diagnostic, 0.16),
            percentile(theta_diagnostic, 0.84)),
        "  normalized Pi12 residual = {:.6f} [{:.6f},{:.6f}]".format(
            percentile(pi12_residual, 0.50), percentile(pi12_residual, 0.16),
            percentile(pi12_residual, 0.84)),
        "  f1 = {:.6f} [{:.6f},{:.6f}] GeV".format(
            percentile(f1_matrix, 0.50), percentile(f1_matrix, 0.16),
            percentile(f1_matrix, 0.84)),
        "  f2 = {:.6f} [{:.6f},{:.6f}] GeV".format(
            percentile(f2_matrix, 0.50), percentile(f2_matrix, 0.16),
            percentile(f2_matrix, 0.84)),
        "",
        "Ds1 legacy comparison only",
        "  method: external f1 anchor + two-point overlap algebra",
        f"  theta = 35.3 deg; fA = {f_a_ds:.6f} GeV; fB = {f_b_ds:.6f} GeV",
        f"  input f1 = {f1_ds:.6f} GeV",
        f"  derived rho_AB = {rho_ds:+.6f}",
        f"  derived f2 = {f2_ds:.6f} GeV",
        f"  anchor clipping used: {clipped}",
        "",
        "Bs1 sector",
        "  method: basis-constant inputs + Pi12=0 diagonal closure",
        "  complete local two-point AA/BB/AB OPE used: no",
        f"  theta = {THETA_BS_DEG:.1f} deg; fA = {f_a_bs:.6f} GeV; fB = {f_b_bs:.6f} GeV",
        "  the central basis inputs alone do not always admit |rho_AB|<=1 at this angle",
        "  invalid closure samples are rejected in the Monte Carlo",
        "  accepted closure f1 = {} [{},{}] GeV".format(
            bs_row["f1_median_GeV"], bs_row["f1_p16_GeV"], bs_row["f1_p84_GeV"]),
        "  accepted closure f2 = {} [{},{}] GeV".format(
            bs_row["f2_median_GeV"], bs_row["f2_p16_GeV"], bs_row["f2_p84_GeV"]),
        "",
        "Conclusion",
        f"  Ds uses external theta={THETA_DS_DEG:.1f} deg and independent matrix-QCDSR f1/f2 projections",
        "  at LO+d3+d5.  The matrix angle is a diagnostic, not the nominal input.",
        "  The Bs normalization remains a basis-input closure at the external angle.",
        "  Local condensates are included only in the ordinary two-point OPE;",
        "  this is independent of their exclusion from the photon transition LCSR.",
    ]
    path = OUT / "decay_constant_provenance_audit.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
