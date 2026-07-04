"""Redo Borel/threshold stability scans using squared s0 notation.

The historical scans in this project often used ``sqrt(s0)`` as the scanned
variable.  For the publication pass we keep the user's preferred notation:
``s0`` is the squared continuum threshold in GeV^2.  Internally, some helper
functions still accept ``sqrt(s0)``; this script converts at the boundary and
writes only squared-threshold labels in the new outputs.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib.pyplot as plt

from bs1_stage2_baseline import FPERP_S_1GEV, evaluate_point as evaluate_bs_basis
from final_stage2_uncertainty_scan import precompute_f1_basis, set_photon_shape_params
from lattice_photon_normalization_comparison import (
    central_inputs_for_scenario,
    evaluate as evaluate_ds_basis,
)
from stage1_tensor_gb_soft2p_estimate import width_keV
from stage2_axial_g1_three_particle import F1_integral
from stage3_bs1_physical_decay_constants import physical_decay_constants
from stage3_ds1_physical_residue_normalization import calibrate_from_f1


THETA_DEG = 35.3
DS_F1_ANCHOR = 0.345
DS_FT = 0.256
BS_FA = 0.305
BS_FT = 0.285


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def relative_span(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    midpoint = 0.5 * (float(arr.max()) + float(arr.min()))
    if abs(midpoint) < 1.0e-14:
        return float("nan")
    return float((float(arr.max()) - float(arr.min())) / abs(midpoint))


def ds_physical_row(M2: float, s0: float, f1_axial: float, f1_basis: dict[str, float]) -> dict[str, float | str]:
    vals = central_inputs_for_scenario("lattice_fperp_s")
    basis = evaluate_ds_basis(vals, M2, math.sqrt(s0), THETA_DEG, f1_axial, f1_basis)
    theta = math.radians(THETA_DEG)
    s = math.sin(theta)
    c = math.cos(theta)
    fA = vals["f_ds1"]
    fB = DS_FT * vals["m_ds1"] / (vals["mc"] + vals["ms"])
    calibrated = calibrate_from_f1(theta, fA, fB, DS_F1_ANCHOR)
    if calibrated is None:
        raise RuntimeError("Ds physical-current calibration failed")
    rho, f1_phys, f2_phys, clipped = calibrated
    g2460 = (s * fA * basis["GA"] + c * fB * basis["GB"]) / f1_phys
    g2536 = (c * fA * basis["GA"] - s * fB * basis["GB"]) / f2_phys
    return {
        "sector": "Ds",
        "state": "Ds1(2460), Ds1(2536)",
        "scenario": "lattice_fperp_s",
        "M2": M2,
        "s0": s0,
        "theta_deg": THETA_DEG,
        "GA_basis": basis["GA"],
        "GB_basis": basis["GB"],
        "f1_GeV": f1_phys,
        "f2_GeV": f2_phys,
        "rho_AB": rho,
        "f1_target_clipped": int(clipped),
        "g_2460": g2460,
        "g_2536": g2536,
        "Gamma_2460_keV": width_keV(vals["m_ds1"], vals["m_ds"], g2460),
        "Gamma_2536_keV": width_keV(vals["m_ds1_2536"], vals["m_ds"], g2536),
        "chi_effective": basis["chi_effective"],
    }


def bs_inputs(m_initial: float) -> dict[str, float]:
    qq = -(0.240) ** 3
    ss = 0.8 * qq
    chi = FPERP_S_1GEV / ss
    return {
        "mc": 4.18,
        "ms": 0.093,
        "m_ds1": m_initial,
        "m_ds": 5.36692,
        "f_ds1": BS_FA,
        "f_ds": 0.2303,
        "fT": BS_FT,
        "ss": ss,
        "chi": chi,
        "f3g": -0.0039,
        "ec": -1.0 / 3.0,
        "es": -1.0 / 3.0,
        "fperp_s_used": FPERP_S_1GEV,
    }


def bs_physical_row(
    target: dict[str, float | str],
    M2: float,
    s0: float,
    f1_axial: float,
    f1_basis: dict[str, float],
) -> dict[str, float | str]:
    vals = bs_inputs(float(target["mass"]))
    set_photon_shape_params(-2.1, 3.8)
    basis = evaluate_bs_basis(vals, M2, math.sqrt(s0), THETA_DEG, f1_axial, f1_basis)
    theta = math.radians(THETA_DEG)
    constants = physical_decay_constants(theta, vals["f_ds1"], basis["fB_effective"])
    if constants is None:
        raise RuntimeError(f"Bs physical-current normalization failed for {target['state']}")
    rho, f1_phys, f2_phys = constants
    s = math.sin(theta)
    c = math.cos(theta)
    g_low = (s * vals["f_ds1"] * basis["GA"] + c * basis["fB_effective"] * basis["GB"]) / f1_phys
    g_high = (c * vals["f_ds1"] * basis["GA"] - s * basis["fB_effective"] * basis["GB"]) / f2_phys
    quoted = str(target["quoted_combo"])
    g_quoted = g_high if quoted == "high" else g_low
    gamma_quoted = width_keV(vals["m_ds1"], vals["m_ds"], g_quoted)
    return {
        "sector": "Bs",
        "state": str(target["state"]),
        "scenario": "lattice_fperp_s",
        "M2": M2,
        "s0": s0,
        "theta_deg": THETA_DEG,
        "quoted_combo": quoted,
        "GA_basis": basis["GA"],
        "GB_basis": basis["GB"],
        "f1_GeV": f1_phys,
        "f2_GeV": f2_phys,
        "rho_AB": rho,
        "g_low": g_low,
        "g_high": g_high,
        "g_quoted": g_quoted,
        "Gamma_low_keV": width_keV(vals["m_ds1"], vals["m_ds"], g_low),
        "Gamma_high_keV": width_keV(vals["m_ds1"], vals["m_ds"], g_high),
        "Gamma_quoted_keV": gamma_quoted,
        "chi_effective": vals["chi"],
    }


def make_ds_plot(rows: list[dict[str, float | str]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.6), sharex=True)
    panels = [
        ("g_2460", r"$g_1$ [GeV$^{-1}$]"),
        ("Gamma_2460_keV", r"$\Gamma[D_{s1}(2460)\to D_s\gamma]$ [keV]"),
        ("g_2536", r"$g_2$ [GeV$^{-1}$]"),
        ("Gamma_2536_keV", r"$\Gamma[D_{s1}(2536)\to D_s\gamma]$ [keV]"),
    ]
    colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    s0_values = sorted({float(r["s0"]) for r in rows})
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for color, s0 in zip(colors, s0_values):
            subset = [r for r in rows if abs(float(r["s0"]) - s0) < 1.0e-10]
            ax.plot(
                [float(r["M2"]) for r in subset],
                [float(r[key]) for r in subset],
                color=color,
                linewidth=1.6,
                label=rf"$s_0={s0:.2f}\,\mathrm{{GeV}}^2$",
            )
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.28, linewidth=0.7)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    axes[0, 1].legend(frameon=False, fontsize=8.0)
    fig.suptitle(r"$D_{s1}\to D_s\gamma$ stability, lattice photon normalization")
    fig.tight_layout()
    fig.savefig(OUT / "redo_ds1_stability_windows.png", dpi=240)
    fig.savefig(OUT / "redo_ds1_stability_windows.pdf")


def make_bs_plot(rows: list[dict[str, float | str]], target_name: str, filename_base: str) -> None:
    subset_all = [r for r in rows if r["state"] == target_name]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.6), sharex=True)
    panels = [
        ("g_quoted", r"$g_{\rm quoted}$ [GeV$^{-1}$]"),
        ("Gamma_quoted_keV", r"$\Gamma_{\rm quoted}$ [keV]"),
        ("GA_basis", r"$g_A^{\rm basis}$ [GeV$^{-1}$]"),
        ("GB_basis", r"$g_B^{\rm basis}$ [GeV$^{-1}$]"),
    ]
    colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    s0_values = sorted({float(r["s0"]) for r in subset_all})
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for color, s0 in zip(colors, s0_values):
            subset = [r for r in subset_all if abs(float(r["s0"]) - s0) < 1.0e-10]
            ax.plot(
                [float(r["M2"]) for r in subset],
                [float(r[key]) for r in subset],
                color=color,
                linewidth=1.6,
                label=rf"$s_0={s0:.2f}\,\mathrm{{GeV}}^2$",
            )
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.28, linewidth=0.7)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    axes[0, 1].legend(frameon=False, fontsize=8.0)
    fig.suptitle(rf"${target_name}\to B_s\gamma$ stability, lattice photon normalization")
    fig.tight_layout()
    fig.savefig(OUT / f"{filename_base}.png", dpi=240)
    fig.savefig(OUT / f"{filename_base}.pdf")


def add_summary_block(lines: list[str], title: str, rows: list[dict[str, float | str]], keys: list[str]) -> None:
    lines.append(title)
    lines.append("-" * len(title))
    lines.append(f"grid points: {len(rows)}")
    lines.append(
        "M2 range: {:.2f}-{:.2f} GeV^2; s0 range: {:.2f}-{:.2f} GeV^2".format(
            min(float(r["M2"]) for r in rows),
            max(float(r["M2"]) for r in rows),
            min(float(r["s0"]) for r in rows),
            max(float(r["s0"]) for r in rows),
        )
    )
    for key in keys:
        values = [float(r[key]) for r in rows]
        lines.append(
            "{}: min {:+.6g}; max {:+.6g}; relative span {:.3f}".format(
                key,
                min(values),
                max(values),
                relative_span(values),
            )
        )
    lines.append("")


def main() -> None:
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()

    ds_s0_values = [2.85**2, 2.95**2, 3.05**2, 3.10**2]
    ds_rows: list[dict[str, float | str]] = []
    for s0 in ds_s0_values:
        for M2 in np.linspace(2.0, 7.0, 31):
            ds_rows.append(ds_physical_row(float(M2), float(s0), f1_axial, f1_basis))
    write_csv(OUT / "redo_ds1_stability_windows.csv", ds_rows)
    make_ds_plot(ds_rows)

    targets = [
        {
            "state": "B_{s1}(5750)",
            "mass": 5.750,
            "quoted_combo": "low",
            "s0_values": [6.05**2, 6.15**2, 6.25**2, 6.35**2],
            "filename": "redo_bs1_low_stability_windows",
        },
        {
            "state": "B_{s1}(5830)",
            "mass": 5.82870,
            "quoted_combo": "high",
            "s0_values": [6.15**2, 6.25**2, 6.35**2, 6.40**2],
            "filename": "redo_bs1_5830_stability_windows",
        },
    ]
    bs_rows: list[dict[str, float | str]] = []
    for target in targets:
        for s0 in target["s0_values"]:
            for M2 in np.linspace(6.0, 18.0, 31):
                bs_rows.append(bs_physical_row(target, float(M2), float(s0), f1_axial, f1_basis))
    write_csv(OUT / "redo_bs1_stability_windows.csv", bs_rows)
    for target in targets:
        make_bs_plot(bs_rows, str(target["state"]), str(target["filename"]))

    lines = [
        "Redo stability-window scan",
        "==========================",
        "Preferred central scenario: lattice f_gamma,s^perp normalization.",
        "Thresholds are written as squared s0 in GeV^2; helper functions receive sqrt(s0) internally.",
        f"theta = {THETA_DEG:.1f} deg.",
        "",
        "Suggested plotting/selection convention:",
        "  Ds broad scan: M2 = 2-7 GeV^2, s0 = 8.12-9.61 GeV^2.",
        "  Bs broad scan: M2 = 6-18 GeV^2, s0 chosen above each pole.",
        "  Final uncertainty windows should be selected from the visible plateau.",
        "",
    ]
    add_summary_block(
        lines,
        "Ds1 scan",
        ds_rows,
        ["g_2460", "Gamma_2460_keV", "g_2536", "Gamma_2536_keV", "GA_basis", "GB_basis"],
    )
    for target in targets:
        subset = [r for r in bs_rows if r["state"] == target["state"]]
        add_summary_block(
            lines,
            f"{target['state']} scan",
            subset,
            ["g_quoted", "Gamma_quoted_keV", "GA_basis", "GB_basis"],
        )
    lines.extend(
        [
            "Wrote:",
            f"  {OUT / 'redo_ds1_stability_windows.csv'}",
            f"  {OUT / 'redo_ds1_stability_windows.pdf'}",
            f"  {OUT / 'redo_bs1_stability_windows.csv'}",
            f"  {OUT / 'redo_bs1_low_stability_windows.pdf'}",
            f"  {OUT / 'redo_bs1_5830_stability_windows.pdf'}",
        ]
    )
    summary_path = OUT / "redo_stability_windows_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
