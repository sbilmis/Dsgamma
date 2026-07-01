"""Additional publication figures for the Ds1 -> Ds gamma analysis."""

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

from final_stage2_uncertainty_scan import matched_f1_integral, precompute_f1_basis
from stage1_axial_g1_baseline import central_inputs, g1_stage1
from stage1_tensor_gb_hard_candidate import hard_tensor_gb
from stage1_tensor_gb_soft2p_estimate import gb_soft_tensor_da, width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from stage2_stability_plots import evaluate as evaluate_full
from stage2_tensor_gb_three_particle_estimate import gb_three_particle_from_integral


COLORS = {
    "blue": "#1f77b4",
    "green": "#2ca02c",
    "red": "#d62728",
    "orange": "#ff7f0e",
    "purple": "#9467bd",
    "gray": "#7f7f7f",
}


def central_setup():
    inputs = central_inputs()
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    f1_tensor = matched_f1_integral(inputs, f1_basis)
    fT = 0.256
    fB = fT * inputs["m_ds1"] / (inputs["mc"] + inputs["ms"])
    return inputs, f1_axial, f1_basis, f1_tensor, fB


def evaluate_components(M2=4.5, s0_root=2.55):
    inputs, f1_axial, _, f1_tensor, fB = central_setup()
    s0 = s0_root * s0_root
    axial1 = g1_stage1(M2, s0, **inputs)
    axial2 = g1_stage2(M2, s0, inputs, f1_axial)
    gb_hard, _ = hard_tensor_gb(M2, s0, inputs, fB)
    gb_soft = gb_soft_tensor_da(axial2, inputs, fB)
    gb_3p = gb_three_particle_from_integral(axial2, inputs, fB, f1_tensor)
    ga_stage1 = axial1["g1_GeV_inv"]
    ga_3p = axial2["F1_delta_g1"]
    return {
        "GA_stage1": ga_stage1,
        "GA_3p": ga_3p,
        "GA_total": ga_stage1 + ga_3p,
        "GB_hard": gb_hard,
        "GB_soft2p": gb_soft,
        "GB_3p": gb_3p,
        "GB_total": gb_hard + gb_soft + gb_3p,
    }


def threshold_plot(f1_axial, f1_basis):
    rows = []
    s0_values = np.linspace(2.50, 2.60, 31)
    M2_values = (3.5, 4.5, 5.5)
    for M2 in M2_values:
        for s0_root in s0_values:
            row = evaluate_full(float(M2), float(s0_root), f1_axial, f1_basis)
            row["plot_M2"] = M2
            rows.append(row)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), sharex=True)
    for M2 in M2_values:
        subset = [r for r in rows if abs(r["plot_M2"] - M2) < 1e-12]
        axes[0].plot(
            [r["s0_root"] for r in subset],
            [r["Gamma2460_keV"] for r in subset],
            marker="o",
            markersize=2.5,
            linewidth=1.6,
            label=rf"$M^2={M2:.1f}$ GeV$^2$",
        )
        axes[1].plot(
            [r["s0_root"] for r in subset],
            [r["Gamma2536_keV"] for r in subset],
            marker="o",
            markersize=2.5,
            linewidth=1.6,
            label=rf"$M^2={M2:.1f}$ GeV$^2$",
        )
    axes[0].set_ylabel(r"$\Gamma_{2460}$ [keV]")
    axes[1].set_ylabel(r"$\Gamma_{2536}$ [keV]")
    for ax in axes:
        ax.set_xlabel(r"$\sqrt{s_0}$ [GeV]")
        ax.grid(True, alpha=0.28, linewidth=0.7)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle(r"Continuum-threshold dependence, $\theta=35.3^\circ$")
    fig.tight_layout()
    fig.savefig(OUT / "stage2_threshold_dependence.png", dpi=220)
    fig.savefig(OUT / "stage2_threshold_dependence.pdf")
    write_csv(OUT / "stage2_threshold_dependence.csv", rows)
    return rows


def mixing_plot(f1_axial, f1_basis):
    rows = []
    theta_values = np.linspace(20.0, 50.0, 151)
    M2 = 4.5
    s0_root = 2.55
    for theta in theta_values:
        row = evaluate_full(M2, s0_root, f1_axial, f1_basis, theta_deg=float(theta))
        row["theta_deg"] = float(theta)
        rows.append(row)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    axes[0].plot(
        [r["theta_deg"] for r in rows],
        [r["Gamma2460_keV"] for r in rows],
        color=COLORS["blue"],
        linewidth=1.8,
        label=r"$D_{s1}(2460)$",
    )
    axes[0].plot(
        [r["theta_deg"] for r in rows],
        [r["Gamma2536_keV"] for r in rows],
        color=COLORS["red"],
        linewidth=1.8,
        label=r"$D_{s1}(2536)$",
    )
    axes[1].plot(
        [r["theta_deg"] for r in rows],
        [r["G2460"] for r in rows],
        color=COLORS["blue"],
        linewidth=1.8,
        label=r"$G_{2460}$",
    )
    axes[1].plot(
        [r["theta_deg"] for r in rows],
        [r["G2536"] for r in rows],
        color=COLORS["red"],
        linewidth=1.8,
        label=r"$G_{2536}$",
    )
    for ax in axes:
        ax.axvline(35.3, color=COLORS["gray"], linestyle="--", linewidth=1.1)
        ax.set_xlabel(r"$\theta$ [deg]")
        ax.grid(True, alpha=0.28, linewidth=0.7)
        ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel(r"$\Gamma$ [keV]")
    axes[1].set_ylabel(r"$G$ [GeV$^{-1}$]")
    fig.suptitle(r"Mixing-angle dependence at $M^2=4.5$ GeV$^2$, $\sqrt{s_0}=2.55$ GeV")
    fig.tight_layout()
    fig.savefig(OUT / "stage2_mixing_angle_dependence.png", dpi=220)
    fig.savefig(OUT / "stage2_mixing_angle_dependence.pdf")
    write_csv(OUT / "stage2_mixing_angle_dependence.csv", rows)
    return rows


def contribution_breakdown_plot():
    comp = evaluate_components()
    theta = math.radians(35.3)
    rows = []
    pieces = [
        ("A Stage 1", comp["GA_stage1"], math.sin(theta) * comp["GA_stage1"], math.cos(theta) * comp["GA_stage1"]),
        ("A 3p", comp["GA_3p"], math.sin(theta) * comp["GA_3p"], math.cos(theta) * comp["GA_3p"]),
        ("B hard", comp["GB_hard"], math.cos(theta) * comp["GB_hard"], -math.sin(theta) * comp["GB_hard"]),
        ("B 2p DA", comp["GB_soft2p"], math.cos(theta) * comp["GB_soft2p"], -math.sin(theta) * comp["GB_soft2p"]),
        ("B 3p DA", comp["GB_3p"], math.cos(theta) * comp["GB_3p"], -math.sin(theta) * comp["GB_3p"]),
    ]
    for label, raw, g2460, g2536 in pieces:
        rows.append({"piece": label, "raw_current_G": raw, "G2460_piece": g2460, "G2536_piece": g2536})

    labels = [r["piece"] for r in rows]
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.2, 4.3))
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.bar(x - width / 2, [r["G2460_piece"] for r in rows], width, label=r"$G_{2460}$")
    ax.bar(x + width / 2, [r["G2536_piece"] for r in rows], width, label=r"$G_{2536}$")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel(r"Contribution to physical $G$ [GeV$^{-1}$]")
    ax.set_title(r"Central contribution breakdown, $M^2=4.5$ GeV$^2$, $\sqrt{s_0}=2.55$ GeV")
    ax.grid(True, axis="y", alpha=0.28, linewidth=0.7)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "stage2_contribution_breakdown.png", dpi=220)
    fig.savefig(OUT / "stage2_contribution_breakdown.pdf")
    write_csv(OUT / "stage2_contribution_breakdown.csv", rows)
    return rows


def literature_comparison_plot():
    rows = [
        {"method": "This work", "low": 27.6, "central": 38.3, "high": 50.7, "kind": "band"},
        {"method": "Colangelo LCSR", "low": 19.0, "central": 24.0, "high": 29.0, "kind": "band"},
        {"method": "VMD", "low": 3.3, "central": 3.3, "high": 3.3, "kind": "point"},
        {"method": "QM [15]", "low": 6.2, "central": 6.2, "high": 6.2, "kind": "point"},
        {"method": "QM [11]", "low": 5.08, "central": 5.08, "high": 5.08, "kind": "point"},
    ]
    y = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    for i, row in enumerate(rows):
        lo = row["central"] - row["low"]
        hi = row["high"] - row["central"]
        ax.errorbar(
            row["central"],
            i,
            xerr=np.array([[lo], [hi]]),
            fmt="o",
            color=COLORS["blue"] if row["method"] == "This work" else COLORS["gray"],
            capsize=4 if row["kind"] == "band" else 0,
            markersize=6,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([r["method"] for r in rows])
    ax.invert_yaxis()
    ax.set_xlabel(r"$\Gamma[D_{s1}(2460)\to D_s\gamma]$ [keV]")
    ax.set_title("Comparison with representative literature estimates")
    ax.grid(True, axis="x", alpha=0.28, linewidth=0.7)
    fig.tight_layout()
    fig.savefig(OUT / "literature_comparison_2460.png", dpi=220)
    fig.savefig(OUT / "literature_comparison_2460.pdf")
    write_csv(OUT / "literature_comparison_2460.csv", rows)
    return rows


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    _, f1_axial, f1_basis, _, _ = central_setup()
    threshold_rows = threshold_plot(f1_axial, f1_basis)
    mixing_rows = mixing_plot(f1_axial, f1_basis)
    breakdown_rows = contribution_breakdown_plot()
    literature_rows = literature_comparison_plot()

    def range_line(rows, key):
        values = np.array([float(r[key]) for r in rows])
        return f"{key}: {values.min():.4g} to {values.max():.4g}"

    lines = [
        "Publication figures",
        "===================",
        "Generated threshold, mixing-angle, contribution-breakdown, and literature-comparison figures.",
        range_line(threshold_rows, "Gamma2460_keV"),
        range_line(threshold_rows, "Gamma2536_keV"),
        range_line(mixing_rows, "Gamma2460_keV"),
        range_line(mixing_rows, "Gamma2536_keV"),
        "Contribution breakdown at M2=4.5 GeV^2, sqrt(s0)=2.55 GeV:",
    ]
    for row in breakdown_rows:
        lines.append(
            f"  {row['piece']}: G2460 {row['G2460_piece']:+.5f}, "
            f"G2536 {row['G2536_piece']:+.5f}"
        )
    lines.append("Literature comparison entries:")
    for row in literature_rows:
        lines.append(f"  {row['method']}: {row['low']} to {row['high']} keV")
    path = OUT / "publication_plots_summary.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
