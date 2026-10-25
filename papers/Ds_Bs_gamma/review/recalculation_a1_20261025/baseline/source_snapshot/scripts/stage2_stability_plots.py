"""Central Stage-2 stability plots in the Borel/threshold window."""

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
from stage1_axial_g1_baseline import central_inputs
from stage1_tensor_gb_hard_candidate import hard_tensor_gb
from stage1_tensor_gb_soft2p_estimate import gb_soft_two_particle, width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from stage2_tensor_gb_three_particle_estimate import gb_em_da, gb_three_particle_from_integral


def evaluate(M2, s0_root, f1_axial, f1_basis, theta_deg=35.3):
    inputs = central_inputs()
    f1_tensor = matched_f1_integral(inputs, f1_basis)
    fT = 0.256
    fB = fT * inputs["m_ds1"] / (inputs["mc"] + inputs["ms"])
    s0 = s0_root * s0_root
    theta = math.radians(theta_deg)

    axial = g1_stage2(M2, s0, inputs, f1_axial)
    GA = axial["g1_stage2_GeV_inv"]
    GB_hard, _ = hard_tensor_gb(M2, s0, inputs, fB)
    GB_soft = gb_soft_two_particle(axial, inputs, fB)
    GB_3p = gb_three_particle_from_integral(axial, inputs, fB, f1_tensor)
    GB_em, _ = gb_em_da(axial, inputs, fB)
    GB = GB_hard + GB_soft + GB_3p + GB_em
    G2460 = math.sin(theta) * GA + math.cos(theta) * GB
    G2536 = math.cos(theta) * GA - math.sin(theta) * GB
    return {
        "M2": M2,
        "s0_root": s0_root,
        "GA": GA,
        "GB": GB,
        "GB_em": GB_em,
        "G2460": G2460,
        "G2536": G2536,
        "Gamma2460_keV": width_keV(inputs["m_ds1"], inputs["m_ds"], G2460),
        "Gamma2536_keV": width_keV(2.53511, inputs["m_ds"], G2536),
    }


def relative_span(values):
    values = np.asarray(values, dtype=float)
    midpoint = 0.5 * (values.max() + values.min())
    if midpoint == 0:
        return float("nan")
    return float((values.max() - values.min()) / abs(midpoint))


def main():
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows = []
    for s0_root in (2.50, 2.55, 2.60):
        for M2 in np.linspace(3.0, 6.0, 25):
            rows.append(evaluate(float(M2), s0_root, f1_axial, f1_basis))

    csv_path = OUT / "stage2_stability_central.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.5), sharex=True)
    panels = [
        ("G2460", r"$G_{2460}$ [GeV$^{-1}$]"),
        ("Gamma2460_keV", r"$\Gamma_{2460}$ [keV]"),
        ("G2536", r"$G_{2536}$ [GeV$^{-1}$]"),
        ("Gamma2536_keV", r"$\Gamma_{2536}$ [keV]"),
    ]
    colors = {2.50: "#1f77b4", 2.55: "#2ca02c", 2.60: "#d62728"}
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for s0_root in (2.50, 2.55, 2.60):
            subset = [r for r in rows if abs(r["s0_root"] - s0_root) < 1e-12]
            ax.plot(
                [r["M2"] for r in subset],
                [r[key] for r in subset],
                marker="o",
                markersize=2.8,
                linewidth=1.5,
                color=colors[s0_root],
                label=rf"$\sqrt{{s_0}}={s0_root:.2f}$ GeV",
            )
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.28, linewidth=0.7)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    axes[0, 1].legend(frameon=False, fontsize=9)
    fig.suptitle(r"Central matched Stage-2 stability, $\theta=35.3^\circ$")
    fig.tight_layout()
    png_path = OUT / "stage2_stability_central.png"
    pdf_path = OUT / "stage2_stability_central.pdf"
    fig.savefig(png_path, dpi=220)
    fig.savefig(pdf_path)

    lines = [
        "Central matched Stage-2 stability",
        "=================================",
        "Inputs fixed to central values; theta = 35.3 deg.",
        "Curves use sqrt(s0) = 2.50, 2.55, 2.60 GeV and M^2 in [3,6] GeV^2.",
    ]
    for key in ("G2460", "Gamma2460_keV", "G2536", "Gamma2536_keV", "GA", "GB"):
        values = np.array([r[key] for r in rows])
        lines.append(
            f"{key}: min {values.min():+.6g}; max {values.max():+.6g}; "
            f"relative span {relative_span(values):.3f}"
        )
    lines.append(f"Wrote table to {csv_path}")
    lines.append(f"Wrote figures to {png_path} and {pdf_path}")
    summary_path = OUT / "stage2_stability_central_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
