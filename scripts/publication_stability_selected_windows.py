"""Publication stability plots for the selected final windows.

The figure shows the two standard stability demonstrations:
  - vary M^2 at selected fixed s0 values;
  - vary s0 at selected fixed M^2 values.

We plot the physical-current form factors because these are the direct LCSR
outputs; widths can be obtained from the final phase-space formula.
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

from redo_stability_windows import bs_physical_row, ds_physical_row
from stage2_axial_g1_three_particle import F1_integral
from final_stage2_uncertainty_scan import precompute_f1_basis


DS_WINDOW = {
    "M2": (3.0, 4.5),
    "s0": (7.84, 8.70),
    "s0_lines": (7.84, 8.27, 8.70),
    "M2_lines": (3.0, 3.75, 4.5),
}

BS5830_WINDOW = {
    "state": "B_{s1}(5830)",
    "mass": 5.82870,
    "quoted_combo": "high",
    "M2": (10.0, 14.0),
    "s0": (37.82, 40.96),
    "s0_lines": (37.82, 39.39, 40.96),
    "M2_lines": (10.0, 12.0, 14.0),
}

BSLOW_WINDOW = {
    "state": "B_{s1}^{low}(5750)",
    "mass": 5.750,
    "quoted_combo": "low",
    "M2": (10.0, 14.0),
    "s0": (36.60, 40.32),
    "s0_lines": (36.60, 38.46, 40.32),
    "M2_lines": (10.0, 12.0, 14.0),
}

COLORS = {
    "low": "#1f77b4",
    "mid": "#2ca02c",
    "high": "#d62728",
    "secondary": "#9467bd",
    "gray": "#7f7f7f",
}


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


def ds_eval(M2: float, s0: float, f1_axial: float, f1_basis: dict[str, float]) -> dict[str, float | str]:
    row = ds_physical_row(M2, s0, f1_axial, f1_basis)
    row["g1_abs"] = abs(float(row["g_2460"]))
    row["g2_abs"] = abs(float(row["g_2536"]))
    return row


def bs_eval(target: dict[str, float | str], M2: float, s0: float, f1_axial: float, f1_basis: dict[str, float]) -> dict[str, float | str]:
    row = bs_physical_row(target, M2, s0, f1_axial, f1_basis)
    row["g_quoted_abs"] = abs(float(row["g_quoted"]))
    return row


def build_rows(f1_axial: float, f1_basis: dict[str, float]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []

    for s0 in DS_WINDOW["s0_lines"]:
        for M2 in np.linspace(*DS_WINDOW["M2"], 31):
            row = ds_eval(float(M2), float(s0), f1_axial, f1_basis)
            row.update({"scan": "M2", "fixed_s0": s0, "fixed_M2": ""})
            rows.append(row)
    for M2 in DS_WINDOW["M2_lines"]:
        for s0 in np.linspace(*DS_WINDOW["s0"], 31):
            row = ds_eval(float(M2), float(s0), f1_axial, f1_basis)
            row.update({"scan": "s0", "fixed_s0": "", "fixed_M2": M2})
            rows.append(row)

    for target in (BSLOW_WINDOW, BS5830_WINDOW):
        for s0 in target["s0_lines"]:
            for M2 in np.linspace(*target["M2"], 31):
                row = bs_eval(target, float(M2), float(s0), f1_axial, f1_basis)
                row.update({"scan": "M2", "fixed_s0": s0, "fixed_M2": ""})
                rows.append(row)
        for M2 in target["M2_lines"]:
            for s0 in np.linspace(*target["s0"], 31):
                row = bs_eval(target, float(M2), float(s0), f1_axial, f1_basis)
                row.update({"scan": "s0", "fixed_s0": "", "fixed_M2": M2})
                rows.append(row)
    return rows


def style_axis(ax) -> None:
    ax.grid(True, alpha=0.26, linewidth=0.7)
    ax.tick_params(direction="in", top=True, right=True)


def plot_publication(rows: list[dict[str, float | str]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2))
    color_cycle = [COLORS["low"], COLORS["mid"], COLORS["high"]]

    # Ds: vary M2 at selected s0.
    ax = axes[0, 0]
    for color, s0 in zip(color_cycle, DS_WINDOW["s0_lines"]):
        subset = [r for r in rows if r["sector"] == "Ds" and r["scan"] == "M2" and abs(float(r["fixed_s0"]) - s0) < 1.0e-10]
        ax.plot([float(r["M2"]) for r in subset], [float(r["g1_abs"]) for r in subset], color=color, linewidth=1.8, label=rf"$s_0={s0:.2f}$")
        ax.plot([float(r["M2"]) for r in subset], [float(r["g2_abs"]) for r in subset], color=color, linewidth=1.4, linestyle="--")
    ax.set_title(r"$D_{s1}\to D_s\gamma$: Borel stability")
    ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    ax.set_ylabel(r"$|g|$ [GeV$^{-1}$]")
    ax.legend(frameon=False, fontsize=8, title=r"$s_0$ [GeV$^2$]")
    ax.text(0.03, 0.08, r"solid: $2460$" "\n" r"dashed: $2536$", transform=ax.transAxes, fontsize=8)
    style_axis(ax)

    # Ds: vary s0 at selected M2.
    ax = axes[0, 1]
    for color, M2 in zip(color_cycle, DS_WINDOW["M2_lines"]):
        subset = [r for r in rows if r["sector"] == "Ds" and r["scan"] == "s0" and abs(float(r["fixed_M2"]) - M2) < 1.0e-10]
        ax.plot([float(r["s0"]) for r in subset], [float(r["g1_abs"]) for r in subset], color=color, linewidth=1.8, label=rf"$M^2={M2:.2f}$")
        ax.plot([float(r["s0"]) for r in subset], [float(r["g2_abs"]) for r in subset], color=color, linewidth=1.4, linestyle="--")
    ax.set_title(r"$D_{s1}\to D_s\gamma$: threshold stability")
    ax.set_xlabel(r"$s_0$ [GeV$^2$]")
    ax.set_ylabel(r"$|g|$ [GeV$^{-1}$]")
    ax.legend(frameon=False, fontsize=8, title=r"$M^2$ [GeV$^2$]")
    style_axis(ax)

    # Bs: vary M2 at selected s0.
    ax = axes[1, 0]
    for target, linestyle, label in ((BSLOW_WINDOW, "-", r"$B_{s1}^{low}$"), (BS5830_WINDOW, "--", r"$B_{s1}(5830)$")):
        for color, s0 in zip(color_cycle, target["s0_lines"]):
            subset = [
                r for r in rows
                if r["sector"] == "Bs" and r["state"] == target["state"] and r["scan"] == "M2"
                and abs(float(r["fixed_s0"]) - s0) < 1.0e-10
            ]
            ax.plot([float(r["M2"]) for r in subset], [float(r["g_quoted_abs"]) for r in subset], color=color, linewidth=1.7, linestyle=linestyle)
    ax.set_title(r"$B_{s1}\to B_s\gamma$: Borel stability")
    ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    ax.set_ylabel(r"$|g_{\rm quoted}|$ [GeV$^{-1}$]")
    ax.text(0.03, 0.78, r"solid: $B_{s1}^{low}$" "\n" r"dashed: $B_{s1}(5830)$", transform=ax.transAxes, fontsize=8)
    handles = [plt.Line2D([0], [0], color=c, lw=1.8) for c in color_cycle]
    labels = [rf"$s_0={s0:.2f}$" for s0 in BS5830_WINDOW["s0_lines"]]
    ax.legend(handles, labels, frameon=False, fontsize=8, title=r"$s_0$ [GeV$^2$]", loc="upper right")
    style_axis(ax)

    # Bs: vary s0 at selected M2.
    ax = axes[1, 1]
    for target, linestyle in ((BSLOW_WINDOW, "-"), (BS5830_WINDOW, "--")):
        for color, M2 in zip(color_cycle, target["M2_lines"]):
            subset = [
                r for r in rows
                if r["sector"] == "Bs" and r["state"] == target["state"] and r["scan"] == "s0"
                and abs(float(r["fixed_M2"]) - M2) < 1.0e-10
            ]
            ax.plot([float(r["s0"]) for r in subset], [float(r["g_quoted_abs"]) for r in subset], color=color, linewidth=1.7, linestyle=linestyle)
    ax.set_title(r"$B_{s1}\to B_s\gamma$: threshold stability")
    ax.set_xlabel(r"$s_0$ [GeV$^2$]")
    ax.set_ylabel(r"$|g_{\rm quoted}|$ [GeV$^{-1}$]")
    handles = [plt.Line2D([0], [0], color=c, lw=1.8) for c in color_cycle]
    labels = [rf"$M^2={M2:.0f}$" for M2 in BS5830_WINDOW["M2_lines"]]
    ax.legend(handles, labels, frameon=False, fontsize=8, title=r"$M^2$ [GeV$^2$]", loc="upper right")
    style_axis(ax)

    fig.tight_layout()
    fig.savefig(OUT / "publication_selected_window_stability.png", dpi=300)
    fig.savefig(OUT / "publication_selected_window_stability.pdf")


def central_norms(rows: list[dict[str, float | str]]) -> dict[str, float]:
    ds_center = min(
        [r for r in rows if r["sector"] == "Ds"],
        key=lambda r: abs(float(r["M2"]) - 3.75) + abs(float(r["s0"]) - 8.27),
    )
    bs_low_center = min(
        [r for r in rows if r["sector"] == "Bs" and r["state"] == BSLOW_WINDOW["state"]],
        key=lambda r: abs(float(r["M2"]) - 12.0) + abs(float(r["s0"]) - 38.46),
    )
    bs_high_center = min(
        [r for r in rows if r["sector"] == "Bs" and r["state"] == BS5830_WINDOW["state"]],
        key=lambda r: abs(float(r["M2"]) - 12.0) + abs(float(r["s0"]) - 39.39),
    )
    return {
        "Ds2460": abs(float(ds_center["g_2460"])),
        "Ds2536": abs(float(ds_center["g_2536"])),
        "Bslow": abs(float(bs_low_center["g_quoted"])),
        "Bs5830": abs(float(bs_high_center["g_quoted"])),
    }


def plot_publication_normalized(rows: list[dict[str, float | str]]) -> None:
    norms = central_norms(rows)
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2))
    color_cycle = [COLORS["low"], COLORS["mid"], COLORS["high"]]

    ax = axes[0, 0]
    for color, s0 in zip(color_cycle, DS_WINDOW["s0_lines"]):
        subset = [r for r in rows if r["sector"] == "Ds" and r["scan"] == "M2" and abs(float(r["fixed_s0"]) - s0) < 1.0e-10]
        ax.plot([float(r["M2"]) for r in subset], [float(r["g1_abs"]) / norms["Ds2460"] for r in subset], color=color, linewidth=1.8, label=rf"$s_0={s0:.2f}$")
        ax.plot([float(r["M2"]) for r in subset], [float(r["g2_abs"]) / norms["Ds2536"] for r in subset], color=color, linewidth=1.4, linestyle="--")
    ax.set_title(r"$D_{s1}\to D_s\gamma$: Borel stability")
    ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    ax.set_ylabel(r"$|g|/|g|_{\rm c}$")
    ax.legend(frameon=False, fontsize=8, title=r"$s_0$ [GeV$^2$]")
    ax.text(0.03, 0.08, r"solid: $2460$" "\n" r"dashed: $2536$", transform=ax.transAxes, fontsize=8)
    style_axis(ax)

    ax = axes[0, 1]
    for color, M2 in zip(color_cycle, DS_WINDOW["M2_lines"]):
        subset = [r for r in rows if r["sector"] == "Ds" and r["scan"] == "s0" and abs(float(r["fixed_M2"]) - M2) < 1.0e-10]
        ax.plot([float(r["s0"]) for r in subset], [float(r["g1_abs"]) / norms["Ds2460"] for r in subset], color=color, linewidth=1.8, label=rf"$M^2={M2:.2f}$")
        ax.plot([float(r["s0"]) for r in subset], [float(r["g2_abs"]) / norms["Ds2536"] for r in subset], color=color, linewidth=1.4, linestyle="--")
    ax.set_title(r"$D_{s1}\to D_s\gamma$: threshold stability")
    ax.set_xlabel(r"$s_0$ [GeV$^2$]")
    ax.set_ylabel(r"$|g|/|g|_{\rm c}$")
    ax.legend(frameon=False, fontsize=8, title=r"$M^2$ [GeV$^2$]")
    style_axis(ax)

    ax = axes[1, 0]
    for target, norm_key, linestyle in ((BSLOW_WINDOW, "Bslow", "-"), (BS5830_WINDOW, "Bs5830", "--")):
        for color, s0 in zip(color_cycle, target["s0_lines"]):
            subset = [
                r for r in rows
                if r["sector"] == "Bs" and r["state"] == target["state"] and r["scan"] == "M2"
                and abs(float(r["fixed_s0"]) - s0) < 1.0e-10
            ]
            ax.plot([float(r["M2"]) for r in subset], [float(r["g_quoted_abs"]) / norms[norm_key] for r in subset], color=color, linewidth=1.7, linestyle=linestyle)
    ax.set_title(r"$B_{s1}\to B_s\gamma$: Borel stability")
    ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    ax.set_ylabel(r"$|g_{\rm quoted}|/|g_{\rm quoted}|_{\rm c}$")
    ax.text(0.03, 0.78, r"solid: $B_{s1}^{low}$" "\n" r"dashed: $B_{s1}(5830)$", transform=ax.transAxes, fontsize=8)
    handles = [plt.Line2D([0], [0], color=c, lw=1.8) for c in color_cycle]
    labels = [rf"$s_0={s0:.2f}$" for s0 in BS5830_WINDOW["s0_lines"]]
    ax.legend(handles, labels, frameon=False, fontsize=8, title=r"$s_0$ [GeV$^2$]", loc="upper right")
    style_axis(ax)

    ax = axes[1, 1]
    for target, norm_key, linestyle in ((BSLOW_WINDOW, "Bslow", "-"), (BS5830_WINDOW, "Bs5830", "--")):
        for color, M2 in zip(color_cycle, target["M2_lines"]):
            subset = [
                r for r in rows
                if r["sector"] == "Bs" and r["state"] == target["state"] and r["scan"] == "s0"
                and abs(float(r["fixed_M2"]) - M2) < 1.0e-10
            ]
            ax.plot([float(r["s0"]) for r in subset], [float(r["g_quoted_abs"]) / norms[norm_key] for r in subset], color=color, linewidth=1.7, linestyle=linestyle)
    ax.set_title(r"$B_{s1}\to B_s\gamma$: threshold stability")
    ax.set_xlabel(r"$s_0$ [GeV$^2$]")
    ax.set_ylabel(r"$|g_{\rm quoted}|/|g_{\rm quoted}|_{\rm c}$")
    handles = [plt.Line2D([0], [0], color=c, lw=1.8) for c in color_cycle]
    labels = [rf"$M^2={M2:.0f}$" for M2 in BS5830_WINDOW["M2_lines"]]
    ax.legend(handles, labels, frameon=False, fontsize=8, title=r"$M^2$ [GeV$^2$]", loc="upper right")
    style_axis(ax)

    for ax in axes.flat:
        ax.axhline(1.0, color=COLORS["gray"], linewidth=0.9, linestyle=":", alpha=0.9)
        ax.set_ylim(0.55, 1.45)
    fig.tight_layout()
    fig.savefig(OUT / "publication_selected_window_stability_normalized.png", dpi=300)
    fig.savefig(OUT / "publication_selected_window_stability_normalized.pdf")


def summarize(rows: list[dict[str, float | str]]) -> str:
    lines = [
        "Publication selected-window stability plot",
        "==========================================",
        "Form factors are plotted as absolute values for visual clarity.",
        "Ds solid/dashed curves denote Ds1(2460)/Ds1(2536).",
        "Bs solid/dashed curves denote lower diagnostic/Bs1(5830).",
    ]
    for sector, key in (("Ds", "g1_abs"), ("Ds", "g2_abs"), ("Bs", "g_quoted_abs")):
        subset = [r for r in rows if r["sector"] == sector]
        if sector == "Ds" and key == "g2_abs":
            label = "Ds g2"
        elif sector == "Ds":
            label = "Ds g1"
        else:
            label = "Bs quoted g"
        vals = np.array([float(r[key]) for r in subset], dtype=float)
        lines.append(f"{label}: {vals.min():.5g} to {vals.max():.5g} GeV^-1")
    lines.extend(
        [
            f"Wrote {OUT / 'publication_selected_window_stability.csv'}",
            f"Wrote {OUT / 'publication_selected_window_stability.pdf'}",
            f"Wrote {OUT / 'publication_selected_window_stability.png'}",
            f"Wrote {OUT / 'publication_selected_window_stability_normalized.pdf'}",
            f"Wrote {OUT / 'publication_selected_window_stability_normalized.png'}",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows = build_rows(f1_axial, f1_basis)
    write_csv(OUT / "publication_selected_window_stability.csv", rows)
    plot_publication(rows)
    plot_publication_normalized(rows)
    summary = summarize(rows)
    (OUT / "publication_selected_window_stability_summary.txt").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
