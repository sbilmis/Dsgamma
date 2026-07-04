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
    "s0_lines": (8.0, 8.3, 8.6),
    "M2_lines": (3.0, 3.5, 4.0),
}

BS5830_WINDOW = {
    "state": "B_{s1}(5830)",
    "mass": 5.82870,
    "quoted_combo": "high",
    "M2": (10.0, 14.0),
    "s0": (37.82, 40.96),
    "s0_lines": (38.0, 39.0, 40.0),
    "M2_lines": (10.0, 12.0, 14.0),
}

BSLOW_WINDOW = {
    "state": "B_{s1}^{low}(5750)",
    "mass": 5.750,
    "quoted_combo": "low",
    "M2": (10.0, 14.0),
    "s0": (36.60, 40.32),
    "s0_lines": (37.0, 38.0, 39.0),
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


def fmt_scan_value(value: float) -> str:
    if abs(value - round(value)) < 1.0e-10:
        return f"{value:.0f}"
    return f"{value:.1f}"


def loosen_ylim(ax, values: list[float], fraction: float = 1.80) -> None:
    """Use a deliberately loose y range so stability is not visually overstated."""
    ymin = min(values)
    ymax = max(values)
    center = 0.5 * (ymin + ymax)
    half = 0.5 * (ymax - ymin)
    if half <= 0.0:
        half = max(abs(center) * 0.15, 1.0e-3)
    half *= 1.0 + fraction
    ax.set_ylim(center - half, center + half)


def add_inside_title(ax, title: str) -> None:
    ax.text(
        0.04,
        0.92,
        title,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
    )


def add_clean_legend(ax, *, title: str) -> None:
    ax.legend(
        frameon=True,
        facecolor="white",
        framealpha=0.92,
        edgecolor="white",
        fontsize=8,
        title=title,
        loc="upper right",
    )


def set_channel_ylim(ax, spec: dict[str, object], values: list[float]) -> None:
    if "ylim" in spec:
        ymin, ymax = spec["ylim"]
        ax.set_ylim(float(ymin), float(ymax))
    else:
        loosen_ylim(ax, values)


def plot_publication(rows: list[dict[str, float | str]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2))
    color_cycle = [COLORS["low"], COLORS["mid"], COLORS["high"]]

    # Ds: vary M2 at selected s0.
    ax = axes[0, 0]
    for color, s0 in zip(color_cycle, DS_WINDOW["s0_lines"]):
        subset = [r for r in rows if r["sector"] == "Ds" and r["scan"] == "M2" and abs(float(r["fixed_s0"]) - s0) < 1.0e-10]
        ax.plot([float(r["M2"]) for r in subset], [float(r["g1_abs"]) for r in subset], color=color, linewidth=1.8, label=rf"$s_0={fmt_scan_value(s0)}$")
        ax.plot([float(r["M2"]) for r in subset], [float(r["g2_abs"]) for r in subset], color=color, linewidth=1.4, linestyle="--")
    ax.set_title(r"$D_{s1}\to D_s\gamma$: Borel stability")
    ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    ax.set_ylabel(r"$|g|$ [GeV$^{-1}$]")
    add_clean_legend(ax, title=r"$s_0$ [GeV$^2$]")
    ax.text(0.03, 0.08, r"solid: $2460$" "\n" r"dashed: $2536$", transform=ax.transAxes, fontsize=8)
    style_axis(ax)

    # Ds: vary s0 at selected M2.
    ax = axes[0, 1]
    for color, M2 in zip(color_cycle, DS_WINDOW["M2_lines"]):
        subset = [r for r in rows if r["sector"] == "Ds" and r["scan"] == "s0" and abs(float(r["fixed_M2"]) - M2) < 1.0e-10]
        ax.plot([float(r["s0"]) for r in subset], [float(r["g1_abs"]) for r in subset], color=color, linewidth=1.8, label=rf"$M^2={fmt_scan_value(M2)}$")
        ax.plot([float(r["s0"]) for r in subset], [float(r["g2_abs"]) for r in subset], color=color, linewidth=1.4, linestyle="--")
    ax.set_title(r"$D_{s1}\to D_s\gamma$: threshold stability")
    ax.set_xlabel(r"$s_0$ [GeV$^2$]")
    ax.set_ylabel(r"$|g|$ [GeV$^{-1}$]")
    add_clean_legend(ax, title=r"$M^2$ [GeV$^2$]")
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
    labels = [rf"$s_0={fmt_scan_value(s0)}$" for s0 in BS5830_WINDOW["s0_lines"]]
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
        key=lambda r: abs(float(r["M2"]) - 3.5) + abs(float(r["s0"]) - 8.3),
    )
    bs_low_center = min(
        [r for r in rows if r["sector"] == "Bs" and r["state"] == BSLOW_WINDOW["state"]],
        key=lambda r: abs(float(r["M2"]) - 12.0) + abs(float(r["s0"]) - 38.0),
    )
    bs_high_center = min(
        [r for r in rows if r["sector"] == "Bs" and r["state"] == BS5830_WINDOW["state"]],
        key=lambda r: abs(float(r["M2"]) - 12.0) + abs(float(r["s0"]) - 39.0),
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
        ax.plot([float(r["M2"]) for r in subset], [float(r["g1_abs"]) / norms["Ds2460"] for r in subset], color=color, linewidth=1.8, label=rf"$s_0={fmt_scan_value(s0)}$")
        ax.plot([float(r["M2"]) for r in subset], [float(r["g2_abs"]) / norms["Ds2536"] for r in subset], color=color, linewidth=1.4, linestyle="--")
    ax.set_title(r"$D_{s1}\to D_s\gamma$: Borel stability")
    ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    ax.set_ylabel(r"$|g|/|g|_{\rm c}$")
    ax.legend(frameon=False, fontsize=8, title=r"$s_0$ [GeV$^2$]", loc="upper right")
    ax.text(0.03, 0.08, r"solid: $2460$" "\n" r"dashed: $2536$", transform=ax.transAxes, fontsize=8)
    style_axis(ax)

    ax = axes[0, 1]
    for color, M2 in zip(color_cycle, DS_WINDOW["M2_lines"]):
        subset = [r for r in rows if r["sector"] == "Ds" and r["scan"] == "s0" and abs(float(r["fixed_M2"]) - M2) < 1.0e-10]
        ax.plot([float(r["s0"]) for r in subset], [float(r["g1_abs"]) / norms["Ds2460"] for r in subset], color=color, linewidth=1.8, label=rf"$M^2={fmt_scan_value(M2)}$")
        ax.plot([float(r["s0"]) for r in subset], [float(r["g2_abs"]) / norms["Ds2536"] for r in subset], color=color, linewidth=1.4, linestyle="--")
    ax.set_title(r"$D_{s1}\to D_s\gamma$: threshold stability")
    ax.set_xlabel(r"$s_0$ [GeV$^2$]")
    ax.set_ylabel(r"$|g|/|g|_{\rm c}$")
    ax.legend(frameon=False, fontsize=8, title=r"$M^2$ [GeV$^2$]", loc="upper right")
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
    labels = [rf"$s_0={fmt_scan_value(s0)}$" for s0 in BS5830_WINDOW["s0_lines"]]
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


CHANNEL_SPECS = [
    {
        "channel_id": "ds1_2460",
        "title": r"$D_{s1}(2460)\to D_s\gamma$",
        "sector": "Ds",
        "value_key": "g1_abs",
        "window": DS_WINDOW,
        "state": "",
        "ylabel": r"$|g_1|$ [GeV$^{-1}$]",
        "ylim": (0.25, 0.60),
    },
    {
        "channel_id": "ds1_2536",
        "title": r"$D_{s1}(2536)\to D_s\gamma$",
        "sector": "Ds",
        "value_key": "g2_abs",
        "window": DS_WINDOW,
        "state": "",
        "ylabel": r"$|g_2|$ [GeV$^{-1}$]",
        "ylim": (0.00, 0.06),
    },
    {
        "channel_id": "bs1_low_5750",
        "title": r"$B_{s1}^{low}(5750)\to B_s\gamma$",
        "sector": "Bs",
        "value_key": "g_quoted_abs",
        "window": BSLOW_WINDOW,
        "state": BSLOW_WINDOW["state"],
        "ylabel": r"$|g_1|$ [GeV$^{-1}$]",
        "ylim": (0.25, 0.50),
    },
    {
        "channel_id": "bs1_5830",
        "title": r"$B_{s1}(5830)\to B_s\gamma$",
        "sector": "Bs",
        "value_key": "g_quoted_abs",
        "window": BS5830_WINDOW,
        "state": BS5830_WINDOW["state"],
        "ylabel": r"$|g_2|$ [GeV$^{-1}$]",
        "ylim": (0.03, 0.05),
    },
]


def channel_subset(rows: list[dict[str, float | str]], spec: dict[str, object], scan: str, fixed_key: str, fixed_value: float) -> list[dict[str, float | str]]:
    subset = [
        r for r in rows
        if r["sector"] == spec["sector"]
        and r["scan"] == scan
        and abs(float(r[fixed_key]) - fixed_value) < 1.0e-10
    ]
    if spec["sector"] == "Bs":
        subset = [r for r in subset if r["state"] == spec["state"]]
    return subset


def plot_channel_absolute(rows: list[dict[str, float | str]], spec: dict[str, object]) -> None:
    window = spec["window"]
    value_key = str(spec["value_key"])
    color_cycle = [COLORS["low"], COLORS["mid"], COLORS["high"]]
    stem = f"publication_stability_{spec['channel_id']}"

    fig, ax = plt.subplots(figsize=(4.8, 3.7))
    y_values: list[float] = []
    for color, s0 in zip(color_cycle, window["s0_lines"]):
        subset = channel_subset(rows, spec, "M2", "fixed_s0", float(s0))
        y = [float(r[value_key]) for r in subset]
        y_values.extend(y)
        ax.plot(
            [float(r["M2"]) for r in subset],
            y,
            color=color,
            linewidth=1.8,
            label=rf"$s_0={fmt_scan_value(float(s0))}$",
        )
    add_inside_title(ax, str(spec["title"]))
    ax.set_xlabel(r"$M^2$ [GeV$^2$]")
    ax.set_ylabel(str(spec["ylabel"]))
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    style_axis(ax)
    set_channel_ylim(ax, spec, y_values)
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}_M2.png", dpi=300)
    fig.savefig(OUT / f"{stem}_M2.pdf")

    fig, ax = plt.subplots(figsize=(4.8, 3.7))
    y_values = []
    for color, M2 in zip(color_cycle, window["M2_lines"]):
        subset = channel_subset(rows, spec, "s0", "fixed_M2", float(M2))
        y = [float(r[value_key]) for r in subset]
        y_values.extend(y)
        ax.plot(
            [float(r["s0"]) for r in subset],
            y,
            color=color,
            linewidth=1.8,
            label=rf"$M^2={fmt_scan_value(float(M2))}$",
        )
    add_inside_title(ax, str(spec["title"]))
    ax.set_xlabel(r"$s_0$ [GeV$^2$]")
    ax.set_ylabel(str(spec["ylabel"]))
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    style_axis(ax)
    set_channel_ylim(ax, spec, y_values)
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}_s0.png", dpi=300)
    fig.savefig(OUT / f"{stem}_s0.pdf")


def plot_channel_figures(rows: list[dict[str, float | str]]) -> None:
    for spec in CHANNEL_SPECS:
        plot_channel_absolute(rows, spec)


def summarize(rows: list[dict[str, float | str]]) -> str:
    lines = [
        "Publication selected-window stability plot",
        "==========================================",
        "The legacy combined plots are kept for comparison.",
        "The channel-separated plots are preferred because they avoid solid/dashed overloading.",
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
            "Wrote channel-separated plots:",
        ]
    )
    for spec in CHANNEL_SPECS:
        lines.append(f"  {OUT / ('publication_stability_' + str(spec['channel_id']) + '_M2.pdf')}")
        lines.append(f"  {OUT / ('publication_stability_' + str(spec['channel_id']) + '_s0.pdf')}")
    return "\n".join(lines) + "\n"


def main() -> None:
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows = build_rows(f1_axial, f1_basis)
    write_csv(OUT / "publication_selected_window_stability.csv", rows)
    plot_publication(rows)
    plot_publication_normalized(rows)
    plot_channel_figures(rows)
    summary = summarize(rows)
    (OUT / "publication_selected_window_stability_summary.txt").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
