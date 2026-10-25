"""Restyle the saved OPE scan and show Pi_2 instead of the normalized g_2.

This only reads the previously checked numerical scan. No OPE, inputs,
residues, or production files are changed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".cache" / "matplotlib"))
sys.dont_write_bytecode = True

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, MultipleLocator
import numpy as np


DATA = HERE / "ds2536_signed_contributions_data.csv"
STEM = "ds2536_signed_invariants"
TERMS = ("pert", "tw2", "tw3", "tw4_total", "total")
LABELS = ("Perturbative", "Twist 2", "Twist 3", "Twist 4", "Total")
COLORS = ("#0047ab", "#cc3300", "#14864b", "#8c4ab0", "#111111")
STYLES = ("-", "--", "-.", ":", "-")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    source_hash = sha256(DATA)
    protected_before = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (ROOT / "final_manuscript").rglob("*") if path.is_file()
    }
    with DATA.open(newline="") as handle:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]
    m2_grid = sorted({r["M2_GeV2"] for r in rows})
    s0_grid = sorted({r["s0_GeV2"] for r in rows})
    by_point = {(r["M2_GeV2"], r["s0_GeV2"]): r for r in rows}
    midpoint = by_point[(3.75, 9.5)]
    c, s = math.cos(math.radians(26.6)), math.sin(math.radians(26.6))
    projection_error = max(
        abs(r[f"physical_numerator_{term}"] - (c*r[f"Pi_A_{term}"] - s*r[f"Pi_B_{term}"]))
        for r in rows for term in TERMS
    )
    component_sum_error = max(
        abs(sum(r[f"physical_numerator_{term}"] for term in TERMS[:-1]) - r["physical_numerator_total"])
        for r in rows
    )
    normalization_error = max(
        abs(r["N2_GeV_minus4"] * r["physical_numerator_total"] - r["g2_total"])
        for r in rows
    )
    assert max(projection_error, component_sum_error, normalization_error) < 1e-12

    # Match the sans-serif math, colours, boxed axes, inward ticks, and
    # ungridded presentation of the supplied manuscript and existing builders.
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "mathtext.fontset": "dejavusans", "font.size": 9,
        "axes.labelsize": 10.5, "axes.titlesize": 11,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.linewidth": 0.8, "svg.fonttype": "path",
        "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.45))
    fig.subplots_adjust(left=0.075, right=0.988, bottom=0.19, top=0.76, wspace=0.42)
    handles = [Line2D([], [], color=color, ls=style, lw=2.0 if term == "total" else 1.8, label=label)
               for term, label, color, style in zip(TERMS, LABELS, COLORS, STYLES)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.53, 0.985),
               ncol=5, frameon=False, fontsize=9.5, handlelength=2.8, columnspacing=1.8)
    specs = (
        ("Pi_A", r"(a) Axial current", r"$\widehat\Pi_A\;[\mathrm{GeV}^3]$"),
        ("Pi_B", r"(b) Tensor current", r"$\widehat\Pi_B\;[\mathrm{GeV}^3]$"),
        ("physical_numerator", r"(c) Rotated current", r"$\widehat\Pi_2\;[\mathrm{GeV}^3]$"),
    )
    bounds = []
    for ax, (prefix, title, ylabel) in zip(axes, specs):
        ymin, ymax = 0.0, 0.0
        ax.axhline(0, color="0.45", lw=0.7, ls=":", zorder=0)
        for term, color, style in zip(TERMS, COLORS, STYLES):
            arr = np.array([[by_point[(m2, s0)][f"{prefix}_{term}"] for m2 in m2_grid] for s0 in s0_grid])
            lower, upper = arr.min(axis=0), arr.max(axis=0)
            central = arr[s0_grid.index(9.5)]
            ax.fill_between(m2_grid, lower, upper, color=color, alpha=0.11, linewidth=0, zorder=1)
            ax.plot(m2_grid, central, color=color, ls=style, lw=2.0 if term == "total" else 1.8,
                    zorder=4 if term == "total" else 2)
            ymin, ymax = min(ymin, float(lower.min())), max(ymax, float(upper.max()))
        bounds.append((ymin, ymax))
        ax.set_title(title, pad=8)
        ax.set_xlabel(r"$M^2\;[\mathrm{GeV}^2]$", labelpad=5)
        ax.set_ylabel(ylabel, labelpad=4)
        ax.set_xlim(3.0, 4.5)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.tick_params(which="major", length=4)
        ax.tick_params(which="minor", length=2)
        ax.grid(False)
    ymin, ymax = min(bounds[0][0], bounds[1][0]), max(bounds[0][1], bounds[1][1])
    pad = 0.09*(ymax-ymin)
    for ax in axes[:2]:
        ax.set_ylim(ymin-pad, ymax+pad)
    ymin, ymax = bounds[2]
    pad = 0.10*(ymax-ymin)
    axes[2].set_ylim(ymin-pad, ymax+pad)
    axes[2].text(0.96, 0.055, r"$\theta=26.6^\circ$", transform=axes[2].transAxes,
                 ha="right", va="bottom", fontsize=8.5, color="0.22")
    fig.savefig(HERE / f"{STEM}.png", dpi=240)
    fig.savefig(HERE / f"{STEM}.svg")
    plt.close(fig)

    table = [
        {"contribution": term, "Pi_A_GeV3": midpoint[f"Pi_A_{term}"],
         "Pi_B_GeV3": midpoint[f"Pi_B_{term}"], "Pi_2_GeV3": midpoint[f"physical_numerator_{term}"]}
        for term in ("pert", "tw2", "tw3", "tw4", "3p_g", "3p_gamma", "tw4_total", "total")
    ]
    with (HERE / f"{STEM}_midpoint.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    assert sha256(DATA) == source_hash
    assert all(sha256(ROOT / name) == digest for name, digest in protected_before.items())
    metadata = {
        "figure": str(HERE / f"{STEM}.png"),
        "source_data": str(DATA), "source_data_sha256": source_hash,
        "plot_script_sha256": sha256(Path(__file__)),
        "third_panel": "Pi_2 = cos(theta) Pi_A - sin(theta) Pi_B; not the normalized g2",
        "all_panel_units": "GeV^3", "central_s0_GeV2": 9.5,
        "shading": "Pointwise envelope over the existing 21-point s0 scan from 9 to 10 GeV^2",
        "style_reference": "/Users/sbilmis/Downloads/r1-2.pdf, Figure 2",
        "max_abs_projection_error_GeV3": projection_error,
        "max_abs_component_sum_error_GeV3": component_sum_error,
        "max_abs_g2_normalization_error_GeV_minus1": normalization_error,
        "source_data_and_manuscript_unchanged": True,
        "midpoint_Pi2_GeV3": midpoint["physical_numerator_total"],
    }
    (HERE / f"{STEM}_provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
