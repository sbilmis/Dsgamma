"""Plot the preserved Ds1(2536) transition OPE for Referee A, major item ii.

Uses exactly the fixed inputs of the manuscript Borel-scan builder. This is
a decomposition of the existing baseline, not a new uncertainty extraction.
All new data, plots, and provenance are written alongside this script.
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
os.environ.setdefault("XDG_CACHE_HOME", str(HERE / ".cache"))
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np

from corrected_transition_central_python import DS_REFERENCE, ds_inputs
from final_stage2_uncertainty_scan import set_photon_shape_params
from rohrwild_transition_exact import (
    physical_couplings,
    precompute_convolutions,
    transition_invariants,
)


STEM = "ds2536_signed_contributions"
M2_GRID = np.linspace(3.0, 4.5, 61)
S0_GRID = np.linspace(9.0, 10.0, 21)
CENTRAL_M2 = 3.75
CENTRAL_S0 = 9.5
PIECES = ("pert", "tw2", "tw3", "tw4_total")
PLOT_PIECES = (*PIECES, "total")
LABELS = {
    "pert": "Perturbative",
    "tw2": "Twist 2",
    "tw3": "Twist 3",
    "tw4_total": "Twist 4 (all retained terms)",
    "total": "Total",
}
COLORS = {
    "pert": "#2463A6",
    "tw2": "#D55E00",
    "tw3": "#14864B",
    "tw4_total": "#8C4AB0",
    "total": "#111111",
}
STYLES = {"pert": "-", "tw2": "--", "tw3": "-.", "tw4_total": ":", "total": "-"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def manuscript_manifest() -> dict[str, str]:
    return {
        str(path.relative_to(ROOT)): sha256(path)
        for path in sorted((ROOT / "final_manuscript").rglob("*"))
        if path.is_file()
    }


def calculate(m2: float, s0: float, inputs: dict, convolutions: dict) -> tuple[dict, dict]:
    inv = transition_invariants(m2, s0, inputs, convolutions=convolutions)
    physical = physical_couplings(
        inv,
        theta_deg=DS_REFERENCE["theta_deg"],
        m_state_1=DS_REFERENCE["m1"],
        m_state_2=DS_REFERENCE["m2"],
        f_1=DS_REFERENCE["f1"],
        f_2=DS_REFERENCE["f2"],
        m_p=inputs["m_ds"],
        f_p=inputs["f_ds"],
        m_q=inputs["mc"],
        m_s=inputs["ms"],
    )
    theta = math.radians(DS_REFERENCE["theta_deg"])
    c, s = math.cos(theta), math.sin(theta)
    row = {"M2_GeV2": m2, "s0_GeV2": s0, "theta_deg": DS_REFERENCE["theta_deg"]}
    raw_pieces = ("pert", "tw2", "tw3", "tw4", "3p_g", "3p_gamma")
    for basis in ("A", "B"):
        for piece in raw_pieces:
            row[f"Pi_{basis}_{piece}"] = float(inv[f"T_{basis}_{piece}"])
        row[f"Pi_{basis}_tw4_total"] = sum(
            row[f"Pi_{basis}_{piece}"] for piece in ("tw4", "3p_g", "3p_gamma")
        )
        row[f"Pi_{basis}_total"] = float(inv[f"T_{basis}"])
    for piece in (*raw_pieces, "tw4_total"):
        projection = c * row[f"Pi_A_{piece}"] - s * row[f"Pi_B_{piece}"]
        row[f"physical_numerator_{piece}"] = projection
        row[f"g2_{piece}"] = float(physical["N_2"]) * projection
    row["physical_numerator_total"] = c * row["Pi_A_total"] - s * row["Pi_B_total"]
    row["g2_total"] = float(physical["g_2"])
    row["g2_from_A"] = float(physical["g_2_A"])
    row["g2_from_B"] = float(physical["g_2_B"])
    row["N2_GeV_minus4"] = float(physical["N_2"])
    row["f2_GeV"] = DS_REFERENCE["f2"]
    return row, {**inv, **physical}


def check_baseline(rows: list[dict], central_values: dict) -> dict:
    by_point = {(round(r["M2_GeV2"], 8), round(r["s0_GeV2"], 8)): r for r in rows}
    old_borel = [
        r for r in read_csv(ROOT / "figures/figure1_borel_form_factors/figure1_borel_form_factors.csv")
        if r["state"] == "Ds1_2536"
    ]
    differences = [
        abs(by_point[(round(float(r["M2"]), 8), round(float(r["s0"]), 8))]["g2_total"] - float(r["g_signed"]))
        for r in old_borel
    ]
    old_central = {
        r["key"].removeprefix("Ds_2536."): float(r["value"])
        for r in read_csv(ROOT / "outputs/corrected_transition_central_python.csv")
        if r["key"].startswith("Ds_2536.")
    }
    central_differences = [abs(float(central_values[k]) - v) for k, v in old_central.items()]
    sum_differences = {
        panel: max(abs(sum(r[f"{panel}_{piece}"] for piece in PIECES) - r[f"{panel}_total"]) for r in rows)
        for panel in ("Pi_A", "Pi_B", "g2")
    }
    assert len(differences) == 93, f"Unexpected Borel reference size: {len(differences)}"
    assert max(differences) < 1e-12
    assert max(central_differences) < 1e-12
    assert max(sum_differences.values()) < 1e-12
    assert all(math.isfinite(float(v)) for r in rows for v in r.values())
    return {
        "saved_borel_points_checked": len(differences),
        "saved_borel_max_abs_g2_difference": max(differences),
        "saved_central_values_checked": len(central_differences),
        "saved_central_max_abs_difference": max(central_differences),
        "max_abs_component_sum_errors": sum_differences,
        "max_abs_axial_twist3": max(abs(r["Pi_A_tw3"]) for r in rows),
    }


def draw(rows: list[dict]) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "mathtext.fontset": "stix",
        "font.size": 12, "axes.labelsize": 15, "axes.titlesize": 14,
        "xtick.labelsize": 11, "ytick.labelsize": 11,
        "axes.linewidth": 1.0, "savefig.facecolor": "white",
        "svg.fonttype": "path",
    })
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 6.2))
    fig.subplots_adjust(left=0.065, right=0.988, bottom=0.25, top=0.745, wspace=0.30)
    fig.suptitle(r"$D_{s1}(2536)\to D_s\gamma$  |  $\theta=26.6^\circ$", y=0.985, fontsize=18)
    fig.text(0.5, 0.905, "Signed contributions across the manuscript Borel window", ha="center", fontsize=12.5)
    handles = [
        Line2D([], [], color=COLORS[piece], linestyle=STYLES[piece],
               linewidth=2.8 if piece == "total" else 2.2, label=LABELS[piece])
        for piece in PLOT_PIECES
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.866),
               ncol=5, frameon=False, handlelength=2.8, columnspacing=1.7, fontsize=11.5)
    panels = (
        ("Pi_A", r"(a) Axial basis", r"$\widehat\Pi_A\;[\mathrm{GeV}^3]$"),
        ("Pi_B", r"(b) Tensor basis", r"$\widehat\Pi_B\;[\mathrm{GeV}^3]$"),
        ("g2", r"(c) Physical higher-state amplitude", r"$g_2\;[\mathrm{GeV}^{-1}]$"),
    )
    by_point = {(round(r["M2_GeV2"], 8), round(r["s0_GeV2"], 8)): r for r in rows}
    bounds = []
    for ax, (prefix, title, ylabel) in zip(axes, panels):
        panel_min, panel_max = 0.0, 0.0
        ax.axhline(0, color="#8A8A8A", linewidth=0.8, zorder=0)
        for piece in PLOT_PIECES:
            values = np.array([
                [by_point[(round(float(m2), 8), round(float(s0), 8))][f"{prefix}_{piece}"]
                 for m2 in M2_GRID]
                for s0 in S0_GRID
            ])
            lower, upper = values.min(axis=0), values.max(axis=0)
            central = values[len(S0_GRID) // 2]
            ax.fill_between(M2_GRID, lower, upper, color=COLORS[piece], alpha=0.13,
                            linewidth=0, zorder=1 if piece != "total" else 3)
            ax.plot(M2_GRID, central, color=COLORS[piece], linestyle=STYLES[piece],
                    linewidth=2.1 if piece != "total" else 2.7,
                    zorder=2 if piece != "total" else 4)
            panel_min, panel_max = min(panel_min, float(lower.min())), max(panel_max, float(upper.max()))
        bounds.append((panel_min, panel_max))
        ax.set_title(title, pad=13)
        ax.set_xlabel(r"$M^2\;[\mathrm{GeV}^2]$", labelpad=9)
        ax.set_ylabel(ylabel, labelpad=8)
        ax.set_xlim(3.0, 4.5)
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.xaxis.set_minor_locator(MultipleLocator(0.125))
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.tick_params(which="major", length=5)
        ax.tick_params(which="minor", length=3)
        ax.grid(axis="y", color="#DADADA", linewidth=0.5, alpha=0.55, zorder=0)
    basis_min = min(bounds[0][0], bounds[1][0])
    basis_max = max(bounds[0][1], bounds[1][1])
    basis_padding = 0.08 * (basis_max - basis_min)
    for ax in axes[:2]:
        ax.set_ylim(basis_min - basis_padding, basis_max + basis_padding)
    physical_padding = 0.08 * (bounds[2][1] - bounds[2][0])
    axes[2].set_ylim(bounds[2][0] - physical_padding, bounds[2][1] + physical_padding)
    midpoint = by_point[(CENTRAL_M2, CENTRAL_S0)]
    axes[2].text(0.04, 0.04,
                 rf"$g_2(3.75,\,9.5)={midpoint['g2_total']:+.4f}\;\mathrm{{GeV}}^{{-1}}$",
                 transform=axes[2].transAxes, fontsize=11,
                 bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=3))
    fig.text(0.5, 0.12,
             r"Lines: $s_0=9.5\,\mathrm{GeV}^2$; bands: $s_0\in[9,10]\,\mathrm{GeV}^2$ (threshold variation only).",
             ha="center", fontsize=11.5)
    fig.text(0.5, 0.065,
             r"Fixed reference inputs, $f_2=166.67\,\mathrm{MeV}$; $M_1^2=M_2^2=2M^2$. Axial twist 3 vanishes at $u_0=1/2$.",
             ha="center", fontsize=11)
    fig.savefig(HERE / f"{STEM}.png", dpi=200)
    fig.savefig(HERE / f"{STEM}.svg")
    plt.close(fig)


def main() -> None:
    protected_before = manuscript_manifest()
    baseline_paths = [
        ROOT / "figures/figure1_borel_form_factors/figure1_borel_form_factors.csv",
        ROOT / "figures/figure1_borel_form_factors/build_figure1_data.py",
        ROOT / "outputs/corrected_transition_central_python.csv",
    ]
    baseline_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in baseline_paths}
    set_photon_shape_params(-2.1, 3.8)
    convolutions = precompute_convolutions()
    inputs = ds_inputs()
    rows = []
    central_values = None
    for s0 in S0_GRID:
        for m2 in M2_GRID:
            row, values = calculate(float(m2), float(s0), inputs, convolutions)
            rows.append(row)
            if abs(m2 - CENTRAL_M2) < 1e-12 and abs(s0 - CENTRAL_S0) < 1e-12:
                central_values = values
    assert central_values is not None
    checks = check_baseline(rows, central_values)
    write_csv(HERE / f"{STEM}_data.csv", rows)
    midpoint = next(r for r in rows if r["M2_GeV2"] == CENTRAL_M2 and r["s0_GeV2"] == CENTRAL_S0)
    summary = [
        {"contribution": piece, "Pi_A_GeV3": midpoint[f"Pi_A_{piece}"],
         "Pi_B_GeV3": midpoint[f"Pi_B_{piece}"], "g2_GeV_minus1": midpoint[f"g2_{piece}"]}
        for piece in ("pert", "tw2", "tw3", "tw4", "3p_g", "3p_gamma", "tw4_total", "total")
    ]
    write_csv(HERE / f"{STEM}_midpoint.csv", summary)
    draw(rows)
    source_files = {Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, "__file__", None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == ".py" and (path.is_relative_to(ROOT / "scripts") or path.is_relative_to(ROOT.parents[1] / "shared")):
                source_files.add(path)
    assert manuscript_manifest() == protected_before, "Protected manuscript changed during the run"
    assert all(sha256(ROOT / name) == value for name, value in baseline_hashes.items())
    checks["protected_manuscript_unchanged"] = True
    checks["reference_source_and_data_unchanged"] = True
    provenance = {
        "purpose": "Review-only decomposition of the preserved Ds1(2536) baseline",
        "reference_manuscript": "/Users/sbilmis/Downloads/r1-2.pdf",
        "channel": "Ds1(2536) -> Ds gamma",
        "theta_deg": DS_REFERENCE["theta_deg"],
        "fixed_inputs": inputs,
        "pole_and_residue_reference": DS_REFERENCE,
        "M2_grid_GeV2": M2_GRID.tolist(),
        "s0_grid_GeV2": S0_GRID.tolist(),
        "central_s0_GeV2": CENTRAL_S0,
        "band_meaning": "Pointwise min/max over the 21 sampled s0 values; not a probability interval",
        "twist4_grouping": "two-particle tw4 + gluonic three-particle + electromagnetic three-particle",
        "baseline_status": "Current production formulas, including existing continuum prescription; no pending correction adopted",
        "fixed_residue_note": "f2 is held at the same reference value as the manuscript Borel scan; transition s0 is not the two-point residue threshold",
        "checks": checks,
        "baseline_sha256": baseline_hashes,
        "source_sha256": {str(path): sha256(path) for path in sorted(source_files)},
        "midpoint": midpoint,
        "g2_range_over_sampled_window": [min(r["g2_total"] for r in rows), max(r["g2_total"] for r in rows)],
    }
    (HERE / f"{STEM}_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"figure": str(HERE / f"{STEM}.png"), "checks": checks,
                      "midpoint_contributions": summary,
                      "g2_range": provenance["g2_range_over_sampled_window"]}, indent=2))


if __name__ == "__main__":
    main()
