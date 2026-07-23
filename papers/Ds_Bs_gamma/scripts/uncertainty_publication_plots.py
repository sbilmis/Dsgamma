"""Publication plots for the Stage-2 uncertainty analysis."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))

import matplotlib.pyplot as plt


COLORS = {
    "blue": "#1f77b4",
    "red": "#d62728",
    "green": "#2ca02c",
    "orange": "#ff7f0e",
    "gray": "#7f7f7f",
}


LABELS = {
    "borel_threshold": r"$M^2,\sqrt{s_0}$",
    "theta": r"$\theta$",
    "quark_and_hadron_masses": "masses",
    "decay_constants": "decay constants",
    "condensates_and_chi": r"condensates, $\chi$",
    "photon_DA_shape": "photon DA shape",
    "all_inputs_fixed_theta": r"all inputs, fixed $\theta$",
    "all_inputs_and_theta": r"all inputs + $\theta$",
}


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def gaussian_pdf(x, mean, sigma):
    return np.exp(-0.5 * ((x - mean) / sigma) ** 2) / (sigma * math.sqrt(2.0 * math.pi))


def mc_histogram_plots():
    rows = read_csv(OUT / "final_stage2_uncertainty_scan.csv")
    summary_rows = []
    for ensemble, suffix, title_extra in (
        ("theta_fixed", "theta_fixed", r"$\theta=35.3^\circ$"),
        ("theta_scan_25_45", "theta_scan", r"$25^\circ\leq\theta\leq45^\circ$"),
    ):
        subset = [r for r in rows if r["ensemble"] == ensemble]
        fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9))
        for ax, key, label, color in (
            (axes[0], "Gamma2460_keV", r"$\Gamma_{2460}$", COLORS["blue"]),
            (axes[1], "Gamma2536_keV", r"$\Gamma_{2536}$", COLORS["red"]),
        ):
            values = np.array([float(r[key]) for r in subset], dtype=float)
            mean = float(np.mean(values))
            sigma = float(np.std(values, ddof=1))
            median = float(np.percentile(values, 50.0))
            p16 = float(np.percentile(values, 16.0))
            p84 = float(np.percentile(values, 84.0))
            spread = float(np.ptp(values))
            x_min = max(0.0, float(values.min()) - 0.08 * spread)
            x_max = float(values.max()) + 0.08 * spread
            x = np.linspace(x_min, x_max, 400)

            ax.hist(values, bins=28, density=True, alpha=0.45, color=color, edgecolor="white")
            ax.plot(x, gaussian_pdf(x, mean, sigma), color="black", linewidth=1.8)
            ax.axvline(mean, color="black", linestyle="-", linewidth=1.0)
            ax.axvline(p16, color=color, linestyle="--", linewidth=1.0)
            ax.axvline(p84, color=color, linestyle="--", linewidth=1.0)
            ax.set_xlabel(r"$\Gamma$ [keV]")
            ax.set_ylabel("density")
            ax.set_title(label)
            ax.grid(True, alpha=0.25, linewidth=0.7)
            ax.text(
                0.97,
                0.93,
                "\n".join(
                    [
                        rf"$\mu={mean:.2f}$ keV",
                        rf"$\sigma={sigma:.2f}$ keV",
                        rf"med. $={median:.2f}$ keV",
                    ]
                ),
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
            )
            summary_rows.append(
                {
                    "ensemble": ensemble,
                    "observable": key,
                    "mean_keV": mean,
                    "sigma_keV": sigma,
                    "median_keV": median,
                    "p16_keV": p16,
                    "p84_keV": p84,
                }
            )
        fig.suptitle(f"Monte Carlo width distribution, {title_extra}")
        fig.tight_layout()
        fig.savefig(OUT / f"mc_width_histograms_{suffix}.png", dpi=220)
        fig.savefig(OUT / f"mc_width_histograms_{suffix}.pdf")
    return summary_rows


def budget_rank_plot():
    rows = read_csv(OUT / "uncertainty_budget_one_at_time_summary.csv")
    component_rows = [
        r for r in rows if r["group"] not in {"all_inputs_fixed_theta", "all_inputs_and_theta"}
    ]
    sorted_2460 = sorted(component_rows, key=lambda r: float(r["Gamma2460_half_width"]), reverse=True)
    sorted_2536 = sorted(component_rows, key=lambda r: float(r["Gamma2536_half_width"]), reverse=True)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    for ax, sorted_rows, key, title, color in (
        (axes[0], sorted_2460, "Gamma2460_half_width", r"$D_{s1}(2460)$", COLORS["blue"]),
        (axes[1], sorted_2536, "Gamma2536_half_width", r"$D_{s1}(2536)$", COLORS["red"]),
    ):
        labels = [LABELS[r["group"]] for r in sorted_rows]
        values = [float(r[key]) for r in sorted_rows]
        y = np.arange(len(labels))
        ax.barh(y, values, color=color, alpha=0.82)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlabel(r"16--84 half-width [keV]")
        ax.set_title(title)
        ax.grid(True, axis="x", alpha=0.25, linewidth=0.7)
    fig.suptitle("Dominant uncertainty sources from one-at-a-time scan")
    fig.tight_layout()
    fig.savefig(OUT / "uncertainty_budget_ranked.png", dpi=220)
    fig.savefig(OUT / "uncertainty_budget_ranked.pdf")

    rank_rows = []
    for decay, sorted_rows, key in (
        ("Gamma2460", sorted_2460, "Gamma2460_half_width"),
        ("Gamma2536", sorted_2536, "Gamma2536_half_width"),
    ):
        for rank, row in enumerate(sorted_rows, start=1):
            rank_rows.append(
                {
                    "decay": decay,
                    "rank": rank,
                    "group": row["group"],
                    "half_width_keV": float(row[key]),
                }
            )
    return rank_rows


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    hist_summary = mc_histogram_plots()
    rank_rows = budget_rank_plot()
    write_csv(OUT / "mc_width_gaussian_fit_summary.csv", hist_summary)
    write_csv(OUT / "uncertainty_budget_ranked.csv", rank_rows)

    lines = [
        "Uncertainty publication plots",
        "=============================",
        "Gaussian overlays use the Monte Carlo sample mean and sample standard deviation.",
        "",
        "Gaussian-fit summaries:",
    ]
    for row in hist_summary:
        lines.append(
            f"  {row['ensemble']} {row['observable']}: "
            f"mean {row['mean_keV']:.3g} keV, sigma {row['sigma_keV']:.3g} keV; "
            f"median {row['median_keV']:.3g} keV"
        )
    lines.extend(["", "Dominant one-at-a-time sources:"])
    for decay in ("Gamma2460", "Gamma2536"):
        top = [r for r in rank_rows if r["decay"] == decay][:3]
        lines.append(
            f"  {decay}: "
            + ", ".join(f"{r['group']} ({r['half_width_keV']:.3g} keV)" for r in top)
        )
    lines.extend(
        [
            "",
            "Figures:",
            f"  {OUT / 'mc_width_histograms_theta_fixed.pdf'}",
            f"  {OUT / 'mc_width_histograms_theta_scan.pdf'}",
            f"  {OUT / 'uncertainty_budget_ranked.pdf'}",
        ]
    )
    summary_path = OUT / "uncertainty_publication_plots_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
