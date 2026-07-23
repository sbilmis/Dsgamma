"""Focused Bs1 Borel-window diagnostics.

This script reuses the diagnostics from ``window_diagnostics.py`` but restricts
the comparison to bottom-strange channels and to candidate Borel intervals.
The purpose is to test whether lowering the Borel window improves the pole
proxy without spoiling stability or the light-cone OPE hierarchy.
"""

from __future__ import annotations

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

from stage2_axial_g1_three_particle import F1_integral
from final_stage2_uncertainty_scan import precompute_f1_basis
from window_diagnostics import sample_window, summarize_window, write_csv


BS_LOW_BASE = {
    "sector": "Bs",
    "state": r"B_{s1}(5750)",
    "mass": 5.750,
    "s0_min": 6.05**2,
    "s0_max": 6.35**2,
    "quoted": "low",
}

BS_5830_BASE = {
    "sector": "Bs",
    "state": r"B_{s1}(5830)",
    "mass": 5.82870,
    "s0_min": 6.15**2,
    "s0_max": 6.40**2,
    "quoted": "high",
}

M2_WINDOWS = [
    (8.0, 12.0),
    (9.0, 13.0),
    (10.0, 14.0),
    (8.0, 14.0),
]


def build_windows() -> list[dict[str, float | str]]:
    windows: list[dict[str, float | str]] = []
    for base, tag in ((BS_LOW_BASE, "low"), (BS_5830_BASE, "5830")):
        for m2_min, m2_max in M2_WINDOWS:
            window = dict(base)
            window.update(
                {
                    "name": f"Bs{tag}_M2_{m2_min:.0f}_{m2_max:.0f}",
                    "M2_min": m2_min,
                    "M2_max": m2_max,
                }
            )
            windows.append(window)
    return windows


def make_plot(summary_rows: list[dict[str, float | str]]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.4), sharex=True)
    labels = [str(r["window"]).replace("Bs", "").replace("_M2_", "\n") for r in summary_rows]
    x = np.arange(len(labels))

    axes[0].bar(
        x - 0.18,
        [float(r["min_spectral_pole_axial"]) for r in summary_rows],
        width=0.36,
        label="axial pole proxy",
        color="#1f77b4",
    )
    axes[0].bar(
        x + 0.18,
        [float(r["min_spectral_pole_tensor"]) for r in summary_rows],
        width=0.36,
        label="tensor pole proxy",
        color="#2ca02c",
    )
    axes[0].axhline(0.5, color="black", linestyle="--", linewidth=1.0, alpha=0.65)
    axes[0].set_ylabel("minimum pole proxy")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(True, axis="y", alpha=0.24)

    axes[1].bar(
        x - 0.18,
        [float(r["Gamma_quoted_keV_span"]) for r in summary_rows],
        width=0.36,
        label=r"$\Gamma$ span",
        color="#d62728",
    )
    axes[1].bar(
        x + 0.18,
        [float(r["max_GA_3p_frac"]) for r in summary_rows],
        width=0.36,
        label="max 3p fraction",
        color="#9467bd",
    )
    axes[1].set_ylabel("relative diagnostic")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(True, axis="y", alpha=0.24)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=0)
    fig.tight_layout()
    fig.savefig(OUT / "bs1_borel_window_candidate_scan.png", dpi=300)
    fig.savefig(OUT / "bs1_borel_window_candidate_scan.pdf")


def main() -> None:
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    all_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []

    for window in build_windows():
        rows = sample_window(window, f1_axial, f1_basis)
        all_rows.extend(rows)
        summary_rows.append(summarize_window(window, rows))

    write_csv(OUT / "bs1_borel_window_candidate_grid.csv", all_rows)
    write_csv(OUT / "bs1_borel_window_candidate_summary.csv", summary_rows)
    make_plot(summary_rows)

    lines = [
        "Bs1 Borel-window candidate scan",
        "================================",
        "All rows use central lattice photon normalization and theta=35.3 deg.",
        "The s0 ranges are kept fixed to isolate the M^2-window dependence.",
        "",
    ]
    for row in summary_rows:
        lines.append(
            "{window}: {state}, M2={M2}, s0={s0} GeV^2; "
            "poleA={pa:.3f}, poleB={pt:.3f}, twist2={tw:.3f}, "
            "Gamma_span={gs:.3f}, max3pA={ga3:.4f}, max3pB={gb3:.4f}".format(
                window=row["window"],
                state=row["state"],
                M2=row["M2_range_GeV2"],
                s0=row["s0_range_GeV2"],
                pa=float(row["min_spectral_pole_axial"]),
                pt=float(row["min_spectral_pole_tensor"]),
                tw=float(row["min_twist2_retention"]),
                gs=float(row["Gamma_quoted_keV_span"]),
                ga3=float(row["max_GA_3p_frac"]),
                gb3=float(row["max_GB_3p_frac"]),
            )
        )

    lines.extend(
        [
            "",
            "Interpretation:",
            "- Lowering M^2 improves the expected pole contribution only if the minimum pole proxy rises.",
            "- It is acceptable only if the three-particle fractions stay small and the width span does not worsen too much.",
            "- The comparison is diagnostic; the quoted window should not be changed merely to satisfy one number.",
            "",
            f"Wrote {OUT / 'bs1_borel_window_candidate_grid.csv'}",
            f"Wrote {OUT / 'bs1_borel_window_candidate_summary.csv'}",
            f"Wrote {OUT / 'bs1_borel_window_candidate_scan.pdf'}",
        ]
    )
    text = "\n".join(lines) + "\n"
    (OUT / "bs1_borel_window_candidate_summary.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
