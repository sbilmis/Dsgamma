"""Python preview fallback for Figure 5."""

from __future__ import annotations

import csv
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(HERE / ".cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DATA = HERE / "figure5_monte_carlo_widths.csv"
OUT_PDF = HERE / "figure5_monte_carlo_widths_python.pdf"
OUT_PNG = HERE / "figure5_monte_carlo_widths_python.png"

PANELS = [
    ("Ds1_2460", r"$D_{s1}(2460)\to D_s\gamma$", (0.0, 78.0)),
    ("Ds1_2536", r"$D_{s1}(2536)\to D_s\gamma$", (0.0, 84.0)),
    ("Bs1_5750", r"$B_{s1}(5750)\to B_s\gamma$", (0.0, 27.0)),
    ("Bs1_5830", r"$B_{s1}(5830)\to B_s\gamma$", (0.0, 43.0)),
]


def load_rows() -> list[dict[str, str]]:
    with DATA.open() as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    rows = load_rows()
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.2), constrained_layout=True)
    for ax, (state_key, title, xlim) in zip(axes.flat, PANELS):
        values = np.asarray([float(row["Gamma_keV"]) for row in rows if row["state_key"] == state_key], dtype=float)
        median = float(np.percentile(values, 50.0))
        p16 = float(np.percentile(values, 16.0))
        p84 = float(np.percentile(values, 84.0))
        ax.hist(values, bins=28, range=xlim, density=True, color="#0047ab", alpha=0.55, edgecolor="white", linewidth=0.5)
        ax.axvline(median, color="#111111", lw=1.4)
        ax.axvline(p16, color="#cc3300", linestyle="--", lw=1.1)
        ax.axvline(p84, color="#cc3300", linestyle="--", lw=1.1)
        ax.set_xlim(*xlim)
        ax.set_title(title, fontsize=10.5)
        ax.set_xlabel(r"$\Gamma$ [keV]")
        ax.set_ylabel("density")
        ax.grid(False)
        ax.tick_params(direction="in", top=True, right=True, labelsize=8.5)
        ax.text(
            0.96,
            0.90,
            rf"${median:.3g}_{{-{median-p16:.2g}}}^{{+{p84-median:.2g}}}$ keV",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.5,
        )
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=240)
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
