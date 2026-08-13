"""Python preview fallback for Figure 4."""

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


DATA = HERE / "figure4_literature_comparison.csv"
OUT_PDF = HERE / "figure4_literature_comparison_python.pdf"
OUT_PNG = HERE / "figure4_literature_comparison_python.png"


def load_rows() -> list[dict[str, str]]:
    with DATA.open() as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    rows = load_rows()
    labels = [
        r"$D_{s1}(2460)$",
        r"$D_{s1}(2536)$",
        r"$B_{s1}(5750)$",
        r"$B_{s1}(5830)$",
    ]
    x = list(range(len(rows)))
    this = [float(row["this_work_median_keV"]) for row in rows]
    p16 = [float(row["this_work_p16_keV"]) for row in rows]
    p84 = [float(row["this_work_p84_keV"]) for row in rows]
    bondar = [float(row["bondar_milstein_2025_keV"]) for row in rows]
    yerr = [[m - lo for m, lo in zip(this, p16)], [hi - m for m, hi in zip(this, p84)]]

    fig, ax = plt.subplots(figsize=(7.1, 3.45), constrained_layout=True)
    ax.errorbar(
        [v - 0.08 for v in x],
        this,
        yerr=yerr,
        fmt="o",
        markersize=4.8,
        capsize=3.0,
        color="#0047ab",
        ecolor="#0047ab",
        elinewidth=1.5,
        label="This work",
    )
    ax.plot(
        [v + 0.08 for v in x],
        bondar,
        linestyle="None",
        marker="s",
        markersize=4.8,
        color="#111111",
        label="Bondar-Milstein",
    )
    ax.set_yscale("log")
    ax.set_ylim(0.02, 500.0)
    ax.set_xlim(-0.45, len(rows) - 0.55)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(r"$\Gamma$ [keV]")
    ax.set_title("Comparison with selected literature estimates", fontsize=11)
    ax.grid(False)
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=9)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=240)
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
