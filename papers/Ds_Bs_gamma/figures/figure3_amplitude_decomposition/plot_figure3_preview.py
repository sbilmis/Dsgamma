"""Python preview fallback for Figure 3."""

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


DATA = HERE / "figure3_amplitude_decomposition.csv"
OUT_PDF = HERE / "figure3_amplitude_decomposition_python.pdf"
OUT_PNG = HERE / "figure3_amplitude_decomposition_python.png"

PANELS = [
    ("Ds", r"$D_{s1}\to D_s\gamma$", ["Ds1_2460", "Ds1_2536"]),
    ("Bs", r"$B_{s1}\to B_s\gamma$", ["Bs1_5750", "Bs1_5830"]),
]
LABELS = {
    "Ds1_2460": r"$D_{s1}(2460)$",
    "Ds1_2536": r"$D_{s1}(2536)$",
    "Bs1_5750": r"$B_{s1}(5750)$",
    "Bs1_5830": r"$B_{s1}(5830)$",
}
COMPONENTS = [
    ("A_component_GeV_inv", r"$A$ piece", "#0047ab", "-"),
    ("B_component_GeV_inv", r"$B$ piece", "#cc3300", "--"),
    ("G_total_GeV_inv", r"$G=A+B$", "#111111", "-"),
]


def load_rows() -> list[dict[str, str]]:
    with DATA.open() as handle:
        return list(csv.DictReader(handle))


def add_legend(ax):
    x0, y0 = 0.60, 0.95
    for idx, (_, label, color, linestyle) in enumerate(COMPONENTS):
        y = y0 - idx * 0.06
        ax.plot([x0, x0 + 0.08], [y, y], transform=ax.transAxes, color=color, linestyle=linestyle, lw=2.2)
        ax.text(x0 + 0.10, y, label, transform=ax.transAxes, ha="left", va="center", fontsize=8.5)


def main() -> None:
    rows = load_rows()
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.35), constrained_layout=True)
    for ax, (sector, title, states) in zip(axes, PANELS):
        panel_rows = [next(row for row in rows if row["state"] == state) for state in states]
        y_base = [1.0, 0.0]
        offsets = [0.22, 0.0, -0.22]
        for key, label, color, linestyle in COMPONENTS:
            values = [float(row[key]) for row in panel_rows]
            for y, value in zip(y_base, values):
                ax.plot([0.0, value], [y + offsets[COMPONENTS.index((key, label, color, linestyle))]] * 2,
                        color=color, linestyle=linestyle, lw=3.0, solid_capstyle="butt")
                ax.plot(value, y + offsets[COMPONENTS.index((key, label, color, linestyle))],
                        marker="o", color=color, markersize=3.8)
        ax.axvline(0.0, color="0.45", lw=1.0)
        ax.set_title(title, fontsize=11)
        ax.set_yticks(y_base)
        ax.set_yticklabels([LABELS[state] for state in states], fontsize=9)
        ax.set_xlabel(r"$G$ contribution [GeV$^{-1}$]")
        ax.set_xlim(-0.90, 0.90)
        ax.set_ylim(-0.55, 1.55)
        ax.grid(False)
        ax.tick_params(direction="in", top=True, right=True, labelsize=9)
        add_legend(ax)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=240)
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
