"""Python preview for the internal s0-stability diagnostic."""

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
from matplotlib.ticker import MaxNLocator


DATA = HERE / "figure1_s0_stability_diagnostic.csv"
OUT_PDF = HERE / "figure1_s0_stability_diagnostic_python.pdf"
OUT_PNG = HERE / "figure1_s0_stability_diagnostic_python.png"


SPECS = [
    ("Ds1_2460", r"$D_{s1}(2460)\to D_s\gamma$", (8.5, 9.5), 0.08),
    ("Ds1_2536", r"$D_{s1}(2536)\to D_s\gamma$", (9.0, 10.0), 0.04),
    ("Bs1_5750", r"$B_{s1}(5750)\to B_s\gamma$", (39.0, 41.0), 0.10),
    ("Bs1_5830", r"$B_{s1}(5830)\to B_s\gamma$", (40.0, 42.0), 0.025),
]

STYLES = [
    {"color": "#0047ab", "linestyle": "-"},
    {"color": "#cc3300", "linestyle": "--"},
    {"color": "#111111", "linestyle": "-."},
]


def load_rows() -> list[dict[str, str]]:
    with DATA.open() as handle:
        return list(csv.DictReader(handle))


def y_limits(curves, min_span):
    vals = [y for curve in curves for _, y in curve]
    ymin, ymax = min(vals), max(vals)
    center = 0.5 * (ymin + ymax)
    half = max(0.5 * (ymax - ymin), 0.5 * min_span, 0.08 * abs(center), 0.002)
    return center - 1.25 * half, center + 1.25 * half


def rows_for(rows, state, m2):
    selected = [
        (float(row["s0"]), float(row["g_abs"]))
        for row in rows
        if row["state"] == state and abs(float(row["M2"]) - m2) < 1e-8
    ]
    return sorted(selected)


def add_legend(ax, m2_values, anchor=(0.58, 0.88)):
    x0, y0 = anchor
    line_len = 0.10
    step = 0.055
    for index, (m2, style) in enumerate(zip(m2_values, STYLES)):
        y = y0 - index * step
        ax.plot([x0, x0 + line_len], [y, y], transform=ax.transAxes,
                color=style["color"], linestyle=style["linestyle"], lw=1.9)
        ax.text(x0 + line_len + 0.018, y, rf"$M^2={m2:g}$",
                transform=ax.transAxes, ha="left", va="center", fontsize=8)


def main() -> None:
    rows = load_rows()
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.4), constrained_layout=True)
    for ax, (state, title, xlim, min_span) in zip(axes.flat, SPECS):
        m2_values = sorted({float(row["M2"]) for row in rows if row["state"] == state})
        curves = [rows_for(rows, state, m2) for m2 in m2_values]
        for style, curve in zip(STYLES, curves):
            xs, ys = zip(*curve)
            ax.plot(xs, ys, color=style["color"], linestyle=style["linestyle"], lw=1.9)
        ax.set_title(title, fontsize=11)
        ax.set_xlim(*xlim)
        ax.set_ylim(*y_limits(curves, min_span))
        ax.set_xlabel(r"$s_0\, [\mathrm{GeV}^2]$")
        ax.set_ylabel(r"$|g|\, [\mathrm{GeV}^{-1}]$")
        ax.grid(False)
        ax.tick_params(direction="in", top=True, right=True, labelsize=9)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        add_legend(ax, m2_values)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=240)
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
