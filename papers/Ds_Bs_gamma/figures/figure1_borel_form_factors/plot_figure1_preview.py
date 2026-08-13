"""Python preview fallback for Figure 1.

The primary source is ``plot_figure1.wl``.  This fallback exists so the figure
can be reviewed while Wolfram/MaTeX startup is unavailable.
"""

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


ROOT = HERE.parents[1]
DATA = HERE / "figure1_borel_form_factors.csv"
OUT_PDF = HERE / "figure1_borel_form_factors_python.pdf"
OUT_PNG = HERE / "figure1_borel_form_factors_python.png"


def f(value: str) -> float:
    return float(value)


def rows_for(rows, *, state, s0):
    selected = []
    for row in rows:
        if row["state"] != state or row["scan"] != "M2":
            continue
        if not row["fixed_s0"] or abs(f(row["fixed_s0"]) - s0) > 1e-8:
            continue
        selected.append((f(row["M2"]), f(row["g_abs"])))
    return sorted(selected)


def y_limits(series, min_span):
    vals = [y for curve in series for _, y in curve]
    ymin, ymax = min(vals), max(vals)
    center = 0.5 * (ymin + ymax)
    half = max(0.5 * (ymax - ymin), 0.5 * min_span, 0.08 * abs(center), 0.002)
    return center - 1.25 * half, center + 1.25 * half


def add_compact_legend(ax, s0_values, styles, anchor):
    x0, y0 = anchor
    line_len = 0.10
    step = 0.055
    for index, (s0, style) in enumerate(zip(s0_values, styles)):
        y = y0 - index * step
        ax.plot(
            [x0, x0 + line_len],
            [y, y],
            transform=ax.transAxes,
            color=style["color"],
            linestyle=style["linestyle"],
            lw=1.9,
            clip_on=False,
        )
        ax.text(
            x0 + line_len + 0.018,
            y,
            rf"$s_0={s0:g}$",
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=8,
            color="0.12",
        )


def main() -> None:
    with DATA.open() as handle:
        rows = list(csv.DictReader(handle))

    specs = [
        {
            "state": "Ds1_2460",
            "s0": [8.5, 9.0, 9.5],
            "xlim": (3.0, 4.5),
            "min_span": 0.08,
            "legend_anchor": (0.08, 0.88),
            "title": r"$D_{s1}(2460)\to D_s\gamma$",
        },
        {
            "state": "Ds1_2536",
            "s0": [9.0, 9.5, 10.0],
            "xlim": (3.0, 4.5),
            "min_span": 0.025,
            "legend_anchor": (0.62, 0.88),
            "title": r"$D_{s1}(2536)\to D_s\gamma$",
        },
        {
            "state": "Bs1_5750",
            "s0": [39.0, 40.0, 41.0],
            "xlim": (10.0, 14.0),
            "min_span": 0.10,
            "legend_anchor": (0.58, 0.88),
            "title": r"$B_{s1}(5750)\to B_s\gamma$",
        },
        {
            "state": "Bs1_5830",
            "s0": [40.0, 41.0, 42.0],
            "xlim": (10.0, 14.0),
            "min_span": 0.025,
            "legend_anchor": (0.58, 0.88),
            "title": r"$B_{s1}(5830)\to B_s\gamma$",
        },
    ]

    styles = [
        {"color": "#0047ab", "linestyle": "-"},
        {"color": "#cc3300", "linestyle": "--"},
        {"color": "#111111", "linestyle": "-."},
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.4), constrained_layout=True)

    for ax, spec in zip(axes.flat, specs):
        series = [
            rows_for(
                rows,
                state=spec["state"],
                s0=s0,
            )
            for s0 in spec["s0"]
        ]
        for style, s0, curve in zip(styles, spec["s0"], series):
            xs, ys = zip(*curve)
            ax.plot(xs, ys, color=style["color"], linestyle=style["linestyle"], lw=1.9)
        ax.set_title(spec["title"], fontsize=11)
        ax.set_xlim(*spec["xlim"])
        ax.set_ylim(*y_limits(series, spec["min_span"]))
        ax.set_xlabel(r"$M^2\, [\mathrm{GeV}^2]$")
        ax.set_ylabel(r"$|g|\, [\mathrm{GeV}^{-1}]$")
        ax.grid(False)
        ax.tick_params(direction="in", top=True, right=True, labelsize=9)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        add_compact_legend(ax, spec["s0"], styles, spec["legend_anchor"])

    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=240)
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
