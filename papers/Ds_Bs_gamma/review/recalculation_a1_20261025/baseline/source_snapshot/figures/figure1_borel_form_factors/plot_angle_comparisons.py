"""Generate Figure-1 comparisons at fixed 35.3 and 90 degree angles.

The transition OPE is the exact final implementation.  For each angle, the
Ds and Bs physical-current residues are reprojected from the corresponding
normalized AA/AB/BB two-point matrices, rather than retaining residues fitted
at a different angle.
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(HERE / ".cache"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from corrected_transition_central_python import bs_inputs, ds_inputs
from ds_angle_benchmark import projected_residues as ds_projected_residues
from final_stage2_uncertainty_scan import set_photon_shape_params
from rohrwild_transition_exact import (
    physical_couplings,
    precompute_convolutions,
    transition_invariants,
)
from twopoint_bs1_matrix_sumrule import (
    central_inputs as bs_twopoint_inputs,
    fitted_threshold as bs_fitted_threshold,
)
from twopoint_ds1_matrix_sumrule import projected_ope


ANGLES = (35.3, 90.0)

SPECS = (
    {
        "sector": "Ds",
        "title": r"$D_{s1}(2460)\to D_s\gamma$",
        "s0": (8.5, 9.0, 9.5),
        "M2": (3.0, 4.5),
        "combo": "low",
        "mass": 2.4595,
        "legend_anchor": (0.08, 0.88),
    },
    {
        "sector": "Ds",
        "title": r"$D_{s1}(2536)\to D_s\gamma$",
        "s0": (9.0, 9.5, 10.0),
        "M2": (3.0, 4.5),
        "combo": "high",
        "mass": 2.53511,
        "legend_anchor": (0.60, 0.88),
    },
    {
        "sector": "Bs",
        "title": r"$B_{s1}(5750)\to B_s\gamma$",
        "s0": (39.0, 40.0, 41.0),
        "M2": (10.0, 14.0),
        "combo": "low",
        "mass": 5.750,
        "legend_anchor": (0.58, 0.88),
    },
    {
        "sector": "Bs",
        "title": r"$B_{s1}(5830)\to B_s\gamma$",
        "s0": (40.0, 41.0, 42.0),
        "M2": (10.0, 14.0),
        "combo": "high",
        "mass": 5.82870,
        "legend_anchor": (0.58, 0.88),
    },
)

STYLES = (
    {"color": "#0047ab", "linestyle": "-"},
    {"color": "#cc3300", "linestyle": "--"},
    {"color": "#111111", "linestyle": "-."},
)


def bs_projected_residues(theta_deg: float) -> tuple[float, float]:
    inp = bs_twopoint_inputs()
    m2 = 7.0
    residues = []
    for channel, mass in enumerate((inp.mass_low, inp.mass_high)):
        s0, _ = bs_fitted_threshold(m2, theta_deg, channel, mass, inp)
        pi, _ = projected_ope(m2, s0, theta_deg, channel, inp)
        if pi <= 0.0:
            raise RuntimeError(f"nonpositive Bs projected residue at theta={theta_deg}")
        residues.append(math.sqrt(pi * math.exp(mass * mass / m2) / (mass * mass)))
    return float(residues[0]), float(residues[1])


def add_legend(ax, s0_values, anchor):
    x0, y0 = anchor
    for index, (s0, style) in enumerate(zip(s0_values, STYLES)):
        y = y0 - 0.055 * index
        ax.plot(
            [x0, x0 + 0.10],
            [y, y],
            transform=ax.transAxes,
            color=style["color"],
            linestyle=style["linestyle"],
            lw=1.9,
            clip_on=False,
        )
        ax.text(
            x0 + 0.118,
            y,
            rf"$s_0={s0:g}$",
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=8,
        )


def evaluate_angle(theta_deg: float):
    set_photon_shape_params(-2.1, 3.8)
    convolutions = precompute_convolutions()
    dvals = ds_inputs()
    ds_f1, ds_f2, _, _ = ds_projected_residues(theta_deg)
    bs_f1, bs_f2 = bs_projected_residues(theta_deg)
    output = []

    for spec in SPECS:
        curves = []
        for s0 in spec["s0"]:
            curve = []
            for m2 in np.linspace(*spec["M2"], 31):
                if spec["sector"] == "Ds":
                    vals = dvals
                    f1, f2 = ds_f1, ds_f2
                    m1, m2_state = 2.4595, 2.53511
                else:
                    vals = bs_inputs(float(spec["mass"]))
                    f1, f2 = bs_f1, bs_f2
                    m1, m2_state = 5.750, 5.82870
                inv = transition_invariants(
                    float(m2), float(s0), vals, convolutions=convolutions
                )
                phys = physical_couplings(
                    inv,
                    theta_deg=theta_deg,
                    m_state_1=m1,
                    m_state_2=m2_state,
                    f_1=f1,
                    f_2=f2,
                    m_p=vals["m_ds"],
                    f_p=vals["f_ds"],
                    m_q=vals["mc"],
                    m_s=vals["ms"],
                )
                key = "1" if spec["combo"] == "low" else "2"
                curve.append(
                    (
                        float(m2),
                        abs(float(phys[f"g_{key}"])),
                        float(phys[f"Gamma_{key}_keV"]),
                    )
                )
            curves.append(curve)
        output.append(curves)
    return output, (ds_f1, ds_f2), (bs_f1, bs_f2)


def plot(theta_deg: float, data, value_index: int, ylabel: str, stem: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.4), constrained_layout=True)
    for ax, spec, curves in zip(axes.flat, SPECS, data):
        for style, curve in zip(STYLES, curves):
            ax.plot(
                [row[0] for row in curve],
                [row[value_index] for row in curve],
                color=style["color"],
                linestyle=style["linestyle"],
                lw=1.9,
            )
        ax.set_title(spec["title"], fontsize=11)
        ax.set_xlim(*spec["M2"])
        ax.set_xlabel(r"$M^2\, [\mathrm{GeV}^2]$")
        ax.set_ylabel(ylabel)
        ax.grid(False)
        ax.tick_params(direction="in", top=True, right=True, labelsize=9)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        add_legend(ax, spec["s0"], spec["legend_anchor"])
        ax.text(
            0.96,
            0.06,
            rf"$\theta={theta_deg:g}^\circ$",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.5,
            color="0.22",
        )
    pdf = HERE / f"{stem}.pdf"
    png = HERE / f"{stem}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=240)
    plt.close(fig)
    print(f"Wrote {pdf}")
    print(f"Wrote {png}")


def main() -> None:
    for theta_deg in ANGLES:
        data, ds_residues, bs_residues = evaluate_angle(theta_deg)
        tag = str(theta_deg).replace(".", "p")
        print(
            f"theta={theta_deg:g}: Ds residues={ds_residues}; "
            f"Bs residues={bs_residues}"
        )
        plot(
            theta_deg,
            data,
            1,
            r"$|g|\, [\mathrm{GeV}^{-1}]$",
            f"figure1_borel_form_factors_theta_{tag}",
        )
        plot(
            theta_deg,
            data,
            2,
            r"$\Gamma\, [\mathrm{keV}]$",
            f"figure1_borel_decay_widths_theta_{tag}",
        )


if __name__ == "__main__":
    main()
