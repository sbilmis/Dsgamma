"""Plot two-point residues and pole fractions versus the mixing angle.

This is an internal diagnostic for the common auxiliary windows requested for
the two physical states in each flavour sector.  Central QCD inputs are used.
The solid/dashed curves use the centres of the M^2 and s0 intervals, while the
bands show the full envelope obtained by scanning those intervals.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from twopoint_bs1_matrix_sumrule import central_inputs as bs_inputs
from twopoint_ds1_matrix_sumrule import (
    CumulativePerturbative,
    Inputs,
    local_components,
)


ANGLES = np.linspace(0.0, 90.0, 361)

SECTORS = {
    "Ds": {
        "title": r"$D_{s1}$ two-point sum rules",
        "inputs": Inputs(),
        "m2": (2.7, 3.3),
        "s0": ((8.8, 9.8), (8.9, 9.9)),
        "theta": 26.6,
        "states": (r"$D_{s1}(2460)$", r"$D_{s1}(2536)$"),
    },
    "Bs": {
        "title": r"$B_{s1}$ two-point sum rules",
        "inputs": bs_inputs(),
        "m2": (7.1, 8.1),
        "s0": ((41.0, 43.0), (42.0, 44.0)),
        "theta": 38.5,
        "states": (r"$B_{s1}(5750)$", r"$B_{s1}(5830)$"),
    },
}


def projected(matrix: np.ndarray, channel: int) -> np.ndarray:
    radians = np.deg2rad(ANGLES)
    if channel == 0:
        x, y = np.sin(radians), np.cos(radians)
    else:
        x, y = np.cos(radians), -np.sin(radians)
    return matrix[0, 0] * x**2 + 2.0 * matrix[0, 1] * x * y + matrix[1, 1] * y**2


def curves(m2: float, s0: float, channel: int, inp: Inputs) -> tuple[np.ndarray, np.ndarray]:
    evaluator = CumulativePerturbative(inp, (m2,), n_grid=12000)
    perturbative = evaluator.matrix(m2, s0)
    full = evaluator.full(m2)
    total = perturbative + sum(
        local_components(m2, inp).values(), np.zeros((2, 2), dtype=float)
    )
    pole_mass = inp.mass_low if channel == 0 else inp.mass_high
    pi_values = projected(total, channel)
    residues = np.where(
        pi_values > 0.0,
        np.sqrt(pi_values * np.exp(pole_mass**2 / m2) / pole_mass**2) * 1000.0,
        np.nan,
    )
    denominator = projected(full, channel)
    pole = np.divide(
        projected(perturbative, channel),
        denominator,
        out=np.full_like(denominator, np.nan),
        where=np.abs(denominator) > 0.0,
    )
    return residues, pole


def central_and_envelope(spec: dict, channel: int) -> tuple[np.ndarray, ...]:
    m2_lo, m2_hi = spec["m2"]
    s0_lo, s0_hi = spec["s0"][channel]
    central = curves(
        0.5 * (m2_lo + m2_hi),
        0.5 * (s0_lo + s0_hi),
        channel,
        spec["inputs"],
    )
    residue_scan: list[np.ndarray] = []
    pole_scan: list[np.ndarray] = []
    for m2 in np.linspace(m2_lo, m2_hi, 9):
        evaluator = CumulativePerturbative(spec["inputs"], (float(m2),), n_grid=12000)
        full = evaluator.full(float(m2))
        locals_sum = sum(
            local_components(float(m2), spec["inputs"]).values(),
            np.zeros((2, 2), dtype=float),
        )
        pole_mass = spec["inputs"].mass_low if channel == 0 else spec["inputs"].mass_high
        denominator = projected(full, channel)
        for s0 in np.linspace(s0_lo, s0_hi, 9):
            perturbative = evaluator.matrix(float(m2), float(s0))
            pi_values = projected(perturbative + locals_sum, channel)
            residue_scan.append(
                np.where(
                    pi_values > 0.0,
                    np.sqrt(pi_values * np.exp(pole_mass**2 / m2) / pole_mass**2) * 1000.0,
                    np.nan,
                )
            )
            pole_scan.append(
                np.divide(
                    projected(perturbative, channel),
                    denominator,
                    out=np.full_like(denominator, np.nan),
                    where=np.abs(denominator) > 0.0,
                )
            )
    residues = np.asarray(residue_scan)
    poles = np.asarray(pole_scan)
    return (
        central[0],
        np.nanmin(residues, axis=0),
        np.nanmax(residues, axis=0),
        central[1],
        np.nanmin(poles, axis=0),
        np.nanmax(poles, axis=0),
    )


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.3), sharex="col")
    styles = (("#1558b0", "-"), ("#111111", "--"))

    for column, sector in enumerate(("Ds", "Bs")):
        spec = SECTORS[sector]
        ax_f = axes[0, column]
        ax_p = axes[1, column]
        for channel, (label, (color, linestyle)) in enumerate(zip(spec["states"], styles)):
            f, f_lo, f_hi, pole, pole_lo, pole_hi = central_and_envelope(spec, channel)
            ax_f.fill_between(ANGLES, f_lo, f_hi, color=color, alpha=0.13, linewidth=0)
            ax_f.plot(ANGLES, f, color=color, linestyle=linestyle, linewidth=1.8, label=label)
            ax_p.fill_between(ANGLES, pole_lo, pole_hi, color=color, alpha=0.13, linewidth=0)
            ax_p.plot(ANGLES, pole, color=color, linestyle=linestyle, linewidth=1.8, label=label)

        for axis in (ax_f, ax_p):
            axis.axvline(spec["theta"], color="#c92828", linestyle="-.", linewidth=1.25)
            axis.axvline(35.3, color="#777777", linestyle=":", linewidth=1.15)
            axis.set_xlim(0.0, 90.0)
            axis.grid(alpha=0.16, linewidth=0.5)
        ax_p.axhline(0.50, color="#777777", linestyle=":", linewidth=1.0)
        ax_f.set_title(spec["title"])
        ax_f.set_ylabel(r"$f_i(\theta)$ [MeV]")
        ax_p.set_ylabel("Pole fraction")
        ax_p.set_xlabel(r"$\theta$ [deg]")
        ax_p.set_ylim(bottom=0.0)
        ax_f.legend(frameon=False, loc="best")
        ax_p.text(
            0.985,
            0.505,
            "0.50",
            transform=ax_p.get_yaxis_transform(),
            ha="right",
            va="bottom",
            color="#666666",
            fontsize=8,
        )
        ax_f.text(
            spec["theta"] + 1.0,
            0.96,
            rf"${spec['theta']:.1f}^\circ$",
            transform=ax_f.get_xaxis_transform(),
            color="#c92828",
            ha="left",
            va="top",
            fontsize=9,
        )

    fig.suptitle(
        r"Two-point residues and pole fractions versus the mixing angle",
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    stem = OUT / "twopoint_residues_pole_fractions_vs_theta_common_windows"
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    print(stem.with_suffix(".png"))
    print(stem.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
