"""Paper-facing LCSR stability plots at the previously determined angles.

The nominal angles are not refitted in the radiative analysis:

    theta_Ds = 26.6 +/- 0.6 deg,
    theta_Bs = 38.5 +/- 0.1 deg.

They are imported from ``mixing_angle_inputs.py`` and originate from the
dedicated analysis arXiv:2408.08014.  The physical-current residues are fixed
to the medians of the accepted final normalization samples.  This deliberately
isolates the Borel/threshold stability of the transition LCSR from a second
scan over the two-point normalization.
"""

from __future__ import annotations

import csv
import math
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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bs1_stage2_baseline import evaluate_point as evaluate_bs_basis
from final_stage2_uncertainty_scan import precompute_f1_basis, set_photon_shape_params
from lattice_photon_normalization_comparison import (
    FPERP_S_1GEV,
    central_inputs_for_scenario,
    evaluate as evaluate_ds_basis,
)
from mixing_angle_inputs import (
    ANGLE_SOURCE,
    THETA_BS_DEG,
    THETA_BS_SIGMA_DEG,
    THETA_DS_DEG,
    THETA_DS_SIGMA_DEG,
)
from stage1_tensor_gb_soft2p_estimate import width_keV
from stage2_axial_g1_three_particle import F1_integral


# Medians of the accepted final-window physical-current normalization samples.
F1_DS = 0.40491088740235476
F2_DS = 0.16827627627950853
F2_BS_5830 = 0.21452584195712590

DS_M2_WINDOW = (3.0, 4.5)
DS_S0_VALUES = (7.5, 8.0, 8.5)
BS_M2_WINDOW = (10.0, 14.0)
BS_S0_VALUES = (38.0, 39.5, 41.0)


def ds_inputs() -> dict[str, float]:
    vals = central_inputs_for_scenario("lattice_fperp_s")
    set_photon_shape_params(vals["omegaA"], vals["omegaV"])
    return vals


def bs_inputs() -> dict[str, float]:
    ss = 0.8 * (-(0.240**3))
    set_photon_shape_params(-2.1, 3.8)
    return {
        "mc": 4.18,
        "ms": 0.093,
        "m_ds1": 5.82870,
        "m_ds": 5.36692,
        "f_ds1": 0.305,
        "f_ds": 0.2303,
        "fT": 0.285,
        "ss": ss,
        "chi": FPERP_S_1GEV / ss,
        "f3g": -0.0039,
        "ec": -1.0 / 3.0,
        "es": -1.0 / 3.0,
        "fperp_s_used": FPERP_S_1GEV,
    }


def ds_widths(
    vals: dict[str, float],
    M2: float,
    s0: float,
    theta_deg: float,
    f1_axial: float,
    f1_basis: dict[str, float],
) -> tuple[float, float]:
    basis = evaluate_ds_basis(
        vals, M2, math.sqrt(s0), theta_deg, f1_axial, f1_basis
    )
    theta = math.radians(theta_deg)
    s, c = math.sin(theta), math.cos(theta)
    f_a = vals["f_ds1"]
    f_b = vals["fT"] * vals["m_ds1"] / (vals["mc"] + vals["ms"])
    g_2460 = (s * f_a * basis["GA"] + c * f_b * basis["GB"]) / F1_DS
    g_2536 = (c * f_a * basis["GA"] - s * f_b * basis["GB"]) / F2_DS
    return (
        width_keV(vals["m_ds1"], vals["m_ds"], g_2460),
        width_keV(vals["m_ds1_2536"], vals["m_ds"], g_2536),
    )


def bs_width(
    vals: dict[str, float],
    M2: float,
    s0: float,
    theta_deg: float,
    f1_axial: float,
    f1_basis: dict[str, float],
) -> float:
    basis = evaluate_bs_basis(
        vals, M2, math.sqrt(s0), theta_deg, f1_axial, f1_basis
    )
    theta = math.radians(theta_deg)
    c, s = math.cos(theta), math.sin(theta)
    g_5830 = (
        c * vals["f_ds1"] * basis["GA"]
        - s * basis["fB_effective"] * basis["GB"]
    ) / F2_BS_5830
    return width_keV(vals["m_ds1"], vals["m_ds"], g_5830)


def angle_envelope(func, central: float, sigma: float) -> tuple[float, float]:
    values = [func(central - sigma), func(central), func(central + sigma)]
    return min(values), max(values)


def make_rows() -> list[dict[str, float | str]]:
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    dvals = ds_inputs()
    bvals = bs_inputs()
    rows: list[dict[str, float | str]] = []

    for s0 in DS_S0_VALUES:
        for M2 in np.linspace(*DS_M2_WINDOW, 31):
            for theta in (
                THETA_DS_DEG - THETA_DS_SIGMA_DEG,
                THETA_DS_DEG,
                THETA_DS_DEG + THETA_DS_SIGMA_DEG,
            ):
                w1, w2 = ds_widths(
                    dvals, float(M2), s0, theta, f1_axial, f1_basis
                )
                rows.extend(
                    [
                        {
                            "sector": "Ds",
                            "state": "Ds1_2460",
                            "scan": "M2",
                            "M2_GeV2": float(M2),
                            "s0_GeV2": s0,
                            "theta_deg": theta,
                            "width_keV": w1,
                        },
                        {
                            "sector": "Ds",
                            "state": "Ds1_2536",
                            "scan": "M2",
                            "M2_GeV2": float(M2),
                            "s0_GeV2": s0,
                            "theta_deg": theta,
                            "width_keV": w2,
                        },
                    ]
                )

    for s0 in BS_S0_VALUES:
        for M2 in np.linspace(*BS_M2_WINDOW, 31):
            for theta in (
                THETA_BS_DEG - THETA_BS_SIGMA_DEG,
                THETA_BS_DEG,
                THETA_BS_DEG + THETA_BS_SIGMA_DEG,
            ):
                rows.append(
                    {
                        "sector": "Bs",
                        "state": "Bs1_5830",
                        "scan": "M2",
                        "M2_GeV2": float(M2),
                        "s0_GeV2": s0,
                        "theta_deg": theta,
                        "width_keV": bs_width(
                            bvals, float(M2), s0, theta, f1_axial, f1_basis
                        ),
                    }
                )

    for s0 in np.linspace(DS_S0_VALUES[0], DS_S0_VALUES[-1], 41):
        for theta in (
            THETA_DS_DEG - THETA_DS_SIGMA_DEG,
            THETA_DS_DEG,
            THETA_DS_DEG + THETA_DS_SIGMA_DEG,
        ):
            w1, w2 = ds_widths(
                dvals, sum(DS_M2_WINDOW) / 2.0, float(s0), theta, f1_axial, f1_basis
            )
            rows.extend(
                [
                    {
                        "sector": "Ds",
                        "state": "Ds1_2460",
                        "scan": "s0",
                        "M2_GeV2": sum(DS_M2_WINDOW) / 2.0,
                        "s0_GeV2": float(s0),
                        "theta_deg": theta,
                        "width_keV": w1,
                    },
                    {
                        "sector": "Ds",
                        "state": "Ds1_2536",
                        "scan": "s0",
                        "M2_GeV2": sum(DS_M2_WINDOW) / 2.0,
                        "s0_GeV2": float(s0),
                        "theta_deg": theta,
                        "width_keV": w2,
                    },
                ]
            )

    for s0 in np.linspace(BS_S0_VALUES[0], BS_S0_VALUES[-1], 41):
        for theta in (
            THETA_BS_DEG - THETA_BS_SIGMA_DEG,
            THETA_BS_DEG,
            THETA_BS_DEG + THETA_BS_SIGMA_DEG,
        ):
            rows.append(
                {
                    "sector": "Bs",
                    "state": "Bs1_5830",
                    "scan": "s0",
                    "M2_GeV2": sum(BS_M2_WINDOW) / 2.0,
                    "s0_GeV2": float(s0),
                    "theta_deg": theta,
                    "width_keV": bs_width(
                        bvals,
                        sum(BS_M2_WINDOW) / 2.0,
                        float(s0),
                        theta,
                        f1_axial,
                        f1_basis,
                    ),
                }
            )
    return rows


def write_rows(rows: list[dict[str, float | str]]) -> None:
    path = OUT / "final_angle_stability.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def subset(
    rows: list[dict[str, float | str]],
    state: str,
    scan: str,
    theta: float,
    s0: float | None = None,
) -> list[dict[str, float | str]]:
    selected = [
        row
        for row in rows
        if row["state"] == state
        and row["scan"] == scan
        and math.isclose(float(row["theta_deg"]), theta)
        and (s0 is None or math.isclose(float(row["s0_GeV2"]), s0))
    ]
    key = "M2_GeV2" if scan == "M2" else "s0_GeV2"
    return sorted(selected, key=lambda row: float(row[key]))


def plot_rows(rows: list[dict[str, float | str]]) -> None:
    colors = ["#2b6cb0", "#2f855a", "#c05621"]
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.3))

    panels = (
        (
            axes[0, 0],
            "Ds1_2460",
            DS_S0_VALUES,
            THETA_DS_DEG,
            THETA_DS_SIGMA_DEG,
            r"$D_{s1}(2460)\to D_s\gamma$",
        ),
        (
            axes[0, 1],
            "Ds1_2536",
            DS_S0_VALUES,
            THETA_DS_DEG,
            THETA_DS_SIGMA_DEG,
            r"$D_{s1}(2536)\to D_s\gamma$",
        ),
        (
            axes[1, 0],
            "Bs1_5830",
            BS_S0_VALUES,
            THETA_BS_DEG,
            THETA_BS_SIGMA_DEG,
            r"$B_{s1}(5830)\to B_s\gamma$",
        ),
    )
    for ax, state, s0_values, theta, sigma, title in panels:
        for color, s0 in zip(colors, s0_values):
            central = subset(rows, state, "M2", theta, s0)
            low = subset(rows, state, "M2", theta - sigma, s0)
            high = subset(rows, state, "M2", theta + sigma, s0)
            x = np.asarray([float(row["M2_GeV2"]) for row in central])
            y = np.asarray([float(row["width_keV"]) for row in central])
            y_min = np.minimum(
                [float(row["width_keV"]) for row in low],
                [float(row["width_keV"]) for row in high],
            )
            y_max = np.maximum(
                [float(row["width_keV"]) for row in low],
                [float(row["width_keV"]) for row in high],
            )
            ax.fill_between(x, y_min, y_max, color=color, alpha=0.12)
            ax.plot(x, y, color=color, lw=1.7, label=rf"$s_0={s0:g}$ GeV$^2$")
        ax.set_title(title + rf", $\theta={theta:.1f}^\circ$")
        ax.set_xlabel(r"$M^2$ [GeV$^2$]")
        ax.set_ylabel(r"$\Gamma$ [keV]")
        ax.grid(alpha=0.25, lw=0.7)
        ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    threshold_specs = (
        ("Ds1_2460", THETA_DS_DEG, THETA_DS_SIGMA_DEG, "#2b6cb0", r"$D_{s1}(2460)$"),
        ("Ds1_2536", THETA_DS_DEG, THETA_DS_SIGMA_DEG, "#805ad5", r"$D_{s1}(2536)$"),
        ("Bs1_5830", THETA_BS_DEG, THETA_BS_SIGMA_DEG, "#c05621", r"$B_{s1}(5830)$"),
    )
    for state, theta, sigma, color, label in threshold_specs:
        central = subset(rows, state, "s0", theta)
        low = subset(rows, state, "s0", theta - sigma)
        high = subset(rows, state, "s0", theta + sigma)
        s0_axis = np.asarray([float(row["s0_GeV2"]) for row in central])
        s0_mid = 0.5 * (s0_axis[0] + s0_axis[-1])
        s0_half_range = 0.5 * (s0_axis[-1] - s0_axis[0])
        x = (s0_axis - s0_mid) / s0_half_range
        y = np.asarray([float(row["width_keV"]) for row in central])
        y0 = y[len(y) // 2]
        y_min = np.minimum(
            [float(row["width_keV"]) for row in low],
            [float(row["width_keV"]) for row in high],
        ) / y0
        y_max = np.maximum(
            [float(row["width_keV"]) for row in low],
            [float(row["width_keV"]) for row in high],
        ) / y0
        ax.fill_between(x, y_min, y_max, color=color, alpha=0.12)
        ax.plot(x, y / y0, color=color, lw=1.7, label=label)
    ax.set_title("Continuum-threshold stability")
    ax.set_xlabel(r"$(s_0-s_{0,\rm mid})/\Delta s_0$")
    ax.set_ylabel(r"$\Gamma(s_0)/\Gamma(s_{0,\rm mid})$")
    ax.grid(alpha=0.25, lw=0.7)
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Final-angle transition-LCSR stability", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "final_angle_stability.pdf", bbox_inches="tight")
    fig.savefig(OUT / "final_angle_stability.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_summary(rows: list[dict[str, float | str]]) -> None:
    lines = [
        "Final-angle transition-LCSR stability",
        "=====================================",
        f"Angle source: {ANGLE_SOURCE}",
        f"theta_Ds = {THETA_DS_DEG:.1f} +/- {THETA_DS_SIGMA_DEG:.1f} deg",
        f"theta_Bs = {THETA_BS_DEG:.1f} +/- {THETA_BS_SIGMA_DEG:.1f} deg",
        "The angles are external inputs and are not refitted here.",
        f"Fixed physical residues: f1_Ds={F1_DS:.6f} GeV, "
        f"f2_Ds={F2_DS:.6f} GeV, f2_Bs5830={F2_BS_5830:.6f} GeV.",
        "The shaded bands show only the quoted angle uncertainty.",
        "Ordinary local condensates are excluded from the transition LCSR.",
        "",
    ]
    for state, theta in (
        ("Ds1_2460", THETA_DS_DEG),
        ("Ds1_2536", THETA_DS_DEG),
        ("Bs1_5830", THETA_BS_DEG),
    ):
        values = [
            float(row["width_keV"])
            for row in rows
            if row["state"] == state
            and row["scan"] == "M2"
            and math.isclose(float(row["theta_deg"]), theta)
        ]
        lines.append(
            f"{state}: central-angle Borel/threshold envelope "
            f"{min(values):.6g}--{max(values):.6g} keV"
        )
    (OUT / "final_angle_stability_summary.txt").write_text("\n".join(lines) + "\n")


def main() -> None:
    rows = make_rows()
    write_rows(rows)
    plot_rows(rows)
    write_summary(rows)
    print(f"Wrote {OUT / 'final_angle_stability.csv'}")
    print(f"Wrote {OUT / 'final_angle_stability.pdf'}")
    print(f"Wrote {OUT / 'final_angle_stability.png'}")
    print(f"Wrote {OUT / 'final_angle_stability_summary.txt'}")


if __name__ == "__main__":
    main()
