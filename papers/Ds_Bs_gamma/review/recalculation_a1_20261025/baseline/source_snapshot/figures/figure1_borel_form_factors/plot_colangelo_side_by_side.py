"""Build an internal side-by-side pure-axial benchmark against Colangelo Fig. 5.

The left panel is the user-supplied image of the published figure. The right
panel evaluates our axial-current implementation at theta=90 degrees using
the Colangelo-local transition prescription, paper-like inputs, and exactly
the three continuum thresholds displayed in the published caption.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(HERE / ".cache"))
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from stage1_axial_g1_baseline import colangelo_inputs
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from validation_checks import photon_params
from corrected_transition_central_python import ds_inputs
from ds_angle_benchmark import projected_residues
from final_stage2_uncertainty_scan import set_photon_shape_params
from rohrwild_transition_exact import (
    physical_couplings,
    precompute_convolutions,
    transition_invariants,
)


THRESHOLDS = (2.50, 2.55, 2.60)
STYLES = (
    {"color": "red", "linestyle": "-", "label": r"$\sqrt{s_0}=2.50$ GeV"},
    {"color": "#22bb22", "linestyle": (0, (7, 6)), "label": r"$\sqrt{s_0}=2.55$ GeV"},
    {"color": "#174cff", "linestyle": (0, (2, 2)), "label": r"$\sqrt{s_0}=2.60$ GeV"},
)


def calculate() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    m2_values = np.linspace(2.0, 7.0, 101)
    # These are the BBK parameters used by the existing Colangelo validation.
    with photon_params(
        kappa=0.2,
        kappa_p=0.0,
        zeta1=0.4,
        zeta1_p=0.0,
        zeta2=0.3,
        zeta2_p=0.0,
    ):
        inputs = colangelo_inputs()
        f1_total, _, _ = F1_integral(u0=0.5)
        for root_s0 in THRESHOLDS:
            for m2 in m2_values:
                result = g1_stage2(
                    float(m2),
                    root_s0**2,
                    inputs,
                    f1_total,
                    transition_scheme="colangelo_local_benchmark",
                )
                rows.append(
                    {
                        "theta_deg": 90.0,
                        "scenario": "colangelo_matched_reconstruction",
                        "M2_GeV2": float(m2),
                        "sqrt_s0_GeV": root_s0,
                        "s0_GeV2": root_s0**2,
                        "g1_GeV_inv": float(result["g1_stage2_GeV_inv"]),
                        "Gamma_keV": float(result["width_stage2_keV"]),
                    }
                )
    return rows


def calculate_current() -> list[dict[str, float]]:
    """Our final pure-A result in the physical transition windows."""
    rows: list[dict[str, float]] = []
    set_photon_shape_params(-2.1, 3.8)
    convolutions = precompute_convolutions()
    inputs = ds_inputs()
    f1, _, _, _ = projected_residues(90.0)
    for s0 in (8.5, 9.0, 9.5):
        for m2 in np.linspace(3.0, 4.5, 61):
            invariants = transition_invariants(
                float(m2), float(s0), inputs, convolutions=convolutions
            )
            result = physical_couplings(
                invariants,
                theta_deg=90.0,
                m_state_1=2.4595,
                m_state_2=2.53511,
                f_1=f1,
                f_2=1.0,
                m_p=inputs["m_ds"],
                f_p=inputs["f_ds"],
                m_q=inputs["mc"],
                m_s=inputs["ms"],
            )
            rows.append(
                {
                    "theta_deg": 90.0,
                    "scenario": "current_final_ope_physical_window",
                    "M2_GeV2": float(m2),
                    "sqrt_s0_GeV": float(np.sqrt(s0)),
                    "s0_GeV2": float(s0),
                    "g1_GeV_inv": float(result["g_1"]),
                    "Gamma_keV": float(result["Gamma_1_keV"]),
                }
            )
    return rows


def write_csv(rows: list[dict[str, float]]) -> Path:
    path = HERE / "colangelo_pure_axial_reconstruction.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def build_figure(
    reference_image: Path,
    benchmark_rows: list[dict[str, float]],
    current_rows: list[dict[str, float]],
) -> tuple[Path, Path]:
    source = mpimg.imread(reference_image)
    # Remove the journal header and printed caption while retaining the axes.
    published_panel = source[25:270, 45:440]

    fig, (ax_ref, ax_match, ax_current) = plt.subplots(
        1, 3, figsize=(14.2, 4.35), constrained_layout=True
    )
    ax_ref.imshow(published_panel)
    ax_ref.axis("off")
    ax_ref.set_title("(a) Colangelo et al., Fig. 5 (published)", fontsize=11)

    for root_s0, style in zip(THRESHOLDS, STYLES):
        subset = [
            row for row in benchmark_rows if row["sqrt_s0_GeV"] == root_s0
        ]
        ax_match.plot(
            [row["M2_GeV2"] for row in subset],
            [row["g1_GeV_inv"] for row in subset],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=style["label"],
        )
    ax_match.set_xlim(2.0, 7.0)
    ax_match.set_ylim(-0.47, -0.15)
    ax_match.set_xlabel(r"$M^2\;[\mathrm{GeV}^2]$")
    ax_match.set_ylabel(r"$g_1\;[\mathrm{GeV}^{-1}]$")
    ax_match.set_title(r"(b) Matched reconstruction, $\theta=90^\circ$", fontsize=11)
    ax_match.tick_params(direction="in", top=True, right=True)
    ax_match.legend(frameon=False, fontsize=7.5, loc="lower left")
    ax_match.text(
        0.97,
        0.95,
        "Colangelo prescription\nand paper-like inputs",
        transform=ax_match.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color="0.25",
    )

    for s0, style in zip((8.5, 9.0, 9.5), STYLES):
        subset = [row for row in current_rows if row["s0_GeV2"] == s0]
        ax_current.plot(
            [row["M2_GeV2"] for row in subset],
            [row["g1_GeV_inv"] for row in subset],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=rf"$s_0={s0:g}\;\mathrm{{GeV}}^2$",
        )
    ax_current.set_xlim(3.0, 4.5)
    ax_current.set_ylim(-0.47, -0.15)
    ax_current.set_xlabel(r"$M^2\;[\mathrm{GeV}^2]$")
    ax_current.set_ylabel(r"$g_1\;[\mathrm{GeV}^{-1}]$")
    ax_current.set_title(r"(c) Our final setup, $\theta=90^\circ$", fontsize=11)
    ax_current.tick_params(direction="in", top=True, right=True)
    ax_current.legend(frameon=False, fontsize=7.5, loc="lower left")
    ax_current.text(
        0.97,
        0.95,
        "Final OPE and\nphysical windows",
        transform=ax_current.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color="0.25",
    )

    pdf = HERE / "internal_colangelo_theta90_side_by_side.pdf"
    png = HERE / "internal_colangelo_theta90_side_by_side.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=240)
    plt.close(fig)
    return pdf, png


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "reference_image",
        type=Path,
        help="User-supplied screenshot containing Colangelo et al. Fig. 5.",
    )
    args = parser.parse_args()
    benchmark_rows = calculate()
    current_rows = calculate_current()
    csv_path = write_csv(benchmark_rows + current_rows)
    pdf, png = build_figure(args.reference_image, benchmark_rows, current_rows)
    print(f"Wrote {csv_path}")
    print(f"Wrote {pdf}")
    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
