"""Build Figure 3 amplitude-decomposition data."""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
FIGURES = HERE.parent
ROOT = HERE.parents[1]
FIGURE2 = FIGURES / "figure2_mixing_angle_widths"
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(FIGURE2))

from build_figure2_data import (
    BS_SPECS,
    DS_SPECS,
    bs_transition_inputs,
    ds_transition_inputs,
    project_residues,
    read_accepted,
    representative_sample,
)
from mixing_angle_inputs import THETA_HQET_DEG
from rohrwild_transition_exact import physical_couplings, precompute_convolutions, transition_invariants


OUT = HERE / "figure3_amplitude_decomposition.csv"


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ds_sample = representative_sample(read_accepted(ROOT / "outputs" / "twopoint_ds1_matrix_mc.csv"))
    bs_sample = representative_sample(read_accepted(ROOT / "outputs" / "twopoint_bs1_matrix_mc.csv"))
    convolutions = precompute_convolutions()
    rows: list[dict[str, float | str]] = []

    theta = THETA_HQET_DEG
    ds_vals = ds_transition_inputs(ds_sample)
    f1_ds, f2_ds = project_residues(ds_sample, theta, "Ds")
    for spec in DS_SPECS:
        invariants = transition_invariants(float(spec["M2"]), float(spec["s0"]), ds_vals, convolutions=convolutions)
        physical = physical_couplings(
            invariants,
            theta_deg=theta,
            m_state_1=ds_vals["m_ds1"],
            m_state_2=ds_vals["m_ds1_2536"],
            f_1=f1_ds,
            f_2=f2_ds,
            m_p=ds_vals["m_ds"],
            f_p=ds_vals["f_ds"],
            m_q=ds_vals["mc"],
            m_s=ds_vals["ms"],
        )
        is_low = spec["combo"] == "low"
        a_piece = physical["g_1_A"] if is_low else physical["g_2_A"]
        b_piece = physical["g_1_B"] if is_low else physical["g_2_B"]
        total = physical["g_1"] if is_low else physical["g_2"]
        rows.append(
            {
                **spec,
                "theta_deg": theta,
                "sqrt_s0": math.sqrt(float(spec["s0"])),
                "A_component_GeV_inv": a_piece,
                "B_component_GeV_inv": b_piece,
                "G_total_GeV_inv": total,
                "Gamma_keV": physical["Gamma_1_keV"] if is_low else physical["Gamma_2_keV"],
                "interference": "constructive" if a_piece * b_piece > 0.0 else "destructive",
            }
        )

    for spec in BS_SPECS:
        vals = bs_transition_inputs(bs_sample, spec)
        f1_bs, f2_bs = project_residues(bs_sample, theta, "Bs")
        invariants = transition_invariants(float(spec["M2"]), float(spec["s0"]), vals, convolutions=convolutions)
        physical = physical_couplings(
            invariants,
            theta_deg=theta,
            m_state_1=float(bs_sample["mass_low_GeV"]),
            m_state_2=float(bs_sample["mass_high_GeV"]),
            f_1=f1_bs,
            f_2=f2_bs,
            m_p=vals["m_ds"],
            f_p=vals["f_ds"],
            m_q=vals["mc"],
            m_s=vals["ms"],
        )
        is_low = spec["combo"] == "low"
        a_piece = physical["g_1_A"] if is_low else physical["g_2_A"]
        b_piece = physical["g_1_B"] if is_low else physical["g_2_B"]
        total = physical["g_1"] if is_low else physical["g_2"]
        rows.append(
            {
                **spec,
                "theta_deg": theta,
                "sqrt_s0": math.sqrt(float(spec["s0"])),
                "A_component_GeV_inv": a_piece,
                "B_component_GeV_inv": b_piece,
                "G_total_GeV_inv": total,
                "Gamma_keV": physical["Gamma_1_keV"] if is_low else physical["Gamma_2_keV"],
                "interference": "constructive" if a_piece * b_piece > 0.0 else "destructive",
            }
        )

    write_csv(OUT, rows)
    for row in rows:
        print(
            "{state}: A={A:+.4f}, B={B:+.4f}, G={G:+.4f} GeV^-1, Gamma={Gamma:.4g} keV ({interference})".format(
                state=row["state"],
                A=float(row["A_component_GeV_inv"]),
                B=float(row["B_component_GeV_inv"]),
                G=float(row["G_total_GeV_inv"]),
                Gamma=float(row["Gamma_keV"]),
                interference=row["interference"],
            )
        )


if __name__ == "__main__":
    main()
