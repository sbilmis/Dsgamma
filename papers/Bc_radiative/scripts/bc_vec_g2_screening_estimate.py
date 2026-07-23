#!/usr/bin/env python3
"""Conservative G2 size screen for Bc1 -> Bc* gamma.

This is not the explicit radiative three-point G2 calculation.  It estimates
how large a missing power correction could be if its relative size is similar
to the completed heavy-heavy two-point G2/Pert ratios in the Bc-mixing study.

For Bc* gamma the coupling denominator contains both the axial Bc1 residue and
the final-state vector Bc* decay constant.  We therefore use

    |delta g/g| <= |delta Pi_3pt/Pi_3pt|
                   + |delta f_i/f_i|
                   + |delta f_Bc*/f_Bc*|

and estimate each residue shift as one half of the moment shift.  Taking the
radiative three-point shift to be bounded by the same maximum two-point G2
ratio gives a conservative envelope

    |delta g/g| <= 2 * max(|Pi_G2/Pi_pert|),
    |delta Gamma/Gamma| <= 4 * max(|Pi_G2/Pi_pert|).
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
MIXING_CONV = Path("/Users/sbilmis/Bc_mixing/BcMixingCoordinateOPEConvergence_WangWindow.csv")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def fnum(row: dict[str, str], key: str) -> float:
    return float(row[key])


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summary(vals: list[float]) -> dict[str, float]:
    vals = sorted(vals)
    n = len(vals)
    return {
        "min": vals[0],
        "median": vals[n // 2],
        "max": vals[-1],
    }


def main() -> None:
    grid = read_csv(OUT / "bc_vec_complete_grid.csv")
    grid = [r for r in grid if r["normalization"] in {"standard_vector_invariant", "reference_f_Bcstar_0p384"}]
    mixing = read_csv(MIXING_CONV)
    mixing_by_point = {(fnum(r, "M2"), fnum(r, "s0")): r for r in mixing}

    rows = []
    for row in grid:
        m2, s0 = fnum(row, "M2"), fnum(row, "s0")
        mix = mixing_by_point[(m2, s0)]
        max_two_point = fnum(mix, "MaxAbsG2OverPert")
        radiative_shift = max_two_point
        axial_residue_shift = 0.5 * max_two_point
        vector_residue_shift = 0.5 * max_two_point
        amp_envelope = radiative_shift + axial_residue_shift + vector_residue_shift
        width_envelope = 2.0 * amp_envelope
        for label, gkey, wkey in [
            ("Bc1(6743)->Bc* gamma", "g1_GeV_inv", "Gamma1_keV"),
            ("Bc1(6750)->Bc* gamma", "g2_GeV_inv", "Gamma2_keV"),
        ]:
            g = fnum(row, gkey)
            width = fnum(row, wkey)
            rows.append(
                {
                    "normalization": row["normalization"],
                    "M2": m2,
                    "s0": s0,
                    "channel": label,
                    "g_pert_GeV_inv": g,
                    "Gamma_pert_keV": width,
                    "two_point_max_abs_G2_over_pert": max_two_point,
                    "estimated_radiative_shift_abs": radiative_shift,
                    "estimated_axial_residue_shift_abs": axial_residue_shift,
                    "estimated_vector_residue_shift_abs": vector_residue_shift,
                    "conservative_abs_delta_g_over_g": amp_envelope,
                    "conservative_abs_delta_Gamma_over_Gamma": width_envelope,
                    "Gamma_shift_envelope_keV": width * width_envelope,
                }
            )

    write_csv(OUT / "bc_vec_g2_screening_estimate.csv", rows)

    lines = [
        "Bc1 -> Bc* gamma G2 screening estimate",
        "======================================",
        "",
        "Purpose:",
        "- This is a conservative size screen, not the final explicit radiative G2 calculation.",
        "- It uses the completed Bc-mixing two-point G2/Pert ratios in the same",
        "  M2={7,8,9} GeV^2 and s0={53,54,55} GeV^2 window.",
        "- The separate contact-support audit finds no ordinary two-channel",
        "  contact term in the vector channel.",
        "",
        "Rule used:",
        "- max two-point moment shift = max(|Pi_G2/Pi_pert|) over AA, AB, BB.",
        "- radiative shift proxy = max two-point moment shift.",
        "- axial residue shift proxy = 0.5 * max two-point moment shift.",
        "- Bc* residue shift proxy = 0.5 * max two-point moment shift.",
        "- conservative |delta g/g| envelope = 2 * max two-point moment shift.",
        "- conservative |delta Gamma/Gamma| envelope = 4 * max two-point moment shift.",
        "",
    ]

    for norm in sorted({r["normalization"] for r in rows}):
        sub_norm = [r for r in rows if r["normalization"] == norm]
        lines.append(f"Normalization: {norm}")
        for channel in sorted({r["channel"] for r in sub_norm}):
            sub = [r for r in sub_norm if r["channel"] == channel]
            amp = summary([r["conservative_abs_delta_g_over_g"] for r in sub])
            width_rel = summary([r["conservative_abs_delta_Gamma_over_Gamma"] for r in sub])
            width_abs = summary([r["Gamma_shift_envelope_keV"] for r in sub])
            lines.append(
                f"- {channel}: |delta g/g| median={amp['median']:.3f}, max={amp['max']:.3f}; "
                f"|delta Gamma/Gamma| median={width_rel['median']:.3f}, max={width_rel['max']:.3f}; "
                f"Gamma envelope median={width_abs['median']:.2f} keV, max={width_abs['max']:.2f} keV."
            )
        lines.append("")

    lines.extend(
        [
            "Interpretation:",
            "- The vector channel is not protected from radiative G2 effects.",
            "- The absolute envelope is especially large for Bc1(6750)->Bc* gamma because",
            "  the perturbative leading-vector width is already large.",
            "- Until the explicit all-positive vector G2 integration is completed,",
            "  the safest wording is leading perturbative vector baseline plus this",
            "  systematic screen.",
        ]
    )
    (OUT / "bc_vec_g2_screening_estimate.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
