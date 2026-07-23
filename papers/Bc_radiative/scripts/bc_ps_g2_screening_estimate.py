#!/usr/bin/env python3
"""Conservative G2 size screen for Bc1 -> Bc gamma.

This is not a substitute for the full radiative three-point G2 reduction.
It answers a narrower question: if the radiative G2 correction is of the same
relative size as the already-computed heavy-heavy two-point G2 correction in
the Bc-mixing study, would neglecting it be harmless?
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


def main() -> None:
    grid = read_csv(OUT / "bc_ps_complete_grid.csv")
    mixing = read_csv(MIXING_CONV)
    mixing_by_point = {(fnum(r, "M2"), fnum(r, "s0")): r for r in mixing}

    rows = []
    for row in grid:
        m2, s0 = fnum(row, "M2"), fnum(row, "s0")
        mix = mixing_by_point[(m2, s0)]
        # This is the strongest two-point G2/Pert ratio in the same window.
        max_two_point = fnum(mix, "MaxAbsG2OverPert")
        # A residue f ~ sqrt(Pi); hence the corresponding normalization shift
        # is roughly half the two-point moment shift.
        residue_shift = 0.5 * max_two_point
        # If the radiative three-point G2 piece is comparable to the two-point
        # G2 ratio, the coupling shift is delta_3pt - delta_residue.  We quote
        # a conservative absolute envelope by adding magnitudes.
        amp_envelope = max_two_point + residue_shift
        width_envelope = 2.0 * amp_envelope
        for label, gkey, wkey in [
            ("Bc1(6743)->Bc gamma", "g1_GeV_inv", "Gamma1_keV"),
            ("Bc1(6750)->Bc gamma", "g2_GeV_inv", "Gamma2_keV"),
        ]:
            g = fnum(row, gkey)
            width = fnum(row, wkey)
            rows.append(
                {
                    "M2": m2,
                    "s0": s0,
                    "channel": label,
                    "g_pert_GeV_inv": g,
                    "Gamma_pert_keV": width,
                    "two_point_max_abs_G2_over_pert": max_two_point,
                    "estimated_residue_shift_abs": residue_shift,
                    "conservative_abs_delta_g_over_g": amp_envelope,
                    "conservative_abs_delta_Gamma_over_Gamma": width_envelope,
                    "Gamma_shift_envelope_keV": width * width_envelope,
                }
            )

    write_csv(OUT / "bc_ps_g2_screening_estimate.csv", rows)

    max_amp = max(r["conservative_abs_delta_g_over_g"] for r in rows)
    max_width = max(r["conservative_abs_delta_Gamma_over_Gamma"] for r in rows)
    median_amp = sorted(r["conservative_abs_delta_g_over_g"] for r in rows)[len(rows) // 2]
    median_width = sorted(r["conservative_abs_delta_Gamma_over_Gamma"] for r in rows)[len(rows) // 2]

    lines = [
        "Bc1 -> Bc gamma G2 screening estimate",
        "======================================",
        "",
        "Purpose:",
        "- This is a conservative size screen, not the final radiative G2 calculation.",
        "- It uses the completed Bc-mixing two-point G2/Pert ratios in the same",
        "  M2={7,8,9} GeV^2 and s0={53,54,55} GeV^2 window.",
        "",
        "Rule used:",
        "- max two-point moment shift = max(|Pi_G2/Pi_pert|) over AA, AB, BB.",
        "- residue shift estimate = 0.5 * max two-point moment shift.",
        "- conservative |delta g/g| envelope = radiative shift + residue shift,",
        "  with the radiative shift provisionally taken as the same max two-point",
        "  ratio.",
        "- conservative |delta Gamma/Gamma| envelope = 2 |delta g/g|.",
        "",
        f"Median conservative |delta g/g| envelope: {median_amp:.3f}",
        f"Maximum conservative |delta g/g| envelope: {max_amp:.3f}",
        f"Median conservative |delta Gamma/Gamma| envelope: {median_width:.3f}",
        f"Maximum conservative |delta Gamma/Gamma| envelope: {max_width:.3f}",
        "",
        "Interpretation:",
        "- The screen is not small enough to justify silently dropping radiative G2",
        "  at publication level.",
        "- If the full three-point G2 reduction gives much smaller numbers, we can",
        "  then quote it as a checked negligible correction.",
        "- Until that explicit reduction is done, the safer statement is that the",
        "  present widths are leading-perturbative hard-QCDSR baselines.",
    ]
    (OUT / "bc_ps_g2_screening_estimate.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
