#!/usr/bin/env python3
"""Summarize perturbative Bc gamma results with the conservative G2 envelope."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def q(vals, x):
    return float(np.quantile(np.asarray(vals, dtype=float), x))


def main() -> None:
    mc_summary = read_csv(OUT / "bc_ps_complete_monte_carlo_summary.csv")
    screen = read_csv(OUT / "bc_ps_g2_screening_estimate.csv")

    rows = []
    for channel in sorted({r["channel"] for r in screen}):
        g_row = next(
            r
            for r in mc_summary
            if r["channel"] == channel and r["observable"].startswith("g")
        )
        w_row = next(
            r
            for r in mc_summary
            if r["channel"] == channel and r["observable"].startswith("Gamma")
        )
        sub = [r for r in screen if r["channel"] == channel]
        dg_env = [f(r, "conservative_abs_delta_g_over_g") for r in sub]
        dw_env = [f(r, "conservative_abs_delta_Gamma_over_Gamma") for r in sub]
        w_shift = [f(r, "Gamma_shift_envelope_keV") for r in sub]
        rows.append(
            {
                "channel": channel,
                "g_pert_median_GeV_inv": f(g_row, "median"),
                "g_pert_q16_GeV_inv": f(g_row, "q16"),
                "g_pert_q84_GeV_inv": f(g_row, "q84"),
                "Gamma_pert_median_keV": f(w_row, "median"),
                "Gamma_pert_q16_keV": f(w_row, "q16"),
                "Gamma_pert_q84_keV": f(w_row, "q84"),
                "median_abs_delta_g_over_g_G2_env": q(dg_env, 0.5),
                "max_abs_delta_g_over_g_G2_env": max(dg_env),
                "median_abs_delta_Gamma_over_Gamma_G2_env": q(dw_env, 0.5),
                "max_abs_delta_Gamma_over_Gamma_G2_env": max(dw_env),
                "median_Gamma_G2_shift_envelope_keV": q(w_shift, 0.5),
                "max_Gamma_G2_shift_envelope_keV": max(w_shift),
            }
        )

    csv_path = OUT / "bc_ps_g2_corrected_summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Bc1 -> Bc gamma perturbative result plus conservative G2 envelope",
        "================================================================",
        "",
        "This is not the final explicit G2 correction.  It is the conservative",
        "systematic envelope inferred from the same-window heavy-heavy G2/Pert",
        "ratios after including the residue-normalization effect.",
        "",
    ]
    for row in rows:
        lines.append(row["channel"])
        lines.append(
            "  g_pert = {g:.4g} [{g16:.4g}, {g84:.4g}] GeV^-1".format(
                g=row["g_pert_median_GeV_inv"],
                g16=row["g_pert_q16_GeV_inv"],
                g84=row["g_pert_q84_GeV_inv"],
            )
        )
        lines.append(
            "  Gamma_pert = {w:.4g} [{w16:.4g}, {w84:.4g}] keV".format(
                w=row["Gamma_pert_median_keV"],
                w16=row["Gamma_pert_q16_keV"],
                w84=row["Gamma_pert_q84_keV"],
            )
        )
        lines.append(
            "  conservative G2 envelope: median |dg/g|={dg:.3f}, max={dgmax:.3f}; "
            "median |dGamma/Gamma|={dw:.3f}, max={dwmax:.3f}".format(
                dg=row["median_abs_delta_g_over_g_G2_env"],
                dgmax=row["max_abs_delta_g_over_g_G2_env"],
                dw=row["median_abs_delta_Gamma_over_Gamma_G2_env"],
                dwmax=row["max_abs_delta_Gamma_over_Gamma_G2_env"],
            )
        )
        lines.append(
            "  Gamma G2 envelope: median +/-{med:.3g} keV, max +/-{mx:.3g} keV".format(
                med=row["median_Gamma_G2_shift_envelope_keV"],
                mx=row["max_Gamma_G2_shift_envelope_keV"],
            )
        )
        lines.append("")

    lines.append("Conclusion:")
    lines.append("- G2 is not proven negligible by this screen.")
    lines.append("- The quoted perturbative widths should carry this as a systematic envelope")
    lines.append("  until the full radiative G2 double-Borel reduction is completed.")
    txt_path = OUT / "bc_ps_g2_corrected_summary.txt"
    txt_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
