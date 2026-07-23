#!/usr/bin/env python3
"""Nonrelativistic E1 sanity check for Bc1 -> Bc^(*) gamma.

This is not the final QCD sum-rule result. It is a controlled benchmark that
captures the expected spin-selection pattern before the heavy-heavy hard-photon
sum-rule calculation is derived.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


ALPHA = 1.0 / 137.035999084
GEV_TO_KEV = 1.0e6


def photon_energy(mi: float, mf: float) -> float:
    return (mi * mi - mf * mf) / (2.0 * mi)


def e1_width_keV(mi: float, mf: float, e_eff: float, r_gevinv: float, spin_weight: float) -> float:
    """Simple E1 estimate in keV.

    The normalization uses Gamma = (4/9) alpha e_eff^2 E_gamma^3 |r|^2 Ef/Mi.
    The `spin_weight` encodes the assumed 1P1-3P1 mixing projection.
    """

    egamma = photon_energy(mi, mf)
    gamma_gev = (4.0 / 9.0) * ALPHA * e_eff * e_eff * egamma**3 * r_gevinv**2
    gamma_gev *= mf / mi
    return gamma_gev * spin_weight * GEV_TO_KEV


def main() -> None:
    out_dir = Path(__file__).resolve().parents[1] / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Bondar-Milstein mass pattern: M(Bc*) - M(Bc) = 61 MeV,
    # M(Bc1(6743)) - M(Bc) = 497 MeV, and M(6750)-M(6743) = 8 MeV.
    m_bc = 6.2749
    m_bc_star = m_bc + 0.061
    states = {
        "Bc1(6743)": 6.743,
        "Bc1(6750)": 6.751,
    }

    # Constituent masses only for this quark-model E1 benchmark.
    m_c_eff = 1.7
    m_b_eff = 4.8
    e_c = 2.0 / 3.0
    e_antib = 1.0 / 3.0
    e_eff = (m_b_eff * e_c - m_c_eff * e_antib) / (m_b_eff + m_c_eff)

    theta_deg = 35.3
    theta = math.radians(theta_deg)
    r_gevinv = 1.5

    rows = []
    for state, mi in states.items():
        if state == "Bc1(6743)":
            w_to_bc = math.cos(theta) ** 2
            w_to_bc_star = math.sin(theta) ** 2
        else:
            w_to_bc = math.sin(theta) ** 2
            w_to_bc_star = math.cos(theta) ** 2

        rows.append(
            {
                "channel": f"{state} -> Bc gamma",
                "M_initial_GeV": mi,
                "M_final_GeV": m_bc,
                "E_gamma_GeV": photon_energy(mi, m_bc),
                "spin_projection": w_to_bc,
                "Gamma_keV": e1_width_keV(mi, m_bc, e_eff, r_gevinv, w_to_bc),
            }
        )
        rows.append(
            {
                "channel": f"{state} -> Bc* gamma",
                "M_initial_GeV": mi,
                "M_final_GeV": m_bc_star,
                "E_gamma_GeV": photon_energy(mi, m_bc_star),
                "spin_projection": w_to_bc_star,
                "Gamma_keV": e1_width_keV(mi, m_bc_star, e_eff, r_gevinv, w_to_bc_star),
            }
        )

    csv_path = out_dir / "stage0_e1_benchmark.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    md_path = out_dir / "stage0_e1_benchmark.md"
    with md_path.open("w") as f:
        f.write("# Stage-0 E1 benchmark\n\n")
        f.write("This is a nonrelativistic sanity check, not the final QCD sum-rule result.\n\n")
        f.write(f"Inputs: theta = {theta_deg:.1f} deg, r = {r_gevinv:.2f} GeV^-1, ")
        f.write(f"e_eff = {e_eff:.4f}.\n\n")
        f.write("| Channel | E_gamma (GeV) | projection | Gamma (keV) |\n")
        f.write("|---|---:|---:|---:|\n")
        for row in rows:
            f.write(
                f"| {row['channel']} | {row['E_gamma_GeV']:.4f} | "
                f"{row['spin_projection']:.3f} | {row['Gamma_keV']:.2f} |\n"
            )

    for row in rows:
        print(f"{row['channel']}: {row['Gamma_keV']:.2f} keV")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()

