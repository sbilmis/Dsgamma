#!/usr/bin/env python3
"""Audit the Bc1 -> Bc* gamma vector projector and width normalization.

The vector-channel projection uses

    V_mnr = p_n eps_{m r p q} - (p.q) eps_{m n r p}

and extracts a coefficient g_V from the three-point sum rule.  This script
checks the polarization sum used to convert that coefficient into a width and
relates g_V to a one-momentum E1-style coupling.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

ALPHA_EM = 1.0 / 137.035999084
GEV_TO_KEV = 1.0e6


def photon_energy(mi: float, mf: float) -> float:
    return (mi * mi - mf * mf) / (2.0 * mi)


def levi_civita4(a, b, c, d):
    return float(np.linalg.det(np.array([a, b, c, d], dtype=float)))


def vector_tensor_polsum(mi: float, mf: float) -> dict[str, float]:
    """Return the brute-force polarization average for the V_mnr basis."""

    eg = photon_energy(mi, mf)
    ef = math.sqrt(mf * mf + eg * eg)
    p_vec = np.array([ef, 0.0, 0.0, -eg])
    q_vec = np.array([eg, 0.0, 0.0, eg])

    pol_a = [
        np.array([0.0, 1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 0.0, 1.0]),
    ]
    pol_g = [
        np.array([0.0, 1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0, 0.0]),
    ]
    pol_v = [
        np.array([0.0, 1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0, 0.0]),
        np.array([eg / mf, 0.0, 0.0, -ef / mf]),
    ]

    def amp(ea, ep, ev):
        total = 0.0
        pdotq = mi * eg
        basis = np.eye(4)
        for mu in range(4):
            for nu in range(4):
                for rho in range(4):
                    eps_mrpq = levi_civita4(basis[mu], basis[rho], p_vec, q_vec)
                    eps_mnrp = levi_civita4(basis[mu], basis[nu], basis[rho], p_vec)
                    v_mnr = p_vec[nu] * eps_mrpq - pdotq * eps_mnrp
                    total += ea[mu] * ep[nu] * ev[rho] * v_mnr
        return total

    summed = 0.0
    for ea in pol_a:
        for ep in pol_g:
            for ev in pol_v:
                summed += amp(ea, ep, ev) ** 2
    avg = summed / 3.0
    analytic = (2.0 / 3.0) * mi * mi * eg * eg * (2.0 * mf * mf + eg * eg)
    width_coeff = (
        eg / (8.0 * math.pi * mi * mi)
        * (4.0 * math.pi * ALPHA_EM)
        * avg
        * GEV_TO_KEV
    )
    e1_coeff = (ALPHA_EM / 3.0) * eg**3 * GEV_TO_KEV
    return {
        "mi_GeV": mi,
        "mf_GeV": mf,
        "Egamma_GeV": eg,
        "polsum_avg_numeric_GeV6": avg,
        "polsum_avg_analytic_GeV6": analytic,
        "numeric_over_analytic": avg / analytic,
        "Gamma_per_gV2_keV": width_coeff,
        "Gamma_per_gE12_keV": e1_coeff,
        "Gamma_coeff_ratio_gV_to_gE1": width_coeff / e1_coeff,
        "gE1_over_gV": math.sqrt(width_coeff / e1_coeff),
        "sqrt_2mf2_plus_E2": math.sqrt(2.0 * mf * mf + eg * eg),
    }


def main() -> None:
    rows = [
        {"channel": "Bc1(6743)->Bc* gamma", **vector_tensor_polsum(6.743, 6.3359)},
        {"channel": "Bc1(6750)->Bc* gamma", **vector_tensor_polsum(6.751, 6.3359)},
    ]
    csv_path = OUT / "bc_vec_projector_width_audit.csv"
    txt_path = OUT / "bc_vec_projector_width_audit.txt"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Bc1 -> Bc* gamma projector/width audit",
        "========================================",
        "",
        "For V_mnr = p_n eps_{m r p q} - (p.q) eps_{m n r p},",
        "the explicit polarization sum gives",
        "  (1/3) sum_pol |eta_A eps_gamma eta_V V|^2",
        "  = (2/3) m_i^2 E_gamma^2 (2 m_f^2 + E_gamma^2).",
        "",
        "Therefore",
        "  Gamma = (alpha_em / 3) g_V^2 (2 m_f^2 + E_gamma^2) E_gamma^3,",
        "or equivalently Gamma = (alpha_em / 3) g_E1^2 E_gamma^3",
        "with g_E1 = sqrt(2 m_f^2 + E_gamma^2) g_V.",
        "",
    ]
    for row in rows:
        lines.append(
            f"{row['channel']}: numeric/analytic polsum = "
            f"{row['numeric_over_analytic']:.12g}, "
            f"g_E1/g_V = {row['gE1_over_gV']:.6g} "
            f"(sqrt(2 mf^2 + E^2) = {row['sqrt_2mf2_plus_E2']:.6g})"
        )
    txt_path.write_text("\n".join(lines) + "\n")
    print(txt_path.read_text())


if __name__ == "__main__":
    main()
