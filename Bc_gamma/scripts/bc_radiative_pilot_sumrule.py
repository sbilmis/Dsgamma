#!/usr/bin/env python3
"""Pilot hard-photon sum-rule numerics for Bc1 -> Bc^(*) gamma.

This script converts the compact triangle-core projections into a first
numerical estimate.  The Bc gamma channel uses the same calibrated spectral
kernel as the Ds hard-photon calculation.  The Bcstar gamma channel is marked
as a diagonal pilot because its full double Borel kernel still needs an
independent audit.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

ALPHA_EM = 1.0 / 137.035999084
GEV_TO_KEV = 1.0e6


def gauss_legendre_integral(func, a: float, b: float, n: int = 360) -> float:
    if b <= a:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (b - a) * x + 0.5 * (b + a)
    return float(0.5 * (b - a) * np.sum(w * func(t)))


def kallen(s, m1, m2):
    return (s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)


def line_kernel(s, mi, mj, a, b):
    """No-charge hard-photon spectral kernel for -4 i [a + b/(p.q)].

    mi is the photon-emitting quark mass and mj is the spectator mass.  a and b
    may be arrays/functions of s through numpy broadcasting.
    """

    s = np.asarray(s, dtype=float)
    lam = np.sqrt(np.maximum(kallen(s, mj, mi), 0.0))
    num = np.clip(s - mj**2 + mi**2 - lam, 1e-300, None)
    den = np.clip(s - mj**2 + mi**2 + lam, 1e-300, None)
    log_term = np.log(num / den)
    return -3.0 / (8.0 * math.pi**2) * (
        2.0 * a * log_term
        + (b / mi**2) * (mj**2 - mi**2 - s) * lam / s**2
    )


@dataclass(frozen=True)
class Inputs:
    mb: float = 4.18
    mc: float = 1.27
    e_c: float = 2.0 / 3.0
    e_bbar: float = 1.0 / 3.0
    theta_deg: float = 43.30
    m_bc: float = 6.2749
    m_bcstar: float = 6.2749 + 0.061
    m_low: float = 6.743
    m_high: float = 6.751
    f_bc: float = 0.371
    f_bcstar: float = 0.384
    f_axial: float = 0.373


def rho_ps_A(s, inp: Inputs):
    mb, mc = inp.mb, inp.mc
    # c-line core: -4i [mc + mc^2(mc-mb)/(p.q)]
    c = line_kernel(s, mc, mb, mc, mc**2 * (mc - mb))
    # anti-b-line core: -4i [mb + mb^2(mc-mb)/(p.q)]
    bbar = line_kernel(s, mb, mc, mb, mb**2 * (mc - mb))
    return inp.e_bbar * bbar - inp.e_c * c


def rho_ps_B(s, inp: Inputs):
    mb, mc = inp.mb, inp.mc
    den = mb + mc
    c_a = (2.0 * mc**2 - mb * mc) / den
    c_b = mc**2 * s / den
    b_a = mb * mc / den
    b_b = mb**2 * s / den
    c = line_kernel(s, mc, mb, c_a, c_b)
    bbar = line_kernel(s, mb, mc, b_a, b_b)
    return inp.e_bbar * bbar - inp.e_c * c


def rho_vec_A_diag(s, inp: Inputs):
    """Diagonal pilot for the epsilon structure in Bc1 -> Bc* gamma."""

    mb, mc = inp.mb, inp.mc
    # Terms without 1/(p.q) are folded into a(s); terms with 1/(p.q) into b(s).
    c_a = (-4.0 * mb * mc + 2.0 * mc**2) / s
    c_b = 4.0 * mc**2
    b_a = (-2.0 * mb**2 + 4.0 * mb * mc) / s
    b_b = 4.0 * mb**2
    c = line_kernel(s, mc, mb, c_a, c_b)
    bbar = line_kernel(s, mb, mc, b_a, b_b)
    return inp.e_bbar * bbar - inp.e_c * c


def rho_vec_B_diag(s, inp: Inputs, pq_phys: float):
    """Diagonal pilot for J_B in Bc1 -> Bc* gamma.

    The pq/p2 term is evaluated with the physical p.q of the channel.  This is
    not the final double-Borel treatment.
    """

    mb, mc = inp.mb, inp.mc
    den = mb + mc
    c_a = (
        2.0 * mc / den
        + (2.0 * mb**2 * mc - 4.0 * mb * mc**2 + 2.0 * mc**3) / (den * s)
        + 2.0 * mc * pq_phys / (den * s)
    )
    c_b = (-4.0 * mb * mc**2 + 4.0 * mc**3) / den
    b_a = (
        2.0 * mb / den
        + (-6.0 * mb**3 + 4.0 * mb**2 * mc + 2.0 * mb * mc**2) / (den * s)
        + 2.0 * mb * pq_phys / (den * s)
    )
    b_b = (-4.0 * mb**3 + 4.0 * mb**2 * mc) / den
    c = line_kernel(s, mc, mb, c_a, c_b)
    bbar = line_kernel(s, mb, mc, b_a, b_b)
    return inp.e_bbar * bbar - inp.e_c * c


def hard_integral(M2: float, s0: float, rho_func, inp: Inputs) -> float:
    lower = (inp.mb + inp.mc) ** 2
    return gauss_legendre_integral(lambda x: np.exp(-x / M2) * rho_func(x), lower, s0)


def pref_ps(m_initial: float, M2: float, inp: Inputs) -> float:
    return (
        math.exp((m_initial**2 + inp.m_bc**2) / (2.0 * M2))
        * (inp.mb + inp.mc)
        / (inp.m_bc**2 * inp.f_bc * m_initial * inp.f_axial)
    )


def pref_vec(m_initial: float, M2: float, inp: Inputs) -> float:
    return (
        math.exp((m_initial**2 + inp.m_bcstar**2) / (2.0 * M2))
        / (inp.m_bcstar * inp.f_bcstar * m_initial * inp.f_axial)
    )


def photon_energy(mi: float, mf: float) -> float:
    return (mi * mi - mf * mf) / (2.0 * mi)


def width_ps_keV(g: float, mi: float, mf: float) -> float:
    eg = photon_energy(mi, mf)
    return ALPHA_EM / 3.0 * g * g * eg**3 * GEV_TO_KEV


def levi_civita4(a, b, c, d):
    mat = np.array([a, b, c, d], dtype=float)
    return float(np.linalg.det(mat))


def width_vec_keV(g: float, mi: float, mf: float) -> float:
    """Width for M=e g eta_A^mu eps_gamma^nu eta_V^rho V_mnr."""

    eg = photon_energy(mi, mf)
    p_vec = np.array([math.sqrt(mf * mf + eg * eg), 0.0, 0.0, -eg])
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
    # Massive vector moving along -z.
    pol_v = [
        np.array([0.0, 1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0, 0.0]),
        np.array([eg / mf, 0.0, 0.0, -p_vec[0] / mf]),
    ]

    def amp(ea, ep, ev):
        total = 0.0
        pdotq = mi * eg
        for mu in range(4):
            for nu in range(4):
                for rho in range(4):
                    basis_mu_rho = np.eye(4)[mu], np.eye(4)[rho]
                    eps_mrpq = levi_civita4(basis_mu_rho[0], basis_mu_rho[1], p_vec, q_vec)
                    eps_mnrp = levi_civita4(np.eye(4)[mu], np.eye(4)[nu], np.eye(4)[rho], p_vec)
                    v = p_vec[nu] * eps_mrpq - pdotq * eps_mnrp
                    total += ea[mu] * ep[nu] * ev[rho] * v
        return total

    summed = 0.0
    for ea in pol_a:
        for ep in pol_g:
            for ev in pol_v:
                summed += amp(ea, ep, ev) ** 2
    avg = summed / 3.0
    return eg / (8.0 * math.pi * mi**2) * (4.0 * math.pi * ALPHA_EM) * g * g * avg * GEV_TO_KEV


def compute_point(M2: float, s0: float, inp: Inputs) -> list[dict[str, float | str]]:
    th = math.radians(inp.theta_deg)
    s, c = math.sin(th), math.cos(th)

    int_ps_A = hard_integral(M2, s0, lambda x: rho_ps_A(x, inp), inp)
    int_ps_B = hard_integral(M2, s0, lambda x: rho_ps_B(x, inp), inp)

    rows = []
    for label, mi, mix_A, mix_B in [
        ("Bc1(6743)", inp.m_low, s, c),
        ("Bc1(6750)", inp.m_high, c, -s),
    ]:
        g = pref_ps(mi, M2, inp) * (mix_A * int_ps_A + mix_B * int_ps_B)
        rows.append(
            {
                "channel": f"{label} -> Bc gamma",
                "status": "hard_QCDSR_baseline",
                "M2": M2,
                "s0": s0,
                "g_GeV_inv": g,
                "Gamma_keV": width_ps_keV(g, mi, inp.m_bc),
            }
        )

    for label, mi, mix_A, mix_B in [
        ("Bc1(6743)", inp.m_low, s, c),
        ("Bc1(6750)", inp.m_high, c, -s),
    ]:
        pq_phys = 0.5 * (mi * mi - inp.m_bcstar * inp.m_bcstar)
        int_v_A = hard_integral(M2, s0, lambda x: rho_vec_A_diag(x, inp), inp)
        int_v_B = hard_integral(M2, s0, lambda x: rho_vec_B_diag(x, inp, pq_phys), inp)
        g = pref_vec(mi, M2, inp) * (mix_A * int_v_A + mix_B * int_v_B)
        rows.append(
            {
                "channel": f"{label} -> Bc* gamma",
                "status": "diagonal_pilot_needs_double_Borel_audit",
                "M2": M2,
                "s0": s0,
                "g_GeV_inv": g,
                "Gamma_keV": width_vec_keV(g, mi, inp.m_bcstar),
            }
        )
    return rows


def summarize(rows: list[dict[str, float | str]]) -> list[dict[str, float | str]]:
    out = []
    channels = sorted(set(str(r["channel"]) for r in rows))
    for ch in channels:
        subset = [r for r in rows if r["channel"] == ch]
        gvals = np.array([float(r["g_GeV_inv"]) for r in subset])
        wvals = np.array([float(r["Gamma_keV"]) for r in subset])
        out.append(
            {
                "channel": ch,
                "status": subset[0]["status"],
                "g_central": float(np.median(gvals)),
                "g_min": float(np.min(gvals)),
                "g_max": float(np.max(gvals)),
                "Gamma_central_keV": float(np.median(wvals)),
                "Gamma_min_keV": float(np.min(wvals)),
                "Gamma_max_keV": float(np.max(wvals)),
            }
        )
    return out


def main() -> None:
    inp = Inputs()
    rows = []
    for M2 in [7.0, 8.0, 9.0]:
        for s0 in [53.0, 54.0, 55.0]:
            rows.extend(compute_point(M2, s0, inp))

    grid_path = OUT / "bc_radiative_pilot_grid.csv"
    with grid_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = summarize(rows)
    summary_path = OUT / "bc_radiative_pilot_summary.csv"
    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)

    txt_path = OUT / "bc_radiative_pilot_summary.txt"
    lines = [
        "Bc radiative pilot hard-photon sum-rule summary",
        "================================================",
        "",
        "Inputs: mb=4.18 GeV, mc=1.27 GeV, theta=43.3 deg,",
        "M2 in {7,8,9} GeV^2 and s0 in {53,54,55} GeV^2.",
        "Decay constants: f_Bc=0.371 GeV, f_Bc*=0.384 GeV, f_Bc1=0.373 GeV.",
        "Bc* gamma channels are diagonal-pilot estimates pending the full double-Borel audit.",
        "",
    ]
    for r in summary:
        lines.append(
            "{channel}: g={g:+.4g} [{gl:+.4g},{gh:+.4g}] GeV^-1, "
            "Gamma={w:.3g} [{wl:.3g},{wh:.3g}] keV ({status})".format(
                channel=r["channel"],
                g=float(r["g_central"]),
                gl=float(r["g_min"]),
                gh=float(r["g_max"]),
                w=float(r["Gamma_central_keV"]),
                wl=float(r["Gamma_min_keV"]),
                wh=float(r["Gamma_max_keV"]),
                status=r["status"],
            )
        )
    txt_path.write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"Wrote {grid_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {txt_path}")


if __name__ == "__main__":
    main()

