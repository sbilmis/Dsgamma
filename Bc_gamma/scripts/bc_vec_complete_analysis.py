#!/usr/bin/env python3
"""Leading Bc1 -> Bc* gamma analysis with explicit vector normalization.

The older pilot table used the correct hard-photon vector traces but divided
by fixed input constants.  This script makes the normalization transparent:

1. f1 and f2 are recomputed from the axial-vector Bc-mixing two-point sum
   rules, exactly as in bc_ps_complete_analysis.py.
2. f_Bc* is computed from a perturbative vector-current two-point sum rule.
   We keep both the direct spin-projector normalization and the standard
   vector-invariant normalization.  The latter is the one appropriate for
   <0|cbar gamma_mu b|Bc*>=m_Bc* f_Bc* eta_mu.
3. A reference-normalization column with f_Bc*=0.384 GeV is also kept to show
   the residual sensitivity to this input.

The three-point OPE here is still the leading hard-photon perturbative
contribution.  Radiative G^2/contact terms are not included in these numbers.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".mplconfig"))
(OUT / ".mplconfig").mkdir(exist_ok=True)

ALPHA_EM = 1.0 / 137.035999084
GEV_TO_KEV = 1.0e6

_LEGENDRE_CACHE: dict[int, tuple[np.ndarray, np.ndarray]] = {}


@dataclass(frozen=True)
class Inputs:
    mb: float = 4.18
    mc: float = 1.27
    e_c: float = 2.0 / 3.0
    e_bbar: float = 1.0 / 3.0
    theta_deg: float = 43.2950875223707
    theta_sigma_deg: float = 0.1513163252823016
    m_bcstar: float = 6.3359
    m_low: float = 6.743
    m_high: float = 6.751
    f_bcstar_ref: float = 0.384
    f_bcstar_ref_sigma: float = 0.038


def gauss_legendre_integral(func, a: float, b: float, n: int = 180) -> float:
    if b <= a:
        return 0.0
    if n not in _LEGENDRE_CACHE:
        _LEGENDRE_CACHE[n] = np.polynomial.legendre.leggauss(n)
    x, w = _LEGENDRE_CACHE[n]
    t = 0.5 * (b - a) * x + 0.5 * (b + a)
    return float(0.5 * (b - a) * np.sum(w * func(t)))


def kallen_mass_sq(s, m1, m2):
    return (s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)


def line_kernel(s, mi, mj, a, b):
    """No-charge hard-photon spectral kernel for -4 i [a + b/(p.q)]."""

    s = np.asarray(s, dtype=float)
    lam = np.sqrt(np.maximum(kallen_mass_sq(s, mj, mi), 0.0))
    num = np.clip(s - mj**2 + mi**2 - lam, 1e-300, None)
    den = np.clip(s - mj**2 + mi**2 + lam, 1e-300, None)
    log_term = np.log(num / den)
    return -3.0 / (8.0 * math.pi**2) * (
        2.0 * a * log_term
        + (b / mi**2) * (mj**2 - mi**2 - s) * lam / s**2
    )


def rho_vec_A_diag(s, inp: Inputs):
    """Leading hard-photon epsilon-structure density for J_A."""

    mb, mc = inp.mb, inp.mc
    c_a = (-4.0 * mb * mc + 2.0 * mc**2) / s
    c_b = 4.0 * mc**2
    b_a = (-2.0 * mb**2 + 4.0 * mb * mc) / s
    b_b = 4.0 * mb**2
    c = line_kernel(s, mc, mb, c_a, c_b)
    bbar = line_kernel(s, mb, mc, b_a, b_b)
    return inp.e_bbar * bbar - inp.e_c * c


def rho_vec_B_diag(s, inp: Inputs, pq_phys: float):
    """Leading hard-photon epsilon-structure density for J_B.

    The explicit pq/p^2 term is evaluated at the physical p.q of the channel.
    This is the same convention used in the earlier vector pilot.
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


def two_point_kallen(s, mb, mc):
    return s**2 + mb**4 + mc**4 - 2 * s * mb**2 - 2 * s * mc**2 - 2 * mb**2 * mc**2


def two_point_num(channel: str, s, mb, mc):
    kp = 0.5 * (s + mc**2 - mb**2)
    k2 = mc**2
    if channel == "AA":
        return -4.0 / s * (2.0 * kp**2 + (3.0 * mb * mc + k2) * s - 3.0 * kp * s)
    if channel == "AB":
        return 12.0 * (-(mb + mc) * kp + mc * s)
    if channel == "BB":
        return -4.0 * (4.0 * kp**2 + (3.0 * mb * mc - k2) * s - 3.0 * kp * s)
    raise ValueError(channel)


def two_point_norm(channel: str, mb, mc):
    if channel == "AA":
        return 1.0
    if channel == "AB":
        return 1.0 / (mb + mc)
    if channel == "BB":
        return 1.0 / (mb + mc) ** 2
    raise ValueError(channel)


def rho_twopoint(channel: str, s, mb, mc):
    lam = np.sqrt(np.maximum(two_point_kallen(s, mb, mc), 0.0))
    return (
        3.0
        / (16.0 * math.pi**2)
        * lam
        / s
        * two_point_num(channel, s, mb, mc)
        * two_point_norm(channel, mb, mc)
    )


def two_point_moment(channel: str, M2: float, s0: float, inp: Inputs) -> float:
    lower = (inp.mb + inp.mc) ** 2
    return gauss_legendre_integral(
        lambda x: np.exp(-x / M2) * rho_twopoint(channel, x, inp.mb, inp.mc),
        lower,
        s0,
    )


def normalize_theta(theta: float) -> float:
    while theta < 0:
        theta += math.pi / 2.0
    while theta >= math.pi / 2.0:
        theta -= math.pi / 2.0
    return theta


def twopoint_residues(M2: float, s0: float, inp: Inputs):
    aa = two_point_moment("AA", M2, s0, inp)
    ab = two_point_moment("AB", M2, s0, inp)
    bb = two_point_moment("BB", M2, s0, inp)
    theta = normalize_theta(0.5 * math.atan2(-2.0 * ab, aa - bb))
    st, ct = math.sin(theta), math.cos(theta)
    pi1 = st * st * aa + 2.0 * st * ct * ab + ct * ct * bb
    pi2 = ct * ct * aa - 2.0 * st * ct * ab + st * st * bb
    f1 = math.sqrt(max(0.0, math.exp(inp.m_low**2 / M2) * pi1 / inp.m_low**2))
    f2 = math.sqrt(max(0.0, math.exp(inp.m_high**2 / M2) * pi2 / inp.m_high**2))
    return {
        "theta_sumrule_deg": math.degrees(theta),
        "Pi_1": pi1,
        "Pi_2": pi2,
        "f1_GeV": f1,
        "f2_GeV": f2,
    }


def rho_bcstar_vector_current(s, inp: Inputs):
    """Perturbative vector-current density for <0|cbar gamma_mu b|Bc*>.

    With the spin-1 invariant convention used here, the equal-mass limit is
    proportional to (s+2m_Q^2) sqrt(1-4m_Q^2/s), as expected for an S-wave
    vector current.
    """

    mb, mc = inp.mb, inp.mc
    s = np.asarray(s, dtype=float)
    kp = 0.5 * (s + mc**2 - mb**2)
    num = -8.0 * kp**2 / s - 4.0 * mc**2 + 12.0 * kp + 12.0 * mb * mc
    lam = np.sqrt(np.maximum(two_point_kallen(s, mb, mc), 0.0))
    return 3.0 / (16.0 * math.pi**2) * lam / s * num


def bcstar_decay_constant(M2v: float, s0v: float, inp: Inputs, convention_scale: float = 1.0) -> float:
    lower = (inp.mb + inp.mc) ** 2
    moment = gauss_legendre_integral(
        lambda x: np.exp(-x / M2v) * convention_scale * rho_bcstar_vector_current(x, inp),
        lower,
        s0v,
    )
    return math.sqrt(max(0.0, math.exp(inp.m_bcstar**2 / M2v) * moment / inp.m_bcstar**2))


def photon_energy(mi: float, mf: float) -> float:
    return (mi * mi - mf * mf) / (2.0 * mi)


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
                    eps_mrpq = levi_civita4(np.eye(4)[mu], np.eye(4)[rho], p_vec, q_vec)
                    eps_mnrp = levi_civita4(np.eye(4)[mu], np.eye(4)[nu], np.eye(4)[rho], p_vec)
                    basis = p_vec[nu] * eps_mrpq - pdotq * eps_mnrp
                    total += ea[mu] * ep[nu] * ev[rho] * basis
        return total

    summed = 0.0
    for ea in pol_a:
        for ep in pol_g:
            for ev in pol_v:
                summed += amp(ea, ep, ev) ** 2
    avg = summed / 3.0
    return eg / (8.0 * math.pi * mi**2) * (4.0 * math.pi * ALPHA_EM) * g * g * avg * GEV_TO_KEV


def vec_couplings(M2: float, s0: float, inp: Inputs, theta_deg: float, f1: float, f2: float, f_bcstar: float):
    th = math.radians(theta_deg)
    st, ct = math.sin(th), math.cos(th)
    pq1 = 0.5 * (inp.m_low**2 - inp.m_bcstar**2)
    pq2 = 0.5 * (inp.m_high**2 - inp.m_bcstar**2)
    int_A = hard_integral(M2, s0, lambda x: rho_vec_A_diag(x, inp), inp)
    int_B_1 = hard_integral(M2, s0, lambda x: rho_vec_B_diag(x, inp, pq1), inp)
    int_B_2 = hard_integral(M2, s0, lambda x: rho_vec_B_diag(x, inp, pq2), inp)
    pref1 = math.exp((inp.m_low**2 + inp.m_bcstar**2) / (2.0 * M2)) / (
        inp.m_bcstar * f_bcstar * inp.m_low * f1
    )
    pref2 = math.exp((inp.m_high**2 + inp.m_bcstar**2) / (2.0 * M2)) / (
        inp.m_bcstar * f_bcstar * inp.m_high * f2
    )
    g1 = pref1 * (st * int_A + ct * int_B_1)
    g2 = pref2 * (ct * int_A - st * int_B_2)
    return {
        "int_A": int_A,
        "int_B_1": int_B_1,
        "int_B_2": int_B_2,
        "g1_GeV_inv": g1,
        "g2_GeV_inv": g2,
        "Gamma1_keV": width_vec_keV(g1, inp.m_low, inp.m_bcstar),
        "Gamma2_keV": width_vec_keV(g2, inp.m_high, inp.m_bcstar),
    }


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summary_stats(vals):
    arr = np.asarray(vals, dtype=float)
    return {
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
        "q16": float(np.quantile(arr, 0.16)),
        "q84": float(np.quantile(arr, 0.84)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def plot_stability(rows: list[dict], norm_label: str, key: str, ylabel: str, path: Path, title: str):
    subset = [r for r in rows if r["normalization"] == norm_label]
    plt.figure(figsize=(5.0, 3.6))
    for s0 in sorted({r["s0"] for r in subset}):
        sub = [r for r in subset if r["s0"] == s0]
        sub.sort(key=lambda r: r["M2"])
        plt.plot([r["M2"] for r in sub], [abs(r[key]) for r in sub], marker="o", label=fr"$s_0={s0:.0f}$")
    plt.xlabel(r"$M^2\,[\mathrm{GeV}^2]$")
    plt.ylabel(ylabel)
    plt.text(0.04, 0.92, title, transform=plt.gca().transAxes)
    plt.grid(alpha=0.25)
    plt.legend(frameon=False, fontsize=9)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def plot_fbcstar_stability(rows: list[dict], key: str, ylabel: str, path: Path):
    plt.figure(figsize=(5.0, 3.6))
    for s0v in sorted({r["s0_vector"] for r in rows}):
        sub = [r for r in rows if r["s0_vector"] == s0v]
        sub.sort(key=lambda r: r["M2_vector"])
        plt.plot(
            [r["M2_vector"] for r in sub],
            [r[key] for r in sub],
            marker="o",
            label=fr"$s_{{0,V}}={s0v:.0f}$",
        )
    plt.xlabel(r"$M_V^2\,[\mathrm{GeV}^2]$")
    plt.ylabel(ylabel)
    plt.text(0.04, 0.92, r"$B_c^\ast$ two-point", transform=plt.gca().transAxes)
    plt.grid(alpha=0.25)
    plt.legend(frameon=False, fontsize=9)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def run_grid(inp: Inputs):
    rows = []
    fstar_rows = []
    for M2v in [5.0, 6.0, 7.0, 8.0]:
        for s0v in [40.0, 42.0, 44.0]:
            fstar_rows.append(
                {
                    "M2_vector": M2v,
                    "s0_vector": s0v,
                    "f_Bcstar_projector_GeV": bcstar_decay_constant(M2v, s0v, inp, 1.0),
                    "f_Bcstar_standard_GeV": bcstar_decay_constant(M2v, s0v, inp, 1.0 / 3.0),
                }
            )
    fstar_projector_central = float(np.median([r["f_Bcstar_projector_GeV"] for r in fstar_rows]))
    fstar_standard_central = float(np.median([r["f_Bcstar_standard_GeV"] for r in fstar_rows]))

    for M2 in [7.0, 8.0, 9.0]:
        for s0 in [53.0, 54.0, 55.0]:
            res = twopoint_residues(M2, s0, inp)
            for norm_label, fstar in [
                ("projector_diagnostic", fstar_projector_central),
                ("standard_vector_invariant", fstar_standard_central),
                ("reference_f_Bcstar_0p384", inp.f_bcstar_ref),
            ]:
                rad = vec_couplings(M2, s0, inp, inp.theta_deg, res["f1_GeV"], res["f2_GeV"], fstar)
                rows.append(
                    {
                        "normalization": norm_label,
                        "M2": M2,
                        "s0": s0,
                        "f_Bcstar_used_GeV": fstar,
                        **res,
                        **rad,
                    }
                )
    return fstar_rows, rows


def run_monte_carlo(n: int = 1000, seed: int = 20260712):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        mb = rng.uniform(4.16, 4.21)
        mc = rng.uniform(1.25, 1.29)
        M2 = rng.uniform(7.0, 9.0)
        s0 = rng.uniform(53.0, 55.0)
        M2v = rng.uniform(5.0, 8.0)
        s0v = rng.uniform(40.0, 44.0)
        theta = rng.normal(43.2950875223707, 0.1513163252823016)
        inp = Inputs(mb=mb, mc=mc, theta_deg=theta)
        if s0 <= (mb + mc) ** 2 or s0v <= (mb + mc) ** 2:
            continue
        res = twopoint_residues(M2, s0, inp)
        fstar_lo = bcstar_decay_constant(M2v, s0v, inp)
        fstar_standard = bcstar_decay_constant(M2v, s0v, inp, 1.0 / 3.0)
        fstar_ref = max(0.1, rng.normal(inp.f_bcstar_ref, inp.f_bcstar_ref_sigma))
        for norm_label, fstar in [
            ("projector_diagnostic", fstar_lo),
            ("standard_vector_invariant", fstar_standard),
            ("reference_f_Bcstar_0p384", fstar_ref),
        ]:
            rad = vec_couplings(M2, s0, inp, theta, res["f1_GeV"], res["f2_GeV"], fstar)
            rows.append(
                {
                    "index": len(rows) + 1,
                    "normalization": norm_label,
                    "mb": mb,
                    "mc": mc,
                    "M2": M2,
                    "s0": s0,
                    "M2_vector": M2v,
                    "s0_vector": s0v,
                    "theta_deg": theta,
                    "f_Bcstar_used_GeV": fstar,
                    **res,
                    **rad,
                }
            )
    return rows


def main():
    inp = Inputs()
    fstar_rows, grid_rows = run_grid(inp)
    write_csv(OUT / "bc_vec_fbcstar_twopoint_grid.csv", fstar_rows)
    write_csv(OUT / "bc_vec_complete_grid.csv", grid_rows)

    fstar_summary = [
        {
            "quantity": "f_Bcstar_projector_GeV",
            **summary_stats([r["f_Bcstar_projector_GeV"] for r in fstar_rows]),
        },
        {
            "quantity": "f_Bcstar_standard_GeV",
            **summary_stats([r["f_Bcstar_standard_GeV"] for r in fstar_rows]),
        },
    ]
    write_csv(OUT / "bc_vec_fbcstar_twopoint_summary.csv", fstar_summary)

    grid_summary = []
    for norm_label in sorted({r["normalization"] for r in grid_rows}):
        subset = [r for r in grid_rows if r["normalization"] == norm_label]
        for channel, gkey, wkey in [
            ("Bc1(6743)->Bc* gamma", "g1_GeV_inv", "Gamma1_keV"),
            ("Bc1(6750)->Bc* gamma", "g2_GeV_inv", "Gamma2_keV"),
        ]:
            grid_summary.append(
                {
                    "normalization": norm_label,
                    "channel": channel,
                    "observable": gkey,
                    **summary_stats([r[gkey] for r in subset]),
                }
            )
            grid_summary.append(
                {
                    "normalization": norm_label,
                    "channel": channel,
                    "observable": wkey,
                    **summary_stats([r[wkey] for r in subset]),
                }
            )
    write_csv(OUT / "bc_vec_complete_grid_summary.csv", grid_summary)

    samples = run_monte_carlo()
    write_csv(OUT / "bc_vec_complete_monte_carlo.csv", samples)

    mc_summary = []
    for norm_label in sorted({r["normalization"] for r in samples}):
        subset = [r for r in samples if r["normalization"] == norm_label]
        for channel, gkey, wkey in [
            ("Bc1(6743)->Bc* gamma", "g1_GeV_inv", "Gamma1_keV"),
            ("Bc1(6750)->Bc* gamma", "g2_GeV_inv", "Gamma2_keV"),
        ]:
            mc_summary.append(
                {
                    "normalization": norm_label,
                    "channel": channel,
                    "observable": gkey,
                    **summary_stats([r[gkey] for r in subset]),
                }
            )
            mc_summary.append(
                {
                    "normalization": norm_label,
                    "channel": channel,
                    "observable": wkey,
                    **summary_stats([r[wkey] for r in subset]),
                }
            )
    write_csv(OUT / "bc_vec_complete_monte_carlo_summary.csv", mc_summary)

    plot_fbcstar_stability(
        fstar_rows,
        "f_Bcstar_standard_GeV",
        r"$f_{B_c^\ast}\,[\mathrm{GeV}]$",
        OUT / "bc_vec_fbcstar_standard_M2_stability.pdf",
    )
    plot_stability(
        grid_rows,
        "standard_vector_invariant",
        "g1_GeV_inv",
        r"$|g_1|\,[\mathrm{GeV}^{-1}]$",
        OUT / "bc_vec_standard_g1_M2_stability.pdf",
        r"$B_{c1}(6743)\to B_c^\ast\gamma$",
    )
    plot_stability(
        grid_rows,
        "standard_vector_invariant",
        "g2_GeV_inv",
        r"$|g_2|\,[\mathrm{GeV}^{-1}]$",
        OUT / "bc_vec_standard_g2_M2_stability.pdf",
        r"$B_{c1}(6750)\to B_c^\ast\gamma$",
    )
    plot_stability(
        grid_rows,
        "standard_vector_invariant",
        "Gamma1_keV",
        r"$\Gamma\,[\mathrm{keV}]$",
        OUT / "bc_vec_standard_gamma1_M2_stability.pdf",
        r"$B_{c1}(6743)\to B_c^\ast\gamma$",
    )
    plot_stability(
        grid_rows,
        "standard_vector_invariant",
        "Gamma2_keV",
        r"$\Gamma\,[\mathrm{keV}]$",
        OUT / "bc_vec_standard_gamma2_M2_stability.pdf",
        r"$B_{c1}(6750)\to B_c^\ast\gamma$",
    )

    lines = ["Leading Bc1 -> Bc* gamma summary", "=================================", ""]
    lines.append("Vector-current two-point normalization:")
    for item in fstar_summary:
        lines.append(
            "{quantity}: median={median:.4g}, 16-84%=[{q16:.4g},{q84:.4g}], "
            "range=[{min:.4g},{max:.4g}]".format(**item)
        )
    lines.append("The standard_vector_invariant value is the physical f_Bcstar convention used in the hadronic matrix element.")
    lines.append("")
    lines.append("Monte Carlo leading hard-photon vector results:")
    for item in mc_summary:
        lines.append(
            "{normalization} {channel} {observable}: median={median:.4g}, "
            "16-84%=[{q16:.4g},{q84:.4g}], range=[{min:.4g},{max:.4g}]".format(**item)
        )
    lines.append("")
    lines.append(
        "Status: perturbative hard-photon three-point OPE. Radiative G^2 is "
        "quoted through the separate conservative screening estimate; the "
        "contact-support audit finds no ordinary two-channel contact term."
    )
    lines.append("Stability plots written for f_Bcstar, g1, g2 and the two vector widths.")
    (OUT / "bc_vec_complete_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
