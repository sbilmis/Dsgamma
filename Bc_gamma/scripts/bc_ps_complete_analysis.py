#!/usr/bin/env python3
"""Controlled Bc1 -> Bc gamma analysis.

This script upgrades the older pilot table in two ways:

1. It extracts the physical two-point residues f1 and f2 from the same
   perturbative heavy-heavy Bc-mixing sum rules used for the mixing-angle
   analysis.
2. It recomputes only the pseudoscalar final-state radiative channel
   Bc1 -> Bc gamma, producing stability tables, publication-style plots and a
   Monte Carlo uncertainty sample.

The three-point radiative OPE in this script is still the leading hard-photon
perturbative contribution.  Gluon-condensate corrections to the radiative
three-point function are not guessed here; they are documented as a missing
power correction to be derived or bounded separately.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".mplconfig"))
(OUT / ".mplconfig").mkdir(exist_ok=True)

import matplotlib.pyplot as plt  # noqa: E402


ALPHA_EM = 1.0 / 137.035999084
GEV_TO_KEV = 1.0e6


@dataclass(frozen=True)
class Inputs:
    mb: float = 4.18
    mc: float = 1.27
    e_c: float = 2.0 / 3.0
    e_bbar: float = 1.0 / 3.0
    theta_deg: float = 43.2950875223707
    theta_sigma_deg: float = 0.1513163252823016
    m_bc: float = 6.2749
    m_low: float = 6.743
    m_high: float = 6.751
    # Adopted pseudoscalar decay constant benchmark.  The numerical uncertainty
    # is intentionally conservative until a dedicated Bc two-point update is
    # added to this folder.
    f_bc: float = 0.371
    f_bc_sigma: float = 0.037


_LEGENDRE_CACHE: dict[int, tuple[np.ndarray, np.ndarray]] = {}


def gauss_legendre_integral(func, a: float, b: float, n: int = 160) -> float:
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


def rho_ps_A(s, inp: Inputs):
    mb, mc = inp.mb, inp.mc
    c = line_kernel(s, mc, mb, mc, mc**2 * (mc - mb))
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


def hard_integral(M2: float, s0: float, rho_func, inp: Inputs) -> float:
    lower = (inp.mb + inp.mc) ** 2
    return gauss_legendre_integral(lambda x: np.exp(-x / M2) * rho_func(x), lower, s0)


def photon_energy(mi: float, mf: float) -> float:
    return (mi * mi - mf * mf) / (2.0 * mi)


def width_ps_keV(g: float, mi: float, mf: float) -> float:
    eg = photon_energy(mi, mf)
    return ALPHA_EM / 3.0 * g * g * eg**3 * GEV_TO_KEV


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
        "theta_deg": math.degrees(theta),
        "Pi_AA": aa,
        "Pi_AB": ab,
        "Pi_BB": bb,
        "Pi_1": pi1,
        "Pi_2": pi2,
        "f1_GeV": f1,
        "f2_GeV": f2,
    }


def ps_couplings(M2: float, s0: float, inp: Inputs, theta_deg: float, f1: float, f2: float):
    int_A = hard_integral(M2, s0, lambda x: rho_ps_A(x, inp), inp)
    int_B = hard_integral(M2, s0, lambda x: rho_ps_B(x, inp), inp)
    th = math.radians(theta_deg)
    st, ct = math.sin(th), math.cos(th)
    pref1 = (
        math.exp((inp.m_low**2 + inp.m_bc**2) / (2.0 * M2))
        * (inp.mb + inp.mc)
        / (inp.m_bc**2 * inp.f_bc * inp.m_low * f1)
    )
    pref2 = (
        math.exp((inp.m_high**2 + inp.m_bc**2) / (2.0 * M2))
        * (inp.mb + inp.mc)
        / (inp.m_bc**2 * inp.f_bc * inp.m_high * f2)
    )
    g1 = pref1 * (st * int_A + ct * int_B)
    g2 = pref2 * (ct * int_A - st * int_B)
    return {
        "int_A": int_A,
        "int_B": int_B,
        "g1_GeV_inv": g1,
        "g2_GeV_inv": g2,
        "Gamma1_keV": width_ps_keV(g1, inp.m_low, inp.m_bc),
        "Gamma2_keV": width_ps_keV(g2, inp.m_high, inp.m_bc),
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


def plot_stability(rows: list[dict], key: str, ylabel: str, path: Path, title: str):
    plt.figure(figsize=(5.0, 3.6))
    for s0 in sorted({r["s0"] for r in rows}):
        sub = [r for r in rows if r["s0"] == s0]
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


def plot_hist(samples: list[dict], key: str, path: Path, title: str):
    vals = np.asarray([r[key] for r in samples], dtype=float)
    stats = summary_stats(vals)
    plt.figure(figsize=(5.0, 3.6))
    plt.hist(vals, bins=36, density=True, alpha=0.72, color="#4C78A8")
    plt.axvline(stats["median"], color="black", lw=1.2, label="median")
    plt.axvline(stats["q16"], color="#D62728", lw=1.0, ls="--", label="16-84%")
    plt.axvline(stats["q84"], color="#D62728", lw=1.0, ls="--")
    plt.xlabel(r"$\Gamma\,[\mathrm{keV}]$")
    plt.ylabel("density")
    plt.text(0.04, 0.92, title, transform=plt.gca().transAxes)
    plt.grid(alpha=0.2)
    plt.legend(frameon=False, fontsize=9)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def run_grid(inp: Inputs):
    rows = []
    residue_rows = []
    for M2 in [7.0, 8.0, 9.0]:
        for s0 in [53.0, 54.0, 55.0]:
            res = twopoint_residues(M2, s0, inp)
            residue_rows.append({"M2": M2, "s0": s0, **res})
            rad = ps_couplings(M2, s0, inp, inp.theta_deg, res["f1_GeV"], res["f2_GeV"])
            rows.append({"M2": M2, "s0": s0, **res, **rad})
    return residue_rows, rows


def run_monte_carlo(n: int = 1000, seed: int = 20260711):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        mb = rng.uniform(4.16, 4.21)
        mc = rng.uniform(1.25, 1.29)
        M2 = rng.uniform(7.0, 9.0)
        s0 = rng.uniform(53.0, 55.0)
        theta = rng.normal(43.2950875223707, 0.1513163252823016)
        f_bc = max(0.1, rng.normal(0.371, 0.037))
        inp = Inputs(mb=mb, mc=mc, theta_deg=theta, f_bc=f_bc)
        if s0 <= (mb + mc) ** 2:
            continue
        res = twopoint_residues(M2, s0, inp)
        rad = ps_couplings(M2, s0, inp, theta, res["f1_GeV"], res["f2_GeV"])
        rows.append(
            {
                "index": len(rows) + 1,
                "mb": mb,
                "mc": mc,
                "M2": M2,
                "s0": s0,
                "theta_deg": theta,
                "f_Bc_GeV": f_bc,
                **res,
                **rad,
            }
        )
    return rows


def main():
    inp = Inputs()
    residue_rows, grid_rows = run_grid(inp)
    write_csv(OUT / "bc1_twopoint_residue_grid.csv", residue_rows)
    write_csv(OUT / "bc_ps_complete_grid.csv", grid_rows)

    residue_summary = []
    for key in ["theta_deg", "f1_GeV", "f2_GeV"]:
        residue_summary.append({"quantity": key, **summary_stats([r[key] for r in residue_rows])})
    write_csv(OUT / "bc1_twopoint_residue_summary.csv", residue_summary)

    grid_summary = []
    for channel, gkey, wkey in [
        ("Bc1(6743)->Bc gamma", "g1_GeV_inv", "Gamma1_keV"),
        ("Bc1(6750)->Bc gamma", "g2_GeV_inv", "Gamma2_keV"),
    ]:
        grid_summary.append({"channel": channel, "observable": gkey, **summary_stats([r[gkey] for r in grid_rows])})
        grid_summary.append({"channel": channel, "observable": wkey, **summary_stats([r[wkey] for r in grid_rows])})
    write_csv(OUT / "bc_ps_complete_grid_summary.csv", grid_summary)

    samples = run_monte_carlo()
    write_csv(OUT / "bc_ps_complete_monte_carlo.csv", samples)
    mc_summary = []
    for channel, gkey, wkey in [
        ("Bc1(6743)->Bc gamma", "g1_GeV_inv", "Gamma1_keV"),
        ("Bc1(6750)->Bc gamma", "g2_GeV_inv", "Gamma2_keV"),
    ]:
        mc_summary.append({"channel": channel, "observable": gkey, **summary_stats([r[gkey] for r in samples])})
        mc_summary.append({"channel": channel, "observable": wkey, **summary_stats([r[wkey] for r in samples])})
    write_csv(OUT / "bc_ps_complete_monte_carlo_summary.csv", mc_summary)

    plot_stability(grid_rows, "g1_GeV_inv", r"$|g_1|\,[\mathrm{GeV}^{-1}]$", OUT / "bc_ps_g1_M2_stability.pdf", r"$B_{c1}(6743)\to B_c\gamma$")
    plot_stability(grid_rows, "g2_GeV_inv", r"$|g_2|\,[\mathrm{GeV}^{-1}]$", OUT / "bc_ps_g2_M2_stability.pdf", r"$B_{c1}(6750)\to B_c\gamma$")
    plot_hist(samples, "Gamma1_keV", OUT / "bc_ps_gamma1_mc_hist.pdf", r"$B_{c1}(6743)\to B_c\gamma$")
    plot_hist(samples, "Gamma2_keV", OUT / "bc_ps_gamma2_mc_hist.pdf", r"$B_{c1}(6750)\to B_c\gamma$")

    lines = ["Controlled Bc1 -> Bc gamma summary", "====================================", ""]
    lines.append("Two-point residue grid (perturbative Bc-mixing moment):")
    for item in residue_summary:
        lines.append(
            "{quantity}: median={median:.4g}, 16-84%=[{q16:.4g},{q84:.4g}], "
            "range=[{min:.4g},{max:.4g}]".format(**item)
        )
    lines.append("")
    lines.append("Monte Carlo radiative results, leading perturbative three-point OPE:")
    for item in mc_summary:
        lines.append(
            "{channel} {observable}: median={median:.4g}, 16-84%=[{q16:.4g},{q84:.4g}], "
            "range=[{min:.4g},{max:.4g}]".format(**item)
        )
    lines.append("")
    lines.append("Radiative condensate status: not included; requires a separate three-point background-field derivation.")
    (OUT / "bc_ps_complete_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
