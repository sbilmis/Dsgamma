"""Stage-1 axial-current benchmark for Ds1 -> Ds gamma.

This script implements the Colangelo-De Fazio-Ozpineci single-axial-current
sum rule for the coupling g1 in D_{s1}(2460) -> D_s gamma, adapted to the
bookkeeping agreed in this project.

Scope:
  - axial current only: J_A = sbar gamma_mu gamma5 c;
  - Stage 1: perturbative hard photon + two-particle photon DAs;
  - all explicit m_s terms in the published formula are retained;
  - the three-particle DA term F1 is set to zero and added later in Stage 2.

The output is a CSV table over M^2 and s0, plus a small text summary.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT))

import photon_da as pda


def gauss_legendre_integral(func, a: float, b: float, n: int = 240) -> float:
    """Numerically integrate func on [a,b] using Gauss-Legendre quadrature."""
    if b <= a:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (b - a) * x + 0.5 * (b + a)
    return float(0.5 * (b - a) * np.sum(w * func(t)))


def kallen(s, mc, ms):
    """Kallen function lambda(s, mc^2, ms^2)."""
    return (s - (mc + ms) ** 2) * (s - (mc - ms) ** 2)


def rho_piece(s, emitting_mass, spectator_mass, charge):
    """One charge-labeled part of rho^P(s).

    This implements the expression inside Colangelo et al. Eq. (g1), before the
    explicit ``-(s <-> c)'' subtraction.  The literal substitution is

      m_s -> emitting_mass, m_c -> spectator_mass, e_s -> charge.
    """
    s = np.asarray(s, dtype=float)
    mi = emitting_mass
    mj = spectator_mass
    lam = np.sqrt(np.maximum(kallen(s, mj, mi), 0.0))
    # The log is finite at threshold by limit; clipping avoids endpoint roundoff.
    num = np.clip(s - mj**2 + mi**2 - lam, 1e-300, None)
    den = np.clip(s - mj**2 + mi**2 + lam, 1e-300, None)
    log_term = np.log(num / den)
    braces = (
        2.0 * mi * log_term
        + (mj - mi) * (mj**2 - mi**2 - s) * lam / s**2
    )
    return -3.0 * charge * braces / (8.0 * math.pi**2)


def rho_p_colangelo_g1(s, mc, ms, ec, es):
    """Perturbative spectral density rho^P(s) for the axial g1 sum rule."""
    return rho_piece(s, ms, mc, es) - rho_piece(s, mc, ms, ec)


def hbar_gamma(u0: float) -> float:
    """bar H_gamma(u0) = int_0^u0 du' int_0^u' dv h_gamma(v)."""
    return gauss_legendre_integral(
        lambda u: (u0 - u) * pda.h_gamma_colangelo(u), 0.0, u0, n=240
    )


def psi_v_integral(u0: float) -> float:
    """Psi^v(u0) = int_0^u0 du psi^v(u)."""
    return gauss_legendre_integral(pda.psi_v_colangelo, 0.0, u0, n=240)


def g1_stage1(
    M2: float,
    s0: float,
    *,
    mc: float,
    ms: float,
    m_ds1: float,
    m_ds: float,
    f_ds1: float,
    f_ds: float,
    ss: float,
    chi: float,
    f3g: float,
    ec: float,
    es: float,
    u0: float = 0.5,
    chi_positive_convention: bool = True,
) -> dict[str, float]:
    """Compute the Stage-1 axial coupling g1 at one Borel point."""
    lower = (mc + ms) ** 2
    pert = gauss_legendre_integral(
        lambda s: np.exp(-s / M2) * rho_p_colangelo_g1(s, mc, ms, ec, es),
        lower,
        s0,
        n=320,
    )

    exp_mc = math.exp(-mc**2 / M2)
    exp_s0 = math.exp(-s0 / M2)

    heavy_local = (
        ec
        * exp_mc
        * ss
        * (
            1.0
            - mc * ms / M2
            + (ms**2 / (2.0 * M2)) * (1.0 + mc**2 / M2)
        )
    )

    # Colangelo et al. print chi = -3.15 GeV^-2 in a convention where the g1
    # sum rule contains - e_s <sbar s> chi phi_gamma.  In this project we store
    # chi as the positive BBK magnitude +3.15.  To preserve the same physical
    # matrix element, the project convention flips this term relative to the
    # literal printed formula.
    chi_sign = 1.0 if chi_positive_convention else -1.0
    twist2 = chi_sign * es * ss * (exp_mc - exp_s0) * M2 * chi * pda.phi_gamma(u0)

    A = float(pda.A_t4_colangelo(u0))
    hbar = hbar_gamma(u0)
    twist4_2p = (
        -es
        * ss
        * exp_mc
        * (
            -0.25 * (A - 8.0 * hbar) * (1.0 + mc**2 / M2)
        )
    )

    twist3_2p = 2.0 * es * f3g * mc * exp_mc * psi_v_integral(u0)

    qcd_side = pert + heavy_local + twist2 + twist4_2p + twist3_2p
    prefactor = (
        math.exp((m_ds1**2 + m_ds**2) / (2.0 * M2))
        * (mc + ms)
        / (m_ds1 * f_ds1 * m_ds**2 * f_ds)
    )
    g1 = prefactor * qcd_side

    qgamma = (m_ds1**2 - m_ds**2) / (2.0 * m_ds1)
    width_gev = (1.0 / 137.036) / 3.0 * g1**2 * qgamma**3
    return {
        "M2": M2,
        "s0": s0,
        "g1_GeV_inv": g1,
        "width_keV": width_gev * 1.0e6,
        "qcd_side": qcd_side,
        "pert": pert,
        "heavy_local": heavy_local,
        "twist2": twist2,
        "twist4_2p": twist4_2p,
        "twist3_2p": twist3_2p,
    }


def central_inputs():
    """Central project inputs for the charm-sector axial baseline."""
    qq = -(0.240) ** 3
    return {
        "mc": 1.27,
        "ms": 0.093,
        "m_ds1": 2.4595,
        "m_ds": 1.96835,
        "f_ds1": 0.225,
        "f_ds": 0.2499,
        "ss": 0.8 * qq,
        "chi": 3.15,
        "f3g": -0.0039,
        "ec": 2.0 / 3.0,
        "es": -1.0 / 3.0,
    }


def colangelo_inputs():
    """Approximate inputs quoted around Colangelo et al. Eq. (g1)."""
    qq = -(0.245) ** 3
    return {
        "mc": 1.35,
        "ms": 0.125,
        "m_ds1": 2.4595,
        "m_ds": 1.96835,
        "f_ds1": 0.225,
        "f_ds": 0.266,
        "ss": 0.8 * qq,
        # User fixed positive chi for our convention.  The paper prints a
        # negative value, so we also run a paper-sign diagnostic below.
        "chi": 3.15,
        "f3g": -0.0039,
        "ec": 2.0 / 3.0,
        "es": -1.0 / 3.0,
    }


def run_grid(label: str, inputs: dict[str, float]):
    rows = []
    for s0_root in (2.50, 2.55, 2.60):
        s0 = s0_root**2
        for M2 in np.linspace(3.0, 6.0, 13):
            row = g1_stage1(float(M2), s0, **inputs)
            row["input_set"] = label
            row["s0_root"] = s0_root
            rows.append(row)
    return rows


def run_grid_printed_colangelo_chi_sign(label: str, inputs: dict[str, float]):
    rows = []
    for s0_root in (2.50, 2.55, 2.60):
        s0 = s0_root**2
        for M2 in np.linspace(3.0, 6.0, 13):
            row = g1_stage1(
                float(M2),
                s0,
                **inputs,
                chi_positive_convention=False,
            )
            row["input_set"] = label
            row["s0_root"] = s0_root
            rows.append(row)
    return rows


def summarize(rows, label):
    gs = np.array([r["g1_GeV_inv"] for r in rows])
    ws = np.array([r["width_keV"] for r in rows])
    return (
        f"{label}: g1 range {gs.min():+.4f} to {gs.max():+.4f} GeV^-1; "
        f"width range {ws.min():.2f} to {ws.max():.2f} keV"
    )


def main():
    all_rows = []
    all_rows.extend(run_grid("project_inputs_chi_positive_stage1", central_inputs()))
    all_rows.extend(run_grid("colangelo_inputs_chi_positive_stage1", colangelo_inputs()))

    all_rows.extend(
        run_grid_printed_colangelo_chi_sign(
            "colangelo_positive_chi_inserted_in_printed_formula_diagnostic",
            colangelo_inputs(),
        )
    )

    csv_path = OUT / "stage1_axial_g1_baseline.csv"
    keys = [
        "input_set",
        "M2",
        "s0_root",
        "s0",
        "g1_GeV_inv",
        "width_keV",
        "qcd_side",
        "pert",
        "heavy_local",
        "twist2",
        "twist4_2p",
        "twist3_2p",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in all_rows:
            writer.writerow({k: row[k] for k in keys})

    summary_lines = [
        "Stage-1 axial-current g1 baseline",
        "=================================",
        "Three-particle DA F1 term is omitted here by construction.",
        summarize(
            [r for r in all_rows if r["input_set"] == "project_inputs_chi_positive_stage1"],
            "Project inputs, chi=+3.15",
        ),
        summarize(
            [r for r in all_rows if r["input_set"] == "colangelo_inputs_chi_positive_stage1"],
            "Colangelo-like inputs, chi=+3.15",
        ),
        summarize(
            [r for r in all_rows if r["input_set"] == "colangelo_positive_chi_inserted_in_printed_formula_diagnostic"],
            "Diagnostic: positive chi inserted in printed negative-chi formula",
        ),
        f"Wrote grid to {csv_path}",
    ]
    summary_path = OUT / "stage1_axial_g1_baseline_summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
