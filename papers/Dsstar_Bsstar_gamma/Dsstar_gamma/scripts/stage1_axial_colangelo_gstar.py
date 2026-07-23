"""Stage-1 axial-current benchmark for Ds1 -> Ds* gamma.

This script implements the Colangelo-De Fazio-Ozpineci sum rule for the
dimensionless epsilon-structure coupling in

    < gamma(q) Ds*(p) | Ds1(P) >
      = i e g_* eps_{alpha beta rho sigma}
          eta^alpha xi^{* beta} eps^{* rho} q^sigma .

Scope:
  - axial basis current only, J_A = sbar gamma_mu gamma5 Q;
  - final vector current, J_nu = sbar gamma_nu Q;
  - full Eq. (5.4) terms that are expressible with the photon DAs already in
    photon_da.py, including the F2 and F3 three-particle integrals;
  - output is contribution-by-contribution so the cancellation can be checked.

The mixed physical-state mapping is added as a controlled axial-only baseline:

    f_1 g_1^* = sin(theta) f_A g_A^* ,
    f_2 g_2^* = cos(theta) f_A g_A^* .

The tensor-current contribution, proportional to f_B g_B^*, is not included in
this Stage-1 script.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Dsstar_gamma" / "outputs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))
sys.path.insert(0, str(ROOT.parents[1] / "shared"))
sys.path.insert(0, str(ROOT / "scripts"))

import photon_da as pda


ALPHA_EM = 1.0 / 137.036


def gl_quad(func, a: float, b: float, n: int = 180) -> float:
    """Gauss-Legendre integration on a finite interval."""
    if b <= a:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (b - a) * x + 0.5 * (a + b)
    try:
        vals = func(t)
    except Exception:
        vals = np.array([func(float(tt)) for tt in t], dtype=float)
    return float(0.5 * (b - a) * np.sum(w * vals))


def kallen(s, m1, m2):
    return (s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)


def rho_gstar(s, mc, ms, ec, es):
    """Perturbative spectral density in Colangelo et al. Eq. (5.5)."""
    s = np.asarray(s, dtype=float)
    lam = np.sqrt(np.maximum(kallen(s, mc, ms), 0.0))

    def log_piece(mj, mi, charge):
        num = np.clip(s - mj**2 + mi**2 - lam, 1e-300, None)
        den = np.clip(s - mj**2 + mi**2 + lam, 1e-300, None)
        return charge * np.log(num / den)

    return (
        3.0
        * ms
        * mc
        / (4.0 * math.pi**2)
        * (log_piece(mc, ms, es) + log_piece(ms, mc, ec))
    )


def H_gamma(u0: float) -> float:
    """H_gamma(u0) = int_0^u0 dv h_gamma(v)."""
    return gl_quad(pda.h_gamma_colangelo, 0.0, u0, n=220)


def Hbar_gamma(u0: float) -> float:
    """Hbar_gamma(u0) = int_0^u0 du int_0^u dv h_gamma(v)."""
    return gl_quad(lambda v: (u0 - v) * pda.h_gamma_colangelo(v), 0.0, u0, n=220)


def Psi_v(u0: float) -> float:
    """Psi^v(u0) = int_0^u0 du psi^v(u)."""
    return gl_quad(pda.psi_v_colangelo, 0.0, u0, n=220)


def psi_a_prime(u0: float) -> float:
    """Numerical derivative of psi^a at u0."""
    h = 1.0e-5
    up = min(1.0 - 1.0e-8, u0 + h)
    um = max(1.0e-8, u0 - h)
    return float((pda.psi_a(up) - pda.psi_a(um)) / (up - um))


def F2(aq, aqb):
    """F2 = S + S_tilde + T1 - T2 - T3 + T4."""
    return (
        pda.S_3p(aq, aqb)
        + pda.St_3p(aq, aqb)
        + pda.T1_3p(aq, aqb)
        - pda.T2_3p(aq, aqb)
        - pda.T3_3p(aq, aqb)
        + pda.T4_3p(aq, aqb)
    )


def F3(aq, aqb):
    """F3 = A + V."""
    return pda.A_3p(aq, aqb) + pda.V_3p(aq, aqb)


def F2_integral(u0: float = 0.5) -> float:
    """The two F2 integrals in Eq. (5.4)."""

    def first_v(v):
        upper = u0 / max(1.0 - v, 1e-15)

        def inner(ag):
            return F2(u0 - (1.0 - v) * ag, 1.0 - u0 - v * ag)

        return gl_quad(inner, 0.0, upper, n=90)

    def second_v(v):
        upper = (1.0 - u0) / max(v, 1e-15)

        def inner(ag):
            return F2(u0 - (1.0 - v) * ag, 1.0 - u0 - v * ag)

        return gl_quad(inner, 0.0, upper, n=90)

    return gl_quad(first_v, 0.0, 1.0 - u0, n=90) + gl_quad(
        second_v, 1.0 - u0, 1.0, n=90
    )


def F3_integral(u0: float = 0.5) -> float:
    """The F3 combination in Eq. (5.4).

    The first double integral contains d alpha_g / alpha_g^2.  The photon DA
    F3 carries alpha_g^2, so the endpoint is finite; we still avoid evaluating
    exactly at zero for numerical stability.
    """

    eps = 1.0e-10

    def first_aqb(aqb):
        lower = max(u0 - aqb, eps)
        upper = 1.0 - aqb

        def inner(ag):
            return F3(1.0 - aqb - ag, aqb) / (ag**2)

        return gl_quad(inner, lower, upper, n=90)

    def second_aqb(aqb):
        den = max(u0 - aqb, eps)
        return F3(1.0 - u0, aqb) / den

    return gl_quad(first_aqb, 0.0, u0, n=90) - gl_quad(
        second_aqb, 0.0, u0, n=90
    )


def width_epsilon_keV(mA: float, mV: float, g: float) -> float:
    """Width for M = i e g epsilon(eta,xi,eps,q)."""
    delta = mA**2 - mV**2
    h33 = delta**2 * (mA**2 + mV**2) / (6.0 * mA**2 * mV**2)
    phase = delta / (16.0 * math.pi * mA**3)
    return 4.0 * math.pi * ALPHA_EM * phase * h33 * g**2 * 1.0e6


def gstar_axial_sum_rule(
    M2: float,
    s0: float,
    *,
    mc: float,
    ms: float,
    mA: float,
    mV: float,
    fA: float,
    fV: float,
    ss: float,
    chi: float,
    f3g: float,
    ec: float,
    es: float,
    u0: float = 0.5,
    chi_positive_convention: bool = True,
    f2_int: float | None = None,
    f3_int: float | None = None,
) -> dict[str, float]:
    """Compute the axial-current epsilon coupling g_A^* at one Borel point."""
    lower = (mc + ms) ** 2
    exp_mc = math.exp(-mc**2 / M2)
    exp_s0 = math.exp(-s0 / M2)
    exp_diff = exp_mc - exp_s0

    pert = gl_quad(
        lambda s: np.exp(-s / M2) * rho_gstar(s, mc, ms, ec, es),
        lower,
        s0,
        n=260,
    )

    heavy_local = (
        ec
        * mc
        * exp_mc
        * ss
        * (1.0 - (ms**2 / M2) * (1.0 - mc**2 / M2))
    )

    # Eq. (5.4) is printed with chi in the older sign convention.  Our project
    # stores chi as the positive magnitude.  The physical sign is therefore
    # flipped for the leading-twist term when chi_positive_convention=True.
    chi_sign = -1.0 if chi_positive_convention else 1.0
    twist2 = chi_sign * es * mc * ss * exp_diff * M2 * chi * pda.phi_gamma(u0)

    A = float(pda.A_t4_colangelo(u0))
    H = H_gamma(u0)
    Hbar = Hbar_gamma(u0)
    twist4_2p = (
        es
        * mc
        * ss
        * exp_mc
        * (
            -0.25 * (mc**2 / M2) * A
            - H * (1.0 - u0)
            - Hbar * (1.0 - 2.0 * mc**2 / M2)
        )
    )

    psiv = float(pda.psi_v_colangelo(u0))
    psia = float(pda.psi_a(u0))
    psia_p = psi_a_prime(u0)
    psiv_int = Psi_v(u0)
    twist3_2p = (
        es
        * f3g
        * M2
        * exp_diff
        * (
            0.25 * (1.0 - u0) * psia_p
            - 0.25 * psia
            - psiv_int * (1.0 + 2.0 * mc**2 / M2)
            + (1.0 - u0) * psiv
        )
    )

    if f2_int is None:
        f2_int = F2_integral(u0)
    if f3_int is None:
        f3_int = F3_integral(u0)

    twist4_3p = mc * es * ss * exp_mc * f2_int
    twist3_3p = -es * f3g * M2 * exp_diff * f3_int

    qcd_side = pert + heavy_local + twist2 + twist4_2p + twist3_2p + twist4_3p + twist3_3p
    prefactor = math.exp((mA**2 + mV**2) / (2.0 * M2)) / (mA * fA * mV * fV)
    gstar = prefactor * qcd_side

    return {
        "M2": M2,
        "s0": s0,
        "gA_star": gstar,
        "width_A_keV": width_epsilon_keV(mA, mV, gstar),
        "qcd_side": qcd_side,
        "pert": pert,
        "heavy_local": heavy_local,
        "twist2_chi_phi": twist2,
        "twist4_2p": twist4_2p,
        "twist3_2p_f3g": twist3_2p,
        "twist4_3p_F2": twist4_3p,
        "twist3_3p_F3": twist3_3p,
        "F2_integral": f2_int,
        "F3_integral": f3_int,
    }


def central_inputs():
    qq = -(0.240) ** 3
    return {
        "mc": 1.27,
        "ms": 0.093,
        "mA": 2.4595,
        "mA2": 2.53511,
        "mV": 2.1122,
        "fA": 0.225,
        # Calculated in Dsstar_gamma/scripts/twopoint_vector_decay_constant.py.
        "fV": 0.227,
        "f1": 0.345,
        "f2": 0.379,
        "ss": 0.8 * qq,
        "chi": 3.15,
        "f3g": -0.0039,
        "ec": 2.0 / 3.0,
        "es": -1.0 / 3.0,
    }


def colangelo_inputs():
    qq = -(0.245) ** 3
    return {
        "mc": 1.35,
        "ms": 0.125,
        "mA": 2.4595,
        "mA2": 2.53511,
        "mV": 2.1122,
        "fA": 0.225,
        "fV": 0.266,
        "f1": 0.225,
        "f2": 0.225,
        "ss": 0.8 * qq,
        "chi": 3.15,
        "f3g": -0.0039,
        "ec": 2.0 / 3.0,
        "es": -1.0 / 3.0,
    }


def run_grid(
    label,
    vals,
    *,
    chi_positive_convention=True,
    M2_values=None,
    s0_values=None,
    window_label="rounded_final",
):
    f2_int = F2_integral(0.5)
    f3_int = F3_integral(0.5)
    rows = []
    if M2_values is None:
        M2_values = np.linspace(3.0, 6.0, 13)
    if s0_values is None:
        s0_values = (7.5, 8.0, 8.5)
    for s0 in s0_values:
        for M2 in M2_values:
            args = {k: vals[k] for k in ("mc", "ms", "mA", "mV", "fA", "fV", "ss", "chi", "f3g", "ec", "es")}
            row = gstar_axial_sum_rule(
                float(M2),
                float(s0),
                **args,
                chi_positive_convention=chi_positive_convention,
                f2_int=f2_int,
                f3_int=f3_int,
            )
            theta = math.radians(35.3)
            g1_phys = math.sin(theta) * vals["fA"] * row["gA_star"] / vals["f1"]
            g2_phys = math.cos(theta) * vals["fA"] * row["gA_star"] / vals["f2"]
            row.update(
                {
                    "input_set": label,
                    "s0_window_label": window_label,
                    "theta_deg": 35.3,
                    "g1_2460_axial_only": g1_phys,
                    "g2_2536_axial_only": g2_phys,
                    "Gamma2460_axial_only_keV": width_epsilon_keV(vals["mA"], vals["mV"], g1_phys),
                    "Gamma2536_axial_only_keV": width_epsilon_keV(vals["mA2"], vals["mV"], g2_phys),
                }
            )
            rows.append(row)
    return rows


def summarize(rows, label):
    gA = np.array([r["gA_star"] for r in rows])
    g1 = np.array([r["g1_2460_axial_only"] for r in rows])
    g2 = np.array([r["g2_2536_axial_only"] for r in rows])
    w1 = np.array([r["Gamma2460_axial_only_keV"] for r in rows])
    w2 = np.array([r["Gamma2536_axial_only_keV"] for r in rows])
    return [
        f"{label}:",
        f"  g_A* range {gA.min():+.4f} to {gA.max():+.4f}",
        f"  axial-only g_1*(2460) range {g1.min():+.4f} to {g1.max():+.4f}; Gamma {w1.min():.3f} to {w1.max():.3f} keV",
        f"  axial-only g_2*(2536) range {g2.min():+.4f} to {g2.max():+.4f}; Gamma {w2.min():.3f} to {w2.max():.3f} keV",
    ]


def main():
    rows = []
    rows.extend(run_grid("project_inputs_positive_chi", central_inputs(), chi_positive_convention=True))
    rows.extend(run_grid("project_inputs_printed_chi_sign_diagnostic", central_inputs(), chi_positive_convention=False))
    rows.extend(run_grid("colangelo_like_positive_chi", colangelo_inputs(), chi_positive_convention=True))
    rows.extend(run_grid("colangelo_like_printed_chi_sign_diagnostic", colangelo_inputs(), chi_positive_convention=False))
    rows.extend(
        run_grid(
            "colangelo_like_original_window_positive_chi",
            colangelo_inputs(),
            chi_positive_convention=True,
            M2_values=np.linspace(4.0, 6.0, 9),
            s0_values=(2.50**2, 2.55**2, 2.60**2),
            window_label="colangelo_original",
        )
    )

    keys = [
        "input_set",
        "s0_window_label",
        "M2",
        "s0",
        "theta_deg",
        "gA_star",
        "width_A_keV",
        "g1_2460_axial_only",
        "Gamma2460_axial_only_keV",
        "g2_2536_axial_only",
        "Gamma2536_axial_only_keV",
        "qcd_side",
        "pert",
        "heavy_local",
        "twist2_chi_phi",
        "twist4_2p",
        "twist3_2p_f3g",
        "twist4_3p_F2",
        "twist3_3p_F3",
        "F2_integral",
        "F3_integral",
    ]
    csv_path = OUT / "stage1_axial_colangelo_gstar.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in keys})

    lines = [
        "Stage-1 axial-current Ds1 -> Ds* gamma benchmark",
        "================================================",
        "The physical rows are axial-only: tensor-current J_B is not included yet.",
        "Windows used here: M^2 = 3.0...6.0 GeV^2, s0 = 7.5, 8.0, 8.5 GeV^2.",
        "",
    ]
    for label in (
        "project_inputs_positive_chi",
        "project_inputs_printed_chi_sign_diagnostic",
        "colangelo_like_positive_chi",
        "colangelo_like_printed_chi_sign_diagnostic",
        "colangelo_like_original_window_positive_chi",
    ):
        lines.extend(summarize([r for r in rows if r["input_set"] == label], label))
        lines.append("")
    lines.append(f"Wrote grid to {csv_path}")

    summary_path = OUT / "stage1_axial_colangelo_gstar_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
