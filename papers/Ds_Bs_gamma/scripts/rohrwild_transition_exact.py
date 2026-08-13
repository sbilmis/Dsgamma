"""Exact post-double-Borel Rohrwild transition invariants.

This module is the numerical implementation of the formulas derived in
``notes/DsBs_gamma_Mathematica_derivation.tex``.  In particular:

* the two virtualities are kept independent until the double Borel transform;
* no physical pole mass is inserted in a QCD-side numerator;
* no standalone local ``<sbar s>`` transition term is included;
* the channel-specific continuum prescription is retained: the leading
  twist-2 term carries ``E_Q-E_0``, whereas the twist-3, twist-4, and
  gluonic three-particle terms carry the raw Borel factor ``E_Q``;
* the separate gauge-completion electromagnetic three-particle term follows
  Rohrwild's ``I_F``-type subtraction and carries ``E_Q-E_0``;
* the heavy-propagator ``x_alpha F^{alpha beta} gamma_beta`` term has no
  ``(slash k + m_Q)`` numerator, so its tensor-current trace vanishes;
* the electromagnetic ``P`` functional is represented by the derivative of
  its support-corrected line density.

The functions return the unnormalised Borel invariants ``T_A`` and ``T_B``.
Physical couplings are formed only afterwards, with the physical residues
``f_1`` and ``f_2`` and the appropriate pole mass for each state.
"""

from __future__ import annotations

import math
import sys
from functools import lru_cache
from pathlib import Path
from typing import Callable, Union

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1] / "shared"))
sys.path.insert(0, str(ROOT / "scripts"))

import photon_da as pda
from stage1_axial_g1_baseline import (
    gauss_legendre_integral,
    rho_p_colangelo_g1,
)
from stage1_tensor_gb_hard_candidate import rho_tensor_hard_candidate


NumericArray = Union[np.ndarray, float]
ArrayFunction4 = Callable[
    [NumericArray, NumericArray, NumericArray, float],
    NumericArray,
]
ArrayFunction3 = Callable[
    [NumericArray, NumericArray, NumericArray],
    NumericArray,
]


def _integrate_scalar(
    func: Callable[[float], float], lo: float, hi: float, n: int
) -> float:
    """Gauss--Legendre integral for a scalar callback."""

    if hi <= lo:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(n)
    points = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
    values = np.asarray([func(float(point)) for point in points], dtype=float)
    return float(0.5 * (hi - lo) * np.sum(w * values))


def line_projection(
    func: ArrayFunction4, z: float = 0.5, *, n: int = 80
) -> float:
    r"""Evaluate the delta-supported three-particle line projection.

    The heavy-line momentum fraction is

    .. math:: a=\alpha_{\bar q}+v\alpha_g.

    The two integration domains are the explicit solution of
    ``delta(z-a)`` on the physical simplex.
    """

    if not 0.0 < z < 1.0:
        raise ValueError("z must lie strictly between zero and one")

    xg, wg = np.polynomial.legendre.leggauss(n)

    def integrate_ag(v: float, upper: float) -> float:
        if upper <= 0.0:
            return 0.0
        ag = 0.5 * upper * (xg + 1.0)
        aq = 1.0 - z - (1.0 - v) * ag
        aqb = z - v * ag
        values = np.asarray(func(aq, aqb, ag, v), dtype=float)
        return float(0.5 * upper * np.sum(wg * values))

    domain_1 = _integrate_scalar(
        lambda v: integrate_ag(v, (1.0 - z) / (1.0 - v)),
        0.0,
        z,
        n,
    )
    domain_2 = _integrate_scalar(
        lambda v: integrate_ag(v, z / v),
        z,
        1.0,
        n,
    )
    return domain_1 + domain_2


def p_line(
    func: ArrayFunction3, z: float = 0.5, *, n: int = 44
) -> float:
    r"""Evaluate the support-corrected line density ``L_z[func]``.

    This is a direct numerical implementation of Eq. (app-P-line-density) in
    the polished manuscript.  The lower limit ``beta=z`` in the second term is
    essential; extending it to zero creates the spurious logarithmic
    divergence in the printed symmetric appendix formula.
    """

    if not 0.0 < z < 1.0:
        raise ValueError("z must lie strictly between zero and one")

    nodes, weights = np.polynomial.legendre.leggauss(n)

    def fixed_quad_vector(
        func_1d: Callable[[np.ndarray], NumericArray],
        lo: float,
        hi: float,
    ) -> float:
        if hi <= lo:
            return 0.0
        points = 0.5 * (hi - lo) * nodes + 0.5 * (hi + lo)
        values = np.asarray(func_1d(points), dtype=float)
        return float(0.5 * (hi - lo) * np.sum(weights * values))

    # The alpha_q' integral in the first line does not occur in the
    # integrand and therefore equals alpha_q.
    def term_1_at_aq(aq: float) -> float:
        def inner(aqb: np.ndarray) -> np.ndarray:
            ag = 1.0 - aq - aqb
            return (
                aq
                * (z - aqb)
                / (ag * ag)
                * np.asarray(func(aq, aqb, ag), dtype=float)
            )

        return fixed_quad_vector(inner, 0.0, z)

    term_1 = _integrate_scalar(term_1_at_aq, 0.0, 1.0 - z, n)

    def term_2_at_aq(aq: float) -> float:
        def beta_integrand(beta: float) -> float:
            inner = fixed_quad_vector(
                lambda aqbp: func(aq, aqbp, 1.0 - aq - aqbp),
                0.0,
                beta,
            )
            return z * inner / (beta * beta)

        return _integrate_scalar(beta_integrand, z, 1.0 - aq, n)

    term_2 = _integrate_scalar(term_2_at_aq, 0.0, 1.0 - z, n)

    def simplex_primitive(aq_upper: float) -> float:
        def at_aqp(aqp: float) -> float:
            return fixed_quad_vector(
                lambda aqb: func(aqp, aqb, 1.0 - aqp - aqb),
                0.0,
                1.0 - aqp,
            )

        return _integrate_scalar(at_aqp, 0.0, aq_upper, n)

    def term_3_at_aq(aq: float) -> float:
        return z * simplex_primitive(aq) / ((1.0 - aq) ** 2)

    term_3 = _integrate_scalar(term_3_at_aq, 0.0, 1.0 - z, n)
    return term_1 - term_2 + term_3


def p_line_prime(
    func: ArrayFunction3,
    z: float = 0.5,
    *,
    n: int = 44,
    step: float = 2.0e-4,
) -> float:
    """Five-point numerical derivative of the support-corrected line density."""

    if not 2.0 * step < z < 1.0 - 2.0 * step:
        raise ValueError("z and step do not leave room for a five-point stencil")
    values = [
        p_line(func, z - 2.0 * step, n=n),
        p_line(func, z - step, n=n),
        p_line(func, z + step, n=n),
        p_line(func, z + 2.0 * step, n=n),
    ]
    return float((values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * step))


def f_g_sigma(aq, aqb, ag, v):
    del ag, v
    return (
        pda.S_3p(aq, aqb)
        + pda.St_3p(aq, aqb)
        - pda.T1_3p(aq, aqb)
        - pda.T2_3p(aq, aqb)
        + pda.T3_3p(aq, aqb)
        + pda.T4_3p(aq, aqb)
    )


def f_g_xgamma(aq, aqb, ag, v):
    del ag
    return 2.0 * v * (
        -pda.S_3p(aq, aqb)
        - pda.T3_3p(aq, aqb)
        + pda.T2_3p(aq, aqb)
    )


def f_g_axial(aq, aqb, ag, v):
    return f_g_sigma(aq, aqb, ag, v) + f_g_xgamma(aq, aqb, ag, v)


def f_em_axial(aq, aqb, ag, v):
    del ag
    return (1.0 - 2.0 * v) * pda.Sg_3p(aq, aqb) - pda.T4g_3p(aq, aqb)


def f_em_tensor_base(aq, aqb, ag, v):
    del ag, v
    return pda.Sg_3p(aq, aqb) - pda.T4g_3p(aq, aqb)


def t4gamma(aq, aqb, ag):
    del ag
    return pda.T4g_3p(aq, aqb)


@lru_cache(maxsize=None)
def precompute_convolutions(
    z: float = 0.5,
    line_order: int = 80,
    p_order: int = 44,
    derivative_step: float = 2.0e-4,
) -> dict[str, float]:
    """Precompute the DA-only convolution constants at a fixed Borel saddle."""

    j_g_sigma = line_projection(f_g_sigma, z, n=line_order)
    j_g_xgamma = line_projection(f_g_xgamma, z, n=line_order)
    j_g_axial = line_projection(f_g_axial, z, n=line_order)
    j_em_axial = line_projection(f_em_axial, z, n=line_order)
    j_em_tensor_base = line_projection(f_em_tensor_base, z, n=line_order)
    l_t4 = p_line(t4gamma, z, n=p_order)
    l_t4_prime = p_line_prime(
        t4gamma,
        z,
        n=p_order,
        step=derivative_step,
    )
    return {
        "z": z,
        "J_g_sigma": j_g_sigma,
        "J_g_xgamma": j_g_xgamma,
        "J_g_axial": j_g_axial,
        "J_g_closure": j_g_axial - j_g_sigma - j_g_xgamma,
        "J_em_axial": j_em_axial,
        "J_em_tensor_base": j_em_tensor_base,
        "L_T4gamma": l_t4,
        "Lprime_T4gamma": l_t4_prime,
    }


def bar_psi_v(a: float) -> float:
    r"""Rohrwild ``bar psi^(V) = 2 Psi^v`` in the paper's dictionary."""

    return float(pda.psi_v(a))


def bar_psi_v_prime(a: float) -> float:
    """Analytic derivative of ``bar_psi_v``."""

    omega_a = pda.val("omegaA")
    omega_v = pda.val("omegaV")
    t = 2.0 * a - 1.0
    # psi^(V) = -20 u(1-u)t
    #           + 15/16 (omegaA-3 omegaV) u(1-u)t(7t^2-3).
    c = 15.0 / 16.0 * (omega_a - 3.0 * omega_v)
    duub = 1.0 - 2.0 * a
    base_prime = duub * t + 2.0 * a * (1.0 - a)
    shape_prime = (
        base_prime * (7.0 * t * t - 3.0)
        + a * (1.0 - a) * t * 28.0 * t
    )
    return float(-20.0 * base_prime + c * shape_prime)


def width_keV(m_initial: float, m_final: float, coupling: float) -> float:
    """E1 width for the coupling convention used in the manuscript."""

    alpha = 1.0 / 137.036
    q_gamma = (m_initial * m_initial - m_final * m_final) / (2.0 * m_initial)
    return alpha * coupling * coupling * q_gamma**3 * 1.0e6 / 3.0


def transition_invariants(
    M2: float,
    s0: float,
    inputs: dict[str, float],
    *,
    a0: float = 0.5,
    convolutions: dict[str, float] | None = None,
) -> dict[str, Union[float, str]]:
    """Return every exact Rohrwild post-Borel contribution to ``T_A,T_B``."""

    if M2 <= 0.0:
        raise ValueError("M2 must be positive")
    m_q = float(inputs["mc"])
    m_s = float(inputs["ms"])
    e_q = float(inputs["ec"])
    e_s = float(inputs["es"])
    ss = float(inputs["ss"])
    d = m_q + m_s
    r_q = m_q / d
    threshold = d * d
    if s0 <= threshold:
        raise ValueError("s0 must exceed the perturbative threshold")

    conv = convolutions or precompute_convolutions(a0)
    e_mass = math.exp(-(m_q * m_q) / M2)
    e_cont = math.exp(-s0 / M2)
    delta_exp = e_mass - e_cont

    t_a_pert = gauss_legendre_integral(
        lambda spectral_s: np.exp(-spectral_s / M2)
        * rho_p_colangelo_g1(spectral_s, m_q, m_s, e_q, e_s),
        threshold,
        s0,
        n=1600,
    )
    t_b_pert = gauss_legendre_integral(
        lambda spectral_s: np.exp(-spectral_s / M2)
        * rho_tensor_hard_candidate(spectral_s, m_q, m_s, e_q, e_s),
        threshold,
        s0,
        n=1600,
    )

    t_a_tw2 = (
        e_s
        * ss
        * delta_exp
        * M2
        * float(inputs["chi"])
        * float(pda.phi_gamma(a0))
    )
    t_a_tw4 = (
        e_s
        * ss
        * e_mass
        * 0.25
        * (float(pda.A_t4(a0)) + 2.0 * float(pda.B_t4(a0)))
        * (1.0 + m_q * m_q / M2)
    )
    bar_psi = bar_psi_v(a0)
    bar_psi_prime = bar_psi_v_prime(a0)
    t_a_tw3 = e_s * float(inputs["f3g"]) * m_q * e_mass * bar_psi

    t_b_tw2 = r_q * t_a_tw2
    t_b_tw4 = r_q * t_a_tw4
    derivative_product = -bar_psi + (1.0 - a0) * bar_psi_prime
    t_b_tw3 = (
        e_s
        * float(inputs["f3g"])
        * e_mass
        / d
        * (m_q * m_q * bar_psi - 0.5 * M2 * derivative_product)
    )

    t_a_3pg = e_s * ss * e_mass * conv["J_g_axial"]
    t_b_3pg = e_s * ss * e_mass * r_q * conv["J_g_sigma"]
    t_a_3pem = (
        e_q
        * ss
        * delta_exp
        * (conv["J_em_axial"] + 2.0 * conv["Lprime_T4gamma"])
    )
    t_b_3pem = (
        e_q
        * ss
        * delta_exp
        * r_q
        * (conv["J_em_tensor_base"] + 2.0 * conv["Lprime_T4gamma"])
    )

    t_a = t_a_pert + t_a_tw2 + t_a_tw3 + t_a_tw4 + t_a_3pg + t_a_3pem
    t_b = t_b_pert + t_b_tw2 + t_b_tw3 + t_b_tw4 + t_b_3pg + t_b_3pem
    return {
        "transition_scheme": "rohrwild_nonlocal_exact_post_borel",
        "ordinary_local_transition_condensate": 0.0,
        "M2": M2,
        "s0": s0,
        "a0": a0,
        "E_Q": e_mass,
        "E_0": e_cont,
        "Delta_E_Q": delta_exp,
        "T_A_pert": t_a_pert,
        "T_A_tw2": t_a_tw2,
        "T_A_tw3": t_a_tw3,
        "T_A_tw4": t_a_tw4,
        "T_A_3p_g": t_a_3pg,
        "T_A_3p_gamma": t_a_3pem,
        "T_A": t_a,
        "T_B_pert": t_b_pert,
        "T_B_tw2": t_b_tw2,
        "T_B_tw3": t_b_tw3,
        "T_B_tw4": t_b_tw4,
        "T_B_3p_g": t_b_3pg,
        "T_B_3p_gamma": t_b_3pem,
        "T_B": t_b,
        "bar_psi_V": bar_psi,
        "bar_psi_V_prime": bar_psi_prime,
        **conv,
    }


def physical_couplings(
    invariants: dict[str, Union[float, str]],
    *,
    theta_deg: float,
    m_state_1: float,
    m_state_2: float,
    f_1: float,
    f_2: float,
    m_p: float,
    f_p: float,
    m_q: float,
    m_s: float,
) -> dict[str, float]:
    """Project ``T_A,T_B`` and normalize with the two physical pole residues."""

    if min(m_state_1, m_state_2, f_1, f_2, m_p, f_p) <= 0.0:
        raise ValueError("masses and decay constants must be positive")
    theta = math.radians(theta_deg)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    t_a = float(invariants["T_A"])
    t_b = float(invariants["T_B"])
    M2 = float(invariants["M2"])
    common = (m_q + m_s) / (m_p * m_p * f_p)
    n_1 = (
        math.exp((m_state_1 * m_state_1 + m_p * m_p) / (2.0 * M2))
        * common
        / (m_state_1 * f_1)
    )
    n_2 = (
        math.exp((m_state_2 * m_state_2 + m_p * m_p) / (2.0 * M2))
        * common
        / (m_state_2 * f_2)
    )
    g_1_a = n_1 * sin_theta * t_a
    g_1_b = n_1 * cos_theta * t_b
    g_2_a = n_2 * cos_theta * t_a
    g_2_b = -n_2 * sin_theta * t_b
    g_1 = g_1_a + g_1_b
    g_2 = g_2_a + g_2_b
    return {
        "N_1": n_1,
        "N_2": n_2,
        "g_1": g_1,
        "g_2": g_2,
        "g_1_A": g_1_a,
        "g_1_B": g_1_b,
        "g_2_A": g_2_a,
        "g_2_B": g_2_b,
        "Gamma_1_keV": width_keV(m_state_1, m_p, g_1),
        "Gamma_2_keV": width_keV(m_state_2, m_p, g_2),
    }
