"""Independent Python implementation of the retained LO QCD sum rules.

Only NumPy and the Python standard library are required.  Every function uses
the conventions in ../conventions.md.  The transition function returns a
term-by-term ledger so that local-condensate diagnostics cannot be silently
mixed into the quoted nonlocal LCSR.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from math import exp, log, pi, sqrt
from typing import Callable

import numpy as np


ALPHA_EM = 0.0072973525693
Q_S = -1.0 / 3.0
Q_C = 2.0 / 3.0
Q_B = -1.0 / 3.0


@dataclass(frozen=True)
class QCDInput:
    mQ: float
    ms: float = 0.093
    eQ: float = Q_C
    es: float = Q_S
    qq: float = -(0.24**3)
    kappa_s: float = 0.8
    m0sq: float = 0.8
    chi: float = -3.15
    f3: float = -0.0039
    omega_v: float = 3.8
    omega_a: float = -2.1
    kappa_da: float = 0.2
    zeta1: float = 0.4
    zeta2: float = 0.3

    @property
    def ss(self) -> float:
        return self.kappa_s * self.qq

    @property
    def mixed_ss(self) -> float:
        return self.m0sq * self.ss


def lambda_kallen(s: np.ndarray | float, a: float, b: float):
    return s * s + a * a + b * b - 2.0 * (s * a + s * b + a * b)


@lru_cache(maxsize=32)
def _legendre(n: int):
    return np.polynomial.legendre.leggauss(n)


def integrate_gauss(
    func: Callable[[np.ndarray], np.ndarray | float],
    lo: float,
    hi: float,
    n: int = 128,
) -> float:
    if hi <= lo:
        return 0.0
    x, w = _legendre(n)
    z = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
    return float(0.5 * (hi - lo) * np.sum(w * np.asarray(func(z), dtype=float)))


def rho_twopoint(s: np.ndarray | float, inp: QCDInput, entry: str) -> np.ndarray:
    """Exact-mass LO transverse AA, AB, BB densities.

    The common normalization is fixed by the standard axial density.  A direct
    cut-trace derivation gives

      rho_AA = sqrt(lambda)/(8 pi^2) (s-S^2)(d^2+2s)/s^2
      rho_AB = sqrt(lambda)/(8 pi^2) 3d(s-S^2)/(S s)
      rho_BB = sqrt(lambda)/(8 pi^2) (s-S^2)(s+2d^2)/(S^2 s).
    """
    s = np.asarray(s, dtype=float)
    m, ms = inp.mQ, inp.ms
    big_s, d = m + ms, m - ms
    lam = np.maximum(lambda_kallen(s, m * m, ms * ms), 0.0)
    root = np.sqrt(lam)
    pref = root / (8.0 * pi**2)
    if entry == "AA":
        kernel = (s - big_s**2) * (d * d + 2.0 * s) / s**2
    elif entry in {"AB", "BA"}:
        kernel = 3.0 * d * (s - big_s**2) / (big_s * s)
    elif entry == "BB":
        kernel = (s - big_s**2) * (s + 2.0 * d * d) / (big_s**2 * s)
    else:
        raise ValueError(entry)
    return pref * kernel


def rho_vector(s: np.ndarray | float, inp: QCDInput) -> np.ndarray:
    """Exact-mass LO transverse vector-current density."""
    s = np.asarray(s, dtype=float)
    m, ms = inp.mQ, inp.ms
    lam = np.maximum(lambda_kallen(s, m * m, ms * ms), 0.0)
    root = np.sqrt(lam)
    bracket = (
        2.0
        - (m * m + ms * ms - 6.0 * m * ms) / s
        - (m * m - ms * ms) ** 2 / s**2
    )
    return root * bracket / (8.0 * pi**2)


def local_twopoint_matrix(m2: float, inp: QCDInput) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the explicit d=3, finite-ms, and d=5 Borel matrices.

    The d=3 entries and the d=5 entries follow from direct rest-frame Dirac
    traces before Borel transformation.  The finite-ms d=3 Wilson correction
    is the EOM-consistent outer-product extension of the standard AA result.
    """
    m, ms, big_s = inp.mQ, inp.ms, inp.mQ + inp.ms
    e = exp(-m * m / m2)
    d3 = inp.ss * e * np.array(
        [[m, m * m / big_s], [m * m / big_s, m**3 / big_s**2]], dtype=float
    )
    corr_aa = inp.ss * e * (
        -(m * m * ms) / (2.0 * m2) + (m**3 * ms * ms) / (2.0 * m2**2)
    )
    v = np.array([1.0, m / big_s])
    d3_ms = corr_aa * np.outer(v, v)
    mixed = inp.mixed_ss
    d5 = e * np.array(
        [
            [
                -mixed * m**3 / (4.0 * m2**2),
                -mixed
                / (4.0 * big_s)
                * (m**4 / m2**2 - m * m / m2 - 1.0),
            ],
            [
                -mixed
                / (4.0 * big_s)
                * (m**4 / m2**2 - m * m / m2 - 1.0),
                mixed
                / big_s**2
                * (-m**5 / (4.0 * m2**2) + m**3 / (2.0 * m2)),
            ],
        ],
        dtype=float,
    )
    return d3, d3_ms, d5


def borel_twopoint_matrix(m2: float, s0: float, inp: QCDInput) -> dict[str, np.ndarray]:
    sth = (inp.mQ + inp.ms) ** 2
    pert = np.empty((2, 2), dtype=float)
    for a, entry in [((0, 0), "AA"), ((0, 1), "AB"), ((1, 0), "AB"), ((1, 1), "BB")]:
        pert[a] = integrate_gauss(
            lambda s, e=entry: np.exp(-s / m2) * rho_twopoint(s, inp, e),
            sth,
            s0,
        )
    d3, d3_ms, d5 = local_twopoint_matrix(m2, inp)
    return {"pert": pert, "d3": d3, "d3_ms": d3_ms, "d5": d5, "total": pert + d3 + d3_ms + d5}


def borel_vector(m2: float, s0: float, inp: QCDInput) -> dict[str, float]:
    sth = (inp.mQ + inp.ms) ** 2
    pert = integrate_gauss(
        lambda s: np.exp(-s / m2) * rho_vector(s, inp), sth, s0
    )
    # Parity flips the condensate Wilson coefficients relative to AA.
    d3, d3_ms, d5 = local_twopoint_matrix(m2, inp)
    local = -(d3[0, 0] + d3_ms[0, 0] + d5[0, 0])
    return {"pert": pert, "local": local, "total": pert + local}


def rotation(theta_deg: float) -> np.ndarray:
    th = np.deg2rad(theta_deg)
    return np.array([[np.sin(th), np.cos(th)], [np.cos(th), -np.sin(th)]])


def projected_twopoint(
    m2: float, s0: float, inp: QCDInput, theta_deg: float
) -> dict[str, np.ndarray]:
    raw = borel_twopoint_matrix(m2, s0, inp)
    r = rotation(theta_deg)
    return {key: r @ value @ r.T for key, value in raw.items()}


def _mass_from_function(func: Callable[[float], float], m2: float) -> float:
    tau = 1.0 / m2
    h = max(1.0e-5, tau * 2.0e-4)
    p_plus = func(1.0 / (tau + h))
    p_minus = func(1.0 / (tau - h))
    p0 = func(m2)
    r1 = -(p_plus - p_minus) / (2.0 * h)
    return sqrt(max(r1 / p0, 0.0)) if p0 > 0.0 else float("nan")


def initial_observables(
    m2: float,
    s0: float,
    inp: QCDInput,
    theta_deg: float,
    state: int,
    physical_mass: float,
) -> dict[str, float]:
    def pi_at(x):
        return float(projected_twopoint(x, s0, inp, theta_deg)["total"][state, state])

    pieces = projected_twopoint(m2, s0, inp, theta_deg)
    total = float(pieces["total"][state, state])
    mass_sr = _mass_from_function(pi_at, m2)
    f = exp(physical_mass**2 / (2.0 * m2)) * sqrt(max(total, 0.0)) / physical_mass
    sth = (inp.mQ + inp.ms) ** 2
    smax = sth + 60.0 * m2
    full_pert = np.empty((2, 2))
    for a, entry in [((0, 0), "AA"), ((0, 1), "AB"), ((1, 0), "AB"), ((1, 1), "BB")]:
        full_pert[a] = integrate_gauss(
            lambda s, e=entry: np.exp(-s / m2) * rho_twopoint(s, inp, e),
            sth,
            smax,
            192,
        )
    r = rotation(theta_deg)
    local = pieces["d3"] + pieces["d3_ms"] + pieces["d5"]
    full = r @ full_pert @ r.T + local
    pc = total / float(full[state, state]) if full[state, state] != 0 else float("nan")
    d5frac = abs(float(pieces["d5"][state, state])) / max(abs(total), 1e-30)
    offdiag = abs(float(pieces["total"][0, 1])) / sqrt(
        max(abs(float(pieces["total"][0, 0] * pieces["total"][1, 1])), 1e-30)
    )
    return {
        "Pi": total,
        "mass_sr": mass_sr,
        "f": f,
        "pc": pc,
        "d5frac": d5frac,
        "offdiag": offdiag,
    }


def vector_observables(
    m2: float, s0: float, inp: QCDInput, physical_mass: float
) -> dict[str, float]:
    def pi_at(x):
        return borel_vector(x, s0, inp)["total"]

    pieces = borel_vector(m2, s0, inp)
    total = pieces["total"]
    mass_sr = _mass_from_function(pi_at, m2)
    f = exp(physical_mass**2 / (2.0 * m2)) * sqrt(max(total, 0.0)) / physical_mass
    sth = (inp.mQ + inp.ms) ** 2
    full_pert = integrate_gauss(
        lambda s: np.exp(-s / m2) * rho_vector(s, inp),
        sth,
        sth + 60.0 * m2,
        192,
    )
    full = full_pert + pieces["local"]
    return {
        "Pi": total,
        "mass_sr": mass_sr,
        "f": f,
        "pc": total / full if full != 0 else float("nan"),
        "localfrac": abs(pieces["local"]) / max(abs(total), 1e-30),
    }


# Photon distribution amplitudes at mu=1 GeV.
def phi_gamma(u):
    u = np.asarray(u)
    return 6.0 * u * (1.0 - u)


def gegenbauer_c2_half(x):
    return 0.5 * (3.0 * x * x - 1.0)


def h_gamma(u, inp: QCDInput):
    x = 2.0 * np.asarray(u) - 1.0
    return -10.0 * gegenbauer_c2_half(x)


def psi_v(u, inp: QCDInput):
    x = 2.0 * np.asarray(u) - 1.0
    return 5.0 * (3.0 * x * x - 1.0) + (3.0 / 64.0) * (
        15.0 * inp.omega_v - 5.0 * inp.omega_a
    ) * (3.0 - 30.0 * x * x + 35.0 * x**4)


def psi_a(u, inp: QCDInput):
    x = 2.0 * np.asarray(u) - 1.0
    return (
        (1.0 - x * x)
        * (5.0 * x * x - 1.0)
        * 2.5
        * (1.0 + 9.0 * inp.omega_v / 16.0 - 3.0 * inp.omega_a / 16.0)
    )


def dpsi_a(u, inp: QCDInput):
    h = 1.0e-5
    return float((psi_a(u + h, inp) - psi_a(u - h, inp)) / (2.0 * h))


def A_twist4(u, inp: QCDInput):
    u = np.asarray(u)
    ub = 1.0 - u
    eps = 1.0e-15
    logu = np.log(np.maximum(u, eps))
    logub = np.log(np.maximum(ub, eps))
    first = 40.0 * u**2 * ub**2 * (3.0 * inp.kappa_da + 1.0)
    second = -24.0 * inp.zeta2 * (
        u * ub * (2.0 + 13.0 * u * ub)
        + 2.0 * u**3 * (10.0 - 15.0 * u + 6.0 * u**2) * logu
        + 2.0 * ub**3 * (10.0 - 15.0 * ub + 6.0 * ub**2) * logub
    )
    return first + second


def cumulative(func: Callable[[np.ndarray], np.ndarray], u: float, n: int = 96) -> float:
    return integrate_gauss(func, 0.0, u, n)


def H_gamma(u: float, inp: QCDInput) -> float:
    return cumulative(lambda x: h_gamma(x, inp), u)


def Hbar_gamma(u: float, inp: QCDInput) -> float:
    # Swap the order of integration: int_0^u (u-t) h(t) dt.
    return integrate_gauss(lambda x: (u - x) * h_gamma(x, inp), 0.0, u, 96)


def Psi_v(u: float, inp: QCDInput) -> float:
    return cumulative(lambda x: psi_v(x, inp), u)


def da_three(name: str, aq, ab, ag, inp: QCDInput):
    k, z1, z2 = inp.kappa_da, inp.zeta1, inp.zeta2
    if name == "V":
        return 540.0 * inp.omega_v * (aq - ab) * aq * ab * ag**2
    if name == "A":
        return 360.0 * aq * ab * ag**2 * (
            1.0 + 0.5 * inp.omega_a * (7.0 * ag - 3.0)
        )
    common = 3.0 * (ab - aq) ** 2 - ag * (1.0 - ag)
    if name == "S":
        return 30.0 * ag**2 * (
            k * (1.0 - ag)
            + z1 * (1.0 - ag) * (1.0 - 2.0 * ag)
            + z2 * common
        )
    if name == "St":
        return -30.0 * ag**2 * (
            k * (1.0 - ag)
            + z1 * (1.0 - ag) * (1.0 - 2.0 * ag)
            + z2 * common
        )
    if name == "T1":
        return -360.0 * z2 * (ab - aq) * aq * ab * ag
    if name == "T2":
        return 30.0 * ag**2 * (ab - aq) * (
            k + z1 * (1.0 - 2.0 * ag) + z2 * (3.0 - 4.0 * ag)
        )
    if name == "T3":
        return -360.0 * z2 * (ab - aq) * aq * ab * ag
    if name == "T4":
        return 30.0 * ag**2 * (ab - aq) * (
            k + z1 * (1.0 - 2.0 * ag) + z2 * (3.0 - 4.0 * ag)
        )
    raise ValueError(name)


def F2(aq, ab, ag, inp: QCDInput):
    return (
        da_three("S", aq, ab, ag, inp)
        + da_three("St", aq, ab, ag, inp)
        + da_three("T1", aq, ab, ag, inp)
        - da_three("T2", aq, ab, ag, inp)
        - da_three("T3", aq, ab, ag, inp)
        + da_three("T4", aq, ab, ag, inp)
    )


def F3(aq, ab, ag, inp: QCDInput):
    return da_three("A", aq, ab, ag, inp) + da_three("V", aq, ab, ag, inp)


def _nested_gauss(
    outer_lo: float,
    outer_hi: float,
    inner_bounds: Callable[[float], tuple[float, float]],
    func: Callable[[float, np.ndarray], np.ndarray],
    n: int = 24,
) -> float:
    xo, wo = _legendre(n)
    vs = 0.5 * (outer_hi - outer_lo) * xo + 0.5 * (outer_hi + outer_lo)
    total = 0.0
    xi, wi = _legendre(n)
    for v, weight in zip(vs, wo):
        lo, hi = inner_bounds(float(v))
        if hi <= lo:
            continue
        ag = 0.5 * (hi - lo) * xi + 0.5 * (hi + lo)
        total += weight * 0.5 * (hi - lo) * float(np.sum(wi * func(float(v), ag)))
    return 0.5 * (outer_hi - outer_lo) * total


def three_particle_integrals(u0: float, inp: QCDInput, n: int = 24) -> tuple[float, float]:
    i_f2a = _nested_gauss(
        0.0,
        1.0 - u0,
        lambda v: (0.0, u0 / (1.0 - v)),
        lambda v, ag: F2(u0 - (1.0 - v) * ag, 1.0 - u0 - v * ag, ag, inp),
        n,
    )
    i_f2b = _nested_gauss(
        1.0 - u0,
        1.0,
        lambda v: (0.0, (1.0 - u0) / v),
        lambda v, ag: F2(u0 - (1.0 - v) * ag, 1.0 - u0 - v * ag, ag, inp),
        n,
    )
    # The apparent endpoint singularities cancel in the bracket.  A mapped
    # Gauss rule avoids evaluating the endpoint itself.
    def f3_for_ab(ab_arr):
        vals = []
        xi, wi = _legendre(n)
        for ab in np.atleast_1d(ab_arr):
            lo, hi = u0 - ab, 1.0 - ab
            ag = 0.5 * (hi - lo) * xi + 0.5 * (hi + lo)
            first = 0.5 * (hi - lo) * np.sum(
                wi * F3(1.0 - ab - ag, ab, ag, inp) / ag**2
            )
            boundary = F3(1.0 - u0, ab, u0 - ab, inp) / (u0 - ab)
            vals.append(first - boundary)
        return np.asarray(vals)

    eps = 1.0e-9
    i_f3 = integrate_gauss(f3_for_ab, 0.0, u0 - eps, n)
    return i_f2a + i_f2b, i_f3


def rho_transition(s: np.ndarray | float, inp: QCDInput) -> np.ndarray:
    """Hard-emission spectral density for the leading epsilon structure."""
    s = np.asarray(s, dtype=float)
    m, ms = inp.mQ, inp.ms
    root = np.sqrt(np.maximum(lambda_kallen(s, m * m, ms * ms), 0.0))
    a_s = s - m * m + ms * ms
    a_q = s - ms * ms + m * m
    tiny = np.finfo(float).tiny
    l_s = np.log(np.maximum((a_s - root) / (a_s + root), tiny))
    l_q = np.log(np.maximum((a_q - root) / (a_q + root), tiny))
    return 3.0 * ms * m * (inp.es * l_s + inp.eQ * l_q) / (4.0 * pi**2)


def transition_terms(
    m2: float,
    s0: float,
    inp: QCDInput,
    u0: float = 0.5,
    three_n: int = 24,
) -> dict[str, float]:
    """Post-double-Borel leading invariant, term by term.

    `local_heavy_diagnostic` reproduces the historical local spectator
    condensate term, but `quoted_total_A` deliberately excludes it.
    """
    m, ms = inp.mQ, inp.ms
    sth = (m + ms) ** 2
    em = exp(-m * m / m2)
    es0 = exp(-s0 / m2)
    diff = em - es0
    hard = integrate_gauss(
        lambda s: np.exp(-s / m2) * rho_transition(s, inp), sth, s0, 160
    )
    local_heavy = (
        inp.eQ
        * m
        * em
        * inp.ss
        * (1.0 - ms * ms / m2 * (1.0 - m * m / m2))
    )
    tw2 = inp.es * m * inp.ss * diff * m2 * inp.chi * float(phi_gamma(u0))
    tw4_2p = inp.es * m * inp.ss * em * (
        -0.25 * m * m / m2 * float(A_twist4(u0, inp))
        - H_gamma(u0, inp) * (1.0 - u0)
        - Hbar_gamma(u0, inp) * (1.0 - 2.0 * m * m / m2)
    )
    tw3_2p = inp.es * inp.f3 * m2 * diff * (
        0.25 * (1.0 - u0) * dpsi_a(u0, inp)
        - 0.25 * float(psi_a(u0, inp))
        - Psi_v(u0, inp) * (1.0 + 2.0 * m * m / m2)
        + (1.0 - u0) * float(psi_v(u0, inp))
    )
    i_f2, i_f3 = three_particle_integrals(u0, inp, three_n)
    tw4_3p = m * inp.es * inp.ss * em * i_f2
    tw3_3p = -inp.es * inp.f3 * m2 * diff * i_f3
    quoted = hard + tw2 + tw4_2p + tw3_2p + tw4_3p + tw3_3p
    r_tensor = m / (m + ms)
    return {
        "hard": hard,
        "tw2": tw2,
        "tw4_2p": tw4_2p,
        "tw3_2p": tw3_2p,
        "tw4_3p": tw4_3p,
        "tw3_3p": tw3_3p,
        "local_heavy_diagnostic": local_heavy,
        "quoted_total_A": quoted,
        "quoted_total_B": r_tensor * quoted,
        "tensor_ratio": r_tensor,
    }


def physical_transition(
    m2: float,
    s0: float,
    inp: QCDInput,
    theta_deg: float,
    state: int,
    mi: float,
    mf: float,
    fi: float,
    fv: float,
    tensor_delta: float = 0.0,
    three_n: int = 24,
) -> dict[str, float]:
    terms = transition_terms(m2, s0, inp, 0.5, three_n)
    ta = terms["quoted_total_A"]
    tb = terms["quoted_total_B"] * (1.0 + tensor_delta)
    r = rotation(theta_deg)
    comp_a = r[state, 0] * ta
    comp_b = r[state, 1] * tb
    tphys = comp_a + comp_b
    pref = exp((mi * mi + mf * mf) / (2.0 * m2)) / (fi * fv * mi * mf)
    g = pref * tphys
    k = (mi * mi - mf * mf) / (2.0 * mi)
    width_gev = (
        ALPHA_EM
        * g
        * g
        * k**3
        * (mi * mi + mf * mf)
        / (3.0 * mi * mi * mf * mf)
    )
    return {
        **terms,
        "component_A": pref * comp_a,
        "component_B": pref * comp_b,
        "g": g,
        "width_keV": width_gev * 1.0e6,
        "k_gamma": k,
    }


def with_updates(inp: QCDInput, **kwargs) -> QCDInput:
    return replace(inp, **kwargs)

