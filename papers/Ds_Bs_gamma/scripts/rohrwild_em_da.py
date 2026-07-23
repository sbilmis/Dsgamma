"""Rohrwild electromagnetic three-particle photon-DA contribution.

This module contains the ``S_gamma,T4_gamma`` sector generated when the
external photon is inserted in the heavy-quark propagator.  It is the
nonlocal background-field replacement for the standalone heavy-charge local
condensate term used in the Colangelo benchmark.  The two prescriptions are
alternatives and must never be summed.

All formulas below use the symmetric double-Borel choice ``u0 = 1/2`` and the
Rohrwild source convention

    S_gamma = 60 alpha_g^2 (alpha_q + alpha_qbar)
              (4 - 7 (alpha_q + alpha_qbar)),
    T4_gamma = 60 alpha_g^2 (alpha_qbar - alpha_q)
               (4 - 7 (alpha_q + alpha_qbar)).

The direct integrals and the P-functional reduce to exact constants.  In the
published appendix the third P-functional integral is printed with a zero
lower limit.  Taken literally that integral is logarithmically divergent.
The double-Borel support condition is alpha_qbar >= u0; with that support the
finite result is

    P[T4_gamma] = log(8) - 1181/320.

This support correction is deliberately exposed as metadata in every result.
It should also be repeated explicitly in the manuscript.
"""

from __future__ import annotations

import math


U0 = 0.5
S_DIRECT = 25.0 / 16.0
T4_DIRECT_VECTOR = 5.0 / 16.0
P_T4_COMPONENTS = (
    -5.0 / 64.0,
    -3.0 / 2.0,
    -135.0 / 64.0 + math.log(8.0),
    -1.0 / 320.0,
)
P_T4_BRACE = sum(P_T4_COMPONENTS)
SUPPORT_NOTE = (
    "P[T4gamma] uses the double-Borel support alpha_qbar>=1/2 in the "
    "third appendix integral; the printed zero lower limit diverges"
)


def continuum_factor(M2: float, s0: float, mQ: float) -> float:
    """Symmetric double-Borel continuum factor E = exp(-mQ^2/M2)-exp(-s0/M2)."""
    return math.exp(-(mQ * mQ) / M2) - math.exp(-s0 / M2)


def axial_em_qcd(M2: float, s0: float, inputs: dict[str, float]) -> dict[str, float | str]:
    """Return the nonlocal electromagnetic-DA OPE term for current J_A.

    The current-specific trace gives ``(1-2u) S_gamma - T4_gamma`` for
    the direct term.  At u0=1/2 both direct integrals vanish by antisymmetry,
    and only the derivative P-functional survives.
    """
    mQ = inputs["mc"]
    mA = inputs["m_ds1"]
    mP = inputs["m_ds"]
    delta_exp = continuum_factor(M2, s0, mQ)
    direct_s = 0.0
    direct_t4 = 0.0
    p_kernel = 2.0 * (mP * mP - mA * mA) / M2 * P_T4_BRACE
    bracket = direct_s + direct_t4 + p_kernel
    qcd = inputs["ec"] * inputs["ss"] * delta_exp * bracket
    return {
        "em_Sgamma_direct_kernel": direct_s,
        "em_T4gamma_direct_kernel": direct_t4,
        "em_T4gamma_P_kernel": p_kernel,
        "em_kernel_total": bracket,
        "em_continuum_factor": delta_exp,
        "em_qcd": qcd,
        "em_support_note": SUPPORT_NOTE,
    }


def tensor_em_qcd(
    M2: float,
    s0: float,
    inputs: dict[str, float],
    *,
    u0: float = U0,
) -> dict[str, float | str]:
    """Return the nonlocal electromagnetic-DA OPE term for current J_B.

    The tensor numerator is reduced with the same physical-residue
    prescription used for the other J_B x-dependent kernels in this project:
    p^2=m_P^2, p.q=(m_A^2-m_P^2)/2 at the u0 Borel saddle.  For u0=1/2 the
    resulting dimensionless kernel is

      (r0-rS) 25/16 - rT 5/16
      + 2 r0 (m_P^2-m_A^2)/M2 P[T4gamma].

    This is a specified numerator prescription, not an independent local
    condensate insertion.
    """
    if not math.isclose(u0, U0, rel_tol=0.0, abs_tol=1.0e-14):
        raise NotImplementedError("tensor_em_qcd is derived only for u0=1/2")

    mQ = inputs["mc"]
    mq = inputs["ms"]
    mA = inputs["m_ds1"]
    mP = inputs["m_ds"]
    pq = (mA * mA - mP * mP) / 2.0
    denom = mQ * (mQ + mq)
    r0 = mQ / (mQ + mq)
    rS = (mQ * mQ + (1.0 - u0) * pq) / denom
    rT = (-mQ * mQ + mP * mP + (1.0 + u0) * pq) / denom

    direct_s = (r0 - rS) * S_DIRECT
    direct_t4 = -rT * T4_DIRECT_VECTOR
    p_kernel = 2.0 * r0 * (mP * mP - mA * mA) / M2 * P_T4_BRACE
    bracket = direct_s + direct_t4 + p_kernel
    delta_exp = continuum_factor(M2, s0, mQ)
    qcd = inputs["ec"] * inputs["ss"] * delta_exp * bracket
    return {
        "em_r0": r0,
        "em_rS": rS,
        "em_rT": rT,
        "em_Sgamma_direct_kernel": direct_s,
        "em_T4gamma_direct_kernel": direct_t4,
        "em_T4gamma_P_kernel": p_kernel,
        "em_kernel_total": bracket,
        "em_continuum_factor": delta_exp,
        "em_qcd": qcd,
        "em_support_note": SUPPORT_NOTE,
        "em_tensor_numerator_prescription": "p2=mP2,pq=(mA2-mP2)/2,u0=1/2",
    }


def self_check() -> None:
    """Exact-constant and finiteness checks used by validation scripts."""
    assert math.isclose(S_DIRECT, 1.5625, rel_tol=0.0, abs_tol=1.0e-15)
    assert math.isclose(T4_DIRECT_VECTOR, 0.3125, rel_tol=0.0, abs_tol=1.0e-15)
    assert math.isfinite(P_T4_BRACE)
    assert P_T4_BRACE < 0.0
    assert math.isclose(
        P_T4_BRACE,
        math.log(8.0) - 1181.0 / 320.0,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    )


if __name__ == "__main__":
    self_check()
    print(f"S direct = {S_DIRECT:+.12f}")
    print(f"T4 direct (vector reference) = {T4_DIRECT_VECTOR:+.12f}")
    print(f"P[T4gamma] brace = {P_T4_BRACE:+.12f}")
    print("P components = " + ", ".join(f"{x:+.12f}" for x in P_T4_COMPONENTS))
    print(SUPPORT_NOTE)
