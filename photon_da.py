"""
photon_da.py  --  Photon distribution amplitudes (DAs) for the LCSR
                  calculation of Ds1 -> Ds gamma (and the Bs counterparts).

WHAT THIS MODULE IS
-------------------
When a quasi-real photon is emitted at long distance, we cannot treat it as a
simple QED vertex on a quark line. Instead its coupling to the soft q-qbar pair
is described by *photon distribution amplitudes* -- universal non-perturbative
functions that encode how the photon's momentum is shared between the quark and
antiquark, organized by TWIST (a measure of the operator's dimension minus spin;
higher twist = more suppressed).

We use the complete set through twist-4 (2- and 3-particle), following:
  [BBK]  Ball, Braun, Kivel, Nucl.Phys.B649(2003)263  [hep-ph/0207307]
  [Rohr] Rohrwild, JHEP0709(2007)073  [arXiv:0708.1405]  -- Appendix A
The definitions and every numerical parameter below are taken from Rohrwild's
Appendix A and Table 1 (renormalization scale mu = 1 GeV).

CONVENTIONS
-----------
  u        = momentum fraction carried by the quark; ubar = 1-u, u in [0,1].
  Two-particle DAs are functions of u.
  Three-particle DAs are functions of alpha = (alpha_q, alpha_qbar, alpha_g)
  with alpha_q + alpha_qbar + alpha_g = 1 (momentum fractions of q, qbar, gluon).
  All DAs are real. Normalizations are fixed so that the leading-twist phi
  integrates to 1:  int_0^1 phi(u) du = 1.

STUDENT CHECK POINTS are marked with '# CHECK:'.
"""

import numpy as np

# ============================================================================
#  NUMERICAL PARAMETERS  (Rohrwild Table 1, mu = 1 GeV)
#  Each has a central value and an uncertainty used later in the error budget.
# ============================================================================
PARAMS = {
    # magnetic susceptibility of the quark condensate (the KEY parameter):
    "chi":     (3.15,   0.30),    # GeV^-2   [also compare Rohrwild's fit 2.85+-0.5]
    # twist-4 two-particle shape parameters:
    "kappa":   (0.15,   0.0),     # dimensionless
    "kappa_p": (-0.05,  0.0),     # "kappa^+"
    "zeta1":   (0.40,   0.0),
    "zeta1_p": (0.0,    0.0),     # "zeta_1^+"
    "zeta2":   (0.30,   0.0),
    "zeta2_p": (0.0,    0.0),     # "zeta_2^+"
    # twist-3 three-particle coupling and shape parameters:
    "f3g":     (-4.0e-3, 2.0e-3), # GeV^2   "f_{3 gamma}"
    "omegaA":  (-2.1,    1.0),    # "omega^A_gamma"
    "omegaV":  (3.8,     1.8),    # "omega^V_gamma"
    # quark condensate (enters as an overall e_q<qbar q> in the matrix elements):
    "qq":      (-(0.240)**3, 3*0.240**2*0.010),  # GeV^3, from -(240+-10 MeV)^3
}

def val(name):
    """Central value of a parameter."""
    return PARAMS[name][0]

# ============================================================================
#  TWIST-2  (leading)  --  phi_gamma(u)
#  Rohrwild Eq.(38): asymptotic form.  Enters with the factor chi*<qbar q>.
# ============================================================================
def phi_gamma(u):
    """Leading-twist photon DA, asymptotic form phi(u)=6 u (1-u)."""
    u = np.asarray(u, dtype=float)
    return 6.0 * u * (1.0 - u)
# CHECK: int_0^1 6u(1-u) du = 6*(1/2 - 1/3) = 1.  (normalization)

# ============================================================================
#  COLANGELO/BBK CONVENTION HELPERS
#  The Colangelo-De Fazio-Ozpineci g1 sum rule uses the Appendix convention
#  printed in hep-ph/0505195.  It differs from the Rohrwild-labelled helpers
#  below for psi^v and A(u), so the numerical Stage-1 benchmark must call these
#  explicitly named functions.
# ============================================================================
def gegenbauer_c2_half(x):
    """C_2^{1/2}(x) = P_2(x), used in h_gamma."""
    x = np.asarray(x, dtype=float)
    return 0.5 * (3.0*x**2 - 1.0)

def A_t4_colangelo(u):
    """Twist-4 two-particle DA mathbb{A}(u), Colangelo Appendix / BBK."""
    u = np.asarray(u, dtype=float)
    ub = 1.0 - u
    k, kp = val("kappa"), val("kappa_p")
    z2, z2p = val("zeta2"), val("zeta2_p")
    lnu = np.log(np.clip(u, 1e-12, 1.0))
    lnub = np.log(np.clip(ub, 1e-12, 1.0))
    bracket = (
        u*ub*(2.0 + 13.0*u*ub)
        + 2.0*u**3*(10.0 - 15.0*u + 6.0*u**2)*lnu
        + 2.0*ub**3*(10.0 - 15.0*ub + 6.0*ub**2)*lnub
    )
    return 40.0*u**2*ub**2*(3.0*k - kp + 1.0) + 8.0*(z2p - 3.0*z2)*bracket

def h_gamma_colangelo(u):
    """Twist-4 h_gamma(u), Colangelo Appendix / BBK."""
    u = np.asarray(u, dtype=float)
    t = 2.0*u - 1.0
    return -10.0*(1.0 + 2.0*val("kappa_p"))*gegenbauer_c2_half(t)

def psi_v_colangelo(u):
    """Twist-3 vector DA psi^v(u), Colangelo Appendix / BBK."""
    u = np.asarray(u, dtype=float)
    t = 2.0*u - 1.0
    wA, wV = val("omegaA"), val("omegaV")
    return (
        5.0*(3.0*t**2 - 1.0)
        + (3.0/64.0)*(15.0*wV - 5.0*wA)
        *(3.0 - 30.0*t**2 + 35.0*t**4)
    )

# ============================================================================
#  TWIST-3  (2-particle)  --  psi^v(u), psi^a(u)
#  Rohrwild Eqs.(63),(64). These come with the coupling f3gamma.
# ============================================================================
def psi_v(u):
    """Twist-3 vector DA psi^(V)(u), Rohrwild Eq.(63)."""
    u = np.asarray(u, dtype=float)
    wA, wV = val("omegaA"), val("omegaV")
    t = 2.0*u - 1.0                      # = u - ubar
    return (-20.0*u*(1-u)*t
            + (15.0/16.0)*(wA - 3.0*wV)*u*(1-u)*t*(7.0*t**2 - 3.0))

def psi_a(u):
    """Twist-3 axial DA psi^(A)(u), Rohrwild Eq.(64)."""
    u = np.asarray(u, dtype=float)
    wA, wV = val("omegaA"), val("omegaV")
    t = 2.0*u - 1.0
    return ((1.0 - t**2)*(5.0*t**2 - 1.0)
            * 2.5 * (1.0 + (9.0/16.0)*wV - (3.0/16.0)*wA))
# CHECK: psi_a is symmetric under u<->1-u (t->-t), psi_v is antisymmetric.

# ============================================================================
#  TWIST-4  (2-particle)  --  A(u), B(u)
#  Rohrwild Eqs.(39),(40).  A(u) enters leading among tw-4; B(u) is a
#  derivative-type structure. Named *_t4 to avoid clashing with the
#  three-particle A(alpha) below.
# ============================================================================
def A_t4(u):
    """Twist-4 two-particle DA  A(u), Rohrwild Eq.(39)."""
    u = np.asarray(u, dtype=float)
    ub = 1.0 - u
    k, kp = val("kappa"), val("kappa_p")
    z2, z2p = val("zeta2"), val("zeta2_p")
    # guard the logs at endpoints
    lnu  = np.log(np.clip(u,  1e-12, 1.0))
    lnub = np.log(np.clip(ub, 1e-12, 1.0))
    bracket = ( u*ub*(2.0 + 13.0*u*ub)
                + u**3*(10.0 - 15.0*u + 6.0*u**2)*lnu
                + ub**3*(10.0 - 15.0*ub + 6.0*ub**2)*lnub )
    return 40.0*u*ub*(3.0*k - kp + 1.0) + 8.0*(z2p - 3.0*z2)*bracket

def B_t4(u):
    """Twist-4 two-particle DA  B(u), Rohrwild Eq.(40) (integral form done
    analytically). B(u) = 40 * int_0^u dalpha (u-alpha)(1+3kappa^+)
                              [ -1/2 + (3/2)(2alpha-1)^2 ]."""
    u = np.asarray(u, dtype=float)
    kp = val("kappa_p")
    # int_0^u (u-a)[-1/2 + (3/2)(2a-1)^2] da  -- closed form:
    #   let f(a) = -1/2 + (3/2)(2a-1)^2 = -1/2 + (3/2)(4a^2-4a+1) = 6a^2-6a+1
    #   int_0^u (u-a)(6a^2-6a+1) da
    #     = u*int_0^u(6a^2-6a+1)da - int_0^u a(6a^2-6a+1)da
    #     = u*(2u^3-3u^2+u) - (1.5u^4-2u^3+0.5u^2)
    inner = u*(2.0*u**3 - 3.0*u**2 + u) - (1.5*u**4 - 2.0*u**3 + 0.5*u**2)
    return 40.0*(1.0 + 3.0*kp)*inner
# CHECK: B(0)=0 by construction (lower limit of the integral).

# ============================================================================
#  TWIST-4  (3-particle)  --  S, S~ , S_gamma, T1..T4, T4^gamma
#  Rohrwild Eqs.(50),(51),(52),(54)-(58).  Functions of (aq, aqb, ag).
# ============================================================================
def _ag(aq, aqb):     # gluon fraction
    return 1.0 - aq - aqb

def S_3p(aq, aqb):
    """Twist-4 3-particle S(alpha), Rohrwild Eq.(50)."""
    ag = _ag(aq, aqb)
    k, kp = val("kappa"), val("kappa_p")
    z1, z1p = val("zeta1"), val("zeta1_p")
    z2 = val("zeta2")
    return 30.0*ag**2*( (k + kp)*(1-ag)
                        + (z1 + z1p)*(1-ag)*(1-2*ag)
                        + z2*(3.0*(aq-aqb)**2 - ag*(1-ag)) )

def St_3p(aq, aqb):
    """Twist-4 3-particle S-tilde(alpha), Rohrwild Eq.(51)."""
    ag = _ag(aq, aqb)
    k, kp = val("kappa"), val("kappa_p")
    z1, z1p = val("zeta1"), val("zeta1_p")
    z2 = val("zeta2")
    return -30.0*ag**2*( (k - kp)*(1-ag)
                         + (z1 - z1p)*(1-ag)*(1-2*ag)
                         + z2*(3.0*(aq-aqb)**2 - ag*(1-ag)) )

def Sg_3p(aq, aqb):
    """Twist-4 3-particle S_gamma(alpha), Rohrwild Eq.(52)."""
    s = aq + aqb
    ag = _ag(aq, aqb)
    return 60.0*ag**2*s*(4.0 - 7.0*s)

def T1_3p(aq, aqb):
    """Rohrwild Eq.(54)."""
    ag = _ag(aq, aqb)
    z2, z2p = val("zeta2"), val("zeta2_p")
    return -120.0*(3.0*z2 + z2p)*(aq-aqb)*aq*aqb*ag

def T2_3p(aq, aqb):
    """Rohrwild Eq.(55)."""
    ag = _ag(aq, aqb)
    k, kp = val("kappa"), val("kappa_p")
    z1, z1p = val("zeta1"), val("zeta1_p")
    z2 = val("zeta2")
    return 30.0*ag**2*(aq-aqb)*( (k - kp)
                                 + (z1 - z1p)*(1-2*ag) + z2*(3.0-4.0*ag) )

def T3_3p(aq, aqb):
    """Rohrwild Eq.(56)."""
    ag = _ag(aq, aqb)
    z2, z2p = val("zeta2"), val("zeta2_p")
    return -120.0*(3.0*z2 - z2p)*(aq-aqb)*aq*aqb*ag

def T4_3p(aq, aqb):
    """Rohrwild Eq.(57)."""
    ag = _ag(aq, aqb)
    k, kp = val("kappa"), val("kappa_p")
    z1, z1p = val("zeta1"), val("zeta1_p")
    z2 = val("zeta2")
    return 30.0*ag**2*(aq-aqb)*( (k + kp)
                                 + (z1 + z1p)*(1-2*ag) + z2*(3.0-4.0*ag) )

def T4g_3p(aq, aqb):
    """T4^gamma(alpha), Rohrwild Eq.(58)."""
    s = aq + aqb
    ag = _ag(aq, aqb)
    return 60.0*ag**2*(aq-aqb)*(4.0 - 7.0*s)

# ============================================================================
#  TWIST-3  (3-particle)  --  V(alpha), A(alpha)
#  Rohrwild Eqs.(65),(66). These come with the coupling f3gamma.
# ============================================================================
def V_3p(aq, aqb):
    """Twist-3 3-particle V(alpha), Rohrwild Eq.(65)."""
    ag = _ag(aq, aqb)
    return 540.0*val("omegaV")*(aq-aqb)*aq*aqb*ag**2

def A_3p(aq, aqb):
    """Twist-3 3-particle A(alpha), Rohrwild Eq.(66)."""
    ag = _ag(aq, aqb)
    return 360.0*aq*aqb*ag**2*(1.0 + 0.5*val("omegaA")*(7.0*ag - 3.0))


# ============================================================================
#  SELF-TESTS  (run:  python photon_da.py)
#  These verify the normalizations and symmetries a student should check.
# ============================================================================
if __name__ == "__main__":
    print("=== Photon DA module self-tests (mu = 1 GeV) ===\n")

    def quad01(f, n=300):
        x, w = np.polynomial.legendre.leggauss(n)
        u = 0.5*(x + 1.0)
        return 0.5*np.sum(w*f(u))

    def simplex_int(f, n=80):
        # Map aq in [0,1], t in [0,1] to aqb=(1-aq)t.
        # Jacobian is (1-aq).
        x, wx = np.polynomial.legendre.leggauss(n)
        aq = 0.5*(x + 1.0)
        wa = 0.5*wx
        t = aq.copy()
        wt = wa.copy()
        total = 0.0
        for ai, wai in zip(aq, wa):
            aqb = (1.0 - ai)*t
            total += wai*np.sum(wt*(1.0 - ai)*f(ai, aqb))
        return float(total)

    # (1) leading-twist normalization
    n_phi = quad01(phi_gamma)
    print(f"int phi_gamma du            = {n_phi:.6f}   (expect 1)")

    # (2) twist-3 symmetry properties
    us = np.linspace(0,1,11)
    sym_a  = np.allclose(psi_a(us), psi_a(1-us))
    anti_v = np.allclose(psi_v(us), -psi_v(1-us))
    print(f"psi_a symmetric u<->1-u    = {sym_a}   (expect True)")
    print(f"psi_v antisymmetric        = {anti_v}   (expect True)")

    # (3) twist-4 B(0)=0
    print(f"B_t4(0)                    = {float(B_t4(0.0)):.3e}   (expect 0)")

    # (4) three-particle support/normalization spot-checks
    # integrate S_3p over the simplex aq,aqb>0, aq+aqb<1
    print(f"int_simplex S_3p           = {simplex_int(S_3p):.4e}")
    print(f"int_simplex A_3p           = {simplex_int(A_3p):.4e}")
    print(f"int_simplex V_3p           = {simplex_int(V_3p):.4e}  (expect ~0, odd in aq-aqb)")

    print("\nParameters used:")
    for k,(v,e) in PARAMS.items():
        print(f"  {k:9s} = {v:+.5g}  +- {e:.3g}")
