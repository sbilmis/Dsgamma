"""Two-point SVZ baseline for the vector decay constant f_Ds*.

We calculate f_{D_s^*} from the transverse vector-current correlator

    Pi_{mu nu}(q) = i int d^4x e^{iqx}
        <0|T{ (sbar gamma_mu c)(x) (cbar gamma_nu s)(0) }|0>

with

    <0|sbar gamma_mu c|D_s^*(q,eps)> = m_Ds* f_Ds* eps_mu.

This is a controlled LO two-point baseline, not a state-of-the-art NNLO
decay-constant extraction.  It is good enough to remove the previous assumed
normalization from the first Ds*gamma LCSR pass.  The benchmark formula follows
the vector-current setup and LO spectral density in Gelhausen et al.,
arXiv:1305.5432, Eq. (10), with a conservative strange-channel window scan.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Dsstar_gamma" / "outputs"
OUT.mkdir(exist_ok=True)


def gl_quad(func, a: float, b: float, n: int = 300) -> float:
    if b <= a:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (b - a) * x + 0.5 * (a + b)
    return float(0.5 * (b - a) * np.sum(w * func(t)))


def rho_T_LO(s, mQ):
    """LO vector spectral density for m_s neglected in the loop.

    rho_T(s) = 1/(8 pi^2) s (1-z)^2 (2+z), z=m_Q^2/s.
    This rho_T already includes the color factor in the convention of
    Gelhausen et al. for the coefficient of -g_{mu nu}.
    """
    s = np.asarray(s, dtype=float)
    z = mQ * mQ / s
    return s * (1.0 - z) ** 2 * (2.0 + z) / (8.0 * math.pi**2)


def f_vector_LO(M2, s0, *, mQ, mV):
    lower = mQ * mQ
    pert = gl_quad(lambda s: np.exp(-s / M2) * rho_T_LO(s, mQ), lower, s0)
    f2 = math.exp(mV * mV / M2) * pert / (mV * mV)
    return math.sqrt(max(f2, 0.0)), pert


def summarize(values):
    arr = np.array(values, dtype=float)
    return {
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "median": float(np.percentile(arr, 50.0)),
        "p16": float(np.percentile(arr, 16.0)),
        "p84": float(np.percentile(arr, 84.0)),
    }


def main():
    inputs = {
        "mQ": 1.27,
        "mV": 2.1122,
    }
    # Rounded thresholds above the vector ground state.  We keep the same
    # simple notation as the main analysis: s0 is quoted directly in GeV^2.
    s0_values = [6.5, 7.0, 7.5]
    M2_values = np.linspace(2.5, 4.0, 7)

    rows = []
    for s0 in s0_values:
        for M2 in M2_values:
            fV, pert = f_vector_LO(float(M2), float(s0), **inputs)
            rows.append(
                {
                    "M2": float(M2),
                    "s0": float(s0),
                    "mQ": inputs["mQ"],
                    "mV": inputs["mV"],
                    "f_Dsstar_GeV": fV,
                    "perturbative_integral": pert,
                }
            )

    csv_path = OUT / "twopoint_f_Dsstar_vector_LO.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    stats = summarize([r["f_Dsstar_GeV"] for r in rows])
    central, central_pert = f_vector_LO(3.25, 7.0, **inputs)
    # Use a conservative uncertainty from the window envelope and missing
    # higher-order/OPE effects; this is intentionally broader than pure
    # numerical stability.
    half_envelope = 0.5 * (stats["max"] - stats["min"])
    systematic = 0.10 * central
    assigned_err = math.sqrt(half_envelope**2 + systematic**2)

    lines = [
        "Two-point vector-current f_Ds* baseline",
        "=======================================",
        "LO transverse vector-current sum rule with rho_T(s)=s(1-z)^2(2+z)/(8 pi^2).",
        "This is a normalization baseline; NNLO/radiative corrections are not included.",
        "Window scanned: M^2=2.5...4.0 GeV^2, s0=6.5,7.0,7.5 GeV^2.",
        f"Window envelope: f_Ds*={stats['min']:.4f}...{stats['max']:.4f} GeV.",
        f"Central point M^2=3.25 GeV^2, s0=7.0 GeV^2: f_Ds*={central:.4f} GeV.",
        f"Assigned normalization input: f_Ds*={central:.4f} +- {assigned_err:.4f} GeV.",
        f"Wrote grid to {csv_path}",
    ]
    summary_path = OUT / "twopoint_f_Dsstar_vector_LO_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
