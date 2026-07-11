"""Two-point SVZ baseline for f_Bs*.

This is the bottom analogue of the f_Ds* calculation.  We use the LO
transverse vector-current spectral density

    rho_T(s) = s (1-z)^2 (2+z)/(8 pi^2), z=m_b^2/s,

and define

    <0|sbar gamma_mu b|Bs*(q,eps)> = m_Bs* f_Bs* eps_mu.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Bsstar_gamma" / "outputs"
OUT.mkdir(exist_ok=True)


def gl_quad(func, a: float, b: float, n: int = 360) -> float:
    if b <= a:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (b - a) * x + 0.5 * (a + b)
    return float(0.5 * (b - a) * np.sum(w * func(t)))


def rho_T_LO(s, mQ):
    s = np.asarray(s, dtype=float)
    z = mQ * mQ / s
    return s * (1.0 - z) ** 2 * (2.0 + z) / (8.0 * math.pi**2)


def f_vector_LO(M2, s0, *, mQ, mV):
    pert = gl_quad(lambda s: np.exp(-s / M2) * rho_T_LO(s, mQ), mQ * mQ, s0)
    return math.sqrt(max(math.exp(mV * mV / M2) * pert / (mV * mV), 0.0)), pert


def summarize(values):
    arr = np.asarray(values, dtype=float)
    return {
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "median": float(np.percentile(arr, 50.0)),
        "p16": float(np.percentile(arr, 16.0)),
        "p84": float(np.percentile(arr, 84.0)),
    }


def main():
    inputs = {"mQ": 4.18, "mV": 5.4154}
    M2_values = np.linspace(5.0, 8.0, 7)
    s0_values = [33.0, 34.0, 35.0, 36.0]
    rows = []
    for s0 in s0_values:
        for M2 in M2_values:
            fV, pert = f_vector_LO(float(M2), float(s0), **inputs)
            rows.append({"M2": float(M2), "s0": float(s0), "mQ": inputs["mQ"], "mV": inputs["mV"], "f_Bsstar_GeV": fV, "perturbative_integral": pert})

    csv_path = OUT / "twopoint_f_Bsstar_vector_LO.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    stats = summarize([r["f_Bsstar_GeV"] for r in rows])
    central, _ = f_vector_LO(6.5, 34.5, **inputs)
    half_envelope = 0.5 * (stats["max"] - stats["min"])
    systematic = 0.10 * central
    assigned_err = math.sqrt(half_envelope**2 + systematic**2)
    lines = [
        "Two-point vector-current f_Bs* baseline",
        "=======================================",
        "LO transverse vector-current sum rule with rho_T(s)=s(1-z)^2(2+z)/(8 pi^2).",
        "Window scanned: M^2=5.0...8.0 GeV^2, s0=33,34,35,36 GeV^2.",
        f"Window envelope: f_Bs*={stats['min']:.4f}...{stats['max']:.4f} GeV.",
        f"Central point M^2=6.5 GeV^2, s0=34.5 GeV^2: f_Bs*={central:.4f} GeV.",
        f"Assigned normalization input: f_Bs*={central:.4f} +- {assigned_err:.4f} GeV.",
        f"Wrote grid to {csv_path}",
    ]
    summary_path = OUT / "twopoint_f_Bsstar_vector_LO_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
