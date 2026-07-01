"""Candidate Stage-1 tensor-current hard contribution for G_B.

This script derives a first controlled perturbative hard-photon estimate for
the tensor current by calibrating the triangle-core numerator against the known
axial-current spectral density used in ``stage1_axial_g1_baseline.py``.

The calibration rule is intentionally explicit:

  - a constant core term ``-4 i a`` maps to ``2 a log``;
  - a core term ``-4 i b/(p.q)`` maps to
    ``(b/m_i^2) (m_j^2-m_i^2-s) sqrt(lambda)/s^2``.

For the tensor current the core contains ``p^2/(p.q)``.  In this first
candidate implementation we use the diagonal double-dispersion identification
``p^2 -> s`` inside that numerator.  The script checks that the same machinery
reconstructs the published axial perturbative spectral density before applying
it to the tensor current.

This is Stage 1 only: hard photon plus two-particle photon DAs.  Stage 2
background-gluon / three-particle DA terms are not included.
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
sys.path.insert(0, str(ROOT / "scripts"))

from stage1_axial_g1_baseline import (
    central_inputs,
    gauss_legendre_integral,
    g1_stage1,
    kallen,
    rho_p_colangelo_g1,
)
from stage1_tensor_gb_soft2p_estimate import gb_soft_tensor_da, width_keV


def line_kernel(s, mi, mj, a, b):
    """No-charge spectral kernel for one photon-emitting line.

    ``mi`` is the mass of the emitting quark and ``mj`` the spectator mass.
    The parameters ``a`` and ``b`` are defined by the triangle core

        -4 i [ a + b/(p.q) ].

    For tensor-current terms with p^2/(p.q), the caller passes b(s).
    """
    s = np.asarray(s, dtype=float)
    lam = np.sqrt(np.maximum(kallen(s, mj, mi), 0.0))
    num = np.clip(s - mj**2 + mi**2 - lam, 1e-300, None)
    den = np.clip(s - mj**2 + mi**2 + lam, 1e-300, None)
    log_term = np.log(num / den)
    return -3.0 / (8.0 * math.pi**2) * (
        2.0 * a * log_term
        + (b / mi**2) * (mj**2 - mi**2 - s) * lam / s**2
    )


def rho_axial_from_calibrated_core(s, mc, ms, ec, es):
    """Reconstruct Colangelo's axial rho^P from triangle-core parameters."""
    strange = line_kernel(s, ms, mc, ms, ms**2 * (mc - ms))
    heavy = line_kernel(s, mc, ms, mc, mc**2 * (ms - mc))
    return es * strange - ec * heavy


def rho_tensor_hard_candidate(s, mc, ms, ec, es):
    """Candidate perturbative spectral density for the tensor current J_B."""
    denom = mc + ms

    # Strange-line tensor core:
    #   -4 i [ mc ms/(mc+ms) + ms^2 p^2/((mc+ms) p.q) ]
    strange_a = mc * ms / denom
    strange_b = ms**2 * s / denom

    # Heavy-line tensor core:
    #   -4 i [ (2 mc^2 - mc ms)/(mc+ms)
    #          + mc^2 p^2/((mc+ms) p.q) ].
    # The heavy term appears as the ``-(s <-> c)`` contribution in the axial
    # calibration, so the line-kernel b parameter has the opposite orientation.
    heavy_a = (2.0 * mc**2 - mc * ms) / denom
    heavy_b = -(mc**2 * s / denom)

    strange = line_kernel(s, ms, mc, strange_a, strange_b)
    heavy = line_kernel(s, mc, ms, heavy_a, heavy_b)
    return es * strange - ec * heavy


def prefactor_b(M2, inputs, fB):
    mc = inputs["mc"]
    ms = inputs["ms"]
    return (
        math.exp((inputs["m_ds1"]**2 + inputs["m_ds"]**2) / (2.0 * M2))
        * (mc + ms)
        / (inputs["m_ds1"] * fB * inputs["m_ds"]**2 * inputs["f_ds"])
    )


def hard_tensor_gb(M2, s0, inputs, fB):
    lower = (inputs["mc"] + inputs["ms"]) ** 2
    pert = gauss_legendre_integral(
        lambda s: np.exp(-s / M2)
        * rho_tensor_hard_candidate(
            s, inputs["mc"], inputs["ms"], inputs["ec"], inputs["es"]
        ),
        lower,
        s0,
        n=360,
    )
    return prefactor_b(M2, inputs, fB) * pert, pert


def axial_calibration_error(inputs):
    lower = (inputs["mc"] + inputs["ms"]) ** 2
    upper = 2.60**2
    grid = np.linspace(lower * 1.000001, upper, 500)
    expected = rho_p_colangelo_g1(
        grid, inputs["mc"], inputs["ms"], inputs["ec"], inputs["es"]
    )
    reconstructed = rho_axial_from_calibrated_core(
        grid, inputs["mc"], inputs["ms"], inputs["ec"], inputs["es"]
    )
    return float(np.max(np.abs(expected - reconstructed)))


def main():
    inputs = central_inputs()
    fT = 0.256
    fB_options = {
        "fB_equal_fT": fT,
        "fB_converted_fT_M_over_msum": fT
        * inputs["m_ds1"]
        / (inputs["mc"] + inputs["ms"]),
    }
    theta = math.radians(35.3)
    calibration_error = axial_calibration_error(inputs)

    rows = []
    for s0_root in (2.50, 2.55, 2.60):
        s0 = s0_root**2
        for M2 in np.linspace(3.0, 6.0, 13):
            axial = g1_stage1(float(M2), s0, **inputs)
            GA = axial["g1_GeV_inv"]
            for label, fB in fB_options.items():
                GB_soft = gb_soft_tensor_da(axial, inputs, fB)
                GB_hard, hard_qcd = hard_tensor_gb(float(M2), s0, inputs, fB)
                GB = GB_soft + GB_hard
                G2460 = math.sin(theta) * GA + math.cos(theta) * GB
                G2536 = math.cos(theta) * GA - math.sin(theta) * GB
                rows.append(
                    {
                        "fB_scheme": label,
                        "M2": float(M2),
                        "s0_root": s0_root,
                        "s0": s0,
                        "GA_stage1_full": GA,
                        "GB_hard_candidate": GB_hard,
                        "GB_soft2p_tensor_DA": GB_soft,
                        "GB_stage1_candidate": GB,
                        "hard_qcd_side": hard_qcd,
                        "G2460_candidate": G2460,
                        "G2536_candidate": G2536,
                        "Gamma2460_keV_candidate": width_keV(
                            inputs["m_ds1"], inputs["m_ds"], G2460
                        ),
                        "Gamma2536_keV_candidate": width_keV(
                            2.53511, inputs["m_ds"], G2536
                        ),
                        "fB": fB,
                    }
                )

    csv_path = OUT / "stage1_tensor_gb_hard_candidate.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Stage-1 tensor-current hard candidate",
        "=====================================",
        f"Axial spectral-density reconstruction max error: {calibration_error:.3e}",
        "Uses p^2 -> s in the tensor p^2/(p.q) hard-core numerator.",
        "Includes hard candidate + tensor-DA soft two-particle contribution.",
        "Does not include vector-DA diagnostic or Stage-2 three-particle DAs.",
    ]
    for scheme in fB_options:
        subset = [r for r in rows if r["fB_scheme"] == scheme]
        gb_h = np.array([r["GB_hard_candidate"] for r in subset])
        gb = np.array([r["GB_stage1_candidate"] for r in subset])
        w2460 = np.array([r["Gamma2460_keV_candidate"] for r in subset])
        w2536 = np.array([r["Gamma2536_keV_candidate"] for r in subset])
        lines.append(
            f"{scheme}: GB_hard {gb_h.min():+.4f} to {gb_h.max():+.4f}; "
            f"GB_total {gb.min():+.4f} to {gb.max():+.4f}; "
            f"Gamma2460 {w2460.min():.2f} to {w2460.max():.2f} keV; "
            f"Gamma2536 {w2536.min():.2f} to {w2536.max():.2f} keV"
        )
    lines.append(f"Wrote grid to {csv_path}")

    summary_path = OUT / "stage1_tensor_gb_hard_candidate_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
