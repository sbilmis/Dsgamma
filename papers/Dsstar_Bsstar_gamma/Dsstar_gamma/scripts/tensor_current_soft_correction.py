"""Soft-DA tensor-current correction for Ds1 -> Ds* gamma.

This script uses the FeynCalc epsilon-projection ratios in
soft_2p_tensor_current_projection.wl to estimate the tensor-current basis
coupling g_B^* from the already computed axial-current QCD-side components.

Scope:
  included:
    - leading tensor photon DA chi phi_gamma;
    - two-particle twist-4 tensor family;
    - F2 three-particle twist-4 tensor family;
    - f3gamma vector/axial family with the vector-projection ratio as a
      controlled approximation.

  not included:
    - hard-photon tensor-current double spectral density.

The output is therefore a "soft tensor-current correction", not yet the final
full tensor-current LCSR.  The notes make this caveat explicit.
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
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "Dsstar_gamma" / "scripts"))

import matplotlib.pyplot as plt

from axial_only_monte_carlo import (
    N_POINTS,
    SEED,
    evaluate_sample,
    f3_integral_linear_coefficients,
    sample_inputs,
    set_photon_shape_params,
    summarize,
    write_csv,
)
from stage1_axial_colangelo_gstar import F2_integral, width_epsilon_keV


def ratio_tensor_da(vals, ubar=0.5):
    """B/A epsilon-projection ratio for tensor DA family at k=p+ubar q."""
    mA = vals["mA"]
    mV = vals["mV"]
    mc = vals["mc"]
    ms = vals["ms"]
    p2 = mV * mV
    pq = 0.5 * (mA * mA - mV * mV)
    kp = p2 + ubar * pq
    kP = p2 + (1.0 + ubar) * pq
    return (11.0 * kp + kP + 7.0 * pq) / (12.0 * mc * (mc + ms))


def ratio_vector_da(vals, ubar=0.5):
    """B/A epsilon-projection ratio for vector DA at k=p+ubar q."""
    mA = vals["mA"]
    mV = vals["mV"]
    mc = vals["mc"]
    ms = vals["ms"]
    p2 = mV * mV
    pq = 0.5 * (mA * mA - mV * mV)
    return mc * (p2 + pq) / ((mc + ms) * (p2 + ubar * pq))


def prefactor_b(vals, M2, fB):
    return math.exp((vals["mA"] ** 2 + vals["mV"] ** 2) / (2.0 * M2)) / (
        vals["mA"] * fB * vals["mV"] * vals["fV"]
    )


def tensor_decay_constant_current(vals):
    """Effective basis-current decay constant for J_B.

    The current contains P^alpha/(m_Q+m_s), so we use the same conversion as in
    the Ds gamma normalization checks: f_B = f_T m_A/(m_Q+m_s).
    """
    fT = vals.get("fT", 0.256)
    return fT * vals["mA"] / (vals["mc"] + vals["ms"])


def add_soft_tensor_correction(axial, vals, M2, theta_deg):
    fB = tensor_decay_constant_current(vals)
    rt = ratio_tensor_da(vals)
    rv = ratio_vector_da(vals)

    qcd_tensor_family = axial["twist2_chi_phi"] + axial["twist4_2p"] + axial["twist4_3p_F2"]
    qcd_f3_family = axial["twist3_2p_f3g"] + axial["twist3_3p_F3"]
    qcd_B_soft = rt * qcd_tensor_family + rv * qcd_f3_family
    gB_soft = prefactor_b(vals, M2, fB) * qcd_B_soft

    theta = math.radians(theta_deg)
    s = math.sin(theta)
    c = math.cos(theta)
    g1_full_soft = (s * vals["fA"] * axial["gA_star"] + c * fB * gB_soft) / vals["f1"]
    g2_full_soft = (c * vals["fA"] * axial["gA_star"] - s * fB * gB_soft) / vals["f2"]

    return {
        "gB_soft": gB_soft,
        "fB_eff_GeV": fB,
        "ratio_tensor_da": rt,
        "ratio_vector_da": rv,
        "qcd_B_soft": qcd_B_soft,
        "qcd_tensor_family_A": qcd_tensor_family,
        "qcd_f3_family_A": qcd_f3_family,
        "g1_2460_softAB": g1_full_soft,
        "g2_2536_softAB": g2_full_soft,
        "Gamma2460_softAB_keV": width_epsilon_keV(vals["mA"], vals["mV"], g1_full_soft),
        "Gamma2536_softAB_keV": width_epsilon_keV(vals["mA2"], vals["mV"], g2_full_soft),
    }


def make_plot(rows):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5), constrained_layout=True)
    colors = {"legacy_chi_condensate": "#6f6f6f", "lattice_fperp_s": "#1f77b4"}
    labels = {
        "legacy_chi_condensate": r"legacy $\chi_s$",
        "lattice_fperp_s": r"lattice $f_{\gamma,s}^{\perp}$",
    }
    for ax, key, title in (
        (axes[0], "Gamma2460_softAB_keV", r"$D_{s1}(2460)\to D_s^\ast\gamma$"),
        (axes[1], "Gamma2536_softAB_keV", r"$D_{s1}(2536)\to D_s^\ast\gamma$"),
    ):
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            vals = [float(r[key]) for r in rows if r["scenario"] == scenario and r["ensemble"] == "theta_fixed"]
            p99 = np.percentile(vals, 99.0)
            vals = [v for v in vals if v <= p99]
            stats = summarize(vals)
            ax.hist(vals, bins=28, density=True, alpha=0.45, color=colors[scenario], label=labels[scenario])
            ax.axvline(stats["median"], color=colors[scenario], lw=1.8)
            ax.axvline(stats["p16"], color=colors[scenario], lw=1.2, ls="--")
            ax.axvline(stats["p84"], color=colors[scenario], lw=1.2, ls="--")
        ax.text(0.97, 0.93, title, transform=ax.transAxes, ha="right", va="top", fontsize=11)
        ax.set_xlabel(r"$\Gamma$ [keV]")
        ax.set_ylabel("density")
        ax.grid(True, alpha=0.25)
    axes[0].legend(frameon=False, fontsize=9)
    path = OUT / "tensor_soft_monte_carlo_histograms.pdf"
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)


def main():
    rng = np.random.default_rng(SEED + 17)
    f2_int = F2_integral(0.5)
    f3_coeffs = f3_integral_linear_coefficients()
    rows = []
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for ensemble in ("theta_fixed", "theta_scan_25_45"):
            for _ in range(N_POINTS):
                vals = sample_inputs(rng, scenario)
                vals["fT"] = max(0.18, min(0.34, float(rng.normal(0.256, 0.017))))
                theta = 35.3 if ensemble == "theta_fixed" else float(rng.uniform(25.0, 45.0))
                M2 = float(rng.uniform(3.0, 4.5))
                s0 = float(rng.uniform(7.5, 8.5))
                axial = evaluate_sample(vals, M2, s0, theta, f2_int, f3_coeffs)
                corr = add_soft_tensor_correction(axial, vals, M2, theta)
                rows.append(
                    {
                        **axial,
                        **corr,
                        "scenario": scenario,
                        "ensemble": ensemble,
                        "M2": M2,
                        "s0": s0,
                    }
                )

    csv_path = OUT / "tensor_soft_monte_carlo.csv"
    write_csv(csv_path, rows)
    make_plot(rows)

    summary_rows = []
    lines = [
        "Soft tensor-current correction for Ds1 -> Ds* gamma",
        "===================================================",
        f"N={N_POINTS} per scenario/ensemble.",
        "Includes soft photon-DA tensor-current correction using FeynCalc projection ratios.",
        "Does not include the hard-photon tensor-current double spectral density.",
        "",
    ]
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for ensemble in ("theta_fixed", "theta_scan_25_45"):
            subset = [r for r in rows if r["scenario"] == scenario and r["ensemble"] == ensemble]
            entry = {"scenario": scenario, "ensemble": ensemble, "n": len(subset)}
            lines.append(f"{scenario}, {ensemble}, n={len(subset)}")
            for key in (
                "gB_soft",
                "ratio_tensor_da",
                "ratio_vector_da",
                "g1_2460_softAB",
                "Gamma2460_softAB_keV",
                "g2_2536_softAB",
                "Gamma2536_softAB_keV",
            ):
                stats = summarize([float(r[key]) for r in subset])
                for stat_key, value in stats.items():
                    entry[f"{key}_{stat_key}"] = value
                lines.append(
                    f"  {key}: median {stats['median']:.4g} "
                    f"[{stats['p16']:.4g},{stats['p84']:.4g}], "
                    f"mean {stats['mean']:.4g}, std {stats['std']:.4g}"
                )
            lines.append("")
            summary_rows.append(entry)

    summary_csv = OUT / "tensor_soft_monte_carlo_summary.csv"
    write_csv(summary_csv, summary_rows)
    summary_txt = OUT / "tensor_soft_monte_carlo_summary.txt"
    summary_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {csv_path}")
    print(f"Wrote {summary_csv}")
    print(f"Wrote {OUT / 'tensor_soft_monte_carlo_histograms.pdf'}")


if __name__ == "__main__":
    main()
