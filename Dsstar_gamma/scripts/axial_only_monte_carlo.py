"""Monte Carlo uncertainty scan for the axial-only Ds1 -> Ds* gamma baseline.

This is not the final mixed-current prediction because the tensor-current
contribution g_B^* has not yet been derived.  It propagates uncertainties for
the completed axial-current LCSR checkpoint, including the calculated
two-point f_Ds* normalization.
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
import photon_da as pda
from lattice_photon_normalization_comparison import FPERP_S_1GEV, FPERP_S_1GEV_ERR
from stage1_axial_colangelo_gstar import (
    F2_integral,
    F3_integral,
    gstar_axial_sum_rule,
    width_epsilon_keV,
)


SEED = 20260711
N_POINTS = 500


def clipped_normal(rng, mean, sigma, lo=None, hi=None):
    x = float(rng.normal(mean, sigma))
    if lo is not None:
        x = max(lo, x)
    if hi is not None:
        x = min(hi, x)
    return x


def summarize(values):
    arr = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
        "p16": float(np.percentile(arr, 16.0)),
        "median": float(np.percentile(arr, 50.0)),
        "p84": float(np.percentile(arr, 84.0)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def set_photon_shape_params(omega_a, omega_v):
    pda.PARAMS["omegaA"] = (omega_a, pda.PARAMS["omegaA"][1])
    pda.PARAMS["omegaV"] = (omega_v, pda.PARAMS["omegaV"][1])


def f3_integral_linear_coefficients():
    """Precompute F3 integral as c0 + cA omegaA + cV omegaV."""
    old_a = pda.PARAMS["omegaA"]
    old_v = pda.PARAMS["omegaV"]
    try:
        set_photon_shape_params(0.0, 0.0)
        c0 = F3_integral(0.5)
        set_photon_shape_params(1.0, 0.0)
        c_a = F3_integral(0.5) - c0
        set_photon_shape_params(0.0, 1.0)
        c_v = F3_integral(0.5) - c0
    finally:
        pda.PARAMS["omegaA"] = old_a
        pda.PARAMS["omegaV"] = old_v
    return c0, c_a, c_v


def sample_inputs(rng, scenario):
    qq_scale = clipped_normal(rng, 0.240, 0.010, 0.200, 0.280)
    ss_ratio = clipped_normal(rng, 0.8, 0.1, 0.4, 1.2)
    ss = ss_ratio * (-(qq_scale**3))

    if scenario == "legacy_chi_condensate":
        chi = clipped_normal(rng, 3.15, 0.30, 2.0, 4.2)
        fperp = chi * ss
    elif scenario == "lattice_fperp_s":
        fperp = clipped_normal(rng, FPERP_S_1GEV, FPERP_S_1GEV_ERR, -0.080, -0.030)
        chi = fperp / ss
    else:
        raise ValueError(f"unknown scenario {scenario}")

    return {
        "mc": clipped_normal(rng, 1.27, 0.02, 1.15, 1.40),
        "ms": clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
        "mA": clipped_normal(rng, 2.4595, 0.0006),
        "mA2": clipped_normal(rng, 2.53511, 0.00006),
        "mV": clipped_normal(rng, 2.1122, 0.0004),
        "fA": clipped_normal(rng, 0.225, 0.025, 0.150, 0.320),
        "fV": clipped_normal(rng, 0.227, 0.028, 0.150, 0.320),
        "f1": clipped_normal(rng, 0.345, 0.017, 0.250, 0.450),
        "f2": clipped_normal(rng, 0.379, 0.042, 0.240, 0.520),
        "ss": ss,
        "chi": chi,
        "fperp_s_used": fperp,
        "f3g": clipped_normal(rng, -0.0039, 0.0020, -0.010, 0.004),
        "omegaA": clipped_normal(rng, -2.1, 1.0, -6.0, 2.0),
        "omegaV": clipped_normal(rng, 3.8, 1.8, -2.0, 10.0),
        "ec": 2.0 / 3.0,
        "es": -1.0 / 3.0,
    }


def evaluate_sample(vals, M2, s0, theta_deg, f2_int, f3_coeffs):
    set_photon_shape_params(vals["omegaA"], vals["omegaV"])
    c0, c_a, c_v = f3_coeffs
    f3_int = c0 + c_a * vals["omegaA"] + c_v * vals["omegaV"]
    row = gstar_axial_sum_rule(
        M2,
        s0,
        mc=vals["mc"],
        ms=vals["ms"],
        mA=vals["mA"],
        mV=vals["mV"],
        fA=vals["fA"],
        fV=vals["fV"],
        ss=vals["ss"],
        chi=vals["chi"],
        f3g=vals["f3g"],
        ec=vals["ec"],
        es=vals["es"],
        f2_int=f2_int,
        f3_int=f3_int,
    )
    theta = math.radians(theta_deg)
    g1 = math.sin(theta) * vals["fA"] * row["gA_star"] / vals["f1"]
    g2 = math.cos(theta) * vals["fA"] * row["gA_star"] / vals["f2"]
    return {
        **row,
        "theta_deg": theta_deg,
        "g1_2460_axial_only": g1,
        "g2_2536_axial_only": g2,
        "Gamma2460_axial_only_keV": width_epsilon_keV(vals["mA"], vals["mV"], g1),
        "Gamma2536_axial_only_keV": width_epsilon_keV(vals["mA2"], vals["mV"], g2),
        "fV_Dsstar_GeV": vals["fV"],
        "fA_GeV": vals["fA"],
        "f1_GeV": vals["f1"],
        "f2_GeV": vals["f2"],
        "chi_effective": vals["chi"],
        "fperp_s_used_GeV": vals["fperp_s_used"],
        "omegaA": vals["omegaA"],
        "omegaV": vals["omegaV"],
    }


def write_csv(path, rows):
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_plot(rows):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5), constrained_layout=True)
    colors = {"legacy_chi_condensate": "#6f6f6f", "lattice_fperp_s": "#1f77b4"}
    labels = {
        "legacy_chi_condensate": r"legacy $\chi_s$",
        "lattice_fperp_s": r"lattice $f_{\gamma,s}^{\perp}$",
    }
    for ax, key, title in (
        (axes[0], "Gamma2460_axial_only_keV", r"$D_{s1}(2460)\to D_s^\ast\gamma$"),
        (axes[1], "Gamma2536_axial_only_keV", r"$D_{s1}(2536)\to D_s^\ast\gamma$"),
    ):
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            vals = [float(r[key]) for r in rows if r["scenario"] == scenario and r["ensemble"] == "theta_fixed"]
            p99 = np.percentile(vals, 99.0)
            vals = [v for v in vals if v <= p99]
            ax.hist(vals, bins=28, density=True, alpha=0.45, color=colors[scenario], label=labels[scenario])
            stats = summarize(vals)
            ax.axvline(stats["median"], color=colors[scenario], lw=1.8)
            ax.axvline(stats["p16"], color=colors[scenario], lw=1.2, ls="--")
            ax.axvline(stats["p84"], color=colors[scenario], lw=1.2, ls="--")
        ax.text(0.97, 0.93, title, transform=ax.transAxes, ha="right", va="top", fontsize=11)
        ax.set_xlabel(r"$\Gamma$ [keV]")
        ax.set_ylabel("density")
        ax.grid(True, alpha=0.25)
    axes[0].legend(frameon=False, fontsize=9)
    path = OUT / "axial_only_monte_carlo_histograms.pdf"
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)


def main():
    rng = np.random.default_rng(SEED)
    f2_int = F2_integral(0.5)
    f3_coeffs = f3_integral_linear_coefficients()
    rows = []
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for ensemble in ("theta_fixed", "theta_scan_25_45"):
            for _ in range(N_POINTS):
                vals = sample_inputs(rng, scenario)
                theta = 35.3 if ensemble == "theta_fixed" else float(rng.uniform(25.0, 45.0))
                M2 = float(rng.uniform(3.0, 4.5))
                s0 = float(rng.uniform(7.5, 8.5))
                out = evaluate_sample(vals, M2, s0, theta, f2_int, f3_coeffs)
                out.update({"scenario": scenario, "ensemble": ensemble, "M2": M2, "s0": s0})
                rows.append(out)

    csv_path = OUT / "axial_only_monte_carlo.csv"
    write_csv(csv_path, rows)
    make_plot(rows)

    summary_rows = []
    lines = [
        "Axial-only Ds1 -> Ds* gamma Monte Carlo",
        "=======================================",
        f"N={N_POINTS} per scenario/ensemble.",
        "This propagates the completed axial-current LCSR only; tensor-current g_B* is not included.",
        "Sampled M^2 in [3.0,4.5] GeV^2 and s0 in [7.5,8.5] GeV^2.",
        "Sampled f_Ds*=0.227+-0.028 GeV from the two-point baseline.",
        "",
    ]
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for ensemble in ("theta_fixed", "theta_scan_25_45"):
            subset = [r for r in rows if r["scenario"] == scenario and r["ensemble"] == ensemble]
            entry = {"scenario": scenario, "ensemble": ensemble, "n": len(subset)}
            lines.append(f"{scenario}, {ensemble}, n={len(subset)}")
            for key in (
                "g1_2460_axial_only",
                "Gamma2460_axial_only_keV",
                "g2_2536_axial_only",
                "Gamma2536_axial_only_keV",
                "gA_star",
                "fV_Dsstar_GeV",
                "chi_effective",
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
    summary_csv = OUT / "axial_only_monte_carlo_summary.csv"
    write_csv(summary_csv, summary_rows)
    summary_txt = OUT / "axial_only_monte_carlo_summary.txt"
    summary_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {csv_path}")
    print(f"Wrote {summary_csv}")
    print(f"Wrote {OUT / 'axial_only_monte_carlo_histograms.pdf'}")


if __name__ == "__main__":
    main()
