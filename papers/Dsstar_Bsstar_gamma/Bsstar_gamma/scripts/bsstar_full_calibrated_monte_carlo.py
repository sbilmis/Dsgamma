"""Full calibrated LCSR scan for Bs1 -> Bs* gamma.

This is the bottom-strange vector-final-state analogue of the completed
Dsstar_gamma/tensor_current_hard_calibrated.py calculation.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Bsstar_gamma" / "outputs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "Dsstar_gamma" / "scripts"))
sys.path.insert(0, str(ROOT.parent / "Ds_Bs_gamma" / "scripts"))

import matplotlib.pyplot as plt

from axial_only_monte_carlo import f3_integral_linear_coefficients, set_photon_shape_params, summarize, write_csv
from final_stage2_uncertainty_scan import clipped_normal
from lattice_photon_normalization_comparison import FPERP_S_1GEV, FPERP_S_1GEV_ERR
from stage1_axial_colangelo_gstar import F2_integral, gstar_axial_sum_rule, width_epsilon_keV
from stage3_bs1_physical_decay_constants import physical_decay_constants
from tensor_current_hard_calibrated import hard_tensor_contribution
from tensor_current_soft_correction import add_soft_tensor_correction


SEED = 20260711 + 104
N_POINTS = 500

TARGETS = [
    {
        "state": "B_{s1}(5750)",
        "state_key": "Bs1_5750",
        "mA": 5.750,
        "mA_sigma": 0.026,
        "combo": "low",
        "M2_min": 10.0,
        "M2_max": 14.0,
        "s0_min": 36.5,
        "s0_max": 40.5,
    },
    {
        "state": "B_{s1}(5830)",
        "state_key": "Bs1_5830",
        "mA": 5.82870,
        "mA_sigma": 0.00020,
        "combo": "high",
        "M2_min": 10.0,
        "M2_max": 14.0,
        "s0_min": 38.0,
        "s0_max": 41.0,
    },
]


def sample_inputs(rng, target, scenario):
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
        raise ValueError(scenario)

    return {
        "mc": clipped_normal(rng, 4.18, 0.03, 4.05, 4.30),
        "ms": clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
        "mA": clipped_normal(rng, target["mA"], target["mA_sigma"]),
        "mV": clipped_normal(rng, 5.4154, 0.0015),
        "fA": clipped_normal(rng, 0.305, 0.027, 0.220, 0.390),
        "fT": clipped_normal(rng, 0.285, 0.027, 0.200, 0.370),
        "fV": clipped_normal(rng, 0.2615, 0.0336, 0.170, 0.360),
        "ss": ss,
        "chi": chi,
        "fperp_s_used": fperp,
        "f3g": clipped_normal(rng, -0.0039, 0.0020, -0.010, 0.004),
        "ec": -1.0 / 3.0,
        "es": -1.0 / 3.0,
    }


def evaluate_target(vals, target, M2, s0, theta_deg, f2_int, f3_coeffs, rng):
    # add_soft_tensor_correction expects these physical-residue placeholders,
    # but gB_soft itself is independent of them.
    vals = dict(vals)
    theta = math.radians(theta_deg)
    fB_eff = vals["fT"] * vals["mA"] / (vals["mc"] + vals["ms"])
    phys = physical_decay_constants(theta, vals["fA"], fB_eff)
    if phys is None:
        return None
    rho, f1, f2 = phys
    vals["f1"] = f1
    vals["f2"] = f2
    vals["mA2"] = vals["mA"]

    c0, c_a, c_v = f3_coeffs
    omegaA = clipped_normal(rng, -2.1, 1.0, -6.0, 2.0)
    omegaV = clipped_normal(rng, 3.8, 1.8, -2.0, 10.0)
    set_photon_shape_params(omegaA, omegaV)
    f3_int = c0 + c_a * omegaA + c_v * omegaV
    axial = gstar_axial_sum_rule(
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
    soft = add_soft_tensor_correction(axial, vals, M2, theta_deg)
    gB_hard, qcd_hard = hard_tensor_contribution(vals, M2, s0)
    gB_total = soft["gB_soft"] + gB_hard
    s = math.sin(theta)
    c = math.cos(theta)
    g_low = (s * vals["fA"] * axial["gA_star"] + c * fB_eff * gB_total) / f1
    g_high = (c * vals["fA"] * axial["gA_star"] - s * fB_eff * gB_total) / f2
    if target["combo"] == "low":
        g_quote = g_low
    else:
        g_quote = g_high
    return {
        **axial,
        **soft,
        "gB_hard": gB_hard,
        "qcd_hard_tensor": qcd_hard,
        "gB_total": gB_total,
        "rho_AB_diag": rho,
        "f1_GeV": f1,
        "f2_GeV": f2,
        "fB_eff_GeV": fB_eff,
        "g_low": g_low,
        "g_high": g_high,
        "g_quoted": g_quote,
        "Gamma_quoted_keV": width_epsilon_keV(vals["mA"], vals["mV"], g_quote),
        "omegaA": omegaA,
        "omegaV": omegaV,
        "fV_Bsstar_GeV": vals["fV"],
        "chi_effective": vals["chi"],
        "fperp_s_used_GeV": vals["fperp_s_used"],
    }


def plot(rows):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5), constrained_layout=True)
    for ax, state_key, title in (
        (axes[0], "Bs1_5750", r"$B_{s1}(5750)\to B_s^\ast\gamma$"),
        (axes[1], "Bs1_5830", r"$B_{s1}(5830)\to B_s^\ast\gamma$"),
    ):
        for scenario, color, label in (
            ("legacy_chi_condensate", "#777777", r"legacy $\chi_s$"),
            ("lattice_fperp_s", "#1f77b4", r"lattice $f_{\gamma,s}^{\perp}$"),
        ):
            vals = [r["Gamma_quoted_keV"] for r in rows if r["state_key"] == state_key and r["scenario"] == scenario and r["ensemble"] == "theta_fixed"]
            vals = [v for v in vals if v <= np.percentile(vals, 99.0)]
            stats = summarize(vals)
            ax.hist(vals, bins=28, density=True, alpha=0.45, color=color, label=label)
            ax.axvline(stats["median"], color=color, lw=1.8)
            ax.axvline(stats["p16"], color=color, lw=1.2, ls="--")
            ax.axvline(stats["p84"], color=color, lw=1.2, ls="--")
        ax.text(0.97, 0.93, title, transform=ax.transAxes, ha="right", va="top", fontsize=11)
        ax.set_xlabel(r"$\Gamma$ [keV]")
        ax.set_ylabel("density")
        ax.grid(True, alpha=0.25)
    axes[0].legend(frameon=False, fontsize=9)
    path = OUT / "bsstar_full_calibrated_histograms.pdf"
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)


def main():
    rng = np.random.default_rng(SEED)
    f2_int = F2_integral(0.5)
    f3_coeffs = f3_integral_linear_coefficients()
    rows = []
    rejected = 0
    for target in TARGETS:
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            for ensemble in ("theta_fixed", "theta_scan_25_45"):
                for _ in range(N_POINTS):
                    vals = sample_inputs(rng, target, scenario)
                    theta = 35.3 if ensemble == "theta_fixed" else float(rng.uniform(25.0, 45.0))
                    M2 = float(rng.uniform(target["M2_min"], target["M2_max"]))
                    s0 = float(rng.uniform(target["s0_min"], target["s0_max"]))
                    out = evaluate_target(vals, target, M2, s0, theta, f2_int, f3_coeffs, rng)
                    if out is None:
                        rejected += 1
                        continue
                    out.update({"state": target["state"], "state_key": target["state_key"], "combo": target["combo"], "scenario": scenario, "ensemble": ensemble, "theta_deg": theta, "M2": M2, "s0": s0})
                    rows.append(out)

    csv_path = OUT / "bsstar_full_calibrated_monte_carlo.csv"
    write_csv(csv_path, rows)
    plot(rows)
    summary_rows = []
    lines = [
        "Full calibrated Bs1 -> Bs* gamma Monte Carlo",
        "============================================",
        f"N target={N_POINTS} per state/scenario/ensemble; accepted={len(rows)}, rejected={rejected}.",
        "Uses f_Bs*=0.2615+-0.0336 GeV from the two-point vector-current baseline.",
        "",
    ]
    for target in TARGETS:
        lines.append(target["state"])
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            for ensemble in ("theta_fixed", "theta_scan_25_45"):
                subset = [r for r in rows if r["state_key"] == target["state_key"] and r["scenario"] == scenario and r["ensemble"] == ensemble]
                if not subset:
                    continue
                entry = {"state": target["state"], "state_key": target["state_key"], "scenario": scenario, "ensemble": ensemble, "n": len(subset)}
                lines.append(f"  {scenario}, {ensemble}, n={len(subset)}")
                for key in ("g_quoted", "Gamma_quoted_keV", "gB_soft", "gB_hard", "gB_total", "fV_Bsstar_GeV", "f1_GeV", "f2_GeV"):
                    stats = summarize([float(r[key]) for r in subset])
                    for stat_key, value in stats.items():
                        entry[f"{key}_{stat_key}"] = value
                    lines.append(f"    {key}: median {stats['median']:.4g} [{stats['p16']:.4g},{stats['p84']:.4g}], mean {stats['mean']:.4g}, std {stats['std']:.4g}")
                summary_rows.append(entry)
        lines.append("")
    summary_csv = OUT / "bsstar_full_calibrated_monte_carlo_summary.csv"
    write_csv(summary_csv, summary_rows)
    summary_txt = OUT / "bsstar_full_calibrated_monte_carlo_summary.txt"
    summary_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {csv_path}")
    print(f"Wrote {summary_csv}")
    print(f"Wrote {OUT / 'bsstar_full_calibrated_histograms.pdf'}")


if __name__ == "__main__":
    main()
