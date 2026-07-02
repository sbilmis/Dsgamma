"""Uncertainty scan and stability plots for Bs1 -> Bs gamma.

This extends the first Bs1 baseline with a transparent Monte Carlo scan.  The
observed Bs1(5830) is treated as the high-combination state.  The lower 1+
partner is kept as a model diagnostic, not as a final prediction.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib.pyplot as plt

from bs1_stage2_baseline import FPERP_S_1GEV, bs1_inputs, evaluate_point
from final_stage2_uncertainty_scan import (
    clipped_normal,
    precompute_f1_basis,
    set_photon_shape_params,
)
from stage2_axial_g1_three_particle import F1_integral


N_POINTS = 500
SEED = 20260704


TARGETS = [
    {
        "state": "Bs1_5830_observed_high_combo",
        "m_initial": 5.82870,
        "m_sigma": 0.00020,
        "quoted_combo": "high",
        "M2_min": 8.0,
        "M2_max": 14.0,
        "s0_min": 6.05,
        "s0_max": 6.25,
    },
    {
        "state": "Bs1_lower_model_5750_low_combo",
        "m_initial": 5.750,
        "m_sigma": 0.026,
        "quoted_combo": "low",
        "M2_min": 8.0,
        "M2_max": 14.0,
        "s0_min": 5.95,
        "s0_max": 6.15,
    },
]


def sample_bs1_inputs(rng, target, scenario):
    vals = bs1_inputs(
        clipped_normal(rng, target["m_initial"], target["m_sigma"]),
        scenario,
    )
    qq_scale = clipped_normal(rng, 0.240, 0.010, 0.200, 0.280)
    ss_ratio = clipped_normal(rng, 0.8, 0.1, 0.4, 1.2)
    ss = ss_ratio * (-(qq_scale**3))
    vals["mc"] = clipped_normal(rng, 4.18, 0.03, 4.05, 4.30)
    vals["ms"] = clipped_normal(rng, 0.093, 0.011, 0.050, 0.140)
    vals["m_ds"] = clipped_normal(rng, 5.36692, 0.00010)
    vals["f_ds"] = clipped_normal(rng, 0.2303, 0.0013, 0.225, 0.236)
    vals["f_ds1"] = clipped_normal(rng, vals["f_ds1"], 0.20 * vals["f_ds1"], 0.070, 0.260)
    vals["fT"] = clipped_normal(rng, vals["fT"], 0.20 * vals["fT"], 0.080, 0.300)
    vals["ss"] = ss
    if scenario == "legacy_chi_condensate":
        vals["chi"] = clipped_normal(rng, 3.15, 0.30, 2.0, 4.2)
    elif scenario == "lattice_fperp_s":
        fperp = clipped_normal(rng, FPERP_S_1GEV, abs(0.0018 * 1.08), -0.080, -0.030)
        vals["chi"] = fperp / vals["ss"]
    vals["fperp_s_used"] = vals["chi"] * vals["ss"]
    vals["f3g"] = clipped_normal(rng, -0.0039, 0.0020, -0.010, 0.004)
    vals["omegaA"] = clipped_normal(rng, -2.1, 1.0, -6.0, 2.0)
    vals["omegaV"] = clipped_normal(rng, 3.8, 1.8, -2.0, 10.0)
    return vals


def evaluate_sample(vals, target, M2, s0_root, theta_deg, f1_axial, f1_basis):
    set_photon_shape_params(vals.pop("omegaA"), vals.pop("omegaV"))
    row = evaluate_point(vals, M2, s0_root, theta_deg, f1_axial, f1_basis)
    row["G_quoted"] = row[f"G_{target['quoted_combo']}_combo"]
    row["Gamma_quoted_keV"] = row[f"Gamma_{target['quoted_combo']}_combo_keV"]
    return row


def summarize(values):
    arr = np.array(values, dtype=float)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
        "p16": float(np.percentile(arr, 16.0)),
        "median": float(np.percentile(arr, 50.0)),
        "p84": float(np.percentile(arr, 84.0)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def write_csv(path, rows):
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_mc():
    rng = np.random.default_rng(SEED)
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows = []
    for target in TARGETS:
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            for ensemble in ("theta_fixed", "theta_scan_25_45"):
                for _ in range(N_POINTS):
                    vals = sample_bs1_inputs(rng, target, scenario)
                    theta_deg = 35.3 if ensemble == "theta_fixed" else float(rng.uniform(25.0, 45.0))
                    row = evaluate_sample(
                        vals,
                        target,
                        float(rng.uniform(target["M2_min"], target["M2_max"])),
                        float(rng.uniform(target["s0_min"], target["s0_max"])),
                        theta_deg,
                        f1_axial,
                        f1_basis,
                    )
                    row.update(
                        {
                            "state": target["state"],
                            "scenario": scenario,
                            "ensemble": ensemble,
                            "quoted_combo": target["quoted_combo"],
                        }
                    )
                    rows.append(row)
    return rows


def write_summary(rows):
    summary_rows = []
    lines = [
        "Bs1 matched Stage-2 Monte Carlo scan",
        "====================================",
        f"{N_POINTS} points per state/scenario/ensemble.",
        "Observed Bs1(5830) is quoted as the high-combination state.",
        "Lower 1+ row is a model-mass diagnostic.",
        "The Bs1 axial and tensor decay constants are assigned 20% widths.",
        "",
    ]
    for state in {r["state"] for r in rows}:
        lines.append(f"State: {state}")
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            for ensemble in ("theta_fixed", "theta_scan_25_45"):
                subset = [
                    r
                    for r in rows
                    if r["state"] == state
                    and r["scenario"] == scenario
                    and r["ensemble"] == ensemble
                ]
                gamma = summarize([r["Gamma_quoted_keV"] for r in subset])
                g = summarize([r["G_quoted"] for r in subset])
                chi = summarize([r["chi_effective"] for r in subset])
                summary_rows.append(
                    {
                        "state": state,
                        "scenario": scenario,
                        "ensemble": ensemble,
                        "G_median": g["median"],
                        "G_p16": g["p16"],
                        "G_p84": g["p84"],
                        "Gamma_median": gamma["median"],
                        "Gamma_p16": gamma["p16"],
                        "Gamma_p84": gamma["p84"],
                        "Gamma_mean": gamma["mean"],
                        "Gamma_std": gamma["std"],
                        "chi_median": chi["median"],
                        "chi_p16": chi["p16"],
                        "chi_p84": chi["p84"],
                    }
                )
                lines.append(
                    f"  {scenario}, {ensemble}: "
                    f"G={g['median']:+.4g} [{g['p16']:+.4g},{g['p84']:+.4g}] GeV^-1; "
                    f"Gamma={gamma['median']:.4g} [{gamma['p16']:.4g},{gamma['p84']:.4g}] keV; "
                    f"chi={chi['median']:.3g} [{chi['p16']:.3g},{chi['p84']:.3g}]"
                )
        lines.append("")
    write_csv(OUT / "bs1_stage2_mc_summary.csv", summary_rows)
    (OUT / "bs1_stage2_mc_summary.txt").write_text("\n".join(lines) + "\n")
    return summary_rows


def plot_histograms(rows):
    state = "Bs1_5830_observed_high_combo"
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=True)
    colors = {"legacy_chi_condensate": "#7f7f7f", "lattice_fperp_s": "#1f77b4"}
    labels = {
        "legacy_chi_condensate": r"legacy $\chi_s$",
        "lattice_fperp_s": r"lattice $f_{\gamma,s}^{\perp}$",
    }
    for ax, ensemble, title in (
        (axes[0], "theta_fixed", r"$\theta=35.3^\circ$"),
        (axes[1], "theta_scan_25_45", r"$25^\circ\leq\theta\leq45^\circ$"),
    ):
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            subset = [
                r["Gamma_quoted_keV"]
                for r in rows
                if r["state"] == state
                and r["ensemble"] == ensemble
                and r["scenario"] == scenario
            ]
            ax.hist(
                subset,
                bins=28,
                density=True,
                histtype="step",
                linewidth=1.8,
                color=colors[scenario],
                label=labels[scenario],
            )
        ax.set_title(title)
        ax.set_xlabel(r"$\Gamma[B_{s1}(5830)\to B_s\gamma]$ [keV]")
        ax.grid(True, alpha=0.25, linewidth=0.7)
        ax.legend(frameon=False, fontsize=8.5)
    axes[0].set_ylabel("density")
    fig.suptitle(r"$B_{s1}(5830)$ Monte Carlo width distributions")
    fig.tight_layout()
    fig.savefig(OUT / "bs1_5830_mc_width_histograms.png", dpi=220)
    fig.savefig(OUT / "bs1_5830_mc_width_histograms.pdf")


def plot_stability():
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    target = TARGETS[0]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.6), sharex=True)
    axes = axes.ravel()
    colors = {"legacy_chi_condensate": "#7f7f7f", "lattice_fperp_s": "#1f77b4"}
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for s0_root in (6.05, 6.15, 6.25):
            rows = []
            for M2 in np.linspace(8.0, 14.0, 25):
                vals = bs1_inputs(target["m_initial"], scenario)
                vals["omegaA"] = -2.1
                vals["omegaV"] = 3.8
                row = evaluate_sample(vals, target, float(M2), s0_root, 35.3, f1_axial, f1_basis)
                rows.append(row)
            label = rf"{scenario.replace('_', ' ')}, $\sqrt{{s_0}}={s0_root:.2f}$"
            style = "-" if scenario == "lattice_fperp_s" else "--"
            axes[0].plot([r["M2"] for r in rows], [r["G_quoted"] for r in rows], style, color=colors[scenario], alpha=0.55, label=label)
            axes[1].plot([r["M2"] for r in rows], [r["Gamma_quoted_keV"] for r in rows], style, color=colors[scenario], alpha=0.55)
            axes[2].plot([r["M2"] for r in rows], [r["GA"] for r in rows], style, color=colors[scenario], alpha=0.55)
            axes[3].plot([r["M2"] for r in rows], [r["GB"] for r in rows], style, color=colors[scenario], alpha=0.55)
    axes[0].set_ylabel(r"$G_{\rm high}$ [GeV$^{-1}$]")
    axes[1].set_ylabel(r"$\Gamma$ [keV]")
    axes[2].set_ylabel(r"$G_A$ [GeV$^{-1}$]")
    axes[3].set_ylabel(r"$G_B$ [GeV$^{-1}$]")
    for ax in axes:
        ax.set_xlabel(r"$M^2$ [GeV$^2$]")
        ax.grid(True, alpha=0.25, linewidth=0.7)
    axes[0].legend(frameon=False, fontsize=6.8)
    fig.suptitle(r"$B_{s1}(5830)\to B_s\gamma$ central stability")
    fig.tight_layout()
    fig.savefig(OUT / "bs1_5830_stability.png", dpi=220)
    fig.savefig(OUT / "bs1_5830_stability.pdf")


def main():
    rows = run_mc()
    write_csv(OUT / "bs1_stage2_mc_scan.csv", rows)
    write_summary(rows)
    plot_histograms(rows)
    plot_stability()
    print((OUT / "bs1_stage2_mc_summary.txt").read_text())


if __name__ == "__main__":
    main()
