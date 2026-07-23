"""Stage-3 Bs1 decay-constant scenario using Pullin-Zwicky inputs.

The older Bs1 rows used HQ-scaled Ds1 decay constants as temporary inputs.
This script keeps the completed Stage-1/2 OPE kernels but replaces those
temporary constants by the heavy-light sum-rule values

    f_B1s  = 0.305 +- 0.027 GeV,
    f_B1s^T = 0.285 +- 0.027 GeV,

treated as a dedicated Bs1 decay-constant scenario.  We do not call this a full
mixed-current two-point solution, since the local AA/BB/AB two-point matrix has
not been derived for the bottom sector.  It is, however, a better Bs1 baseline
than the naive HQ-scaled constants.
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

from bs1_stage2_baseline import FPERP_S_1GEV, evaluate_point
from final_stage2_uncertainty_scan import clipped_normal, precompute_f1_basis, set_photon_shape_params
from stage2_axial_g1_three_particle import F1_integral


N_POINTS = 500
SEED = 20260706

F_B1S_A_PZ = 0.305
F_B1S_A_PZ_ERR = 0.027
F_B1S_T_PZ = 0.285
F_B1S_T_PZ_ERR = 0.027


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


def sample_inputs(rng, target, scenario):
    qq_scale = clipped_normal(rng, 0.240, 0.010, 0.200, 0.280)
    ss_ratio = clipped_normal(rng, 0.8, 0.1, 0.4, 1.2)
    ss = ss_ratio * (-(qq_scale**3))
    if scenario == "legacy_chi_condensate":
        chi = clipped_normal(rng, 3.15, 0.30, 2.0, 4.2)
    elif scenario == "lattice_fperp_s":
        fperp = clipped_normal(rng, FPERP_S_1GEV, abs(0.0018 * 1.08), -0.080, -0.030)
        chi = fperp / ss
    else:
        raise ValueError(f"unknown scenario: {scenario}")

    return {
        "mc": clipped_normal(rng, 4.18, 0.03, 4.05, 4.30),  # m_b
        "ms": clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
        "m_ds1": clipped_normal(rng, target["m_initial"], target["m_sigma"]),
        "m_ds": clipped_normal(rng, 5.36692, 0.00010),
        "f_ds1": clipped_normal(rng, F_B1S_A_PZ, F_B1S_A_PZ_ERR, 0.220, 0.390),
        "f_ds": clipped_normal(rng, 0.2303, 0.0013, 0.225, 0.236),
        "fT": clipped_normal(rng, F_B1S_T_PZ, F_B1S_T_PZ_ERR, 0.200, 0.370),
        "ss": ss,
        "chi": chi,
        "f3g": clipped_normal(rng, -0.0039, 0.0020, -0.010, 0.004),
        "ec": -1.0 / 3.0,
        "es": -1.0 / 3.0,
        "fperp_s_used": chi * ss,
        "omegaA": clipped_normal(rng, -2.1, 1.0, -6.0, 2.0),
        "omegaV": clipped_normal(rng, 3.8, 1.8, -2.0, 10.0),
    }


def evaluate_sample(vals, target, M2, s0_root, theta_deg, f1_axial, f1_basis):
    set_photon_shape_params(vals.pop("omegaA"), vals.pop("omegaV"))
    row = evaluate_point(vals, M2, s0_root, theta_deg, f1_axial, f1_basis)
    row["G_quoted"] = row[f"G_{target['quoted_combo']}_combo"]
    row["Gamma_quoted_keV"] = row[f"Gamma_{target['quoted_combo']}_combo_keV"]
    return row


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


def plot_observed_histograms(rows):
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
    fig.suptitle(r"$B_{s1}(5830)$ PZ decay-constant scenario")
    fig.tight_layout()
    fig.savefig(OUT / "stage3_bs1_pz_width_histograms.png", dpi=220)
    fig.savefig(OUT / "stage3_bs1_pz_width_histograms.pdf")


def main():
    rng = np.random.default_rng(SEED)
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows = []
    for target in TARGETS:
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            for ensemble in ("theta_fixed", "theta_scan_25_45"):
                for _ in range(N_POINTS):
                    vals = sample_inputs(rng, target, scenario)
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
                            "f_B1s_A_PZ_input": vals["f_ds1"],
                            "f_B1s_T_PZ_input": vals["fT"],
                        }
                    )
                    rows.append(row)

    write_csv(OUT / "stage3_bs1_pz_decay_constant_scan.csv", rows)
    plot_observed_histograms(rows)

    summary_rows = []
    lines = [
        "Stage-3 Bs1 Pullin-Zwicky decay-constant scenario",
        "==================================================",
        f"f_B1s={F_B1S_A_PZ:.3f}+-{F_B1S_A_PZ_ERR:.3f} GeV, "
        f"f_B1s^T={F_B1S_T_PZ:.3f}+-{F_B1S_T_PZ_ERR:.3f} GeV.",
        "These replace the temporary HQ-scaled Bs1 constants in the Stage-1/2 kernels.",
        "",
    ]
    for state in (t["state"] for t in TARGETS):
        lines.append(f"State: {state}")
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            for ensemble in ("theta_fixed", "theta_scan_25_45"):
                subset = [
                    r for r in rows
                    if r["state"] == state and r["scenario"] == scenario and r["ensemble"] == ensemble
                ]
                g = summarize([r["G_quoted"] for r in subset])
                gamma = summarize([r["Gamma_quoted_keV"] for r in subset])
                fbeff = summarize([r["fB_effective"] for r in subset])
                summary = {
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
                    "fB_effective_median": fbeff["median"],
                    "fB_effective_p16": fbeff["p16"],
                    "fB_effective_p84": fbeff["p84"],
                }
                summary_rows.append(summary)
                lines.append(
                    f"  {scenario}, {ensemble}: "
                    f"G={g['median']:+.4g} [{g['p16']:+.4g},{g['p84']:+.4g}] GeV^-1; "
                    f"Gamma={gamma['median']:.4g} [{gamma['p16']:.4g},{gamma['p84']:.4g}] keV; "
                    f"fB_eff={fbeff['median']:.3g} [{fbeff['p16']:.3g},{fbeff['p84']:.3g}] GeV"
                )
        lines.append("")

    write_csv(OUT / "stage3_bs1_pz_decay_constant_summary.csv", summary_rows)
    (OUT / "stage3_bs1_pz_decay_constant_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
