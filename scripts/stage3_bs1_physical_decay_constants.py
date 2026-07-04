"""Stage-3 Bs1 physical-current decay constants from the mixed-current model.

This is the bottom-sector analogue of the Ds1 Stage-3 normalization step, but
with a different closure condition.  For Ds1 we calibrated the off-diagonal
two-point overlap from the working f1 anchor.  For Bs1 we do not have such an
anchor, so we determine the overlap by requiring the physical currents

    J1 = sin(theta) JA + cos(theta) JB,
    J2 = cos(theta) JA - sin(theta) JB

to be diagonal in the two-point matrix, Pi12 = 0.  This is not a new local
two-point OPE derivation; it is the self-contained mixed-current normalization
model that uses the same basis-current inputs already used in the Bs1 scan.
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
from stage1_tensor_gb_soft2p_estimate import width_keV
from stage2_axial_g1_three_particle import F1_integral


N_POINTS = 500
SEED = 20260707

F_B1S_A_BASIS = 0.305
F_B1S_A_BASIS_ERR = 0.027
F_B1S_T_BASIS = 0.285
F_B1S_T_BASIS_ERR = 0.027


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


def diagonal_overlap(theta, fA, fB):
    """Return rho_AB that makes the mixed physical two-point correlator diagonal."""
    s = math.sin(theta)
    c = math.cos(theta)
    denom = (c * c - s * s) * fA * fB
    if abs(denom) < 1.0e-12:
        return None
    rho = -s * c * (fA * fA - fB * fB) / denom
    if rho < -1.0 or rho > 1.0:
        return None
    return rho


def physical_decay_constants(theta, fA, fB):
    rho = diagonal_overlap(theta, fA, fB)
    if rho is None:
        return None
    s = math.sin(theta)
    c = math.cos(theta)
    f1_sq = s * s * fA * fA + c * c * fB * fB + 2.0 * s * c * rho * fA * fB
    f2_sq = c * c * fA * fA + s * s * fB * fB - 2.0 * s * c * rho * fA * fB
    if f1_sq <= 0.0 or f2_sq <= 0.0:
        return None
    return rho, math.sqrt(f1_sq), math.sqrt(f2_sq)


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
        "mc": clipped_normal(rng, 4.18, 0.03, 4.05, 4.30),
        "ms": clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
        "m_ds1": clipped_normal(rng, target["m_initial"], target["m_sigma"]),
        "m_ds": clipped_normal(rng, 5.36692, 0.00010),
        "f_ds1": clipped_normal(rng, F_B1S_A_BASIS, F_B1S_A_BASIS_ERR, 0.220, 0.390),
        "f_ds": clipped_normal(rng, 0.2303, 0.0013, 0.225, 0.236),
        "fT": clipped_normal(rng, F_B1S_T_BASIS, F_B1S_T_BASIS_ERR, 0.200, 0.370),
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
    omegaA = vals.pop("omegaA")
    omegaV = vals.pop("omegaV")
    set_photon_shape_params(omegaA, omegaV)
    row = evaluate_point(vals, M2, s0_root, theta_deg, f1_axial, f1_basis)

    theta = math.radians(theta_deg)
    s = math.sin(theta)
    c = math.cos(theta)
    fA = vals["f_ds1"]
    fB = row["fB_effective"]
    decay_constants = physical_decay_constants(theta, fA, fB)
    if decay_constants is None:
        return None
    rho, f1_phys, f2_phys = decay_constants

    g_low_phys = (s * fA * row["GA"] + c * fB * row["GB"]) / f1_phys
    g_high_phys = (c * fA * row["GA"] - s * fB * row["GB"]) / f2_phys
    row["rho_AB_diag"] = rho
    row["f1_phys_GeV"] = f1_phys
    row["f2_phys_GeV"] = f2_phys
    row["G_low_phys"] = g_low_phys
    row["G_high_phys"] = g_high_phys
    row["Gamma_low_phys_keV"] = width_keV(vals["m_ds1"], vals["m_ds"], g_low_phys)
    row["Gamma_high_phys_keV"] = width_keV(vals["m_ds1"], vals["m_ds"], g_high_phys)
    row["G_quoted_phys"] = row[f"G_{target['quoted_combo']}_phys"]
    row["Gamma_quoted_phys_keV"] = row[f"Gamma_{target['quoted_combo']}_phys_keV"]
    row["G_quoted_basis"] = row[f"G_{target['quoted_combo']}_combo"]
    row["Gamma_quoted_basis_keV"] = row[f"Gamma_{target['quoted_combo']}_combo_keV"]
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
                r["Gamma_quoted_phys_keV"]
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
    fig.suptitle(r"$B_{s1}(5830)$ physical-current decay constants")
    fig.tight_layout()
    fig.savefig(OUT / "stage3_bs1_physical_width_histograms.png", dpi=220)
    fig.savefig(OUT / "stage3_bs1_physical_width_histograms.pdf")


def main():
    rng = np.random.default_rng(SEED)
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows = []
    rejected = 0
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
                    if row is None:
                        rejected += 1
                        continue
                    row.update(
                        {
                            "state": target["state"],
                            "scenario": scenario,
                            "ensemble": ensemble,
                            "quoted_combo": target["quoted_combo"],
                            "f_B1s_A_basis_input": vals["f_ds1"],
                            "f_B1s_T_basis_input": vals["fT"],
                        }
                    )
                    rows.append(row)

    write_csv(OUT / "stage3_bs1_physical_decay_constant_scan.csv", rows)
    plot_observed_histograms(rows)

    summary_rows = []
    lines = [
        "Stage-3 Bs1 physical-current decay constants",
        "============================================",
        "Closure: diagonalize the mixed-current two-point matrix, Pi12=0.",
        "Basis-current input distributions are the same bottom-sector hadronic inputs",
        f"used in the previous PZ scenario: fA={F_B1S_A_BASIS:.3f}+-{F_B1S_A_BASIS_ERR:.3f} GeV, "
        f"fT={F_B1S_T_BASIS:.3f}+-{F_B1S_T_BASIS_ERR:.3f} GeV.",
        f"Rejected samples with |rho_AB|>1 or negative f_i^2: {rejected}.",
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
                if not subset:
                    continue
                stats = {
                    key: summarize([r[key] for r in subset])
                    for key in (
                        "f1_phys_GeV",
                        "f2_phys_GeV",
                        "rho_AB_diag",
                        "G_quoted_phys",
                        "Gamma_quoted_phys_keV",
                        "G_quoted_basis",
                        "Gamma_quoted_basis_keV",
                    )
                }
                summary = {
                    "state": state,
                    "scenario": scenario,
                    "ensemble": ensemble,
                    "n": len(subset),
                }
                for key, st in stats.items():
                    summary[f"{key}_median"] = st["median"]
                    summary[f"{key}_p16"] = st["p16"]
                    summary[f"{key}_p84"] = st["p84"]
                    summary[f"{key}_mean"] = st["mean"]
                    summary[f"{key}_std"] = st["std"]
                summary_rows.append(summary)
                lines.append(
                    "  {scenario}, {ensemble}, n={n}: f1={f1:.3f} [{f1l:.3f},{f1h:.3f}] GeV; "
                    "f2={f2:.3f} [{f2l:.3f},{f2h:.3f}] GeV; rho={rho:+.3f} [{rhol:+.3f},{rhoh:+.3f}]".format(
                        scenario=scenario,
                        ensemble=ensemble,
                        n=len(subset),
                        f1=stats["f1_phys_GeV"]["median"],
                        f1l=stats["f1_phys_GeV"]["p16"],
                        f1h=stats["f1_phys_GeV"]["p84"],
                        f2=stats["f2_phys_GeV"]["median"],
                        f2l=stats["f2_phys_GeV"]["p16"],
                        f2h=stats["f2_phys_GeV"]["p84"],
                        rho=stats["rho_AB_diag"]["median"],
                        rhol=stats["rho_AB_diag"]["p16"],
                        rhoh=stats["rho_AB_diag"]["p84"],
                    )
                )
                lines.append(
                    "    physical: G={g:+.4g} [{gl:+.4g},{gh:+.4g}] GeV^-1; "
                    "Gamma={ga:.4g} [{gal:.4g},{gah:.4g}] keV".format(
                        g=stats["G_quoted_phys"]["median"],
                        gl=stats["G_quoted_phys"]["p16"],
                        gh=stats["G_quoted_phys"]["p84"],
                        ga=stats["Gamma_quoted_phys_keV"]["median"],
                        gal=stats["Gamma_quoted_phys_keV"]["p16"],
                        gah=stats["Gamma_quoted_phys_keV"]["p84"],
                    )
                )
                lines.append(
                    "    basis/PZ comparison: G={g:+.4g} [{gl:+.4g},{gh:+.4g}] GeV^-1; "
                    "Gamma={ga:.4g} [{gal:.4g},{gah:.4g}] keV".format(
                        g=stats["G_quoted_basis"]["median"],
                        gl=stats["G_quoted_basis"]["p16"],
                        gh=stats["G_quoted_basis"]["p84"],
                        ga=stats["Gamma_quoted_basis_keV"]["median"],
                        gal=stats["Gamma_quoted_basis_keV"]["p16"],
                        gah=stats["Gamma_quoted_basis_keV"]["p84"],
                    )
                )
        lines.append("")

    write_csv(OUT / "stage3_bs1_physical_decay_constant_summary.csv", summary_rows)
    (OUT / "stage3_bs1_physical_decay_constant_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
