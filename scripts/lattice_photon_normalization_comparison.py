"""Compare legacy chi*<sbar s> with lattice-normalized f_gamma,s^perp.

The legacy setup samples chi and <sbar s> independently.  The lattice scenario
uses the recent direct normalization of the leading-twist strange photon DA,

    f_{gamma,s}^perp(2 GeV) = -51.0(1.8) MeV,

and converts it to the working 1 GeV scale with a simple leading-order tensor
current running factor.  Only the leading-twist chi*<sbar s> normalization is
replaced.  Condensate-dependent higher-twist and three-particle terms keep the
explicit <sbar s> input.
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

from final_stage2_uncertainty_scan import (
    clipped_normal,
    matched_f1_integral,
    precompute_f1_basis,
    set_photon_shape_params,
)
from stage1_axial_g1_baseline import central_inputs
from stage1_tensor_gb_hard_candidate import hard_tensor_gb
from stage1_tensor_gb_soft2p_estimate import gb_soft_tensor_da, width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from stage2_tensor_gb_three_particle_estimate import gb_three_particle_from_integral


N_POINTS = 500
SEED = 20260703

FPERP_S_2GEV = -0.0510
FPERP_S_2GEV_ERR = 0.0018
# LO tensor-current running estimate from 2 GeV down to 1 GeV.
# The factor is kept explicit so the comparison can be updated easily.
TENSOR_RUNNING_2_TO_1 = 1.08
FPERP_S_1GEV = FPERP_S_2GEV * TENSOR_RUNNING_2_TO_1
FPERP_S_1GEV_ERR = FPERP_S_2GEV_ERR * TENSOR_RUNNING_2_TO_1


def sample_base_inputs(rng):
    qq_scale = clipped_normal(rng, 0.240, 0.010, 0.200, 0.280)
    ss_ratio = clipped_normal(rng, 0.8, 0.1, 0.4, 1.2)
    qq = -(qq_scale**3)
    return {
        "mc": clipped_normal(rng, 1.27, 0.02, 1.15, 1.40),
        "ms": clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
        "m_ds1": clipped_normal(rng, 2.4595, 0.0006),
        "m_ds1_2536": clipped_normal(rng, 2.53511, 0.00006),
        "m_ds": clipped_normal(rng, 1.96835, 0.00007),
        "f_ds1": clipped_normal(rng, 0.225, 0.025, 0.150, 0.320),
        "f_ds": clipped_normal(rng, 0.2499, 0.0005, 0.245, 0.255),
        "fT": clipped_normal(rng, 0.256, 0.017, 0.190, 0.330),
        "ss": ss_ratio * qq,
        "f3g": clipped_normal(rng, -0.0039, 0.0020, -0.010, 0.004),
        "ec": 2.0 / 3.0,
        "es": -1.0 / 3.0,
        "omegaA": clipped_normal(rng, -2.1, 1.0, -6.0, 2.0),
        "omegaV": clipped_normal(rng, 3.8, 1.8, -2.0, 10.0),
    }


def sample_inputs(rng, scenario):
    vals = sample_base_inputs(rng)
    if scenario == "legacy_chi_condensate":
        vals["chi"] = clipped_normal(rng, 3.15, 0.30, 2.0, 4.2)
        vals["fperp_s_used"] = vals["chi"] * vals["ss"]
    elif scenario == "lattice_fperp_s":
        fperp = clipped_normal(rng, FPERP_S_1GEV, FPERP_S_1GEV_ERR, -0.080, -0.030)
        vals["chi"] = fperp / vals["ss"]
        vals["fperp_s_used"] = fperp
    else:
        raise ValueError(f"Unknown scenario: {scenario}")
    return vals


def central_inputs_for_scenario(scenario):
    vals = central_inputs()
    vals["m_ds1_2536"] = 2.53511
    vals["fT"] = 0.256
    vals["omegaA"] = -2.1
    vals["omegaV"] = 3.8
    if scenario == "legacy_chi_condensate":
        vals["fperp_s_used"] = vals["chi"] * vals["ss"]
    elif scenario == "lattice_fperp_s":
        vals["chi"] = FPERP_S_1GEV / vals["ss"]
        vals["fperp_s_used"] = FPERP_S_1GEV
    else:
        raise ValueError(f"Unknown scenario: {scenario}")
    return vals


def evaluate(vals, M2, s0_root, theta_deg, f1_axial_total, f1_basis):
    set_photon_shape_params(vals["omegaA"], vals["omegaV"])
    inputs = dict(vals)
    m2536 = inputs.pop("m_ds1_2536")
    fT = inputs.pop("fT")
    inputs.pop("omegaA")
    inputs.pop("omegaV")
    fperp = inputs.pop("fperp_s_used")

    s0 = s0_root * s0_root
    theta = math.radians(theta_deg)
    axial = g1_stage2(M2, s0, inputs, f1_axial_total)
    GA = axial["g1_stage2_GeV_inv"]
    fB = fT * inputs["m_ds1"] / (inputs["mc"] + inputs["ms"])
    GB_hard, _ = hard_tensor_gb(M2, s0, inputs, fB)
    GB_soft = gb_soft_tensor_da(axial, inputs, fB)
    f1_tensor = matched_f1_integral(inputs, f1_basis)
    GB_3p = gb_three_particle_from_integral(axial, inputs, fB, f1_tensor)
    GB = GB_hard + GB_soft + GB_3p
    G2460 = math.sin(theta) * GA + math.cos(theta) * GB
    G2536 = math.cos(theta) * GA - math.sin(theta) * GB
    return {
        "M2": M2,
        "s0_root": s0_root,
        "theta_deg": theta_deg,
        "GA": GA,
        "GB": GB,
        "GB_hard": GB_hard,
        "GB_soft2p": GB_soft,
        "GB_3p": GB_3p,
        "G2460": G2460,
        "G2536": G2536,
        "Gamma2460_keV": width_keV(inputs["m_ds1"], inputs["m_ds"], G2460),
        "Gamma2536_keV": width_keV(m2536, inputs["m_ds"], G2536),
        "chi_effective": inputs["chi"],
        "fperp_s_used_GeV": fperp,
    }


def summarize(values):
    arr = np.array(values, dtype=float)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
        "median": float(np.percentile(arr, 50.0)),
        "p16": float(np.percentile(arr, 16.0)),
        "p84": float(np.percentile(arr, 84.0)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def run_scan():
    rng = np.random.default_rng(SEED)
    f1_axial_total, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows = []
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        central = central_inputs_for_scenario(scenario)
        rows.append(
            {
                "scenario": scenario,
                "ensemble": "central",
                **evaluate(central, 4.5, 2.55, 35.3, f1_axial_total, f1_basis),
            }
        )
        for ensemble in ("theta_fixed", "theta_scan_25_45"):
            for _ in range(N_POINTS):
                vals = sample_inputs(rng, scenario)
                theta_deg = 35.3 if ensemble == "theta_fixed" else float(rng.uniform(25.0, 45.0))
                rows.append(
                    {
                        "scenario": scenario,
                        "ensemble": ensemble,
                        **evaluate(
                            vals,
                            float(rng.uniform(3.0, 6.0)),
                            float(rng.uniform(2.50, 2.60)),
                            theta_deg,
                            f1_axial_total,
                            f1_basis,
                        ),
                    }
                )
    return rows


def write_csv(path, rows):
    fieldnames = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def comparison_plot(rows):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9))
    colors = {"legacy_chi_condensate": "#7f7f7f", "lattice_fperp_s": "#1f77b4"}
    labels = {
        "legacy_chi_condensate": r"legacy $\chi\langle\bar{s}s\rangle$",
        "lattice_fperp_s": r"lattice $f_{\gamma,s}^{\perp}$",
    }
    for ax, key, title in (
        (axes[0], "Gamma2460_keV", r"$D_{s1}(2460)$"),
        (axes[1], "Gamma2536_keV", r"$D_{s1}(2536)$"),
    ):
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            subset = [
                float(r[key])
                for r in rows
                if r["scenario"] == scenario and r["ensemble"] == "theta_fixed"
            ]
            ax.hist(
                subset,
                bins=26,
                density=True,
                histtype="step",
                linewidth=1.8,
                color=colors[scenario],
                label=labels[scenario],
            )
        ax.set_xlabel(r"$\Gamma$ [keV]")
        ax.set_ylabel("density")
        ax.set_title(title + r", $\theta=35.3^\circ$")
        ax.grid(True, alpha=0.25, linewidth=0.7)
        ax.legend(frameon=False, fontsize=8.5)
    fig.suptitle("Photon-DA normalization comparison")
    fig.tight_layout()
    fig.savefig(OUT / "lattice_fperp_comparison_histograms.png", dpi=220)
    fig.savefig(OUT / "lattice_fperp_comparison_histograms.pdf")


def main():
    rows = run_scan()
    write_csv(OUT / "lattice_fperp_comparison_scan.csv", rows)
    comparison_plot(rows)

    summary_rows = []
    lines = [
        "Legacy chi versus lattice f_gamma,s^perp comparison",
        "===================================================",
        f"Legacy central chi*<sbar s> = {central_inputs_for_scenario('legacy_chi_condensate')['fperp_s_used']:+.5f} GeV.",
        f"Lattice input f_gamma,s^perp(2 GeV) = {FPERP_S_2GEV:+.5f} +/- {FPERP_S_2GEV_ERR:.5f} GeV.",
        f"Using LO running factor 2->1 GeV = {TENSOR_RUNNING_2_TO_1:.3f}.",
        f"Lattice scenario f_gamma,s^perp(1 GeV) = {FPERP_S_1GEV:+.5f} +/- {FPERP_S_1GEV_ERR:.5f} GeV.",
        "",
    ]
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for ensemble in ("central", "theta_fixed", "theta_scan_25_45"):
            subset = [r for r in rows if r["scenario"] == scenario and r["ensemble"] == ensemble]
            if ensemble == "central":
                r = subset[0]
                lines.append(
                    f"{scenario}, central: Gamma2460={r['Gamma2460_keV']:.3f} keV, "
                    f"Gamma2536={r['Gamma2536_keV']:.3f} keV, "
                    f"chi_eff={r['chi_effective']:.3f} GeV^-2"
                )
                summary_rows.append(
                    {
                        "scenario": scenario,
                        "ensemble": ensemble,
                        "observable": "central",
                        "median": "",
                        "p16": "",
                        "p84": "",
                        "mean": "",
                        "std": "",
                        "value_Gamma2460_keV": r["Gamma2460_keV"],
                        "value_Gamma2536_keV": r["Gamma2536_keV"],
                    }
                )
                continue
            lines.append(f"{scenario}, {ensemble}:")
            for key in ("Gamma2460_keV", "Gamma2536_keV", "G2460", "G2536"):
                stats = summarize([r[key] for r in subset])
                lines.append(
                    f"  {key}: median {stats['median']:.4g}; "
                    f"68% [{stats['p16']:.4g}, {stats['p84']:.4g}]; "
                    f"mean {stats['mean']:.4g} +/- {stats['std']:.4g}"
                )
                summary_rows.append(
                    {
                        "scenario": scenario,
                        "ensemble": ensemble,
                        "observable": key,
                        **stats,
                        "value_Gamma2460_keV": "",
                        "value_Gamma2536_keV": "",
                    }
                )
            lines.append("")

    write_csv(OUT / "lattice_fperp_comparison_summary.csv", summary_rows)
    summary_path = OUT / "lattice_fperp_comparison_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
