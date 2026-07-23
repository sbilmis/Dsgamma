"""Final-window Monte Carlo uncertainty scan.

This script is the publication-pass uncertainty scan.  It samples the final
working Borel and continuum-threshold windows directly in squared-s0 notation,
uses the preferred lattice photon normalization as the main scenario, and keeps
the legacy chi input as a comparison scenario.

The Ds physical-current normalization is supplied by accepted AA/AB/BB
two-point QCD-sum-rule samples projected with the previous-study angle.  Thus
the external theta input and the projected f1/f2 residues remain correlated;
the former external-f1/overlap closure is not used.
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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bs1_stage2_baseline import evaluate_point as evaluate_bs_basis
from final_stage2_uncertainty_scan import clipped_normal, precompute_f1_basis, set_photon_shape_params
from lattice_photon_normalization_comparison import (
    FPERP_S_1GEV,
    FPERP_S_1GEV_ERR,
    evaluate as evaluate_ds_basis,
    sample_inputs as sample_ds_inputs,
)
from mixing_angle_inputs import (
    THETA_BS_DEG,
    THETA_BS_SIGMA_DEG,
    THETA_HQET_DEG,
    THETA_SENSITIVITY_MAX_DEG,
    THETA_SENSITIVITY_MIN_DEG,
    sample_gaussian_angle,
)
from stage1_tensor_gb_soft2p_estimate import width_keV
from stage2_axial_g1_three_particle import F1_integral
from stage3_bs1_physical_decay_constants import physical_decay_constants


SEED = 20260705
N_POINTS = 500

DS_WINDOW = {
    "M2_min": 3.0,
    "M2_max": 4.5,
    "s0_min": 7.5,
    "s0_max": 8.5,
}

BS_WINDOWS = [
    {
        "window_id": "central_10_14",
        "state": "B_{s1}(5750)",
        "state_key": "Bs1_5750",
        "m_initial": 5.750,
        "m_sigma": 0.026,
        "quoted_combo": "low",
        "M2_min": 10.0,
        "M2_max": 14.0,
        "s0_min": 36.5,
        "s0_max": 40.5,
    },
    {
        "window_id": "central_10_14",
        "state": "B_{s1}(5830)",
        "state_key": "Bs1_5830",
        "m_initial": 5.82870,
        "m_sigma": 0.00020,
        "quoted_combo": "high",
        "M2_min": 10.0,
        "M2_max": 14.0,
        "s0_min": 38.0,
        "s0_max": 41.0,
    },
    {
        "window_id": "crosscheck_9_13",
        "state": "B_{s1}(5750)",
        "state_key": "Bs1_5750",
        "m_initial": 5.750,
        "m_sigma": 0.026,
        "quoted_combo": "low",
        "M2_min": 9.0,
        "M2_max": 13.0,
        "s0_min": 36.5,
        "s0_max": 40.5,
    },
    {
        "window_id": "crosscheck_9_13",
        "state": "B_{s1}(5830)",
        "state_key": "Bs1_5830",
        "m_initial": 5.82870,
        "m_sigma": 0.00020,
        "quoted_combo": "high",
        "M2_min": 9.0,
        "M2_max": 13.0,
        "s0_min": 38.0,
        "s0_max": 41.0,
    },
]


def summarize(values: list[float]) -> dict[str, float]:
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


def sample_bs_inputs(rng: np.random.Generator, target: dict[str, float | str], scenario: str) -> dict[str, float]:
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
        raise ValueError(f"unknown scenario: {scenario}")

    return {
        "mc": clipped_normal(rng, 4.18, 0.03, 4.05, 4.30),
        "ms": clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
        "m_ds1": clipped_normal(rng, float(target["m_initial"]), float(target["m_sigma"])),
        "m_ds": clipped_normal(rng, 5.36692, 0.00010),
        "f_ds1": clipped_normal(rng, 0.305, 0.027, 0.220, 0.390),
        "f_ds": clipped_normal(rng, 0.2303, 0.0013, 0.225, 0.236),
        "fT": clipped_normal(rng, 0.285, 0.027, 0.200, 0.370),
        "ss": ss,
        "chi": chi,
        "f3g": clipped_normal(rng, -0.0039, 0.0020, -0.010, 0.004),
        "ec": -1.0 / 3.0,
        "es": -1.0 / 3.0,
        "fperp_s_used": fperp,
        "omegaA": clipped_normal(rng, -2.1, 1.0, -6.0, 2.0),
        "omegaV": clipped_normal(rng, 3.8, 1.8, -2.0, 10.0),
    }


def evaluate_ds_physical(
    vals: dict[str, float],
    M2: float,
    s0: float,
    residue_sample: dict[str, str],
    scenario: str,
    ensemble: str,
    f1_axial: float,
    f1_basis: dict[str, float],
) -> list[dict[str, float | str]]:
    theta_deg = float(residue_sample["theta_deg"])
    basis = evaluate_ds_basis(vals, M2, math.sqrt(s0), theta_deg, f1_axial, f1_basis)
    theta = math.radians(theta_deg)
    s = math.sin(theta)
    c = math.cos(theta)
    fA = vals["f_ds1"]
    fB = vals["fT"] * vals["m_ds1"] / (vals["mc"] + vals["ms"])
    f1_phys = float(residue_sample["f1_GeV"])
    f2_phys = float(residue_sample["f2_GeV"])
    if f1_phys <= 0.0 or f2_phys <= 0.0:
        return []
    g1 = (s * fA * basis["GA"] + c * fB * basis["GB"]) / f1_phys
    g2 = (c * fA * basis["GA"] - s * fB * basis["GB"]) / f2_phys
    g1_a = s * fA * basis["GA"] / f1_phys
    g1_b = c * fB * basis["GB"] / f1_phys
    g2_a = c * fA * basis["GA"] / f2_phys
    g2_b = -s * fB * basis["GB"] / f2_phys
    common = {
        "sector": "Ds",
        "scenario": scenario,
        "ensemble": ensemble,
        "window_id": "central",
        "M2": M2,
        "s0": s0,
        "theta_deg": theta_deg,
        "f1_GeV": f1_phys,
        "f2_GeV": f2_phys,
        "normalization_scheme": "AA_AB_BB_two_point_LO_d3_d5",
        "two_point_M2_GeV2": float(residue_sample["M2_GeV2"]),
        "two_point_s0_mix_GeV2": float(residue_sample["s0_mix_GeV2"]),
        "two_point_s01_GeV2": float(residue_sample["s01_GeV2"]),
        "two_point_s02_GeV2": float(residue_sample["s02_GeV2"]),
        "GA_basis": basis["GA"],
        "GB_basis": basis["GB"],
        "chi_effective": basis["chi_effective"],
        "fperp_s_used_GeV": basis["fperp_s_used_GeV"],
    }
    return [
        {
            **common,
            "state": "D_{s1}(2460)",
            "state_key": "Ds1_2460",
            "g_GeV_inv": g1,
            "g_A_component_GeV_inv": g1_a,
            "g_B_component_GeV_inv": g1_b,
            "interference": "constructive" if g1_a * g1_b > 0.0 else "destructive",
            "Gamma_keV": width_keV(vals["m_ds1"], vals["m_ds"], g1),
        },
        {
            **common,
            "state": "D_{s1}(2536)",
            "state_key": "Ds1_2536",
            "g_GeV_inv": g2,
            "g_A_component_GeV_inv": g2_a,
            "g_B_component_GeV_inv": g2_b,
            "interference": "constructive" if g2_a * g2_b > 0.0 else "destructive",
            "Gamma_keV": width_keV(vals["m_ds1_2536"], vals["m_ds"], g2),
        },
    ]


def evaluate_bs_physical(
    vals: dict[str, float],
    target: dict[str, float | str],
    M2: float,
    s0: float,
    theta_deg: float,
    scenario: str,
    ensemble: str,
    f1_axial: float,
    f1_basis: dict[str, float],
) -> dict[str, float | str] | None:
    omega_a = vals.pop("omegaA")
    omega_v = vals.pop("omegaV")
    set_photon_shape_params(omega_a, omega_v)
    basis = evaluate_bs_basis(vals, M2, math.sqrt(s0), theta_deg, f1_axial, f1_basis)
    theta = math.radians(theta_deg)
    constants = physical_decay_constants(theta, vals["f_ds1"], basis["fB_effective"])
    if constants is None:
        return None
    rho, f1_phys, f2_phys = constants
    s = math.sin(theta)
    c = math.cos(theta)
    g_low = (s * vals["f_ds1"] * basis["GA"] + c * basis["fB_effective"] * basis["GB"]) / f1_phys
    g_high = (c * vals["f_ds1"] * basis["GA"] - s * basis["fB_effective"] * basis["GB"]) / f2_phys
    g_low_a = s * vals["f_ds1"] * basis["GA"] / f1_phys
    g_low_b = c * basis["fB_effective"] * basis["GB"] / f1_phys
    g_high_a = c * vals["f_ds1"] * basis["GA"] / f2_phys
    g_high_b = -s * basis["fB_effective"] * basis["GB"] / f2_phys
    quoted = str(target["quoted_combo"])
    g = g_low if quoted == "low" else g_high
    return {
        "sector": "Bs",
        "scenario": scenario,
        "ensemble": ensemble,
        "window_id": str(target["window_id"]),
        "state": str(target["state"]),
        "state_key": str(target["state_key"]),
        "M2": M2,
        "s0": s0,
        "theta_deg": theta_deg,
        "f1_GeV": f1_phys,
        "f2_GeV": f2_phys,
        "rho_AB": rho,
        "GA_basis": basis["GA"],
        "GB_basis": basis["GB"],
        "g_low_GeV_inv": g_low,
        "g_high_GeV_inv": g_high,
        "g_GeV_inv": g,
        "g_A_component_GeV_inv": g_low_a if quoted == "low" else g_high_a,
        "g_B_component_GeV_inv": g_low_b if quoted == "low" else g_high_b,
        "interference": (
            "constructive"
            if (g_low_a * g_low_b if quoted == "low" else g_high_a * g_high_b) > 0.0
            else "destructive"
        ),
        "Gamma_keV": width_keV(vals["m_ds1"], vals["m_ds"], g),
        "quoted_combo": quoted,
        "chi_effective": vals["chi"],
        "fperp_s_used_GeV": vals["fperp_s_used"],
    }


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_accepted_twopoint_samples() -> list[dict[str, str]]:
    path = OUT / "twopoint_ds1_matrix_mc.csv"
    with path.open() as f:
        rows = [row for row in csv.DictReader(f) if int(row["accepted"]) == 1]
    if not rows:
        raise RuntimeError(f"no accepted two-point samples in {path}")
    return rows


def bs_rng() -> np.random.Generator:
    """Use a reproducible stream independent of changes in the Ds loop."""

    return np.random.default_rng(SEED + 1)


def write_summary_report(
    summary_rows: list[dict[str, float | str]],
    accepted_row_count: int,
    rejected_bs: int,
) -> None:
    lines = [
        "Final-window Monte Carlo uncertainty scan",
        "=========================================",
        f"Seed: {SEED}",
        f"Accepted rows: {accepted_row_count}",
        f"Rejected bottom-sector physical-normalization samples: {rejected_bs}",
        f"Samples per ensemble/channel/scenario/window: {N_POINTS}",
        "Ds normalization: theta_Ds=26.6+-0.6 deg external input; f1/f2 projected from the AA/AB/BB two-point sum rule.",
        "Bs nominal angle: theta_Bs=38.5+-0.1 deg external input.",
        "HQET 35.3 deg and 25--45 deg rows are benchmark/sensitivity outputs only.",
        "Ds window: M^2=[3.0,4.5] GeV^2, s0=[7.5,8.5] GeV^2.",
        "Bs central window: M^2=[10,14] GeV^2 with state-specific s0 windows.",
        "Bs cross-check window: M^2=[9,13] GeV^2 with the same s0 windows.",
        "",
        "Preferred lattice-photon central-window rows:",
    ]
    for row in summary_rows:
        if row["scenario"] == "lattice_fperp_s" and (
            (row["sector"] == "Ds" and row["ensemble"] == "theta_prior_gaussian")
            or (row["sector"] == "Bs" and row["ensemble"] == "theta_prior_gaussian")
        ) and row["window_id"] in ("central", "central_10_14"):
            lines.append(
                "{state}: g={g:+.4g} [{gl:+.4g},{gh:+.4g}] GeV^-1; "
                "Gamma={ga:.4g} [{gal:.4g},{gah:.4g}] keV; n={n}".format(
                    state=row["state"],
                    g=float(row["g_median_GeV_inv"]),
                    gl=float(row["g_p16_GeV_inv"]),
                    gh=float(row["g_p84_GeV_inv"]),
                    ga=float(row["Gamma_median_keV"]),
                    gal=float(row["Gamma_p16_keV"]),
                    gah=float(row["Gamma_p84_keV"]),
                    n=int(row["n"]),
                )
            )
    lines.extend(
        [
            "",
            "Caveat: Bs rows are conditional on accepted |rho_AB|<=1 closure samples;",
            "the nominal lattice central-window acceptance is 245/500 (5750) and 236/500 (5830).",
            "A direct bottom-sector AA/AB/BB two-point matrix is still required.",
            "",
            "Outputs:",
            f"  {OUT / 'final_window_mc_scan.csv'}",
            f"  {OUT / 'final_window_mc_summary.csv'}",
            f"  {OUT / 'final_window_mc_width_histograms.pdf'}",
            f"  {OUT / 'final_window_bs_window_crosscheck.pdf'}",
        ]
    )
    (OUT / "final_window_mc_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def make_histograms(rows: list[dict[str, float | str]]) -> None:
    plot_rows = [
        r for r in rows
        if r["scenario"] == "lattice_fperp_s"
        and (
            (r["sector"] == "Ds" and r["ensemble"] == "theta_prior_gaussian")
            or (r["sector"] == "Bs" and r["ensemble"] == "theta_prior_gaussian")
        )
        and r["window_id"] in ("central", "central_10_14")
    ]
    panels = [
        ("Ds1_2460", r"$D_{s1}(2460)\to D_s\gamma$"),
        ("Ds1_2536", r"$D_{s1}(2536)\to D_s\gamma$"),
        ("Bs1_5750", r"$B_{s1}(5750)\to B_s\gamma$"),
        ("Bs1_5830", r"$B_{s1}(5830)\to B_s\gamma$"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.4))
    for ax, (state_key, title) in zip(axes.flat, panels):
        values = np.asarray([float(r["Gamma_keV"]) for r in plot_rows if r["state_key"] == state_key], dtype=float)
        stats = summarize(values.tolist())
        x_max = float(np.percentile(values, 97.5))
        if x_max <= stats["p84"]:
            x_max = 1.25 * stats["p84"]
        hist_values = values[values <= x_max]
        ax.hist(hist_values, bins=24, range=(0.0, x_max), density=True, color="#4c78a8", alpha=0.58, edgecolor="white")
        ax.axvline(stats["median"], color="black", linewidth=1.2)
        ax.axvline(stats["p16"], color="#d62728", linestyle="--", linewidth=1.0)
        ax.axvline(stats["p84"], color="#d62728", linestyle="--", linewidth=1.0)
        ax.set_xlim(0.0, 1.08 * x_max)
        ax.set_title(title)
        ax.set_xlabel(r"$\Gamma$ [keV]")
        ax.set_ylabel("density")
        ax.grid(True, alpha=0.24, linewidth=0.7)
        ax.text(
            0.97,
            0.93,
            rf"${stats['median']:.3g}_{{-{stats['median'] - stats['p16']:.2g}}}^{{+{stats['p84'] - stats['median']:.2g}}}$ keV",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
        )
    fig.tight_layout()
    fig.savefig(OUT / "final_window_mc_width_histograms.png", dpi=240)
    fig.savefig(OUT / "final_window_mc_width_histograms.pdf")


def make_bs_window_check(rows: list[dict[str, float | str]]) -> None:
    subset = [
        r for r in rows
        if r["sector"] == "Bs" and r["scenario"] == "lattice_fperp_s" and r["ensemble"] == "theta_prior_gaussian"
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=False)
    colors = {"central_10_14": "#1f77b4", "crosscheck_9_13": "#d62728"}
    labels = {"central_10_14": r"$10\leq M^2\leq14$", "crosscheck_9_13": r"$9\leq M^2\leq13$"}
    for ax, state_key, title in (
        (axes[0], "Bs1_5750", r"$B_{s1}(5750)$"),
        (axes[1], "Bs1_5830", r"$B_{s1}(5830)$"),
    ):
        for window_id in ("central_10_14", "crosscheck_9_13"):
            values = [float(r["Gamma_keV"]) for r in subset if r["state_key"] == state_key and r["window_id"] == window_id]
            ax.hist(values, bins=24, density=True, histtype="step", linewidth=1.8, color=colors[window_id], label=labels[window_id])
        ax.set_title(title)
        ax.set_xlabel(r"$\Gamma$ [keV]")
        ax.grid(True, alpha=0.24, linewidth=0.7)
        ax.legend(frameon=False, fontsize=8.5)
    axes[0].set_ylabel("density")
    fig.suptitle(r"Bottom-sector Borel-window cross-check")
    fig.tight_layout()
    fig.savefig(OUT / "final_window_bs_window_crosscheck.png", dpi=240)
    fig.savefig(OUT / "final_window_bs_window_crosscheck.pdf")


def main() -> None:
    rng_ds = np.random.default_rng(SEED)
    rng_bs = bs_rng()
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    twopoint_samples = read_accepted_twopoint_samples()
    rows: list[dict[str, float | str]] = []
    rejected_bs = 0

    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        ensemble = "theta_prior_gaussian"
        for _ in range(N_POINTS):
            vals = sample_ds_inputs(rng_ds, scenario)
            residue_sample = twopoint_samples[
                int(rng_ds.integers(0, len(twopoint_samples)))
            ]
            rows.extend(
                evaluate_ds_physical(
                    vals,
                    float(rng_ds.uniform(DS_WINDOW["M2_min"], DS_WINDOW["M2_max"])),
                    float(rng_ds.uniform(DS_WINDOW["s0_min"], DS_WINDOW["s0_max"])),
                    residue_sample,
                    scenario,
                    ensemble,
                    f1_axial,
                    f1_basis,
                )
            )

    for target in BS_WINDOWS:
        for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
            for ensemble in (
                "theta_prior_gaussian",
                "theta_hqet_35p3",
                "theta_sensitivity_25_45",
            ):
                for _ in range(N_POINTS):
                    vals = sample_bs_inputs(rng_bs, target, scenario)
                    if ensemble == "theta_prior_gaussian":
                        theta_deg = sample_gaussian_angle(
                            rng_bs, THETA_BS_DEG, THETA_BS_SIGMA_DEG
                        )
                    elif ensemble == "theta_hqet_35p3":
                        theta_deg = THETA_HQET_DEG
                    else:
                        theta_deg = float(
                            rng_bs.uniform(
                                THETA_SENSITIVITY_MIN_DEG,
                                THETA_SENSITIVITY_MAX_DEG,
                            )
                        )
                    row = evaluate_bs_physical(
                        vals,
                        target,
                        float(rng_bs.uniform(float(target["M2_min"]), float(target["M2_max"]))),
                        float(rng_bs.uniform(float(target["s0_min"]), float(target["s0_max"]))),
                        theta_deg,
                        scenario,
                        ensemble,
                        f1_axial,
                        f1_basis,
                    )
                    if row is None:
                        rejected_bs += 1
                    else:
                        rows.append(row)

    write_csv(OUT / "final_window_mc_scan.csv", rows)

    summary_rows: list[dict[str, float | str]] = []
    group_keys = sorted({(r["sector"], r["state"], r["state_key"], r["scenario"], r["ensemble"], r["window_id"]) for r in rows})
    for sector, state, state_key, scenario, ensemble, window_id in group_keys:
        subset = [
            r for r in rows
            if (r["sector"], r["state"], r["state_key"], r["scenario"], r["ensemble"], r["window_id"])
            == (sector, state, state_key, scenario, ensemble, window_id)
        ]
        g_stats = summarize([float(r["g_GeV_inv"]) for r in subset])
        gamma_stats = summarize([float(r["Gamma_keV"]) for r in subset])
        f1_stats = summarize([float(r["f1_GeV"]) for r in subset])
        f2_stats = summarize([float(r["f2_GeV"]) for r in subset])
        summary_rows.append(
            {
                "sector": sector,
                "state": state,
                "scenario": scenario,
                "ensemble": ensemble,
                "window_id": window_id,
                "n": len(subset),
                "g_median_GeV_inv": g_stats["median"],
                "g_p16_GeV_inv": g_stats["p16"],
                "g_p84_GeV_inv": g_stats["p84"],
                "Gamma_median_keV": gamma_stats["median"],
                "Gamma_p16_keV": gamma_stats["p16"],
                "Gamma_p84_keV": gamma_stats["p84"],
                "f1_median_GeV": f1_stats["median"],
                "f1_p16_GeV": f1_stats["p16"],
                "f1_p84_GeV": f1_stats["p84"],
                "f2_median_GeV": f2_stats["median"],
                "f2_p16_GeV": f2_stats["p16"],
                "f2_p84_GeV": f2_stats["p84"],
            }
        )
    write_csv(OUT / "final_window_mc_summary.csv", summary_rows)
    write_summary_report(summary_rows, len(rows), rejected_bs)
    make_histograms(rows)
    make_bs_window_check(rows)


if __name__ == "__main__":
    main()
