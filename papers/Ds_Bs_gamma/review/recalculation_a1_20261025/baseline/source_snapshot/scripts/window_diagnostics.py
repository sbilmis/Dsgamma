"""Window diagnostics for Borel and continuum-threshold choices.

This script compares candidate windows using three transparent diagnostics:

1. perturbative pole-fraction proxy for the explicit spectral integrals;
2. photon-DA continuum-retention proxy for the leading two-particle term;
3. numerical stability and contribution hierarchy inside each window.

The pole-fraction proxy is intentionally named as a proxy because the photon-DA
contact terms do not all have the same explicit continuum integral as the hard
spectral terms.
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

from bs1_stage2_baseline import FPERP_S_1GEV, evaluate_point as evaluate_bs_basis
from final_stage2_uncertainty_scan import (
    matched_f1_integral,
    precompute_f1_basis,
    set_photon_shape_params,
)
from lattice_photon_normalization_comparison import (
    central_inputs_for_scenario,
    evaluate as evaluate_ds_basis,
)
from stage1_axial_g1_baseline import (
    gauss_legendre_integral,
    rho_p_colangelo_g1,
)
from stage1_tensor_gb_hard_candidate import rho_tensor_hard_candidate
from stage1_tensor_gb_soft2p_estimate import width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from stage2_tensor_gb_three_particle_estimate import gb_three_particle_from_integral
from stage3_bs1_physical_decay_constants import physical_decay_constants
from stage3_ds1_physical_residue_normalization import calibrate_from_f1


THETA_DEG = 35.3
DS_F1_ANCHOR = 0.345
DS_FT = 0.256
BS_FA = 0.305
BS_FT = 0.285


WINDOWS = [
    {
        "name": "Ds_old_near_threshold",
        "sector": "Ds",
        "state": "Ds1(2460), Ds1(2536)",
        "M2_min": 3.0,
        "M2_max": 6.0,
        "s0_min": 2.50**2,
        "s0_max": 2.60**2,
        "quoted": "both",
    },
    {
        "name": "Ds_above_pole_broad",
        "sector": "Ds",
        "state": "Ds1(2460), Ds1(2536)",
        "M2_min": 2.0,
        "M2_max": 7.0,
        "s0_min": 2.85**2,
        "s0_max": 3.10**2,
        "quoted": "both",
    },
    {
        "name": "Ds_preferred_candidate",
        "sector": "Ds",
        "state": "Ds1(2460), Ds1(2536)",
        "M2_min": 3.0,
        "M2_max": 4.5,
        "s0_min": 7.5,
        "s0_max": 8.5,
        "quoted": "both",
    },
    {
        "name": "Ds_above_pole_upper_candidate",
        "sector": "Ds",
        "state": "Ds1(2460), Ds1(2536)",
        "M2_min": 3.5,
        "M2_max": 5.5,
        "s0_min": 2.85**2,
        "s0_max": 3.00**2,
        "quoted": "both",
    },
    {
        "name": "Bs_low_old_diagnostic",
        "sector": "Bs",
        "state": "B_{s1}(5750)",
        "mass": 5.750,
        "M2_min": 8.0,
        "M2_max": 14.0,
        "s0_min": 6.05**2,
        "s0_max": 6.35**2,
        "quoted": "low",
    },
    {
        "name": "Bs_low_candidate",
        "sector": "Bs",
        "state": "B_{s1}(5750)",
        "mass": 5.750,
        "M2_min": 10.0,
        "M2_max": 14.0,
        "s0_min": 36.5,
        "s0_max": 40.5,
        "quoted": "low",
    },
    {
        "name": "Bs5830_old_diagnostic",
        "sector": "Bs",
        "state": "B_{s1}(5830)",
        "mass": 5.82870,
        "M2_min": 8.0,
        "M2_max": 14.0,
        "s0_min": 6.15**2,
        "s0_max": 6.40**2,
        "quoted": "high",
    },
    {
        "name": "Bs5830_candidate",
        "sector": "Bs",
        "state": "B_{s1}(5830)",
        "mass": 5.82870,
        "M2_min": 10.0,
        "M2_max": 14.0,
        "s0_min": 38.0,
        "s0_max": 41.0,
        "quoted": "high",
    },
]


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


def span(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    mid = 0.5 * (float(arr.max()) + float(arr.min()))
    if abs(mid) < 1.0e-14:
        return float("nan")
    return float((float(arr.max()) - float(arr.min())) / abs(mid))


def pole_fraction_spectral(inputs: dict[str, float], M2: float, s0: float, kind: str) -> float:
    """Absolute spectral-weight pole proxy for the hard perturbative term.

    Some charge combinations make the signed spectral density change sign over
    the integration domain.  A signed ratio can then become very large or
    negative because of cancellations in the denominator, which is not a useful
    measure of continuum leakage.  The absolute-weight ratio is more stable and
    answers the diagnostic question: how much of the hard spectral weight lies
    below the chosen threshold?
    """
    lower = (inputs["mc"] + inputs["ms"]) ** 2
    upper = 80.0 if inputs["mc"] < 2.0 else 180.0
    if kind == "axial":
        rho = lambda x: rho_p_colangelo_g1(x, inputs["mc"], inputs["ms"], inputs["ec"], inputs["es"])
    elif kind == "tensor":
        rho = lambda x: rho_tensor_hard_candidate(x, inputs["mc"], inputs["ms"], inputs["ec"], inputs["es"])
    else:
        raise ValueError(kind)
    num = gauss_legendre_integral(lambda x: np.exp(-x / M2) * np.abs(rho(x)), lower, s0, n=300)
    den = gauss_legendre_integral(lambda x: np.exp(-x / M2) * np.abs(rho(x)), lower, upper, n=420)
    if abs(den) < 1.0e-16:
        return float("nan")
    return float(num / den)


def da_retention_proxy(inputs: dict[str, float], M2: float, s0: float) -> float:
    """Continuum-retained fraction of the leading exp(mc)-exp(s0) DA term."""
    mc2 = inputs["mc"] ** 2
    if s0 <= mc2:
        return 0.0
    return float(1.0 - math.exp(-(s0 - mc2) / M2))


def ds_inputs() -> dict[str, float]:
    vals = central_inputs_for_scenario("lattice_fperp_s")
    vals["fT"] = DS_FT
    return vals


def bs_inputs(m_initial: float) -> dict[str, float]:
    qq = -(0.240) ** 3
    ss = 0.8 * qq
    chi = FPERP_S_1GEV / ss
    return {
        "mc": 4.18,
        "ms": 0.093,
        "m_ds1": m_initial,
        "m_ds": 5.36692,
        "f_ds1": BS_FA,
        "f_ds": 0.2303,
        "fT": BS_FT,
        "ss": ss,
        "chi": chi,
        "f3g": -0.0039,
        "ec": -1.0 / 3.0,
        "es": -1.0 / 3.0,
        "fperp_s_used": FPERP_S_1GEV,
    }


def evaluate_ds(M2: float, s0: float, f1_axial: float, f1_basis: dict[str, float]) -> dict[str, float | str]:
    vals = ds_inputs()
    basis = evaluate_ds_basis(vals, M2, math.sqrt(s0), THETA_DEG, f1_axial, f1_basis)
    axial = g1_stage2(M2, s0, {k: vals[k] for k in ("mc", "ms", "m_ds1", "m_ds", "f_ds1", "f_ds", "ss", "chi", "f3g", "ec", "es")}, f1_axial)
    theta = math.radians(THETA_DEG)
    s = math.sin(theta)
    c = math.cos(theta)
    fA = vals["f_ds1"]
    fB = DS_FT * vals["m_ds1"] / (vals["mc"] + vals["ms"])
    calibrated = calibrate_from_f1(theta, fA, fB, DS_F1_ANCHOR)
    if calibrated is None:
        raise RuntimeError("Ds calibration failed")
    rho, f1_phys, f2_phys, _ = calibrated
    g2460 = (s * fA * basis["GA"] + c * fB * basis["GB"]) / f1_phys
    g2536 = (c * fA * basis["GA"] - s * fB * basis["GB"]) / f2_phys
    GB_abs = abs(basis["GB"]) if abs(basis["GB"]) > 1.0e-14 else 1.0
    GA_abs = abs(basis["GA"]) if abs(basis["GA"]) > 1.0e-14 else 1.0
    return {
        "sector": "Ds",
        "M2": M2,
        "s0": s0,
        "g_2460": g2460,
        "g_2536": g2536,
        "Gamma_2460_keV": width_keV(vals["m_ds1"], vals["m_ds"], g2460),
        "Gamma_2536_keV": width_keV(vals["m_ds1_2536"], vals["m_ds"], g2536),
        "GA_basis": basis["GA"],
        "GB_basis": basis["GB"],
        "GB_hard_frac": abs(basis["GB_hard"]) / GB_abs,
        "GB_soft2p_frac": abs(basis["GB_soft2p"]) / GB_abs,
        "GB_3p_frac": abs(basis["GB_3p"]) / GB_abs,
        "GA_3p_frac": abs(axial["F1_delta_g1"]) / GA_abs,
        "spectral_pole_axial": pole_fraction_spectral(vals, M2, s0, "axial"),
        "spectral_pole_tensor": pole_fraction_spectral(vals, M2, s0, "tensor"),
        "twist2_retention": da_retention_proxy(vals, M2, s0),
        "rho_AB": rho,
    }


def evaluate_bs(
    target: dict[str, float | str],
    M2: float,
    s0: float,
    f1_axial: float,
    f1_basis: dict[str, float],
) -> dict[str, float | str]:
    vals = bs_inputs(float(target["mass"]))
    set_photon_shape_params(-2.1, 3.8)
    basis = evaluate_bs_basis(vals, M2, math.sqrt(s0), THETA_DEG, f1_axial, f1_basis)
    theta = math.radians(THETA_DEG)
    constants = physical_decay_constants(theta, vals["f_ds1"], basis["fB_effective"])
    if constants is None:
        raise RuntimeError("Bs physical constants failed")
    rho, f1_phys, f2_phys = constants
    s = math.sin(theta)
    c = math.cos(theta)
    g_low = (s * vals["f_ds1"] * basis["GA"] + c * basis["fB_effective"] * basis["GB"]) / f1_phys
    g_high = (c * vals["f_ds1"] * basis["GA"] - s * basis["fB_effective"] * basis["GB"]) / f2_phys
    quoted = str(target["quoted"])
    g_quoted = g_high if quoted == "high" else g_low
    GB_abs = abs(basis["GB"]) if abs(basis["GB"]) > 1.0e-14 else 1.0
    GA_abs = abs(basis["GA"]) if abs(basis["GA"]) > 1.0e-14 else 1.0
    return {
        "sector": "Bs",
        "M2": M2,
        "s0": s0,
        "state": str(target["state"]),
        "g_low": g_low,
        "g_high": g_high,
        "g_quoted": g_quoted,
        "Gamma_low_keV": width_keV(vals["m_ds1"], vals["m_ds"], g_low),
        "Gamma_high_keV": width_keV(vals["m_ds1"], vals["m_ds"], g_high),
        "Gamma_quoted_keV": width_keV(vals["m_ds1"], vals["m_ds"], g_quoted),
        "GA_basis": basis["GA"],
        "GB_basis": basis["GB"],
        "GB_hard_frac": abs(basis["GB_hard"]) / GB_abs,
        "GB_soft2p_frac": abs(basis["GB_soft2p"]) / GB_abs,
        "GB_3p_frac": abs(basis["GB_3p"]) / GB_abs,
        "GA_3p_frac": abs(basis["GB_3p"]) / GA_abs,
        "spectral_pole_axial": pole_fraction_spectral(vals, M2, s0, "axial"),
        "spectral_pole_tensor": pole_fraction_spectral(vals, M2, s0, "tensor"),
        "twist2_retention": da_retention_proxy(vals, M2, s0),
        "rho_AB": rho,
    }


def sample_window(window: dict[str, float | str], f1_axial: float, f1_basis: dict[str, float]) -> list[dict[str, float | str]]:
    m2_values = np.linspace(float(window["M2_min"]), float(window["M2_max"]), 11)
    s0_values = np.linspace(float(window["s0_min"]), float(window["s0_max"]), 7)
    rows: list[dict[str, float | str]] = []
    for M2 in m2_values:
        for s0 in s0_values:
            if window["sector"] == "Ds":
                row = evaluate_ds(float(M2), float(s0), f1_axial, f1_basis)
            else:
                row = evaluate_bs(window, float(M2), float(s0), f1_axial, f1_basis)
            row["window"] = str(window["name"])
            row["window_state"] = str(window["state"])
            row["s0_sqrt"] = math.sqrt(float(s0))
            rows.append(row)
    return rows


def summarize_window(window: dict[str, float | str], rows: list[dict[str, float | str]]) -> dict[str, float | str]:
    if window["sector"] == "Ds":
        stability_keys = ["g_2460", "Gamma_2460_keV", "g_2536", "Gamma_2536_keV"]
    else:
        stability_keys = ["g_quoted", "Gamma_quoted_keV"]
    out: dict[str, float | str] = {
        "window": str(window["name"]),
        "sector": str(window["sector"]),
        "state": str(window["state"]),
        "M2_range_GeV2": f"{float(window['M2_min']):.2f}-{float(window['M2_max']):.2f}",
        "s0_range_GeV2": f"{float(window['s0_min']):.2f}-{float(window['s0_max']):.2f}",
        "sqrt_s0_range_GeV": f"{math.sqrt(float(window['s0_min'])):.2f}-{math.sqrt(float(window['s0_max'])):.2f}",
        "n_grid": len(rows),
        "min_spectral_pole_axial": min(float(r["spectral_pole_axial"]) for r in rows),
        "min_spectral_pole_tensor": min(float(r["spectral_pole_tensor"]) for r in rows),
        "min_twist2_retention": min(float(r["twist2_retention"]) for r in rows),
        "max_GB_3p_frac": max(float(r["GB_3p_frac"]) for r in rows),
        "max_GA_3p_frac": max(float(r["GA_3p_frac"]) for r in rows),
        "max_GB_hard_frac": max(float(r["GB_hard_frac"]) for r in rows),
        "max_GB_soft2p_frac": max(float(r["GB_soft2p_frac"]) for r in rows),
    }
    for key in stability_keys:
        values = [float(r[key]) for r in rows]
        out[f"{key}_min"] = min(values)
        out[f"{key}_max"] = max(values)
        out[f"{key}_span"] = span(values)
    return out


def make_summary_plot(summary_rows: list[dict[str, float | str]]) -> None:
    labels = [str(r["window"]) for r in summary_rows]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.0), sharex=True)
    axes[0].bar(x - 0.2, [float(r["min_spectral_pole_axial"]) for r in summary_rows], width=0.4, label="axial spectral")
    axes[0].bar(x + 0.2, [float(r["min_spectral_pole_tensor"]) for r in summary_rows], width=0.4, label="tensor spectral")
    axes[0].axhline(0.5, color="black", linewidth=1.0, linestyle="--", alpha=0.55)
    axes[0].set_ylabel("minimum pole-fraction proxy")
    axes[0].legend(frameon=False, fontsize=8.5)
    axes[0].grid(True, axis="y", alpha=0.25)

    stability = []
    for row in summary_rows:
        if str(row["sector"]) == "Ds":
            stability.append(float(row["Gamma_2460_keV_span"]))
        else:
            stability.append(float(row["Gamma_quoted_keV_span"]))
    axes[1].bar(x, stability, color="#1f77b4")
    axes[1].set_ylabel("representative width span")
    axes[1].grid(True, axis="y", alpha=0.25)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=35, ha="right")
    fig.suptitle("Window diagnostics: pole proxy and stability")
    fig.tight_layout()
    fig.savefig(OUT / "window_diagnostics_summary.png", dpi=240)
    fig.savefig(OUT / "window_diagnostics_summary.pdf")


def main() -> None:
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    all_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []
    for window in WINDOWS:
        rows = sample_window(window, f1_axial, f1_basis)
        all_rows.extend(rows)
        summary_rows.append(summarize_window(window, rows))

    write_csv(OUT / "window_diagnostics_grid.csv", all_rows)
    write_csv(OUT / "window_diagnostics_summary.csv", summary_rows)
    make_summary_plot(summary_rows)

    lines = [
        "Window diagnostics",
        "==================",
        "Scenario: central lattice f_gamma,s^perp normalization, theta=35.3 deg.",
        "Pole dominance is reported as an absolute perturbative spectral-weight proxy.",
        "The two-particle photon DA continuum entry is a retention proxy.",
        "",
    ]
    for row in summary_rows:
        lines.append(
            "{window}: M2={M2}, s0={s0} GeV^2, sqrt(s0)={rs}; "
            "pole_axial_min={pa:.3f}, pole_tensor_min={pt:.3f}, "
            "twist2_retention_min={tr:.3f}".format(
                window=row["window"],
                M2=row["M2_range_GeV2"],
                s0=row["s0_range_GeV2"],
                rs=row["sqrt_s0_range_GeV"],
                pa=float(row["min_spectral_pole_axial"]),
                pt=float(row["min_spectral_pole_tensor"]),
                tr=float(row["min_twist2_retention"]),
            )
        )
        if row["sector"] == "Ds":
            lines.append(
                "  spans: Gamma2460={:.3f}, Gamma2536={:.3f}, "
                "max GA_3p={:.3f}, max GB_3p={:.3f}".format(
                    float(row["Gamma_2460_keV_span"]),
                    float(row["Gamma_2536_keV_span"]),
                    float(row["max_GA_3p_frac"]),
                    float(row["max_GB_3p_frac"]),
                )
            )
        else:
            lines.append(
                "  spans: Gamma_quoted={:.3f}, max GA_3p={:.3f}, max GB_3p={:.3f}".format(
                    float(row["Gamma_quoted_keV_span"]),
                    float(row["max_GA_3p_frac"]),
                    float(row["max_GB_3p_frac"]),
                )
            )
        lines.append("")
    lines.extend(
        [
            "Interpretation:",
            "- A larger pole-fraction proxy is better; values below about 0.5 are a warning.",
            "- Smaller stability spans are better, but near-zero cancellation channels can have large relative spans.",
            "- Large three-particle fractions indicate that the chosen M2 range is pushing OPE hierarchy.",
            "- Preferred redo windows from this pass:",
            "  Ds1: M2 = 3.0-4.5 GeV^2, s0 = 7.5-8.5 GeV^2.",
            "  Bs1 lower diagnostic: M2 = 10-14 GeV^2, s0 = 36.5-40.5 GeV^2.",
            "  Bs1(5830): M2 = 10-14 GeV^2, s0 = 38.0-41.0 GeV^2.",
            "",
            f"Wrote {OUT / 'window_diagnostics_grid.csv'}",
            f"Wrote {OUT / 'window_diagnostics_summary.csv'}",
            f"Wrote {OUT / 'window_diagnostics_summary.pdf'}",
        ]
    )
    (OUT / "window_diagnostics_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
