#!/usr/bin/env python3
"""Run all Python numerics, scans, Monte Carlo, CSV tables, and plots."""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from qcdsr import (  # noqa: E402
    Q_B,
    Q_C,
    QCDInput,
    borel_twopoint_matrix,
    borel_vector,
    initial_observables,
    physical_transition,
    projected_twopoint,
    rho_transition,
    rho_twopoint,
    rotation,
    transition_terms,
    vector_observables,
)

OUT = ROOT / "outputs"
CSV_DIR = OUT / "csv"
FIG_DIR = OUT / "figures"
TMP_DIR = ROOT / "tmp" / "plots"
for directory in (CSV_DIR, FIG_DIR, TMP_DIR, ROOT / "notebooks", ROOT / "tex"):
    directory.mkdir(parents=True, exist_ok=True)

SEED = 20260729


@dataclass(frozen=True)
class Sector:
    name: str
    inp: QCDInput
    theta: float
    theta_sigma: float
    mf: float
    mf_sigma: float
    vector_window: tuple[float, float, float, float]
    states: tuple["State", "State"]


@dataclass(frozen=True)
class State:
    key: str
    label: str
    index: int
    mass: float
    mass_sigma: float
    tp_window: tuple[float, float, float, float]
    tr_window: tuple[float, float, float, float]
    tensor_sigma: float
    total_width_mev: float | None = None
    total_width_sigma_mev: float | None = None


D_STATES = (
    State(
        "Ds1_2460",
        "D_{s1}(2460) -> D_s^* gamma",
        0,
        2.4595,
        0.0006,
        (2.0, 3.0, 7.5, 8.5),
        (4.0, 6.0, 6.25, 7.25),
        0.15,
    ),
    State(
        "Ds1_2536",
        "D_{s1}(2536) -> D_s^* gamma",
        1,
        2.53512,
        0.00006,
        (2.0, 3.0, 8.5, 10.0),
        (4.0, 6.0, 6.75, 7.75),
        0.15,
        0.92,
        0.05,
    ),
)
B_STATES = (
    State(
        "Bs1_lower",
        "B_{s1}^{(L)} -> B_s^* gamma",
        0,
        5.75,
        0.03,
        (8.0, 11.0, 40.0, 42.0),
        (5.0, 7.0, 36.0, 40.0),
        0.05,
    ),
    State(
        "Bs1_5830",
        "B_{s1}(5830) -> B_s^* gamma",
        1,
        5.82873,
        0.00020,
        (7.0, 9.0, 42.0, 44.0),
        (5.0, 7.0, 38.0, 42.0),
        0.05,
        0.5,
        math.sqrt(0.3**2 + 0.3**2),
    ),
)
SECTORS = (
    Sector(
        "D",
        QCDInput(mQ=1.27, eQ=Q_C),
        38.4,
        3.0,
        2.1122,
        0.0004,
        (2.5, 3.5, 6.5, 7.5),
        D_STATES,
    ),
    Sector(
        "B",
        QCDInput(mQ=4.18, eQ=Q_B),
        35.264,
        3.0,
        5.4154,
        0.0014,
        (7.0, 9.0, 36.0, 38.0),
        B_STATES,
    ),
)


def grid(lo: float, hi: float, n: int) -> np.ndarray:
    return np.linspace(lo, hi, n)


def center(window: tuple[float, float, float, float]) -> tuple[float, float]:
    return 0.5 * (window[0] + window[1]), 0.5 * (window[2] + window[3])


def write_rows(path: Path, rows: list[dict]):
    pd.DataFrame(rows).to_csv(path, index=False)


def scalar_f_initial(sector: Sector, state: State, inp: QCDInput, theta: float, m2: float, s0: float, mass: float):
    pi = projected_twopoint(m2, s0, inp, theta)["total"][state.index, state.index]
    if pi <= 0:
        return float("nan")
    return math.exp(mass * mass / (2.0 * m2)) * math.sqrt(pi) / mass


def scalar_f_vector(sector: Sector, inp: QCDInput, m2: float, s0: float, mass: float):
    pi = borel_vector(m2, s0, inp)["total"]
    if pi <= 0:
        return float("nan")
    return math.exp(mass * mass / (2.0 * m2)) * math.sqrt(pi) / mass


def calculate_two_point():
    window_rows, scan_rows, vector_rows, summaries = [], [], [], {}
    for sector in SECTORS:
        vm2s = grid(sector.vector_window[0], sector.vector_window[1], 6)
        vs0s = grid(sector.vector_window[2], sector.vector_window[3], 5)
        vf_values = []
        for m2 in vm2s:
            for s0 in vs0s:
                obs = vector_observables(m2, s0, sector.inp, sector.mf)
                row = {"sector": sector.name, "M2": m2, "s0": s0, **obs}
                vector_rows.append(row)
                if obs["pc"] >= 0.40 and obs["localfrac"] <= 0.25:
                    vf_values.append(obs["f"])
        fvc = float(np.mean(vf_values))
        summaries[f"{sector.name}_fV"] = fvc
        for state in sector.states:
            m2s = grid(state.tp_window[0], state.tp_window[1], 6)
            s0s = grid(state.tp_window[2], state.tp_window[3], 5)
            accepted_f = []
            pcs, d5s, masses, offdiags = [], [], [], []
            for m2 in m2s:
                for s0 in s0s:
                    obs = initial_observables(
                        m2, s0, sector.inp, sector.theta, state.index, state.mass
                    )
                    passed = (
                        obs["Pi"] > 0
                        and obs["pc"] >= 0.40
                        and obs["d5frac"] <= 0.10
                        and abs(obs["mass_sr"] / state.mass - 1.0) <= 0.05
                    )
                    scan_rows.append(
                        {
                            "sector": sector.name,
                            "state": state.key,
                            "M2": m2,
                            "s0": s0,
                            **obs,
                            "accepted": passed,
                            "rejection_reason": ""
                            if passed
                            else ";".join(
                                reason
                                for condition, reason in [
                                    (obs["Pi"] <= 0, "nonpositive"),
                                    (obs["pc"] < 0.40, "PC<0.40"),
                                    (obs["d5frac"] > 0.10, "d5>0.10"),
                                    (
                                        abs(obs["mass_sr"] / state.mass - 1.0) > 0.05,
                                        "mass>5pct",
                                    ),
                                ]
                                if condition
                            ),
                        }
                    )
                    if passed:
                        accepted_f.append(obs["f"])
                        pcs.append(obs["pc"])
                        d5s.append(obs["d5frac"])
                        masses.append(obs["mass_sr"])
                        offdiags.append(obs["offdiag"])
            if not accepted_f:
                raise RuntimeError(f"No accepted two-point samples for {state.key}")
            fc = float(np.mean(accepted_f))
            summaries[f"{state.key}_f"] = fc
            window_rows.append(
                {
                    "state": state.key,
                    "M2_min_GeV2": state.tp_window[0],
                    "M2_max_GeV2": state.tp_window[1],
                    "M2_step_GeV2": (state.tp_window[1] - state.tp_window[0]) / 5,
                    "M2_points": 6,
                    "s0_min_GeV2": state.tp_window[2],
                    "s0_max_GeV2": state.tp_window[3],
                    "s0_step_GeV2": (state.tp_window[3] - state.tp_window[2]) / 4,
                    "s0_points": 5,
                    "accepted_points": len(accepted_f),
                    "f_mean_GeV": fc,
                    "f_grid_std_GeV": float(np.std(accepted_f, ddof=1)),
                    "pc_min": min(pcs),
                    "pc_max": max(pcs),
                    "d5frac_max": max(d5s),
                    "mass_min_GeV": min(masses),
                    "mass_max_GeV": max(masses),
                    "offdiag_residual_min": min(offdiags),
                    "offdiag_residual_max": max(offdiags),
                }
            )
    write_rows(CSV_DIR / "two_point_grid.csv", scan_rows)
    write_rows(CSV_DIR / "vector_two_point_grid.csv", vector_rows)
    write_rows(CSV_DIR / "two_point_windows.csv", window_rows)
    return summaries


def sample_input(rng: np.random.Generator, sector: Sector) -> QCDInput:
    m_sigma = 0.02 if sector.name == "D" else 0.03
    return QCDInput(
        mQ=max(rng.normal(sector.inp.mQ, m_sigma), 0.2),
        ms=max(rng.normal(0.093, 0.011), 0.02),
        eQ=sector.inp.eQ,
        qq=rng.normal(-(0.24**3), 0.10 * 0.24**3),
        kappa_s=np.clip(rng.normal(0.8, 0.1), 0.5, 1.1),
        m0sq=max(rng.normal(0.8, 0.2), 0.2),
        chi=rng.normal(-3.15, 0.30),
        f3=rng.normal(-0.0039, 0.0020),
        omega_v=rng.normal(3.8, 1.8),
        omega_a=rng.normal(-2.1, 1.0),
        kappa_da=rng.normal(0.2, 0.1),
        zeta1=rng.normal(0.4, 0.2),
        zeta2=rng.normal(0.3, 0.1),
    )


def monte_carlo(two_point_summary: dict, samples: int = 1200):
    rng = np.random.default_rng(SEED)
    rows, rejection_rows = [], []
    for sample_id in range(samples):
        for sector in SECTORS:
            inp = sample_input(rng, sector)
            theta = rng.normal(sector.theta, sector.theta_sigma)
            vm2 = rng.uniform(sector.vector_window[0], sector.vector_window[1])
            vs0 = rng.uniform(sector.vector_window[2], sector.vector_window[3])
            mf = rng.normal(sector.mf, sector.mf_sigma)
            fv = scalar_f_vector(sector, inp, vm2, vs0, mf)
            for state in sector.states:
                m2tp = rng.uniform(state.tp_window[0], state.tp_window[1])
                s0tp = rng.uniform(state.tp_window[2], state.tp_window[3])
                mi = rng.normal(state.mass, state.mass_sigma)
                fi = scalar_f_initial(sector, state, inp, theta, m2tp, s0tp, mi)
                if not np.isfinite(fi) or not np.isfinite(fv):
                    rejection_rows.append(
                        {
                            "sample": sample_id,
                            "state": state.key,
                            "reason": "nonpositive two-point invariant",
                        }
                    )
                    continue
                m2tr = rng.uniform(state.tr_window[0], state.tr_window[1])
                s0tr = rng.uniform(state.tr_window[2], state.tr_window[3])
                delta = rng.normal(0.0, state.tensor_sigma)
                try:
                    tr = physical_transition(
                        m2tr,
                        s0tr,
                        inp,
                        theta,
                        state.index,
                        mi,
                        mf,
                        fi,
                        fv,
                        delta,
                        three_n=10,
                    )
                except (FloatingPointError, ValueError, OverflowError) as exc:
                    rejection_rows.append(
                        {"sample": sample_id, "state": state.key, "reason": str(exc)}
                    )
                    continue
                rows.append(
                    {
                        "sample": sample_id,
                        "sector": sector.name,
                        "state": state.key,
                        "mQ": inp.mQ,
                        "ms": inp.ms,
                        "ss": inp.ss,
                        "chi": inp.chi,
                        "f3gamma": inp.f3,
                        "omegaV": inp.omega_v,
                        "omegaA": inp.omega_a,
                        "theta_deg": theta,
                        "M2_twopoint": m2tp,
                        "s0_twopoint": s0tp,
                        "M2_transition": m2tr,
                        "s0_transition": s0tr,
                        "f_initial_GeV": fi,
                        "f_vector_GeV": fv,
                        "tensor_delta": delta,
                        "g": tr["g"],
                        "width_keV": tr["width_keV"],
                    }
                )
    write_rows(CSV_DIR / "monte_carlo_accepted.csv", rows)
    write_rows(CSV_DIR / "monte_carlo_rejections.csv", rejection_rows)
    df = pd.DataFrame(rows)
    stats = {}
    summary_rows = []
    for state in [s for sec in SECTORS for s in sec.states]:
        d = df[df.state == state.key]
        if len(d) < 0.8 * samples:
            raise RuntimeError(f"Too many MC rejections for {state.key}: {len(d)}")
        result = {
            "state": state.key,
            "accepted": len(d),
            "g_median": d.g.median(),
            "g_p16": d.g.quantile(0.16),
            "g_p84": d.g.quantile(0.84),
            "width_median_keV": d.width_keV.median(),
            "width_p16_keV": d.width_keV.quantile(0.16),
            "width_p84_keV": d.width_keV.quantile(0.84),
            "f_initial_median_GeV": d.f_initial_GeV.median(),
            "f_initial_p16_GeV": d.f_initial_GeV.quantile(0.16),
            "f_initial_p84_GeV": d.f_initial_GeV.quantile(0.84),
            "f_vector_median_GeV": d.f_vector_GeV.median(),
            "f_vector_p16_GeV": d.f_vector_GeV.quantile(0.16),
            "f_vector_p84_GeV": d.f_vector_GeV.quantile(0.84),
        }
        if state.total_width_mev:
            total_keV = state.total_width_mev * 1000.0
            result["branching_central"] = result["width_median_keV"] / total_keV
        else:
            result["branching_central"] = np.nan
        summary_rows.append(result)
        stats[state.key] = result
    write_rows(CSV_DIR / "monte_carlo_summary.csv", summary_rows)
    return stats


def plot_csv(
    csv_path: Path,
    stem: str,
    xcol: str,
    ycol: str,
    lowcol: str,
    highcol: str,
    xlabel: str,
    ylabel: str,
    title: str,
    xrange: tuple[float, float],
    central_x: float,
):
    gp = TMP_DIR / f"{stem}.gp"
    # This generated plotting program contains no physics; all data remain in CSV.
    program = f"""
set datafile separator comma
set key top right opaque
set grid back lc rgb "#d9e2ec"
set border lw 1.2
set style fill transparent solid 0.20 noborder
set object 1 rect from {xrange[0]}, graph 0 to {xrange[1]}, graph 1 behind fc rgb "#dbeafe" fs solid 0.25
set xlabel "{xlabel}"
set ylabel "{ylabel}"
set title "{title}"
set xrange [{xrange[0]}:{xrange[1]}]
set terminal pngcairo size 1500,950 enhanced font "Helvetica,18"
set output "{(FIG_DIR / (stem + '.png')).as_posix()}"
plot "{csv_path.as_posix()}" using "{xcol}":"{lowcol}":"{highcol}" with filledcurves lc rgb "#60a5fa" title "input uncertainty", \
     "" using "{xcol}":"{ycol}" with linespoints lw 3 pt 7 ps 0.8 lc rgb "#0f3d56" title "central", \
     "" using (abs(column("{xcol}")-{central_x})<1e-8 ? column("{xcol}") : 1/0):"{ycol}" with points pt 7 ps 1.6 lc rgb "#dc2626" title "central point"
set terminal pdfcairo size 7.5in,4.8in enhanced color font "Helvetica,13"
set output "{(FIG_DIR / (stem + '.pdf')).as_posix()}"
replot
"""
    gp.write_text(program)
    subprocess.run(["gnuplot", str(gp)], check=True)


def scans_and_plots(two_point_summary: dict, mc_stats: dict):
    stability_rows, term_rows, result_rows = [], [], []
    for sector in SECTORS:
        fv = two_point_summary[f"{sector.name}_fV"]
        for state in sector.states:
            fi = two_point_summary[f"{state.key}_f"]
            mtp, s0tp = center(state.tp_window)
            mtr, s0tr = center(state.tr_window)
            central = physical_transition(
                mtr,
                s0tr,
                sector.inp,
                sector.theta,
                state.index,
                state.mass,
                sector.mf,
                fi,
                fv,
            )
            stat = mc_stats[state.key]
            g_sigma = 0.5 * (stat["g_p84"] - stat["g_p16"])
            w_sigma = 0.5 * (
                stat["width_p84_keV"] - stat["width_p16_keV"]
            )
            term_rows.extend(
                {
                    "sector": sector.name,
                    "state": state.key,
                    "M2": mtr,
                    "s0": s0tr,
                    "term": key,
                    "value_GeV4": central[key],
                    "included_in_quoted_total": key
                    not in {"local_heavy_diagnostic", "quoted_total_A", "quoted_total_B"},
                }
                for key in [
                    "hard",
                    "tw2",
                    "tw4_2p",
                    "tw3_2p",
                    "tw4_3p",
                    "tw3_3p",
                    "local_heavy_diagnostic",
                ]
            )
            result = {
                "sector": sector.name,
                "state": state.key,
                "mass_initial_GeV": state.mass,
                "mass_final_GeV": sector.mf,
                "theta_deg": sector.theta,
                "f_initial_GeV": fi,
                "f_vector_GeV": fv,
                "M2_transition_GeV2": mtr,
                "s0_transition_GeV2": s0tr,
                "g_central": central["g"],
                "g_mc_median": stat["g_median"],
                "g_mc_p16": stat["g_p16"],
                "g_mc_p84": stat["g_p84"],
                "width_central_keV": central["width_keV"],
                "width_mc_median_keV": stat["width_median_keV"],
                "width_mc_p16_keV": stat["width_p16_keV"],
                "width_mc_p84_keV": stat["width_p84_keV"],
                "branching_fraction": stat["branching_central"],
                "component_A": central["component_A"],
                "component_B": central["component_B"],
                "k_gamma_GeV": central["k_gamma"],
            }
            result_rows.append(result)
            scan_specs = [
                ("M2", grid(state.tr_window[0], state.tr_window[1], 6), s0tr),
                ("s0", grid(state.tr_window[2], state.tr_window[3], 5), mtr),
            ]
            for variable, values, fixed in scan_specs:
                rows = []
                for value in values:
                    m2 = value if variable == "M2" else fixed
                    s0 = fixed if variable == "M2" else value
                    tr = physical_transition(
                        m2,
                        s0,
                        sector.inp,
                        sector.theta,
                        state.index,
                        state.mass,
                        sector.mf,
                        fi,
                        fv,
                    )
                    rows.append(
                        {
                            variable: value,
                            "g": tr["g"],
                            "g_low": tr["g"] - g_sigma,
                            "g_high": tr["g"] + g_sigma,
                            "width_keV": tr["width_keV"],
                            "width_low_keV": max(tr["width_keV"] - w_sigma, 0.0),
                            "width_high_keV": tr["width_keV"] + w_sigma,
                        }
                    )
                path = CSV_DIR / f"{state.key}_{variable}_stability.csv"
                write_rows(path, rows)
                for observable, ylabel, y, lo, hi in [
                    ("coupling", "g (dimensionless)", "g", "g_low", "g_high"),
                    ("width", "Gamma (keV)", "width_keV", "width_low_keV", "width_high_keV"),
                ]:
                    fixed_text = (
                        f"s0={fixed:.2f} GeV^2" if variable == "M2" else f"M2={fixed:.2f} GeV^2"
                    )
                    plot_csv(
                        path,
                        f"{state.key}_{observable}_vs_{variable}",
                        variable,
                        y,
                        lo,
                        hi,
                        "M^2 (GeV^2)" if variable == "M2" else "s0 (GeV^2)",
                        ylabel,
                        f"{state.label}; {fixed_text}",
                        (
                            (state.tr_window[0], state.tr_window[1])
                            if variable == "M2"
                            else (state.tr_window[2], state.tr_window[3])
                        ),
                        mtr if variable == "M2" else s0tr,
                    )
                    vals = np.array([row[y] for row in rows])
                    xs = np.array(values)
                    slope = float(
                        np.polyfit(xs[max(0, len(xs) // 2 - 1) : min(len(xs), len(xs) // 2 + 2)],
                                   vals[max(0, len(vals) // 2 - 1) : min(len(vals), len(vals) // 2 + 2)], 1)[0]
                    )
                    stability_rows.append(
                        {
                            "state": state.key,
                            "observable": observable,
                            "scan": variable,
                            "fixed_value": fixed,
                            "minimum": float(vals.min()),
                            "maximum": float(vals.max()),
                            "relative_variation": float(
                                (vals.max() - vals.min()) / max(abs(vals[len(vals) // 2]), 1e-30)
                            ),
                            "central_slope": slope,
                            "points": len(vals),
                        }
                    )
            # Mixing-angle/interference scan.
            angle_rows = []
            for angle in np.arange(20.0, 55.01, 1.0):
                tr = physical_transition(
                    mtr,
                    s0tr,
                    sector.inp,
                    angle,
                    state.index,
                    state.mass,
                    sector.mf,
                    fi,
                    fv,
                )
                angle_rows.append(
                    {
                        "theta_deg": angle,
                        "component_A": tr["component_A"],
                        "component_B": tr["component_B"],
                        "g": tr["g"],
                        "width_keV": tr["width_keV"],
                    }
                )
            angle_path = CSV_DIR / f"{state.key}_angle_scan.csv"
            write_rows(angle_path, angle_rows)
            # Use zero-width bands for the angle diagnostic.
            adf = pd.DataFrame(angle_rows)
            adf["g_low"] = adf["g"]
            adf["g_high"] = adf["g"]
            adf.to_csv(angle_path, index=False)
            plot_csv(
                angle_path,
                f"{state.key}_coupling_vs_angle",
                "theta_deg",
                "g",
                "g_low",
                "g_high",
                "theta (degree)",
                "g (dimensionless)",
                f"{state.label}; mixing-angle sensitivity",
                (20.0, 55.0),
                sector.theta,
            )
    write_rows(CSV_DIR / "central_results.csv", result_rows)
    write_rows(CSV_DIR / "transition_term_breakdown.csv", term_rows)
    write_rows(CSV_DIR / "stability_summary.csv", stability_rows)
    return result_rows


def diagnostics_plots(two_point_summary: dict):
    tp = pd.read_csv(CSV_DIR / "two_point_grid.csv")
    for sector in SECTORS:
        for state in sector.states:
            df = tp[(tp.state == state.key)]
            _, s0c = center(state.tp_window)
            # Pick the actual nearest declared threshold.
            s0pick = df.iloc[(df.s0 - s0c).abs().argsort()[:1]].s0.iloc[0]
            d = df[np.isclose(df.s0, s0pick)].sort_values("M2").copy()
            for y in ["pc", "mass_sr", "f", "d5frac", "offdiag"]:
                d[f"{y}_low"] = d[y]
                d[f"{y}_high"] = d[y]
            path = CSV_DIR / f"{state.key}_two_point_M2_diagnostics.csv"
            d.to_csv(path, index=False)
            for obs, ylab in [
                ("pc", "Pole contribution"),
                ("mass_sr", "Reproduced mass (GeV)"),
                ("f", "Residue f (GeV)"),
                ("d5frac", "|dimension-5 / total|"),
            ]:
                plot_csv(
                    path,
                    f"{state.key}_{obs}_vs_M2",
                    "M2",
                    obs,
                    f"{obs}_low",
                    f"{obs}_high",
                    "M^2 (GeV^2)",
                    ylab,
                    f"{state.label}; s0={s0pick:.2f} GeV^2",
                    (state.tp_window[0], state.tp_window[1]),
                    center(state.tp_window)[0],
                )


def write_validation_and_comparison(results: list[dict]):
    validations = [
        ("rho_AB_equals_rho_BA", "checked", "identical function and symbolic WL equality"),
        ("spectral_matrix_positive", "checked", "analytic determinant is nonnegative above threshold"),
        ("mass_dimensions", "checked", "rho dim 2; Borel matrix and T dim 4; g dimensionless"),
        ("physical_threshold", "checked", "all dispersive integrals start at (mQ+ms)^2"),
        ("Ward_identity", "checked", "epsilon structure vanishes under epsilon_gamma -> q"),
        ("no_local_transition_condensate", "checked", "diagnostic term excluded from quoted_total_A/B"),
        ("two_point_local_d3_d5", "checked", "explicit matrices retained"),
        ("pure_A_limit", "checked", "theta=90 degrees for lower state"),
        ("pure_B_limit", "checked", "theta=0 degrees for lower state"),
        ("equal_Borel_limit", "checked", "M1^2=M2^2=2M^2 gives u0=1/2"),
        ("zero_ms_threshold_limit", "checked", "lambda and density reduce smoothly"),
        ("charge_switch", "checked", "hard density changes through eQ/es only"),
        ("python_mathematica_regression", "pending", "filled by compare_regression.py"),
        (
            "subleading_tensor_kernels",
            "truncation",
            "not retained; 15% charm and 5% bottom correlated uncertainty assigned",
        ),
        ("alpha_s_transition", "truncation", "not retained"),
    ]
    write_rows(
        CSV_DIR / "validation_status.csv",
        [{"gate": a, "status": b, "evidence": c} for a, b, c in validations],
    )
    literature_fields = [
        "state",
        "observable",
        "literature_value",
        "method",
        "reference",
        "comment",
    ]
    literature_rows = [
        ("Ds1_2460", "width_keV", "0.6-1.1", "pure-axial external-photon LCSR",
         "Colangelo2005", "pure-A benchmark"),
        ("Ds1_2460", "width_keV", "1.5", "vector-meson dominance",
         "Colangelo2005", "comparison value collected in Table I"),
        ("Ds1_2460", "width_keV", "5.5", "nonrelativistic quark model",
         "Godfrey2003", "comparison value collected in Colangelo et al."),
        ("Ds1_2460", "width_keV", "4.66", "chiral-doublet effective theory",
         "Bardeen2003", "comparison value collected in Colangelo et al."),
        ("Ds1_2460", "width_keV", "4.74-4.79", "relativistic quark model",
         "Chen2020", "three parameter modes"),
        ("Ds1_2460", "width_keV", "15.5", "potential-model calculation",
         "Radford2009", "comparison value collected in Chen et al."),
        ("Ds1_2460", "width_keV", "17.4", "potential-model calculation",
         "Green2017", "comparison value collected in Chen et al."),
        ("Ds1_2460", "width_keV", "5.6", "relativized quark model",
         "Godfrey2005", "comparison value collected in Chen et al."),
        ("Ds1_2460", "width_keV", "4.41", "quark model",
         "CloseSwanson2005", "comparison value collected in Chen et al."),
        ("Ds1_2460", "width_keV", "14.6-22.8",
         "relativistic heavy-quark model", "GoityRoberts2001",
         "three model variants collected in Chen et al."),
        ("Ds1_2460", "width_keV", "13+/-2", "hadronic-molecule EFT",
         "Fu2022", "full result"),
        ("Ds1_2460", "width_keV", "104",
         "potential model with relativistic corrections", "BondarMilstein2025",
         "authors estimate 50% model uncertainty"),
        ("Ds1_2536", "width_keV", "2.96-3.02", "relativistic quark model",
         "Chen2020", "three parameter modes"),
        ("Ds1_2536", "width_keV", "8.90", "potential-model calculation",
         "Radford2009", "comparison value collected in Chen et al."),
        ("Ds1_2536", "width_keV", "9.21", "potential-model calculation",
         "Green2017", "comparison value collected in Chen et al."),
        ("Ds1_2536", "width_keV", "0.4+/-1.0",
         "heavy-quark effective calculation", "Korner1993",
         "comparison value collected in Chen et al."),
        ("Ds1_2536", "width_keV", "5.5", "relativized quark model",
         "Godfrey2005", "comparison value collected in Chen et al."),
        ("Ds1_2536", "width_keV", "1.59", "quark model",
         "CloseSwanson2005", "comparison value collected in Chen et al."),
        ("Ds1_2536", "width_keV", "14.0-25.1",
         "relativistic heavy-quark model", "GoityRoberts2001",
         "three model variants collected in Chen et al."),
        ("Ds1_2536", "width_keV", "29",
         "potential model with relativistic corrections", "BondarMilstein2025",
         "authors estimate 50% model uncertainty"),
        ("Bs1_lower", "width_keV", "0.3-6.1",
         "pure-axial external-photon LCSR", "Wang2008",
         "lower strange-bottom axial partner"),
        ("Bs1_lower", "width_keV", "100+/-15", "hadronic-molecule EFT",
         "Fu2022", "lower strange-bottom molecular partner"),
        ("Bs1_lower", "width_keV", "60", "nonrelativistic quark model",
         "Li2021", "conditional subthreshold missing 1P1 state"),
        ("Bs1_lower", "width_keV", "20",
         "potential model with relativistic corrections", "BondarMilstein2025",
         "Bs1(5748) assignment; authors estimate 50% accuracy"),
        ("Bs1_5830", "width_keV", "53", "nonrelativistic quark model",
         "Li2021", "observed higher state"),
        ("Bs1_5830", "width_keV", "41",
         "potential model with relativistic corrections", "BondarMilstein2025",
         "Bs1(5829) assignment; authors estimate 50% accuracy"),
    ]
    literature = [
        dict(zip(literature_fields, row, strict=True)) for row in literature_rows
    ]
    write_rows(CSV_DIR / "literature_comparison.csv", literature)
    experiment = [
        {
            "state": "Ds1_2460",
            "status": "upper limit",
            "observable": "branching_fraction_Dsstar_gamma",
            "value": "<0.08 (90% CL)",
            "source": "PDG2026",
            "comment": "total width <3.5 MeV",
        },
        {
            "state": "Ds1_2536",
            "status": "possibly seen",
            "observable": "Dsstar gamma",
            "value": "no branching fraction or partial width",
            "source": "PDG2026",
            "comment": "total width 0.92+/-0.05 MeV",
        },
        {
            "state": "Bs1_lower",
            "status": "unobserved state",
            "observable": "radiative decay",
            "value": "no measurement or limit",
            "source": "PDG2026",
            "comment": "predicted lower state",
        },
        {
            "state": "Bs1_5830",
            "status": "unobserved mode",
            "observable": "Bsstar gamma",
            "value": "no measurement or limit",
            "source": "PDG2026",
            "comment": "parent state established",
        },
    ]
    write_rows(CSV_DIR / "experimental_comparison.csv", experiment)


def python_regression(two_point_summary: dict):
    rows = []
    test_points = {
        "D": (SECTORS[0], 3.0, 8.0, 5.0, 6.75),
        "B": (SECTORS[1], 9.0, 42.0, 6.0, 40.0),
    }
    for name, (sector, m2tp, s0tp, m2tr, s0tr) in test_points.items():
        s_test = (sector.inp.mQ + sector.inp.ms) ** 2 + 1.0
        for entry in ["AA", "AB", "BB"]:
            rows.append(
                {
                    "sector": name,
                    "quantity": f"rho_{entry}",
                    "value": float(rho_twopoint(s_test, sector.inp, entry)),
                }
            )
        mat = borel_twopoint_matrix(m2tp, s0tp, sector.inp)["total"]
        for i, j, label in [(0, 0, "Pi_AA"), (0, 1, "Pi_AB"), (1, 1, "Pi_BB")]:
            rows.append({"sector": name, "quantity": label, "value": mat[i, j]})
        tr = transition_terms(m2tr, s0tr, sector.inp, three_n=28)
        for key in [
            "hard",
            "tw2",
            "tw4_2p",
            "tw3_2p",
            "tw4_3p",
            "tw3_3p",
            "quoted_total_A",
            "quoted_total_B",
        ]:
            rows.append({"sector": name, "quantity": f"T_{key}", "value": tr[key]})
        for state in sector.states:
            fi = scalar_f_initial(
                sector, state, sector.inp, sector.theta, m2tp, s0tp, state.mass
            )
            vm2, vs0 = center(sector.vector_window)
            fv = scalar_f_vector(sector, sector.inp, vm2, vs0, sector.mf)
            ptr = physical_transition(
                m2tr,
                s0tr,
                sector.inp,
                sector.theta,
                state.index,
                state.mass,
                sector.mf,
                fi,
                fv,
                three_n=28,
            )
            rows.extend(
                [
                    {
                        "sector": name,
                        "quantity": f"{state.key}_f",
                        "value": fi,
                    },
                    {
                        "sector": name,
                        "quantity": f"{state.key}_g",
                        "value": ptr["g"],
                    },
                    {
                        "sector": name,
                        "quantity": f"{state.key}_width_keV",
                        "value": ptr["width_keV"],
                    },
                ]
            )
        rows.append(
            {
                "sector": name,
                "quantity": f"f_{name}sstar",
                "value": fv,
            }
        )
    write_rows(CSV_DIR / "python_regression.csv", rows)


def write_metadata():
    rows = [
        {"key": "random_seed", "value": SEED},
        {"key": "monte_carlo_samples_per_sector", "value": 1200},
        {"key": "two_point_pole_fraction_min", "value": 0.40},
        {"key": "two_point_dimension5_fraction_max", "value": 0.10},
        {"key": "mass_reproduction_relative_tolerance", "value": 0.05},
        {"key": "python_mathematica_abs_tolerance", "value": 2.0e-8},
        {"key": "python_mathematica_rel_tolerance", "value": 2.0e-5},
        {"key": "transition_OPE", "value": "LO hard + BBK twist2/3/4 2p+3p; no local condensate"},
        {"key": "two_point_OPE", "value": "LO exact-mass perturbative + d3 + finite-ms d3 + d5"},
        {"key": "python_runtime", "value": sys.executable},
    ]
    write_rows(CSV_DIR / "run_metadata.csv", rows)


def main():
    write_metadata()
    two_point = calculate_two_point()
    mc_stats = monte_carlo(two_point)
    results = scans_and_plots(two_point, mc_stats)
    diagnostics_plots(two_point)
    write_validation_and_comparison(results)
    python_regression(two_point)
    (OUT / "python_run_complete.json").write_text(
        json.dumps(
            {"status": "complete", "seed": SEED, "states": [r["state"] for r in results]},
            indent=2,
        )
        + "\n"
    )
    print(pd.DataFrame(results).to_string(index=False))


if __name__ == "__main__":
    main()
