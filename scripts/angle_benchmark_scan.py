"""Angle-only benchmark scans for radiative-width comparisons.

The purpose is diagnostic: keep the sampled QCD inputs fixed, reconstruct the
two unmixed reduced amplitudes from the mixed physical amplitudes, and rotate
them with a trial mixing angle.  This tests whether discrepancies with selected
literature values can be removed by changing theta alone.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
DSSTAR_OUT = ROOT / "Dsstar_gamma" / "outputs"
BSSTAR_OUT = ROOT / "Bsstar_gamma" / "outputs"

ALPHA_EM = 1.0 / 137.036


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def width_p_to_p_keV(m_initial: float, m_final: float, g: np.ndarray) -> np.ndarray:
    qgamma = (m_initial**2 - m_final**2) / (2.0 * m_initial)
    return ALPHA_EM / 3.0 * g**2 * qgamma**3 * 1.0e6


def width_p_to_v_keV(m_initial: float, m_final: float, g: np.ndarray) -> np.ndarray:
    delta = m_initial**2 - m_final**2
    h33 = delta**2 * (m_initial**2 + m_final**2) / (6.0 * m_initial**2 * m_final**2)
    phase = delta / (16.0 * math.pi * m_initial**3)
    return 4.0 * math.pi * ALPHA_EM * phase * h33 * g**2 * 1.0e6


def stats(values: np.ndarray) -> dict[str, float]:
    return {
        "median": float(np.percentile(values, 50.0)),
        "p16": float(np.percentile(values, 16.0)),
        "p84": float(np.percentile(values, 84.0)),
    }


def reconstruct_ab(theta_deg: np.ndarray | float, f1: np.ndarray, f2: np.ndarray, g1: np.ndarray, g2: np.ndarray):
    theta = np.deg2rad(theta_deg)
    s = np.sin(theta)
    c = np.cos(theta)
    y1 = f1 * g1
    y2 = f2 * g2
    amp_a = s * y1 + c * y2
    amp_b = c * y1 - s * y2
    return amp_a, amp_b


def rotate(amp_a: np.ndarray, amp_b: np.ndarray, f1: np.ndarray, f2: np.ndarray, theta_deg: float):
    theta = math.radians(theta_deg)
    s = math.sin(theta)
    c = math.cos(theta)
    g1 = (s * amp_a + c * amp_b) / f1
    g2 = (c * amp_a - s * amp_b) / f2
    return g1, g2


def best_single(grid: np.ndarray, medians: np.ndarray, target: float) -> tuple[float, float]:
    idx = int(np.argmin(np.abs(medians - target)))
    return float(grid[idx]), float(medians[idx])


def best_pair(grid: np.ndarray, med1: np.ndarray, target1: float, med2: np.ndarray, target2: float) -> tuple[float, float, float]:
    eps = 1.0e-12
    score = (np.log((med1 + eps) / target1)) ** 2 + (np.log((med2 + eps) / target2)) ** 2
    idx = int(np.argmin(score))
    return float(grid[idx]), float(med1[idx]), float(med2[idx])


def ds_pseudoscalar_dataset(scenario: str = "lattice_fperp_s"):
    rows = [
        r for r in read_csv(OUT / "final_window_mc_scan.csv")
        if r["sector"] == "Ds" and r["scenario"] == scenario and r["ensemble"] == "theta_fixed"
    ]
    groups: dict[tuple[str, str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    for r in rows:
        key = (r["M2"], r["s0"], r["scenario"])
        groups[key][r["state_key"]] = r
    f1 = []
    f2 = []
    g1 = []
    g2 = []
    theta = []
    for pair in groups.values():
        if "Ds1_2460" not in pair or "Ds1_2536" not in pair:
            continue
        r1 = pair["Ds1_2460"]
        r2 = pair["Ds1_2536"]
        f1.append(float(r1["f1_GeV"]))
        f2.append(float(r1["f2_GeV"]))
        g1.append(float(r1["g_GeV_inv"]))
        g2.append(float(r2["g_GeV_inv"]))
        theta.append(float(r1["theta_deg"]))
    amp_a, amp_b = reconstruct_ab(np.asarray(theta), np.asarray(f1), np.asarray(f2), np.asarray(g1), np.asarray(g2))
    return amp_a, amp_b, np.asarray(f1), np.asarray(f2)


def ds_vector_dataset(scenario: str = "lattice_fperp_s"):
    rows = [
        r for r in read_csv(DSSTAR_OUT / "tensor_full_calibrated_monte_carlo.csv")
        if r["scenario"] == scenario and r["ensemble"] == "theta_fixed"
    ]
    f1 = np.asarray([float(r["f1_GeV"]) for r in rows])
    f2 = np.asarray([float(r["f2_GeV"]) for r in rows])
    g1 = np.asarray([float(r["g1_2460_fullAB"]) for r in rows])
    g2 = np.asarray([float(r["g2_2536_fullAB"]) for r in rows])
    theta = np.asarray([float(r["theta_deg"]) for r in rows])
    amp_a, amp_b = reconstruct_ab(theta, f1, f2, g1, g2)
    return amp_a, amp_b, f1, f2


def bs_vector_dataset(state_key: str, scenario: str = "lattice_fperp_s"):
    rows = [
        r for r in read_csv(BSSTAR_OUT / "bsstar_full_calibrated_monte_carlo.csv")
        if r["state_key"] == state_key and r["scenario"] == scenario and r["ensemble"] == "theta_fixed"
    ]
    f1 = np.asarray([float(r["f1_GeV"]) for r in rows])
    f2 = np.asarray([float(r["f2_GeV"]) for r in rows])
    glow = np.asarray([float(r["g_low"]) for r in rows])
    ghigh = np.asarray([float(r["g_high"]) for r in rows])
    theta = np.asarray([float(r["theta_deg"]) for r in rows])
    amp_a, amp_b = reconstruct_ab(theta, f1, f2, glow, ghigh)
    return amp_a, amp_b, f1, f2


def scan_ds_pseudoscalar(grid: np.ndarray) -> list[dict[str, float | str]]:
    out = []
    for scenario in ("lattice_fperp_s", "legacy_chi_condensate"):
        amp_a, amp_b, f1, f2 = ds_pseudoscalar_dataset(scenario)
        med2460 = []
        med2536 = []
        for th in grid:
            g1, g2 = rotate(amp_a, amp_b, f1, f2, float(th))
            med2460.append(stats(width_p_to_p_keV(2.4595, 1.96834, g1))["median"])
            med2536.append(stats(width_p_to_p_keV(2.53511, 1.96834, g2))["median"])
        med2460 = np.asarray(med2460)
        med2536 = np.asarray(med2536)
        th_pair, v2460, v2536 = best_pair(grid, med2460, 297.0, med2536, 12.0)
        th_col, v_col = best_single(grid, med2460, 24.0)
        out.extend(
            [
                {"family": "Ds->Ds gamma", "scenario": scenario, "target": "Bondar pair 2460=297,2536=12", "theta_deg": th_pair, "Gamma2460_keV": v2460, "Gamma2536_keV": v2536},
                {"family": "Ds->Ds gamma", "scenario": scenario, "target": "Colangelo 2460 central 24", "theta_deg": th_col, "Gamma2460_keV": v_col, "Gamma2536_keV": float(med2536[np.argmin(np.abs(med2460 - 24.0))])},
            ]
        )
    return out


def scan_ds_vector(grid: np.ndarray) -> list[dict[str, float | str]]:
    out = []
    for scenario in ("lattice_fperp_s", "legacy_chi_condensate"):
        amp_a, amp_b, f1, f2 = ds_vector_dataset(scenario)
        med2460 = []
        med2536 = []
        for th in grid:
            g1, g2 = rotate(amp_a, amp_b, f1, f2, float(th))
            med2460.append(stats(width_p_to_v_keV(2.4595, 2.1122, g1))["median"])
            med2536.append(stats(width_p_to_v_keV(2.53511, 2.1122, g2))["median"])
        med2460 = np.asarray(med2460)
        med2536 = np.asarray(med2536)
        th_pair, v2460, v2536 = best_pair(grid, med2460, 104.0, med2536, 29.0)
        th_2460, v_2460 = best_single(grid, med2460, 104.0)
        th_2536, v_2536 = best_single(grid, med2536, 29.0)
        out.extend(
            [
                {"family": "Ds->Ds* gamma", "scenario": scenario, "target": "Bondar pair 2460=104,2536=29", "theta_deg": th_pair, "Gamma2460_keV": v2460, "Gamma2536_keV": v2536},
                {"family": "Ds->Ds* gamma", "scenario": scenario, "target": "Bondar 2460 only 104", "theta_deg": th_2460, "Gamma2460_keV": v_2460, "Gamma2536_keV": float(med2536[np.argmin(np.abs(med2460 - 104.0))])},
                {"family": "Ds->Ds* gamma", "scenario": scenario, "target": "Bondar 2536 only 29", "theta_deg": th_2536, "Gamma2460_keV": float(med2460[np.argmin(np.abs(med2536 - 29.0))]), "Gamma2536_keV": v_2536},
            ]
        )
    return out


def scan_bs_vector(grid: np.ndarray) -> list[dict[str, float | str]]:
    specs = [
        ("Bs1_5750", "B_{s1}(5750)->Bs* gamma", 5.750, 20.0, "low"),
        ("Bs1_5830", "B_{s1}(5830)->Bs* gamma", 5.82870, 41.0, "high"),
    ]
    out = []
    for scenario in ("lattice_fperp_s", "legacy_chi_condensate"):
        for state_key, label, mass, target, which in specs:
            amp_a, amp_b, f1, f2 = bs_vector_dataset(state_key, scenario)
            med = []
            companion = []
            for th in grid:
                glow, ghigh = rotate(amp_a, amp_b, f1, f2, float(th))
                if which == "low":
                    med.append(stats(width_p_to_v_keV(mass, 5.4154, glow))["median"])
                    companion.append(stats(width_p_to_v_keV(5.82870, 5.4154, ghigh))["median"])
                else:
                    med.append(stats(width_p_to_v_keV(mass, 5.4154, ghigh))["median"])
                    companion.append(stats(width_p_to_v_keV(5.750, 5.4154, glow))["median"])
            med = np.asarray(med)
            th, val = best_single(grid, med, target)
            out.append({"family": label, "scenario": scenario, "target": f"Bondar {target:g}", "theta_deg": th, "Gamma_keV": val})
    return out


def main() -> None:
    grid = np.linspace(0.0, 90.0, 901)
    rows = scan_ds_pseudoscalar(grid) + scan_ds_vector(grid) + scan_bs_vector(grid)
    path = OUT / "angle_benchmark_scan.csv"
    with path.open("w", newline="") as f:
        fieldnames: list[str] = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {path}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
