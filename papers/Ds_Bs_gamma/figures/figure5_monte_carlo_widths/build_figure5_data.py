"""Build Figure 5 Monte Carlo width distributions.

This figure uses the preferred lattice-normalized photon input and varies
transition inputs, two-point residue samples, Borel parameters, continuum
thresholds, and the mixing angle over 25--45 degrees.
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPTS = ROOT / "scripts"
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPTS))

from final_window_monte_carlo import (
    BS_WINDOWS,
    DS_WINDOW,
    evaluate_bs_physical,
    evaluate_ds_physical,
    sample_bs_inputs,
)
from lattice_photon_normalization_comparison import sample_inputs as sample_ds_inputs
from mixing_angle_inputs import THETA_SENSITIVITY_MAX_DEG, THETA_SENSITIVITY_MIN_DEG
from rohrwild_transition_exact import precompute_convolutions
from twopoint_ds1_matrix_sumrule import Inputs as TwoPointInputs
from twopoint_ds1_matrix_sumrule import projected_ope as projected_twopoint_ope


SEED = 20260812
N_POINTS = 1000
SCENARIO = "lattice_fperp_s"
ENSEMBLE = "theta_sensitivity_25_45"
OUT_CSV = HERE / "figure5_monte_carlo_widths.csv"
SUMMARY_CSV = HERE / "figure5_monte_carlo_summary.csv"


def read_accepted(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = [row for row in csv.DictReader(handle) if int(row.get("accepted", "1")) == 1]
    if not rows:
        raise RuntimeError(f"No accepted rows in {path}")
    return rows


def project_ds_residue_sample(sample: dict[str, str], theta_deg: float) -> dict[str, str]:
    """Return a residue sample projected at theta_deg for Ds."""
    inp = TwoPointInputs(
        mc=float(sample["mc_GeV"]),
        ms=float(sample["ms_GeV"]),
        qq=float(sample["ss_GeV3"]),
        kappa_s=1.0,
        m0_sq=float(sample["m0_sq_GeV2"]),
        mass_low=2.4595,
        mass_high=2.53511,
    )
    m2 = float(sample["M2_GeV2"])
    pi1, _ = projected_twopoint_ope(m2, float(sample["s01_GeV2"]), theta_deg, 0, inp)
    pi2, _ = projected_twopoint_ope(m2, float(sample["s02_GeV2"]), theta_deg, 1, inp)
    if pi1 <= 0.0 or pi2 <= 0.0:
        raise RuntimeError("Nonpositive Ds projected two-point OPE")
    out = dict(sample)
    out["theta_deg"] = str(theta_deg)
    out["f1_GeV"] = str(math.sqrt(pi1 * math.exp(inp.mass_low**2 / m2) / inp.mass_low**2))
    out["f2_GeV"] = str(math.sqrt(pi2 * math.exp(inp.mass_high**2 / m2) / inp.mass_high**2))
    return out


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize(values: list[float]) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    return {
        "median": float(np.percentile(arr, 50.0)),
        "p16": float(np.percentile(arr, 16.0)),
        "p84": float(np.percentile(arr, 84.0)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def main() -> None:
    rng = np.random.default_rng(SEED)
    convolutions = precompute_convolutions()
    ds_samples = read_accepted(ROOT / "outputs" / "twopoint_ds1_matrix_mc.csv")
    bs_samples = read_accepted(ROOT / "outputs" / "twopoint_bs1_matrix_mc.csv")
    rows: list[dict[str, float | str]] = []
    rejected = 0

    for _ in range(N_POINTS):
        theta_deg = float(rng.uniform(THETA_SENSITIVITY_MIN_DEG, THETA_SENSITIVITY_MAX_DEG))
        vals = sample_ds_inputs(rng, SCENARIO)
        residue_sample = ds_samples[int(rng.integers(0, len(ds_samples)))]
        residue_sample = project_ds_residue_sample(residue_sample, theta_deg)
        rows.extend(
            evaluate_ds_physical(
                vals,
                float(rng.uniform(DS_WINDOW["M2_min"], DS_WINDOW["M2_max"])),
                float(rng.uniform(DS_WINDOW["s0_min"], DS_WINDOW["s0_max"])),
                residue_sample,
                SCENARIO,
                ENSEMBLE,
                convolutions,
            )
        )

    for target in [w for w in BS_WINDOWS if w["window_id"] == "central_10_14"]:
        for _ in range(N_POINTS):
            theta_deg = float(rng.uniform(THETA_SENSITIVITY_MIN_DEG, THETA_SENSITIVITY_MAX_DEG))
            vals = sample_bs_inputs(rng, target, SCENARIO)
            residue_sample = bs_samples[int(rng.integers(0, len(bs_samples)))]
            row = evaluate_bs_physical(
                vals,
                target,
                residue_sample,
                float(rng.uniform(float(target["M2_min"]), float(target["M2_max"]))),
                float(rng.uniform(float(target["s0_min"]), float(target["s0_max"]))),
                theta_deg,
                SCENARIO,
                ENSEMBLE,
                convolutions,
            )
            if row is None:
                rejected += 1
            else:
                rows.append(row)

    write_csv(OUT_CSV, rows)
    summary_rows: list[dict[str, float | str]] = []
    for state_key in ("Ds1_2460", "Ds1_2536", "Bs1_5750", "Bs1_5830"):
        subset = [row for row in rows if row["state_key"] == state_key]
        gamma = summarize([float(row["Gamma_keV"]) for row in subset])
        g = summarize([float(row["g_GeV_inv"]) for row in subset])
        summary_rows.append(
            {
                "state_key": state_key,
                "n": len(subset),
                "g_median_GeV_inv": g["median"],
                "g_p16_GeV_inv": g["p16"],
                "g_p84_GeV_inv": g["p84"],
                "Gamma_median_keV": gamma["median"],
                "Gamma_p16_keV": gamma["p16"],
                "Gamma_p84_keV": gamma["p84"],
                "Gamma_mean_keV": gamma["mean"],
                "Gamma_std_keV": gamma["std"],
                "Gamma_min_keV": gamma["min"],
                "Gamma_max_keV": gamma["max"],
            }
        )
    write_csv(SUMMARY_CSV, summary_rows)
    print(f"Wrote {OUT_CSV}")
    print(f"Wrote {SUMMARY_CSV}")
    print(f"Rejected bottom samples: {rejected}")
    for row in summary_rows:
        print(
            "{state_key}: Gamma={m:.4g} [{lo:.4g},{hi:.4g}] keV; mean={mean:.4g}, std={std:.4g}, n={n}".format(
                state_key=row["state_key"],
                m=float(row["Gamma_median_keV"]),
                lo=float(row["Gamma_p16_keV"]),
                hi=float(row["Gamma_p84_keV"]),
                mean=float(row["Gamma_mean_keV"]),
                std=float(row["Gamma_std_keV"]),
                n=int(row["n"]),
            )
        )


if __name__ == "__main__":
    main()
