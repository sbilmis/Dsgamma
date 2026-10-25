"""Build the preferred-prior Monte Carlo width distributions.

This figure uses the preferred lattice-normalized photon input and the same
sector-specific Gaussian mixing-angle priors as the final results:
theta_Ds = 26.6 +/- 0.6 degrees and theta_Bs = 38.5 +/- 0.1 degrees.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCENARIO = "lattice_fperp_s"
ENSEMBLE = "theta_prior_gaussian"
SOURCE_CSV = ROOT / "outputs" / "final_window_mc_scan.csv"
OUT_CSV = HERE / "figure5_monte_carlo_widths.csv"
SUMMARY_CSV = HERE / "figure5_monte_carlo_summary.csv"


def read_accepted(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = [row for row in csv.DictReader(handle) if int(row.get("accepted", "1")) == 1]
    if not rows:
        raise RuntimeError(f"No accepted rows in {path}")
    return rows


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
    source_rows = read_accepted(SOURCE_CSV)
    rows = [
        row
        for row in source_rows
        if row["scenario"] == SCENARIO
        and row["ensemble"] == ENSEMBLE
        and (
            (row["sector"] == "Ds" and row["window_id"] == "central")
            or (row["sector"] == "Bs" and row["window_id"] == "central_10_14")
        )
    ]
    if not rows:
        raise RuntimeError("No preferred-prior rows found in the final scan")

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
