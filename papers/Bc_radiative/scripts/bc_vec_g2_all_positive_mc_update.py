#!/usr/bin/env python3
"""Fold the explicit all-positive vector G2 correction into MC samples."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def stats(vals):
    arr = np.asarray(vals, dtype=float)
    return {
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
        "q16": float(np.quantile(arr, 0.16)),
        "q84": float(np.quantile(arr, 0.84)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def main() -> None:
    rng = np.random.default_rng(20260712)
    mc = [
        r for r in read_csv(OUT / "bc_vec_complete_monte_carlo.csv")
        if r["normalization"] == "standard_vector_invariant"
    ]
    g2_grid = read_csv(OUT / "bc_vec_g2_all_positive_direct_borel_grid.csv")
    ratios = {
        ch: np.asarray([float(r["delta_g_over_g"]) for r in g2_grid if r["channel"] == ch], dtype=float)
        for ch in sorted({r["channel"] for r in g2_grid})
    }

    rows = []
    for r in mc:
        for channel, gkey, wkey in [
            ("Bc1(6743)->Bc* gamma", "g1_GeV_inv", "Gamma1_keV"),
            ("Bc1(6750)->Bc* gamma", "g2_GeV_inv", "Gamma2_keV"),
        ]:
            delta = float(rng.choice(ratios[channel]))
            g_pert = float(r[gkey])
            gamma_pert = float(r[wkey])
            g_total = g_pert * (1.0 + delta)
            gamma_total = gamma_pert * (1.0 + delta) ** 2
            rows.append(
                {
                    "source_index": r["index"],
                    "channel": channel,
                    "delta_g_over_g_G2ap": delta,
                    "g_pert_GeV_inv": g_pert,
                    "Gamma_pert_keV": gamma_pert,
                    "g_pert_plus_G2ap_GeV_inv": g_total,
                    "Gamma_pert_plus_G2ap_keV": gamma_total,
                }
            )

    write_csv(OUT / "bc_vec_g2_all_positive_monte_carlo.csv", rows)
    summary_rows = []
    for channel in sorted({r["channel"] for r in rows}):
        sub = [r for r in rows if r["channel"] == channel]
        for obs in [
            "delta_g_over_g_G2ap",
            "g_pert_plus_G2ap_GeV_inv",
            "Gamma_pert_plus_G2ap_keV",
        ]:
            summary_rows.append({"channel": channel, "observable": obs, **stats([float(r[obs]) for r in sub])})
    write_csv(OUT / "bc_vec_g2_all_positive_monte_carlo_summary.csv", summary_rows)

    lines = [
        "Bc1 -> Bc* gamma MC updated with explicit all-positive G2",
        "==========================================================",
        "",
        "The nine-point direct-Borel all-positive G2 delta g/g values are sampled",
        "uniformly and folded into the existing standard-vector MC samples.",
        "",
    ]
    for row in summary_rows:
        if row["observable"] == "Gamma_pert_plus_G2ap_keV":
            lines.append(
                "{channel}: Gamma(pert + G2_ap) = {median:.4g} "
                "[{q16:.4g}, {q84:.4g}] keV".format(**row)
            )
        if row["observable"] == "delta_g_over_g_G2ap":
            lines.append(
                "{channel}: delta g/g = {median:.4g} [{q16:.4g}, {q84:.4g}]".format(**row)
            )
    (OUT / "bc_vec_g2_all_positive_monte_carlo_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
