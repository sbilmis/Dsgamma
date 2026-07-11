"""Uncertainty diagnostics for Bs1 -> Bs* gamma.

Uses rank correlations in the final Monte Carlo sample.  This is a diagnostic
for which sampled/derived quantities control the spread; it is not a refit.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
BSSTAR_OUT = ROOT / "Bsstar_gamma" / "outputs"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def rankdata(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and x[order[j]] == x[order[i]]:
            j += 1
        if j - i > 1:
            ranks[order[i:j]] = 0.5 * (i + j - 1)
        i = j
    return ranks


def corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 4 or np.std(x) == 0.0 or np.std(y) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return corr(rankdata(x), rankdata(y))


def numeric_array(rows: list[dict[str, str]], key: str) -> np.ndarray:
    vals = []
    for r in rows:
        try:
            vals.append(float(r[key]))
        except (KeyError, ValueError):
            vals.append(float("nan"))
    return np.asarray(vals, dtype=float)


def summarize_group(rows: list[dict[str, str]], state_key: str, scenario: str, ensemble: str) -> list[dict[str, str | float]]:
    subset = [
        r for r in rows
        if r["state_key"] == state_key and r["scenario"] == scenario and r["ensemble"] == ensemble
    ]
    y = np.log(np.maximum(numeric_array(subset, "Gamma_quoted_keV"), 1.0e-12))
    variables = [
        ("M2", "input"),
        ("s0", "input"),
        ("theta_deg", "input"),
        ("fV_Bsstar_GeV", "input"),
        ("f1_GeV", "normalization"),
        ("f2_GeV", "normalization"),
        ("chi_effective", "photon input"),
        ("fperp_s_used_GeV", "photon input"),
        ("omegaA", "photon shape"),
        ("omegaV", "photon shape"),
        ("gB_soft", "derived amplitude"),
        ("gB_hard", "derived amplitude"),
        ("gB_total", "derived amplitude"),
        ("gA_star", "derived amplitude"),
    ]
    out = []
    for var, kind in variables:
        x = numeric_array(subset, var)
        mask = np.isfinite(x) & np.isfinite(y)
        if int(np.sum(mask)) < 8:
            continue
        rho = spearman(x[mask], y[mask])
        pear = corr(x[mask], y[mask])
        if not math.isfinite(rho):
            continue
        out.append(
            {
                "state_key": state_key,
                "scenario": scenario,
                "ensemble": ensemble,
                "quantity": var,
                "kind": kind,
                "spearman_abs": abs(rho),
                "spearman": rho,
                "pearson_logGamma": pear,
                "n": int(np.sum(mask)),
            }
        )
    out.sort(key=lambda r: float(r["spearman_abs"]), reverse=True)
    return out


def main() -> None:
    rows = read_csv(BSSTAR_OUT / "bsstar_full_calibrated_monte_carlo.csv")
    all_rows = []
    for state_key in ("Bs1_5750", "Bs1_5830"):
        for scenario in ("lattice_fperp_s", "legacy_chi_condensate"):
            all_rows.extend(summarize_group(rows, state_key, scenario, "theta_fixed"))

    csv_path = OUT / "bsstar_uncertainty_diagnostics.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)

    lines = ["# Bs* Final-State Uncertainty Diagnostics", ""]
    for state_key in ("Bs1_5750", "Bs1_5830"):
        for scenario in ("lattice_fperp_s", "legacy_chi_condensate"):
            subset = [
                r for r in all_rows
                if r["state_key"] == state_key and r["scenario"] == scenario
            ][:8]
            lines.append(f"## {state_key}, {scenario}, theta fixed")
            lines.append("")
            lines.append("| rank | quantity | class | |rho_S| | rho_S |")
            lines.append("|---:|---|---|---:|---:|")
            for i, row in enumerate(subset, 1):
                lines.append(
                    f"| {i} | `{row['quantity']}` | {row['kind']} | "
                    f"{float(row['spearman_abs']):.3f} | {float(row['spearman']):+.3f} |"
                )
            lines.append("")
    lines.extend(
        [
            "Interpretation: for the preferred lattice-photon scenario the strange-photon normalization is no longer the dominant uncertainty by itself; the spread is shared with the physical-current normalization constants and the Borel/threshold window.  In the legacy scenario, the effective susceptibility/condensate normalization remains visibly important.",
            "",
            "The lattice input used here is the strange photon coupling \(f_{\\gamma,s}^{\\perp}(1\\,\\mathrm{GeV})\).  It is not a bottom-specific lattice calculation; it is the same strange-quark photon-DA normalization that enters both charm-strange and bottom-strange channels.",
        ]
    )
    note_path = ROOT / "notes" / "bsstar_uncertainty_diagnostics.md"
    note_path.write_text("\n".join(lines) + "\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {note_path}")


if __name__ == "__main__":
    main()
