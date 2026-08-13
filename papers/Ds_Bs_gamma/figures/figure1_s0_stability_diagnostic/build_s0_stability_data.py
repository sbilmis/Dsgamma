"""Build an internal s0-stability diagnostic for Figure 1.

This is not intended as a default manuscript figure.  It complements the
published Borel-stability plot by scanning the continuum threshold at fixed
representative Borel values.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from final_stage2_uncertainty_scan import precompute_f1_basis
from publication_stability_selected_windows import BS5830_WINDOW, BSLOW_WINDOW
from redo_stability_windows import bs_physical_row, ds_physical_row
from stage2_axial_g1_three_particle import F1_integral


OUT = HERE / "figure1_s0_stability_diagnostic.csv"

SPECS = [
    {
        "sector": "Ds",
        "state": "Ds1_2460",
        "state_label": "D_{s1}(2460)",
        "s0_range": (8.5, 9.5),
        "M2_values": [3.0, 3.75, 4.5],
        "value_key": "g_2460",
    },
    {
        "sector": "Ds",
        "state": "Ds1_2536",
        "state_label": "D_{s1}(2536)",
        "s0_range": (9.0, 10.0),
        "M2_values": [3.0, 3.75, 4.5],
        "value_key": "g_2536",
    },
    {
        "sector": "Bs",
        "state": "Bs1_5750",
        "state_label": "B_{s1}(5750)",
        "target": {**BSLOW_WINDOW, "s0_lines": (39.0, 40.0, 41.0)},
        "s0_range": (39.0, 41.0),
        "M2_values": [10.0, 12.0, 14.0],
        "value_key": "g_quoted",
    },
    {
        "sector": "Bs",
        "state": "Bs1_5830",
        "state_label": "B_{s1}(5830)",
        "target": {**BS5830_WINDOW, "s0_lines": (40.0, 41.0, 42.0)},
        "s0_range": (40.0, 42.0),
        "M2_values": [10.0, 12.0, 14.0],
        "value_key": "g_quoted",
    },
]


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


def main() -> None:
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_basis = precompute_f1_basis()
    rows: list[dict[str, float | str]] = []

    for spec in SPECS:
        s0_values = np.linspace(float(spec["s0_range"][0]), float(spec["s0_range"][1]), 41)
        for m2 in spec["M2_values"]:
            for s0 in s0_values:
                if spec["sector"] == "Ds":
                    row = ds_physical_row(float(m2), float(s0), f1_axial, f1_basis)
                else:
                    row = bs_physical_row(spec["target"], float(m2), float(s0), f1_axial, f1_basis)
                value = row[str(spec["value_key"])]
                rows.append(
                    {
                        "sector": spec["sector"],
                        "state": spec["state"],
                        "state_label": spec["state_label"],
                        "scan": "s0",
                        "M2": float(m2),
                        "s0": float(s0),
                        "g_abs": abs(float(value)),
                    }
                )

    write_csv(OUT, rows)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
