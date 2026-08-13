"""Build the Figure 1 Borel-stability data with panel-specific s0 values."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from publication_stability_selected_windows import BS5830_WINDOW, BSLOW_WINDOW
from redo_stability_windows import bs_physical_row, ds_physical_row
from stage2_axial_g1_three_particle import F1_integral
from final_stage2_uncertainty_scan import precompute_f1_basis


OUT = HERE / "figure1_borel_form_factors.csv"

SPECS = [
    {
        "sector": "Ds",
        "state": "Ds1_2460",
        "state_label": "D_{s1}(2460)",
        "s0_values": [8.5, 9.0, 9.5],
        "M2_range": (3.0, 4.5),
        "value_key": "g_2460",
    },
    {
        "sector": "Ds",
        "state": "Ds1_2536",
        "state_label": "D_{s1}(2536)",
        "s0_values": [9.0, 9.5, 10.0],
        "M2_range": (3.0, 4.5),
        "value_key": "g_2536",
    },
    {
        "sector": "Bs",
        "state": "Bs1_5750",
        "state_label": "B_{s1}(5750)",
        "target": {**BSLOW_WINDOW, "s0_lines": (39.0, 40.0, 41.0)},
        "s0_values": [39.0, 40.0, 41.0],
        "M2_range": (10.0, 14.0),
        "value_key": "g_quoted",
    },
    {
        "sector": "Bs",
        "state": "Bs1_5830",
        "state_label": "B_{s1}(5830)",
        "target": {**BS5830_WINDOW, "s0_lines": (40.0, 41.0, 42.0)},
        "s0_values": [40.0, 41.0, 42.0],
        "M2_range": (10.0, 14.0),
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
        for s0 in spec["s0_values"]:
            for m2 in np.linspace(*spec["M2_range"], 31):
                if spec["sector"] == "Ds":
                    row = ds_physical_row(float(m2), float(s0), f1_axial, f1_basis)
                    value = row[str(spec["value_key"])]
                else:
                    row = bs_physical_row(
                        spec["target"], float(m2), float(s0), f1_axial, f1_basis
                    )
                    value = row[str(spec["value_key"])]
                rows.append(
                    {
                        "sector": spec["sector"],
                        "state": spec["state"],
                        "state_label": spec["state_label"],
                        "scan": "M2",
                        "M2": float(m2),
                        "s0": float(s0),
                        "fixed_s0": float(s0),
                        "g_abs": abs(float(value)),
                    }
                )

    write_csv(OUT, rows)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
