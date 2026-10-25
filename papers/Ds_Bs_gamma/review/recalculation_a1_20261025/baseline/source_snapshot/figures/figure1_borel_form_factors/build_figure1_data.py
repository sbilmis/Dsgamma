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

from corrected_transition_central_python import (
    BS_TWOPOINT_REFERENCE,
    DS_REFERENCE,
    bs_direct_decay_constants,
    bs_inputs,
    ds_inputs,
)
from final_stage2_uncertainty_scan import set_photon_shape_params
from rohrwild_transition_exact import (
    physical_couplings,
    precompute_convolutions,
    transition_invariants,
)


OUT = HERE / "figure1_borel_form_factors.csv"

SPECS = [
    {
        "sector": "Ds",
        "state": "Ds1_2460",
        "state_label": "D_{s1}(2460)",
        "s0_values": [8.5, 9.0, 9.5],
        "M2_range": (3.0, 4.5),
        "value_key": "g_1",
    },
    {
        "sector": "Ds",
        "state": "Ds1_2536",
        "state_label": "D_{s1}(2536)",
        "s0_values": [9.0, 9.5, 10.0],
        "M2_range": (3.0, 4.5),
        "value_key": "g_2",
    },
    {
        "sector": "Bs",
        "state": "Bs1_5750",
        "state_label": "B_{s1}(5750)",
        "s0_values": [39.0, 40.0, 41.0],
        "M2_range": (10.0, 14.0),
        "combo": "low",
        "mass": 5.750,
    },
    {
        "sector": "Bs",
        "state": "Bs1_5830",
        "state_label": "B_{s1}(5830)",
        "s0_values": [40.0, 41.0, 42.0],
        "M2_range": (10.0, 14.0),
        "combo": "high",
        "mass": 5.82870,
    },
]


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    set_photon_shape_params(-2.1, 3.8)
    convolutions = precompute_convolutions()
    ds_vals = ds_inputs()
    bs_f1, bs_f2 = bs_direct_decay_constants()
    rows: list[dict[str, float | str]] = []

    for spec in SPECS:
        for s0 in spec["s0_values"]:
            for m2 in np.linspace(*spec["M2_range"], 31):
                if spec["sector"] == "Ds":
                    invariants = transition_invariants(
                        float(m2), float(s0), ds_vals, convolutions=convolutions
                    )
                    physical = physical_couplings(
                        invariants,
                        theta_deg=float(DS_REFERENCE["theta_deg"]),
                        m_state_1=float(DS_REFERENCE["m1"]),
                        m_state_2=float(DS_REFERENCE["m2"]),
                        f_1=float(DS_REFERENCE["f1"]),
                        f_2=float(DS_REFERENCE["f2"]),
                        m_p=ds_vals["m_ds"],
                        f_p=ds_vals["f_ds"],
                        m_q=ds_vals["mc"],
                        m_s=ds_vals["ms"],
                    )
                    value = physical[str(spec["value_key"])]
                    width = physical[
                        "Gamma_1_keV" if spec["value_key"] == "g_1" else "Gamma_2_keV"
                    ]
                else:
                    vals = bs_inputs(float(spec["mass"]))
                    invariants = transition_invariants(
                        float(m2), float(s0), vals, convolutions=convolutions
                    )
                    physical = physical_couplings(
                        invariants,
                        theta_deg=float(BS_TWOPOINT_REFERENCE["theta_deg"]),
                        m_state_1=float(spec["mass"]),
                        m_state_2=float(spec["mass"]),
                        f_1=bs_f1,
                        f_2=bs_f2,
                        m_p=vals["m_ds"],
                        f_p=vals["f_ds"],
                        m_q=vals["mc"],
                        m_s=vals["ms"],
                    )
                    value = physical["g_1" if spec["combo"] == "low" else "g_2"]
                    width = physical[
                        "Gamma_1_keV" if spec["combo"] == "low" else "Gamma_2_keV"
                    ]
                rows.append(
                    {
                        "sector": spec["sector"],
                        "state": spec["state"],
                        "state_label": spec["state_label"],
                        "scan": "M2",
                        "M2": float(m2),
                        "s0": float(s0),
                        "fixed_s0": float(s0),
                        "theta_deg": float(
                            DS_REFERENCE["theta_deg"]
                            if spec["sector"] == "Ds"
                            else BS_TWOPOINT_REFERENCE["theta_deg"]
                        ),
                        "g_signed": float(value),
                        "g_abs": abs(float(value)),
                        "Gamma_keV": float(width),
                    }
                )

    write_csv(OUT, rows)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
