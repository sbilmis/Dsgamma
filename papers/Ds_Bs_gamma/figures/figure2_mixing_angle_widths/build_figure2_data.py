"""Build Figure 2 angle-dependence data.

The plotted curves use central lattice-normalized photon inputs and the updated
central Borel/threshold choices.  For each value of the mixing angle, the
physical-current residues are projected from the accepted AA/AB/BB two-point
matrix samples, matching the normalization convention used in the final-window
uncertainty scan.
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
OUT = ROOT / "outputs"
SCRIPTS = ROOT / "scripts"
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPTS))

from final_window_monte_carlo import sample_bs_inputs
from lattice_photon_normalization_comparison import central_inputs_for_scenario
from mixing_angle_inputs import THETA_HQET_DEG
from rohrwild_transition_exact import physical_couplings, precompute_convolutions, transition_invariants
from twopoint_ds1_matrix_sumrule import Inputs as TwoPointInputs
from twopoint_ds1_matrix_sumrule import projected_ope as projected_twopoint_ope


THETA_MIN = 25.0
THETA_MAX = 45.0
N_THETA = 161

DS_SPECS = [
    {
        "sector": "Ds",
        "state": "Ds1_2460",
        "label": "D_s1(2460)",
        "M2": 3.75,
        "s0": 9.0,
        "combo": "low",
    },
    {
        "sector": "Ds",
        "state": "Ds1_2536",
        "label": "D_s1(2536)",
        "M2": 3.75,
        "s0": 9.5,
        "combo": "high",
    },
]

BS_SPECS = [
    {
        "sector": "Bs",
        "state": "Bs1_5750",
        "label": "B_s1(5750)",
        "M2": 12.0,
        "s0": 40.0,
        "combo": "low",
        "m_initial": 5.750,
        "m_sigma": 0.0,
        "quoted_combo": "low",
        "window_id": "central_10_14",
    },
    {
        "sector": "Bs",
        "state": "Bs1_5830",
        "label": "B_s1(5830)",
        "M2": 12.0,
        "s0": 41.0,
        "combo": "high",
        "m_initial": 5.82870,
        "m_sigma": 0.0,
        "quoted_combo": "high",
        "window_id": "central_10_14",
    },
]


def read_accepted(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = [row for row in csv.DictReader(handle) if int(row.get("accepted", "1")) == 1]
    if not rows:
        raise RuntimeError(f"No accepted two-point samples in {path}")
    return rows


def representative_sample(rows: list[dict[str, str]], theta_key: str = "theta_deg") -> dict[str, str]:
    """Choose a central accepted row instead of mixing unrelated medians."""
    theta_values = np.asarray([float(row[theta_key]) for row in rows], dtype=float)
    target_theta = float(np.median(theta_values))
    return min(rows, key=lambda row: abs(float(row[theta_key]) - target_theta))


def fperp_to_chi(fperp: float, ss: float) -> float:
    return fperp / ss


def project_residues(sample: dict[str, str], theta_deg: float, sector: str) -> tuple[float, float]:
    if sector == "Ds":
        inp = TwoPointInputs(
            mc=float(sample["mc_GeV"]),
            ms=float(sample["ms_GeV"]),
            qq=float(sample["ss_GeV3"]),
            kappa_s=1.0,
            m0_sq=float(sample["m0_sq_GeV2"]),
            mass_low=2.4595,
            mass_high=2.53511,
        )
        mass_low = 2.4595
        mass_high = 2.53511
    else:
        inp = TwoPointInputs(
            mc=float(sample["mb_GeV"]),
            ms=float(sample["ms_GeV"]),
            qq=float(sample["ss_GeV3"]),
            kappa_s=1.0,
            m0_sq=float(sample["m0_sq_GeV2"]),
            mass_low=float(sample["mass_low_GeV"]),
            mass_high=float(sample["mass_high_GeV"]),
        )
        mass_low = float(sample["mass_low_GeV"])
        mass_high = float(sample["mass_high_GeV"])

    m2 = float(sample["M2_GeV2"])
    pi1, _ = projected_twopoint_ope(m2, float(sample["s01_GeV2"]), theta_deg, 0, inp)
    pi2, _ = projected_twopoint_ope(m2, float(sample["s02_GeV2"]), theta_deg, 1, inp)
    if pi1 <= 0.0 or pi2 <= 0.0:
        raise RuntimeError(f"Nonpositive projected two-point OPE at theta={theta_deg:.3f} deg")
    f1 = math.sqrt(pi1 * math.exp(mass_low**2 / m2) / mass_low**2)
    f2 = math.sqrt(pi2 * math.exp(mass_high**2 / m2) / mass_high**2)
    return f1, f2


def ds_transition_inputs(two_point_sample: dict[str, str]) -> dict[str, float]:
    vals = central_inputs_for_scenario("lattice_fperp_s")
    vals["mc"] = float(two_point_sample["mc_GeV"])
    vals["ms"] = float(two_point_sample["ms_GeV"])
    vals["ss"] = float(two_point_sample["ss_GeV3"])
    vals["chi"] = fperp_to_chi(vals["fperp_s_used"], vals["ss"])
    vals["omegaA"] = -2.1
    vals["omegaV"] = 3.8
    return vals


def bs_transition_inputs(two_point_sample: dict[str, str], spec: dict[str, float | str]) -> dict[str, float]:
    rng = np.random.default_rng(20260812)
    vals = sample_bs_inputs(rng, spec, "lattice_fperp_s")
    vals["mc"] = float(two_point_sample["mb_GeV"])
    vals["ms"] = float(two_point_sample["ms_GeV"])
    vals["ss"] = float(two_point_sample["ss_GeV3"])
    vals["chi"] = fperp_to_chi(vals["fperp_s_used"], vals["ss"])
    vals["m_ds1"] = float(spec["m_initial"])
    vals["m_ds"] = 5.36692
    vals["f_ds"] = 0.2303
    vals["f3g"] = -0.0039
    vals["omegaA"] = -2.1
    vals["omegaV"] = 3.8
    return vals


def build_rows() -> list[dict[str, float | str]]:
    ds_sample = representative_sample(read_accepted(OUT / "twopoint_ds1_matrix_mc.csv"))
    bs_sample = representative_sample(read_accepted(OUT / "twopoint_bs1_matrix_mc.csv"))
    theta_values = np.linspace(THETA_MIN, THETA_MAX, N_THETA)
    convolutions = precompute_convolutions()
    rows: list[dict[str, float | str]] = []

    ds_vals = ds_transition_inputs(ds_sample)
    for theta_deg in theta_values:
        f1, f2 = project_residues(ds_sample, float(theta_deg), "Ds")
        for spec in DS_SPECS:
            invariants = transition_invariants(float(spec["M2"]), float(spec["s0"]), ds_vals, convolutions=convolutions)
            physical = physical_couplings(
                invariants,
                theta_deg=float(theta_deg),
                m_state_1=ds_vals["m_ds1"],
                m_state_2=ds_vals["m_ds1_2536"],
                f_1=f1,
                f_2=f2,
                m_p=ds_vals["m_ds"],
                f_p=ds_vals["f_ds"],
                m_q=ds_vals["mc"],
                m_s=ds_vals["ms"],
            )
            is_low = spec["combo"] == "low"
            rows.append(
                {
                    **spec,
                    "theta_deg": float(theta_deg),
                    "theta_reference_deg": THETA_HQET_DEG,
                    "sqrt_s0": math.sqrt(float(spec["s0"])),
                    "f1_GeV": f1,
                    "f2_GeV": f2,
                    "g_GeV_inv": physical["g_1"] if is_low else physical["g_2"],
                    "g_A_component_GeV_inv": physical["g_1_A"] if is_low else physical["g_2_A"],
                    "g_B_component_GeV_inv": physical["g_1_B"] if is_low else physical["g_2_B"],
                    "Gamma_keV": physical["Gamma_1_keV"] if is_low else physical["Gamma_2_keV"],
                    "two_point_source_theta_deg": float(ds_sample["theta_deg"]),
                    "normalization_scheme": "AA_AB_BB_two_point_projection",
                }
            )

    for spec in BS_SPECS:
        vals = bs_transition_inputs(bs_sample, spec)
        mass_low = float(bs_sample["mass_low_GeV"])
        mass_high = float(bs_sample["mass_high_GeV"])
        for theta_deg in theta_values:
            f1, f2 = project_residues(bs_sample, float(theta_deg), "Bs")
            invariants = transition_invariants(float(spec["M2"]), float(spec["s0"]), vals, convolutions=convolutions)
            physical = physical_couplings(
                invariants,
                theta_deg=float(theta_deg),
                m_state_1=mass_low,
                m_state_2=mass_high,
                f_1=f1,
                f_2=f2,
                m_p=vals["m_ds"],
                f_p=vals["f_ds"],
                m_q=vals["mc"],
                m_s=vals["ms"],
            )
            is_low = spec["combo"] == "low"
            rows.append(
                {
                    **spec,
                    "theta_deg": float(theta_deg),
                    "theta_reference_deg": THETA_HQET_DEG,
                    "sqrt_s0": math.sqrt(float(spec["s0"])),
                    "f1_GeV": f1,
                    "f2_GeV": f2,
                    "g_GeV_inv": physical["g_1"] if is_low else physical["g_2"],
                    "g_A_component_GeV_inv": physical["g_1_A"] if is_low else physical["g_2_A"],
                    "g_B_component_GeV_inv": physical["g_1_B"] if is_low else physical["g_2_B"],
                    "Gamma_keV": physical["Gamma_1_keV"] if is_low else physical["Gamma_2_keV"],
                    "two_point_source_theta_deg": float(bs_sample["theta_deg"]),
                    "normalization_scheme": "AA_AB_BB_two_point_projection",
                }
            )
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


def main() -> None:
    rows = build_rows()
    write_csv(HERE / "figure2_mixing_angle_widths.csv", rows)
    for sector in ("Ds", "Bs"):
        print(sector)
        for state in sorted({str(r["state"]) for r in rows if r["sector"] == sector}):
            values = [float(r["Gamma_keV"]) for r in rows if r["state"] == state]
            print(f"  {state}: min={min(values):.4g} max={max(values):.4g}")


if __name__ == "__main__":
    main()
