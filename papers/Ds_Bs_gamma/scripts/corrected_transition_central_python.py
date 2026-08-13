"""Central numerical reference for the corrected Rohrwild transition LCSR.

The output is intentionally a flat key/value table so that the independent
Mathematica implementation can be compared term by term without any
format-conversion ambiguity.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from bs1_stage2_baseline import FPERP_S_1GEV
from final_stage2_uncertainty_scan import set_photon_shape_params
from rohrwild_transition_exact import (
    physical_couplings,
    precompute_convolutions,
    transition_invariants,
)
from stage1_axial_g1_baseline import central_inputs
from twopoint_ds1_matrix_sumrule import (
    Inputs as TwoPointInputs,
    projected_ope as projected_twopoint_ope,
)


DS_REFERENCE = {
    "M2": 3.75,
    "s0": 8.0,
    "theta_deg": 26.6,
    "m1": 2.4595,
    "m2": 2.53511,
    "f1": 0.40058930493635564,
    "f2": 0.16667036765413243,
}

BS_REFERENCES = (
    {
        "label": "Bs_low",
        "M2": 12.0,
        "s0": 38.5,
        "theta_deg": 38.5,
        "m1": 5.750,
        "m2": 5.750,
        "quoted": "g_1",
    },
    {
        "label": "Bs_high",
        "M2": 12.0,
        "s0": 39.5,
        "theta_deg": 38.5,
        "m1": 5.82870,
        "m2": 5.82870,
        "quoted": "g_2",
    },
)

BS_TWOPOINT_REFERENCE = {
    "M2": 7.0,
    "s01": 44.35,
    "s02": 43.025,
    "theta_deg": 38.5,
    "m1": 5.750,
    "m2": 5.82870,
}


def ds_inputs() -> dict[str, float]:
    values = central_inputs()
    values.update(
        {
            "m_ds1_2536": DS_REFERENCE["m2"],
            "fT": 0.256,
            "omegaA": -2.1,
            "omegaV": 3.8,
        }
    )
    values["fperp_s_used"] = FPERP_S_1GEV
    values["chi"] = FPERP_S_1GEV / values["ss"]
    return values


def bs_inputs(m_initial: float) -> dict[str, float]:
    ss = 0.8 * (-(0.240) ** 3)
    return {
        "mc": 4.18,
        "ms": 0.093,
        "m_ds1": m_initial,
        "m_ds": 5.36692,
        "f_ds1": 0.305,
        "f_ds": 0.2303,
        # An accepted reference point for the provisional bottom-sector
        # diagonal-overlap closure.  The nominal fT=0.285 central value gives
        # |rho_AB|>1 at theta=38.5 degrees and is therefore not a valid
        # physical-normalization point.
        "fT": 0.240,
        "ss": ss,
        "chi": FPERP_S_1GEV / ss,
        "f3g": -0.0039,
        "ec": -1.0 / 3.0,
        "es": -1.0 / 3.0,
        "omegaA": -2.1,
        "omegaV": 3.8,
        "fperp_s_used": FPERP_S_1GEV,
    }


def flatten(prefix: str, values: dict[str, object], rows: list[dict[str, object]]) -> None:
    for key, value in values.items():
        if isinstance(value, (int, float)):
            rows.append({"key": f"{prefix}.{key}", "value": float(value)})


def bs_direct_decay_constants() -> tuple[float, float]:
    inp = TwoPointInputs(
        mc=4.18,
        ms=0.093,
        qq=-(0.240**3),
        kappa_s=0.8,
        m0_sq=0.8,
        mass_low=float(BS_TWOPOINT_REFERENCE["m1"]),
        mass_high=float(BS_TWOPOINT_REFERENCE["m2"]),
    )
    pi1, _ = projected_twopoint_ope(
        float(BS_TWOPOINT_REFERENCE["M2"]),
        float(BS_TWOPOINT_REFERENCE["s01"]),
        float(BS_TWOPOINT_REFERENCE["theta_deg"]),
        0,
        inp,
    )
    pi2, _ = projected_twopoint_ope(
        float(BS_TWOPOINT_REFERENCE["M2"]),
        float(BS_TWOPOINT_REFERENCE["s02"]),
        float(BS_TWOPOINT_REFERENCE["theta_deg"]),
        1,
        inp,
    )
    f1 = (
        pi1
        * math.exp(
            float(BS_TWOPOINT_REFERENCE["m1"]) ** 2
            / float(BS_TWOPOINT_REFERENCE["M2"])
        )
        / float(BS_TWOPOINT_REFERENCE["m1"]) ** 2
    ) ** 0.5
    f2 = (
        pi2
        * math.exp(
            float(BS_TWOPOINT_REFERENCE["m2"]) ** 2
            / float(BS_TWOPOINT_REFERENCE["M2"])
        )
        / float(BS_TWOPOINT_REFERENCE["m2"]) ** 2
    ) ** 0.5
    return f1, f2


def main() -> None:
    set_photon_shape_params(-2.1, 3.8)
    convolutions = precompute_convolutions()
    rows: list[dict[str, object]] = []
    flatten("convolution", convolutions, rows)

    ds = ds_inputs()
    ds_inv = transition_invariants(
        DS_REFERENCE["M2"],
        DS_REFERENCE["s0"],
        ds,
        convolutions=convolutions,
    )
    ds_phys = physical_couplings(
        ds_inv,
        theta_deg=DS_REFERENCE["theta_deg"],
        m_state_1=DS_REFERENCE["m1"],
        m_state_2=DS_REFERENCE["m2"],
        f_1=DS_REFERENCE["f1"],
        f_2=DS_REFERENCE["f2"],
        m_p=ds["m_ds"],
        f_p=ds["f_ds"],
        m_q=ds["mc"],
        m_s=ds["ms"],
    )
    flatten("Ds", ds_inv, rows)
    flatten("Ds", ds_phys, rows)

    summaries = [
        "Corrected Rohrwild central numerical reference (Python)",
        "=======================================================",
        "Transition local condensate: excluded.",
        "Tensor/EM numerator prescription: exact off-shell double Borel.",
        (
            "Ds: g1={:+.10f} GeV^-1, Gamma1={:.10f} keV; "
            "g2={:+.10f} GeV^-1, Gamma2={:.10f} keV"
        ).format(
            ds_phys["g_1"],
            ds_phys["Gamma_1_keV"],
            ds_phys["g_2"],
            ds_phys["Gamma_2_keV"],
        ),
    ]

    bs_f1, bs_f2 = bs_direct_decay_constants()
    for reference in BS_REFERENCES:
        values = bs_inputs(float(reference["m1"]))
        inv = transition_invariants(
            float(reference["M2"]),
            float(reference["s0"]),
            values,
            convolutions=convolutions,
        )
        phys = physical_couplings(
            inv,
            theta_deg=float(reference["theta_deg"]),
            m_state_1=float(BS_TWOPOINT_REFERENCE["m1"]),
            m_state_2=float(BS_TWOPOINT_REFERENCE["m2"]),
            f_1=bs_f1,
            f_2=bs_f2,
            m_p=values["m_ds"],
            f_p=values["f_ds"],
            m_q=values["mc"],
            m_s=values["ms"],
        )
        label = str(reference["label"])
        flatten(label, inv, rows)
        flatten(
            label,
            {
                **phys,
                "f1": bs_f1,
                "f2": bs_f2,
                "two_point_M2": BS_TWOPOINT_REFERENCE["M2"],
                "two_point_s01": BS_TWOPOINT_REFERENCE["s01"],
                "two_point_s02": BS_TWOPOINT_REFERENCE["s02"],
            },
            rows,
        )
        quoted = str(reference["quoted"])
        gamma_key = "Gamma_1_keV" if quoted == "g_1" else "Gamma_2_keV"
        summaries.append(
            "{}: {}={:+.10f} GeV^-1, {}={:.10f} keV; "
            "f1={:.10f}, f2={:.10f}".format(
                label,
                quoted,
                phys[quoted],
                gamma_key,
                phys[gamma_key],
                bs_f1,
                bs_f2,
            )
        )

    csv_path = OUT / "corrected_transition_central_python.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("key", "value"))
        writer.writeheader()
        writer.writerows(rows)

    summary_path = OUT / "corrected_transition_central_python.txt"
    summaries.extend(
        [
            "",
            "DA convolutions:",
            *[
                f"  {key} = {value:+.12e}"
                for key, value in convolutions.items()
                if key != "z"
            ],
            "",
            f"Wrote {csv_path}",
        ]
    )
    summary_path.write_text("\n".join(summaries) + "\n")
    print("\n".join(summaries))


if __name__ == "__main__":
    main()
