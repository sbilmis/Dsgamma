"""Project the mixed two-pole normalization onto Wang's pure-AA channel.

Wang's current is the pure local axial current, so the direct comparison is
with the AA correlator, not with either mixed-current residue f1 or f2.  The
script also reports the axial-current pole overlaps reconstructed from the
mixed residues at the external Ds mixing angle.

The current arXiv v4 value is f_Ds1 = 0.345 +- 0.017 GeV.  Earlier versions
printed 0.245 GeV; arXiv v4 explicitly records that number as an error.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from twopoint_ds1_matrix_sumrule import Inputs, mc_sample, ope_matrix


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

SEED = 20260723
N_SAMPLES = 2000

WANG_F_A_GEV = 0.345
WANG_F_A_ERR_GEV = 0.017
WANG_M2_RANGE = (2.7, 3.3)
WANG_S0_RANGE = (8.8, 9.8)  # 9.3 +- 0.5 GeV^2


def stats(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(array)),
        "std": float(np.std(array, ddof=1)),
        "p16": float(np.percentile(array, 16.0)),
        "median": float(np.percentile(array, 50.0)),
        "p84": float(np.percentile(array, 84.0)),
    }


def direct_aa_samples() -> list[float]:
    rng = np.random.default_rng(SEED)
    values: list[float] = []
    for _ in range(N_SAMPLES):
        inp = mc_sample(rng)
        m2 = float(rng.uniform(*WANG_M2_RANGE))
        s0 = float(rng.uniform(*WANG_S0_RANGE))
        pi_aa = float(ope_matrix(m2, s0, inp)[0][0, 0])
        if pi_aa <= 0.0:
            continue
        values.append(
            math.sqrt(
                pi_aa * math.exp(inp.mass_low**2 / m2) / inp.mass_low**2
            )
        )
    return values


def mixed_projection_samples() -> dict[str, list[float]]:
    with (OUT / "twopoint_ds1_matrix_mc.csv").open() as handle:
        rows = [row for row in csv.DictReader(handle) if int(row["accepted"]) == 1]

    result = {
        "f1": [],
        "f2": [],
        "fA_low_pole": [],
        "fA_high_pole": [],
        "fA_two_pole_effective": [],
        "Pi12_normalized": [],
    }
    m1 = Inputs().mass_low
    m2_state = Inputs().mass_high
    for row in rows:
        theta = math.radians(float(row["theta_deg"]))
        f1 = float(row["f1_GeV"])
        f2 = float(row["f2_GeV"])
        borel = float(row["M2_GeV2"])
        f_a_low = math.sin(theta) * f1
        f_a_high = math.cos(theta) * f2
        # Convert the two physical poles to the equivalent low-mass single-pole
        # residue at the sample's Borel mass.  This is only a diagnostic because
        # Wang's analysis uses a one-pole ansatz.
        f_a_effective = math.sqrt(
            (f_a_low * m1) ** 2
            + (f_a_high * m2_state) ** 2
            * math.exp(-(m2_state**2 - m1**2) / borel)
        ) / m1
        result["f1"].append(f1)
        result["f2"].append(f2)
        result["fA_low_pole"].append(f_a_low)
        result["fA_high_pole"].append(f_a_high)
        result["fA_two_pole_effective"].append(f_a_effective)
        result["Pi12_normalized"].append(float(row["Pi12_normalized"]))
    return result


def main() -> None:
    direct = stats(direct_aa_samples())
    projected = {key: stats(values) for key, values in mixed_projection_samples().items()}
    combined_sigma = math.sqrt(direct["std"] ** 2 + WANG_F_A_ERR_GEV**2)
    pull = (direct["mean"] - WANG_F_A_GEV) / combined_sigma

    rows = [
        {
            "quantity": "Wang_v4_pure_JA",
            "current": "J_A",
            "median_GeV": WANG_F_A_GEV,
            "p16_GeV": WANG_F_A_GEV - WANG_F_A_ERR_GEV,
            "p84_GeV": WANG_F_A_GEV + WANG_F_A_ERR_GEV,
            "comparison_scope": "external pure-current benchmark",
        },
        {
            "quantity": "our_direct_AA_Wang_window",
            "current": "J_A",
            "median_GeV": direct["median"],
            "p16_GeV": direct["p16"],
            "p84_GeV": direct["p84"],
        "comparison_scope": "direct truncated-AA OPE benchmark",
        },
    ]
    for key, current, scope in (
        ("f1", "J_1", "physical mixed-current residue; not directly Wang-comparable"),
        ("f2", "J_2", "physical mixed-current residue; not directly Wang-comparable"),
        ("fA_low_pole", "J_A -> Ds1(2460)", "mixed-basis projection diagnostic"),
        ("fA_high_pole", "J_A -> Ds1(2536)", "mixed-basis projection diagnostic"),
        ("fA_two_pole_effective", "J_A, two poles", "one-pole-equivalent diagnostic"),
    ):
        value = projected[key]
        rows.append(
            {
                "quantity": key,
                "current": current,
                "median_GeV": value["median"],
                "p16_GeV": value["p16"],
                "p84_GeV": value["p84"],
                "comparison_scope": scope,
            }
        )

    csv_path = OUT / "wang_pure_axial_residue_comparison_table.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (OUT / "stage3_residue_benchmark_checks.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "Wang pure-axial decay-constant projection benchmark",
        "====================================================",
        "Wang arXiv:1506.01993v4 uses J_A=Qbar gamma_mu gamma5 q and",
        "defines <0|J_A|A>=f_A m_A epsilon_mu.",
        "The corrected Table-1 result is f_Ds1=0.345+-0.017 GeV;",
        "0.245+-0.017 GeV is the explicitly corrected obsolete value.",
        "",
        "Direct truncated-AA OPE comparison",
        "  our AA result in Wang's M^2=[2.7,3.3] GeV^2 and",
        "  s0=[8.8,9.8] GeV^2 window:",
        "  f_A={:.6f} [{:.6f},{:.6f}] GeV".format(
            direct["median"], direct["p16"], direct["p84"]
        ),
        f"  difference using sample std and Wang error: {pull:+.2f} sigma",
        "  conclusion: the present LO+d3+d5 AA truncation does not reproduce",
        "  Wang's absolute normalization within the quoted uncertainties;",
        "  this does not validate or invalidate the mixed two-pole residues.",
        "",
        "Mixed-current projection at theta_Ds=26.6+-0.6 deg",
        "  f1={:.6f} [{:.6f},{:.6f}] GeV".format(
            projected["f1"]["median"], projected["f1"]["p16"], projected["f1"]["p84"]
        ),
        "  f2={:.6f} [{:.6f},{:.6f}] GeV".format(
            projected["f2"]["median"], projected["f2"]["p16"], projected["f2"]["p84"]
        ),
        "  sin(theta) f1={:.6f} [{:.6f},{:.6f}] GeV".format(
            projected["fA_low_pole"]["median"], projected["fA_low_pole"]["p16"], projected["fA_low_pole"]["p84"]
        ),
        "  cos(theta) f2={:.6f} [{:.6f},{:.6f}] GeV".format(
            projected["fA_high_pole"]["median"], projected["fA_high_pole"]["p16"], projected["fA_high_pole"]["p84"]
        ),
        "  two-pole one-pole-equivalent f_A={:.6f} [{:.6f},{:.6f}] GeV".format(
            projected["fA_two_pole_effective"]["median"], projected["fA_two_pole_effective"]["p16"], projected["fA_two_pole_effective"]["p84"]
        ),
        "  normalized Pi12 residual={:.4f} [{:.4f},{:.4f}]".format(
            projected["Pi12_normalized"]["median"], projected["Pi12_normalized"]["p16"], projected["Pi12_normalized"]["p84"]
        ),
        "",
        "Interpretation",
        "  Wang benchmarks only the pure-J_A correlator scale.  Neither f1 nor",
        "  f2 should be set equal to Wang's result at the physical angle.",
        "  The direct-AA and projected two-pole quantities use different pole",
        "  models and must not be advertised as an exact reproduction test.",
    ]
    text = "\n".join(lines) + "\n"
    (OUT / "wang_pure_axial_validation.txt").write_text(text)
    # Retain the old output name as an explicit redirect for reproducibility.
    (OUT / "stage3_residue_benchmark_checks.txt").write_text(
        "This legacy overlap-fit benchmark has been superseded.\n\n" + text
    )
    print(text)
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
