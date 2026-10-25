"""Check two-point residues in the transition-LCSR auxiliary windows.

This diagnostic intentionally replaces the mass-fitted two-point thresholds
and dedicated two-point Borel windows by the transition-LCSR windows used in
the manuscript.  It reports the resulting residues together with the
two-point pole, dimension-five, and mass diagnostics.  It does not overwrite
the accepted two-point Monte Carlo samples used by the final transition scan.
"""

from __future__ import annotations

import math

import numpy as np

from mixing_angle_inputs import (
    THETA_BS_DEG,
    THETA_BS_SIGMA_DEG,
    THETA_DS_DEG,
    THETA_DS_SIGMA_DEG,
    sample_gaussian_angle,
)
from twopoint_bs1_matrix_sumrule import sample_inputs as sample_bs_inputs
from twopoint_ds1_matrix_sumrule import (
    CumulativePerturbative,
    local_components,
    mc_sample as sample_ds_inputs,
    projected_value,
)


SEED = 20260826
N_SAMPLES = 4000
DELTA_TAU = 2.0e-4

SECTORS = {
    "Ds": {
        "m2": (3.0, 4.5),
        "s0": ((8.5, 9.5), (9.0, 10.0)),
        "theta": (THETA_DS_DEG, THETA_DS_SIGMA_DEG),
        "sample_inputs": sample_ds_inputs,
        "states": ("Ds1(2460)", "Ds1(2536)"),
    },
    "Bs": {
        "m2": (10.0, 14.0),
        "s0": ((39.0, 41.0), (40.0, 42.0)),
        "theta": (THETA_BS_DEG, THETA_BS_SIGMA_DEG),
        "sample_inputs": sample_bs_inputs,
        "states": ("Bs1(5750)", "Bs1(5830)"),
    },
}


def evaluate_sample(
    rng: np.random.Generator,
    sector: str,
) -> list[dict[str, float | str | int]]:
    spec = SECTORS[sector]
    inp = spec["sample_inputs"](rng)
    m2 = float(rng.uniform(*spec["m2"]))
    theta = sample_gaussian_angle(rng, *spec["theta"])
    tau = 1.0 / m2
    m2_plus = 1.0 / (tau + DELTA_TAU)
    m2_minus = 1.0 / (tau - DELTA_TAU)
    evaluator = CumulativePerturbative(
        inp, (m2, m2_plus, m2_minus), n_grid=5000
    )
    masses = (inp.mass_low, inp.mass_high)
    rows: list[dict[str, float | str | int]] = []

    for channel, (state, target_mass, s0_window) in enumerate(
        zip(spec["states"], masses, spec["s0"])
    ):
        s0 = float(rng.uniform(*s0_window))
        locals_now = local_components(m2, inp)
        total = evaluator.matrix(m2, s0) + sum(
            locals_now.values(), np.zeros((2, 2), dtype=float)
        )
        pi_value = projected_value(total, theta, channel)
        perturbative = projected_value(evaluator.matrix(m2, s0), theta, channel)
        perturbative_full = projected_value(evaluator.full(m2), theta, channel)
        pole = perturbative / perturbative_full
        d5 = abs(
            projected_value(locals_now["d5_mixed"], theta, channel)
        ) / abs(pi_value)

        plus_total = evaluator.matrix(m2_plus, s0) + sum(
            local_components(m2_plus, inp).values(),
            np.zeros((2, 2), dtype=float),
        )
        minus_total = evaluator.matrix(m2_minus, s0) + sum(
            local_components(m2_minus, inp).values(),
            np.zeros((2, 2), dtype=float),
        )
        plus = projected_value(plus_total, theta, channel)
        minus = projected_value(minus_total, theta, channel)
        if plus > 0.0 and minus > 0.0:
            mass_sq = -(math.log(plus) - math.log(minus)) / (2.0 * DELTA_TAU)
            mass = math.sqrt(mass_sq) if mass_sq > 0.0 else math.nan
        else:
            mass = math.nan
        residue = (
            math.sqrt(
                pi_value * math.exp(target_mass**2 / m2) / target_mass**2
            )
            if pi_value > 0.0
            else math.nan
        )
        accepted = (
            math.isfinite(residue)
            and math.isfinite(mass)
            and pole >= 0.50
            and d5 <= 0.15
            and abs(mass - target_mass) <= 0.03
        )
        rows.append(
            {
                "sector": sector,
                "state": state,
                "channel": channel,
                "M2_GeV2": m2,
                "s0_GeV2": s0,
                "theta_deg": theta,
                "target_mass_GeV": target_mass,
                "effective_mass_GeV": mass,
                "mass_shift_MeV": 1000.0 * (mass - target_mass),
                "pole_fraction": pole,
                "d5_fraction": d5,
                "f_GeV": residue,
                "accepted": int(accepted),
            }
        )
    return rows


def quantiles(rows: list[dict[str, float | str | int]], key: str) -> tuple[float, float, float]:
    values = np.asarray([float(row[key]) for row in rows], dtype=float)
    values = values[np.isfinite(values)]
    return tuple(float(value) for value in np.percentile(values, (16.0, 50.0, 84.0)))


def format_quantiles(values: tuple[float, float, float], scale: float = 1.0) -> str:
    low, median, high = (scale * value for value in values)
    return f"{median:.6g} [{low:.6g}, {high:.6g}]"


def report(rows: list[dict[str, float | str | int]]) -> str:
    lines = [
        "Two-point residues evaluated in the transition-LCSR windows",
        "===========================================================",
        f"Samples per sector: {N_SAMPLES}; seed: {SEED}",
        "Strict acceptance: pole >= 0.50, d5 <= 0.15, |m_eff-m_pole| <= 30 MeV.",
        "",
    ]
    for sector in ("Ds", "Bs"):
        sector_rows = [row for row in rows if row["sector"] == sector]
        by_sample: dict[int, list[dict[str, float | str | int]]] = {}
        for index, row in enumerate(sector_rows):
            by_sample.setdefault(index // 2, []).append(row)
        joint = sum(
            int(len(sample_rows) == 2 and all(int(row["accepted"]) for row in sample_rows))
            for sample_rows in by_sample.values()
        )
        lines.append(
            f"{sector}: joint two-state acceptance {joint}/{len(by_sample)} "
            f"({100.0 * joint / len(by_sample):.2f}%)"
        )
        for state in SECTORS[sector]["states"]:
            state_rows = [row for row in sector_rows if row["state"] == state]
            accepted = [row for row in state_rows if int(row["accepted"])]
            lines.extend(
                [
                    f"  {state}: accepted {len(accepted)}/{len(state_rows)} "
                    f"({100.0 * len(accepted) / len(state_rows):.2f}%)",
                    f"    theta [deg]: {format_quantiles(quantiles(state_rows, 'theta_deg'))}",
                    f"    f all [MeV]: {format_quantiles(quantiles(state_rows, 'f_GeV'), 1000.0)}",
                    f"    m_eff-m_pole all [MeV]: {format_quantiles(quantiles(state_rows, 'mass_shift_MeV'))}",
                    f"    pole all: {format_quantiles(quantiles(state_rows, 'pole_fraction'))}",
                    f"    d5 all: {format_quantiles(quantiles(state_rows, 'd5_fraction'))}",
                ]
            )
            if accepted:
                lines.append(
                    f"    f accepted [MeV]: {format_quantiles(quantiles(accepted, 'f_GeV'), 1000.0)}"
                )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> None:
    rng = np.random.default_rng(SEED)
    rows: list[dict[str, float | str | int]] = []
    for sector in ("Ds", "Bs"):
        for _ in range(N_SAMPLES):
            rows.extend(evaluate_sample(rng, sector))
    print(report(rows))


if __name__ == "__main__":
    main()
