"""Direct normalized-current AA/AB/BB two-point QCD sum rule for Bs1.

The analytic OPE is the heavy-flavour continuation of the Ds1 calculation in
``twopoint_ds1_matrix_sumrule.py``:

    J_A^mu = sbar gamma^mu gamma5 b,
    J_B^mu = i/(mb+ms) sbar sigma^(mu nu) p_nu gamma5 b.

It contains the exact-mass leading-order perturbative matrix, the local
strange-quark condensate through O(ms^2), and the local mixed condensate.  The
local terms are appropriate here because this is an ordinary two-point SVZ
sum rule; they remain excluded from the external-photon transition LCSR.

The physical residues are projected at the external angle
theta_Bs = 38.5 +/- 0.1 degrees.  The angle obtained by diagonalising the
truncated matrix is retained only as a diagnostic.  Separate thresholds are
fitted to the model lower pole and the observed Bs1(5830) mass.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from mixing_angle_inputs import (
    ANGLE_SOURCE,
    THETA_BS_DEG,
    THETA_BS_SIGMA_DEG,
    sample_gaussian_angle,
)
from twopoint_ds1_matrix_sumrule import (
    CumulativePerturbative,
    Inputs,
    local_components,
    ope_matrix,
    projected_effective_mass,
    projected_ope,
    projected_pole_fraction,
    projected_value,
    rotation_from_matrix,
    rotation_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

SEED = 20260726
N_SAMPLES = 2000

M2_MIN = 6.50
M2_MAX = 7.30
S0_FIT_MIN = 40.0
S0_FIT_MAX = 50.0
S0_MIX_MIN = 42.0
S0_MIX_MAX = 48.0


def central_inputs() -> Inputs:
    return Inputs(
        mc=4.18,
        ms=0.093,
        qq=-(0.240**3),
        kappa_s=0.8,
        m0_sq=0.8,
        mass_low=5.750,
        mass_high=5.82870,
    )


def fitted_threshold(
    m2: float,
    theta_deg: float,
    channel: int,
    target_mass: float,
    inp: Inputs,
    candidates: np.ndarray | None = None,
) -> tuple[float, float]:
    if candidates is None:
        candidates = np.linspace(S0_FIT_MIN, S0_FIT_MAX, 401)
    best_s0 = math.nan
    best_mass = math.nan
    best_distance = math.inf
    for s0 in candidates:
        mass = projected_effective_mass(
            m2, float(s0), theta_deg, channel, inp
        )
        if not math.isfinite(mass):
            continue
        distance = abs(mass - target_mass)
        if distance < best_distance:
            best_distance = distance
            best_s0 = float(s0)
            best_mass = mass
    return best_s0, best_mass


def physical_point(m2: float, inp: Inputs) -> dict[str, float]:
    theta = THETA_BS_DEG
    s01, meff1 = fitted_threshold(m2, theta, 0, inp.mass_low, inp)
    s02, meff2 = fitted_threshold(m2, theta, 1, inp.mass_high, inp)
    pi1, pieces1 = projected_ope(m2, s01, theta, 0, inp)
    pi2, pieces2 = projected_ope(m2, s02, theta, 1, inp)
    f1 = math.sqrt(pi1 * math.exp(inp.mass_low**2 / m2) / inp.mass_low**2)
    f2 = math.sqrt(pi2 * math.exp(inp.mass_high**2 / m2) / inp.mass_high**2)
    pole1 = projected_pole_fraction(m2, s01, theta, 0, inp)
    pole2 = projected_pole_fraction(m2, s02, theta, 1, inp)
    d51 = abs(pieces1["d5_mixed"]) / abs(pi1)
    d52 = abs(pieces2["d5_mixed"]) / abs(pi2)

    angle_samples: list[float] = []
    for s0_mix in np.arange(S0_MIX_MIN, S0_MIX_MAX + 0.001, 0.25):
        matrix, _ = ope_matrix(m2, float(s0_mix), inp)
        matrix_angle, _ = rotation_from_matrix(matrix)
        angle_samples.append(matrix_angle)
    common_matrix, _ = ope_matrix(m2, 45.0, inp)
    rotation = rotation_matrix(theta)
    rotated = rotation @ common_matrix @ rotation.T
    offdiag = abs(float(rotated[0, 1])) / math.sqrt(
        abs(float(rotated[0, 0] * rotated[1, 1]))
    )

    return {
        "M2_GeV2": m2,
        "theta_deg": theta,
        "theta_matrix_median_deg": float(np.median(angle_samples)),
        "theta_matrix_min_deg": float(min(angle_samples)),
        "theta_matrix_max_deg": float(max(angle_samples)),
        "Pi12_normalized_at_s0mix45": offdiag,
        "s01_GeV2": s01,
        "s02_GeV2": s02,
        "sqrt_s01_minus_m1_GeV": math.sqrt(s01) - inp.mass_low,
        "sqrt_s02_minus_m2_GeV": math.sqrt(s02) - inp.mass_high,
        "m_eff1_GeV": meff1,
        "m_eff2_GeV": meff2,
        "Pi1_GeV4": pi1,
        "Pi2_GeV4": pi2,
        "f1_GeV": f1,
        "f2_GeV": f2,
        "pole_fraction1": pole1,
        "pole_fraction2": pole2,
        "d5_fraction1": d51,
        "d5_fraction2": d52,
    }


def central_scan(inp: Inputs) -> list[dict[str, float]]:
    return [
        physical_point(float(m2), inp)
        for m2 in np.arange(M2_MIN, M2_MAX + 0.001, 0.05)
    ]


def accepted_point(row: dict[str, float]) -> bool:
    return (
        row["pole_fraction1"] >= 0.50
        and row["pole_fraction2"] >= 0.50
        and row["d5_fraction1"] <= 0.15
        and row["d5_fraction2"] <= 0.15
        and abs(row["m_eff1_GeV"] - row["mass_low_GeV"]) <= 0.03
        and abs(row["m_eff2_GeV"] - row["mass_high_GeV"]) <= 0.03
        and 0.45 <= row["sqrt_s01_minus_m1_GeV"] <= 1.05
        and 0.45 <= row["sqrt_s02_minus_m2_GeV"] <= 1.05
        and math.isfinite(row["f1_GeV"])
        and math.isfinite(row["f2_GeV"])
    )


def clipped_normal(
    rng: np.random.Generator,
    mean: float,
    sigma: float,
    lo: float,
    hi: float,
) -> float:
    return float(np.clip(rng.normal(mean, sigma), lo, hi))


def sample_inputs(rng: np.random.Generator) -> Inputs:
    condensate_scale = clipped_normal(rng, 0.240, 0.010, 0.205, 0.275)
    return Inputs(
        mc=clipped_normal(rng, 4.18, 0.03, 4.05, 4.30),
        ms=clipped_normal(rng, 0.093, 0.011, 0.050, 0.140),
        qq=-(condensate_scale**3),
        kappa_s=clipped_normal(rng, 0.80, 0.10, 0.45, 1.15),
        m0_sq=clipped_normal(rng, 0.80, 0.20, 0.30, 1.30),
        mass_low=clipped_normal(rng, 5.750, 0.026, 5.670, 5.830),
        mass_high=clipped_normal(rng, 5.82870, 0.00020, 5.8280, 5.8294),
    )


def monte_carlo_scan(n_samples: int = N_SAMPLES) -> list[dict[str, float]]:
    rng = np.random.default_rng(SEED)
    rows: list[dict[str, float]] = []
    delta_tau = 2.0e-4
    candidates = np.linspace(S0_FIT_MIN, S0_FIT_MAX, 401)

    for _ in range(n_samples):
        inp = sample_inputs(rng)
        m2 = float(rng.uniform(M2_MIN, M2_MAX))
        theta = sample_gaussian_angle(rng, THETA_BS_DEG, THETA_BS_SIGMA_DEG)
        s0_mix = float(rng.uniform(S0_MIX_MIN, S0_MIX_MAX))
        tau = 1.0 / m2
        m2_plus = 1.0 / (tau + delta_tau)
        m2_minus = 1.0 / (tau - delta_tau)
        evaluator = CumulativePerturbative(
            inp, (m2, m2_plus, m2_minus), n_grid=5000
        )

        mix_matrix = evaluator.matrix(m2, s0_mix) + sum(
            local_components(m2, inp).values(),
            np.zeros((2, 2), dtype=float),
        )
        theta_matrix, _ = rotation_from_matrix(mix_matrix)
        rotation = rotation_matrix(theta)
        rotated_mix = rotation @ mix_matrix @ rotation.T
        offdiag = abs(float(rotated_mix[0, 1])) / math.sqrt(
            abs(float(rotated_mix[0, 0] * rotated_mix[1, 1]))
        )

        fitted: list[tuple[float, float, float, float, float, float]] = []
        failed = False
        for channel, target_mass in enumerate((inp.mass_low, inp.mass_high)):
            best: tuple[float, float, float] | None = None
            for s0 in candidates:
                plus_matrix = evaluator.matrix(m2_plus, float(s0)) + sum(
                    local_components(m2_plus, inp).values(),
                    np.zeros((2, 2), dtype=float),
                )
                minus_matrix = evaluator.matrix(m2_minus, float(s0)) + sum(
                    local_components(m2_minus, inp).values(),
                    np.zeros((2, 2), dtype=float),
                )
                plus = projected_value(plus_matrix, theta, channel)
                minus = projected_value(minus_matrix, theta, channel)
                if plus <= 0.0 or minus <= 0.0:
                    continue
                mass_sq = -(math.log(plus) - math.log(minus)) / (
                    2.0 * delta_tau
                )
                if mass_sq <= 0.0:
                    continue
                mass = math.sqrt(mass_sq)
                distance = abs(mass - target_mass)
                if best is None or distance < best[0]:
                    best = (distance, float(s0), mass)
            if best is None:
                failed = True
                break

            _, s0, meff = best
            perturbative = evaluator.matrix(m2, s0)
            locals_now = local_components(m2, inp)
            total = perturbative + sum(
                locals_now.values(), np.zeros((2, 2), dtype=float)
            )
            pi_value = projected_value(total, theta, channel)
            pole = projected_value(perturbative, theta, channel) / projected_value(
                evaluator.full(m2), theta, channel
            )
            d5 = abs(
                projected_value(locals_now["d5_mixed"], theta, channel)
            ) / abs(pi_value)
            residue = (
                math.sqrt(
                    pi_value
                    * math.exp(target_mass**2 / m2)
                    / target_mass**2
                )
                if pi_value > 0.0
                else math.nan
            )
            fitted.append((s0, meff, pi_value, pole, d5, residue))

        if failed or len(fitted) != 2:
            continue
        one, two = fitted
        row = {
            "M2_GeV2": m2,
            "s0_mix_GeV2": s0_mix,
            "mb_GeV": inp.mc,
            "ms_GeV": inp.ms,
            "ss_GeV3": inp.ss,
            "m0_sq_GeV2": inp.m0_sq,
            "mass_low_GeV": inp.mass_low,
            "mass_high_GeV": inp.mass_high,
            "theta_deg": theta,
            "theta_matrix_deg": theta_matrix,
            "Pi12_normalized": offdiag,
            "s01_GeV2": one[0],
            "s02_GeV2": two[0],
            "sqrt_s01_minus_m1_GeV": math.sqrt(one[0]) - inp.mass_low,
            "sqrt_s02_minus_m2_GeV": math.sqrt(two[0]) - inp.mass_high,
            "m_eff1_GeV": one[1],
            "m_eff2_GeV": two[1],
            "f1_GeV": one[5],
            "f2_GeV": two[5],
            "pole_fraction1": one[3],
            "pole_fraction2": two[3],
            "d5_fraction1": one[4],
            "d5_fraction2": two[4],
        }
        row["accepted"] = int(accepted_point(row))
        rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        raise RuntimeError(f"no rows to write to {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def quantile_line(
    rows: list[dict[str, float]],
    key: str,
    label: str,
    unit: str = "",
) -> str:
    values = np.asarray([row[key] for row in rows], dtype=float)
    q16, q50, q84 = np.percentile(values, [16.0, 50.0, 84.0])
    return f"{label}: {q50:.6g} [{q16:.6g}, {q84:.6g}] {unit}".rstrip()


def central_summary(rows: list[dict[str, float]]) -> str:
    central = min(rows, key=lambda row: abs(row["M2_GeV2"] - 7.0))
    lines = [
        "Direct Bs1 AA/AB/BB two-point QCD sum rule",
        "==========================================",
        "OPE: exact-mass LO perturbative + local d=3 through ms^2 + local d=5.",
        "No external fA, fT, f1, f2, or overlap parameter is used.",
        f"theta_Bs={THETA_BS_DEG:.1f}+-{THETA_BS_SIGMA_DEG:.1f} deg is external.",
        f"Angle source: {ANGLE_SOURCE}.",
        f"Selected two-point window: M^2=[{M2_MIN:.2f},{M2_MAX:.2f}] GeV^2.",
        "",
        "Central point M^2=7.00 GeV^2:",
    ]
    for key, label, unit in (
        ("s01_GeV2", "s01", "GeV^2"),
        ("s02_GeV2", "s02", "GeV^2"),
        ("m_eff1_GeV", "m_eff1", "GeV"),
        ("m_eff2_GeV", "m_eff2", "GeV"),
        ("f1_GeV", "f1", "GeV"),
        ("f2_GeV", "f2", "GeV"),
        ("pole_fraction1", "pole fraction 1", ""),
        ("pole_fraction2", "pole fraction 2", ""),
        ("d5_fraction1", "d5 fraction 1", ""),
        ("d5_fraction2", "d5 fraction 2", ""),
        ("theta_matrix_median_deg", "matrix-angle diagnostic", "deg"),
        ("Pi12_normalized_at_s0mix45", "normalized off-diagonal residual", ""),
    ):
        lines.append(f"{label}: {central[key]:.9g} {unit}".rstrip())
    return "\n".join(lines) + "\n"


def mc_summary(rows: list[dict[str, float]]) -> str:
    accepted = [row for row in rows if int(row["accepted"]) == 1]
    lines = [
        "Direct Bs1 AA/AB/BB two-point Monte Carlo",
        "==========================================",
        f"Seed: {SEED}; generated valid fits: {len(rows)}; accepted: {len(accepted)}.",
        "Inputs varied: mb, ms, condensates, both pole masses, M^2, and theta_Bs.",
        "Acceptance: pole fractions >=0.50, d5 fractions <=0.15, mass fits within 30 MeV,",
        "and continuum gaps 0.45--1.05 GeV.",
        "",
    ]
    if not accepted:
        lines.append("No samples passed all diagnostics.")
        return "\n".join(lines) + "\n"
    for key, label, unit in (
        ("theta_deg", "theta input", "deg"),
        ("theta_matrix_deg", "matrix-angle diagnostic", "deg"),
        ("Pi12_normalized", "normalized off-diagonal residual", ""),
        ("M2_GeV2", "M^2", "GeV^2"),
        ("s01_GeV2", "s01", "GeV^2"),
        ("s02_GeV2", "s02", "GeV^2"),
        ("f1_GeV", "f1", "GeV"),
        ("f2_GeV", "f2", "GeV"),
        ("pole_fraction1", "pole fraction 1", ""),
        ("pole_fraction2", "pole fraction 2", ""),
    ):
        lines.append(quantile_line(accepted, key, label, unit))
    return "\n".join(lines) + "\n"


def main() -> None:
    central_rows = central_scan(central_inputs())
    central_csv = OUT / "twopoint_bs1_physical_residue_grid.csv"
    central_txt = OUT / "twopoint_bs1_physical_residue_summary.txt"
    write_csv(central_csv, central_rows)
    summary = central_summary(central_rows)
    central_txt.write_text(summary)
    print(summary)

    mc_rows = monte_carlo_scan()
    mc_csv = OUT / "twopoint_bs1_matrix_mc.csv"
    mc_txt = OUT / "twopoint_bs1_matrix_mc_summary.txt"
    write_csv(mc_csv, mc_rows)
    summary_mc = mc_summary(mc_rows)
    mc_txt.write_text(summary_mc)
    print(summary_mc)
    print(f"Wrote {central_csv}")
    print(f"Wrote {central_txt}")
    print(f"Wrote {mc_csv}")
    print(f"Wrote {mc_txt}")


if __name__ == "__main__":
    main()
