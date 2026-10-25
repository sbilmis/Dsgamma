"""Normalized AA/AB/BB two-point QCD sum rule for the Ds1 system.

This script replaces the former ``external f1 + overlap-derived f2`` closure.
It uses the basis currents

    J_A^mu = sbar gamma^mu gamma5 c,
    J_B^mu = i/(mc+ms) sbar sigma^(mu nu) p_nu gamma5 c,

and evaluates the full 2 x 2 transverse correlator matrix at a common,
explicit OPE truncation:

  * exact-mass leading-order perturbative spectral densities;
  * the local strange-quark condensate through O(ms^2);
  * the local mixed condensate <sbar g sigma.G s> (dimension five).

The AA entry reproduces the Colangelo two-point sum rule analytically.  The
AB and BB local terms are generated from the same coordinate-space propagator
and were independently traced in ``twopoint_local_condensate_audit.wl``.

The ordinary local condensates in this file belong to a two-point SVZ sum
rule.  Their presence does not conflict with their exclusion from the
external-photon Rohrwild transition LCSR.

The calculation is deliberately labelled LO+d3+d5: gluon-condensate and NLO
corrections are not silently modelled.  This gives a complete correlation
matrix at a controlled common truncation, rather than mixing unrelated
single-current decay constants.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mixing_angle_inputs import (
    ANGLE_SOURCE,
    THETA_DS_DEG,
    THETA_DS_SIGMA_DEG,
    sample_gaussian_angle,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


@dataclass(frozen=True)
class Inputs:
    mc: float = 1.27
    ms: float = 0.093
    qq: float = -(0.24**3)
    kappa_s: float = 0.8
    m0_sq: float = 0.8
    mass_low: float = 2.4595
    mass_high: float = 2.53511

    @property
    def ss(self) -> float:
        return self.kappa_s * self.qq

    @property
    def mixed_ss(self) -> float:
        return self.m0_sq * self.ss

    @property
    def threshold(self) -> float:
        return (self.mc + self.ms) ** 2


class CumulativePerturbative:
    """Fast continuum integrals for one input point and several Borel masses."""

    def __init__(self, inp: Inputs, m2_values: tuple[float, ...], n_grid: int = 12000):
        self.s = np.linspace(inp.threshold, 80.0, n_grid)
        self.cumulative: dict[float, np.ndarray] = {}
        rho = spectral_matrix(self.s, inp)
        ds = np.diff(self.s)
        for m2 in m2_values:
            integrand = rho * np.exp(-self.s / m2)[:, None, None]
            increments = 0.5 * (integrand[1:] + integrand[:-1]) * ds[:, None, None]
            cumulative = np.empty_like(integrand)
            cumulative[0] = 0.0
            cumulative[1:] = np.cumsum(increments, axis=0)
            self.cumulative[m2] = cumulative

    def matrix(self, m2: float, s0: float) -> np.ndarray:
        values = self.cumulative[m2]
        return np.array(
            [
                [
                    np.interp(s0, self.s, values[:, i, j])
                    for j in range(2)
                ]
                for i in range(2)
            ],
            dtype=float,
        )

    def full(self, m2: float) -> np.ndarray:
        return self.cumulative[m2][-1].copy()


def kallen(s: np.ndarray, m_heavy: float, m_light: float) -> np.ndarray:
    return (s - (m_heavy + m_light) ** 2) * (
        s - (m_heavy - m_light) ** 2
    )


def spectral_matrix(s: np.ndarray, inp: Inputs) -> np.ndarray:
    """Exact-mass LO spectral matrix for the normalized currents.

    The AA entry is algebraically identical to Eq. (two-point f_Ds1) of
    Colangelo-De Fazio-Ozpineci.  Nc=3 is already included.
    """

    mc, ms = inp.mc, inp.ms
    den = mc + ms
    lam_sqrt = np.sqrt(np.clip(kallen(s, mc, ms), 0.0, None))
    common = lam_sqrt / (8.0 * math.pi**2)

    rho_aa = common * (
        2.0
        - (mc * mc + ms * ms + 6.0 * mc * ms) / s
        - (mc * mc - ms * ms) ** 2 / (s * s)
    )
    rho_ab = common * (
        3.0
        * (ms - mc)
        * ((mc + ms) ** 2 - s)
        / (den * s)
    )
    rho_bb = common * (
        -((mc + ms) ** 2 - s)
        * (2.0 * ms * ms - 4.0 * ms * mc + 2.0 * mc * mc + s)
        / (den * den * s)
    )

    out = np.empty((s.size, 2, 2), dtype=float)
    out[:, 0, 0] = rho_aa
    out[:, 0, 1] = rho_ab
    out[:, 1, 0] = rho_ab
    out[:, 1, 1] = rho_bb
    return out


def local_components(m2: float, inp: Inputs) -> dict[str, np.ndarray]:
    """Borel-transformed d=3 and d=5 local contributions.

    The light-propagator expansion used here is

      -<ss>/12 [1 - i ms xslash/4 - ms^2 x^2/8]
      -x^2 <s g sigma.G s>/192,

    with the signs fixed by the independent AA Colangelo expression.
    """

    mc, ms = inp.mc, inp.ms
    den = mc + ms
    ss = inp.ss
    mixed = inp.mixed_ss
    exp_factor = math.exp(-(mc * mc) / m2)

    # Constant quark-condensate term.
    q0 = ss * np.array(
        [
            [mc, mc * mc / den],
            [mc * mc / den, mc**3 / den**2],
        ],
        dtype=float,
    )

    # O(ms) slash-x term.
    q1 = ss * np.array(
        [
            [-ms * mc * mc / (2.0 * m2),
             ms * mc * (1.0 - mc * mc / m2) / (2.0 * den)],
            [ms * mc * (1.0 - mc * mc / m2) / (2.0 * den),
             ms * mc * mc * (1.0 - mc * mc / (2.0 * m2)) / den**2],
        ],
        dtype=float,
    )

    # O(ms^2) x^2 term in the quark-condensate expansion.
    q2 = ss * np.array(
        [
            [ms * ms * mc**3 / (2.0 * m2**2),
             ms * ms * (
                 mc**4 / m2**2 - mc * mc / m2 - 1.0
             ) / (2.0 * den)],
            [ms * ms * (
                 mc**4 / m2**2 - mc * mc / m2 - 1.0
             ) / (2.0 * den),
             ms * ms * (
                 mc**5 / (2.0 * m2**2) - mc**3 / m2
             ) / den**2],
        ],
        dtype=float,
    )

    # Dimension-five mixed condensate from the x^2 propagator term.
    d5 = mixed * np.array(
        [
            [-mc**3 / (4.0 * m2**2),
             -(
                 mc**4 / m2**2 - mc * mc / m2 - 1.0
             ) / (4.0 * den)],
            [-(
                 mc**4 / m2**2 - mc * mc / m2 - 1.0
             ) / (4.0 * den),
             -(
                 mc**5 / (4.0 * m2**2) - mc**3 / (2.0 * m2)
             ) / den**2],
        ],
        dtype=float,
    )

    return {
        "d3_ms0": exp_factor * q0,
        "d3_ms1": exp_factor * q1,
        "d3_ms2": exp_factor * q2,
        "d5_mixed": exp_factor * d5,
    }


def perturbative_matrix(m2: float, s0: float, inp: Inputs, n_grid: int = 5000) -> np.ndarray:
    s = np.linspace(inp.threshold, s0, n_grid)
    rho = spectral_matrix(s, inp)
    weight = np.exp(-s / m2)
    return np.trapezoid(rho * weight[:, None, None], s, axis=0)


def ope_matrix(m2: float, s0: float, inp: Inputs) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    pieces = {"perturbative": perturbative_matrix(m2, s0, inp)}
    pieces.update(local_components(m2, inp))
    total = sum(pieces.values(), np.zeros((2, 2), dtype=float))
    return total, pieces


def rotation_from_matrix(matrix: np.ndarray) -> tuple[float, np.ndarray]:
    """Return theta in [0,90 deg) and diag(R Pi R^T).

    R has rows (sin theta, cos theta) and (cos theta, -sin theta), exactly as
    in the paper's current convention.
    """

    aa, ab, bb = matrix[0, 0], matrix[0, 1], matrix[1, 1]
    theta = 0.5 * math.atan2(-2.0 * ab, aa - bb)
    while theta < 0.0:
        theta += math.pi / 2.0
    while theta >= math.pi / 2.0:
        theta -= math.pi / 2.0
    s, c = math.sin(theta), math.cos(theta)
    rotation = np.array([[s, c], [c, -s]], dtype=float)
    diagonal = rotation @ matrix @ rotation.T
    return math.degrees(theta), diagonal


def physical_residues(m2: float, diagonal: np.ndarray, inp: Inputs) -> tuple[float, float]:
    lam1, lam2 = diagonal[0, 0], diagonal[1, 1]
    if lam1 <= 0.0 or lam2 <= 0.0:
        return math.nan, math.nan
    f1 = math.sqrt(lam1 * math.exp(inp.mass_low**2 / m2) / inp.mass_low**2)
    f2 = math.sqrt(lam2 * math.exp(inp.mass_high**2 / m2) / inp.mass_high**2)
    return f1, f2


def rotation_matrix(theta_deg: float) -> np.ndarray:
    theta = math.radians(theta_deg)
    s, c = math.sin(theta), math.cos(theta)
    return np.array([[s, c], [c, -s]], dtype=float)


def projected_value(matrix: np.ndarray, theta_deg: float, channel: int) -> float:
    vector = rotation_matrix(theta_deg)[channel]
    return float(vector @ matrix @ vector)


def projected_ope(m2: float, s0: float, theta_deg: float, channel: int, inp: Inputs) -> tuple[float, dict[str, float]]:
    matrix, pieces = ope_matrix(m2, s0, inp)
    projected = {
        key: projected_value(value, theta_deg, channel)
        for key, value in pieces.items()
    }
    return projected_value(matrix, theta_deg, channel), projected


def projected_effective_mass(
    m2: float,
    s0: float,
    theta_deg: float,
    channel: int,
    inp: Inputs,
    delta_tau: float = 2.0e-4,
) -> float:
    tau = 1.0 / m2
    plus, _ = projected_ope(1.0 / (tau + delta_tau), s0, theta_deg, channel, inp)
    minus, _ = projected_ope(1.0 / (tau - delta_tau), s0, theta_deg, channel, inp)
    if plus <= 0.0 or minus <= 0.0:
        return math.nan
    mass_sq = -(math.log(plus) - math.log(minus)) / (2.0 * delta_tau)
    return math.sqrt(mass_sq) if mass_sq > 0.0 else math.nan


def fitted_threshold(
    m2: float,
    theta_deg: float,
    channel: int,
    target_mass: float,
    inp: Inputs,
) -> tuple[float, float]:
    candidates = np.linspace(6.55, 10.50, 159)
    best_s0 = math.nan
    best_mass = math.nan
    best_distance = math.inf
    for s0 in candidates:
        mass = projected_effective_mass(m2, float(s0), theta_deg, channel, inp)
        if not math.isfinite(mass):
            continue
        distance = abs(mass - target_mass)
        if distance < best_distance:
            best_distance = distance
            best_s0 = float(s0)
            best_mass = mass
    return best_s0, best_mass


def projected_pole_fraction(m2: float, s0: float, theta_deg: float, channel: int, inp: Inputs) -> float:
    finite = perturbative_matrix(m2, s0, inp)
    full = perturbative_matrix(m2, 80.0, inp, n_grid=20000)
    numerator = projected_value(finite, theta_deg, channel)
    denominator = projected_value(full, theta_deg, channel)
    return numerator / denominator if denominator else math.nan


def physical_threshold_scan(inp: Inputs) -> list[dict[str, float]]:
    """Project at the external Ds angle, then fit one threshold per pole."""

    rows: list[dict[str, float]] = []
    # The mass-fitted pole channels support a lower Borel window than the
    # common-threshold angle scan.  We determine it from the pole and d=5
    # diagnostics instead of inheriting the transition-LCSR window.
    for m2 in np.arange(2.0, 4.5001, 0.05):
        matrix_angle_samples = []
        for s0_mix in np.arange(9.0, 11.0001, 0.1):
            matrix, _ = ope_matrix(float(m2), float(s0_mix), inp)
            matrix_angle, _ = rotation_from_matrix(matrix)
            matrix_angle_samples.append(matrix_angle)
        theta = THETA_DS_DEG

        s01, meff1 = fitted_threshold(
            float(m2), theta, 0, inp.mass_low, inp
        )
        s02, meff2 = fitted_threshold(
            float(m2), theta, 1, inp.mass_high, inp
        )
        if not (math.isfinite(s01) and math.isfinite(s02)):
            continue

        pi1, pieces1 = projected_ope(float(m2), s01, theta, 0, inp)
        pi2, pieces2 = projected_ope(float(m2), s02, theta, 1, inp)
        common_matrix, _ = ope_matrix(float(m2), 10.0, inp)
        rotation = rotation_matrix(theta)
        rotated_common = rotation @ common_matrix @ rotation.T
        offdiag_norm = abs(float(rotated_common[0, 1])) / math.sqrt(
            abs(float(rotated_common[0, 0] * rotated_common[1, 1]))
        )
        f1 = (
            math.sqrt(pi1 * math.exp(inp.mass_low**2 / m2) / inp.mass_low**2)
            if pi1 > 0.0
            else math.nan
        )
        f2 = (
            math.sqrt(pi2 * math.exp(inp.mass_high**2 / m2) / inp.mass_high**2)
            if pi2 > 0.0
            else math.nan
        )
        pole1 = projected_pole_fraction(float(m2), s01, theta, 0, inp)
        pole2 = projected_pole_fraction(float(m2), s02, theta, 1, inp)
        d51 = abs(pieces1["d5_mixed"]) / abs(pi1) if pi1 else math.nan
        d52 = abs(pieces2["d5_mixed"]) / abs(pi2) if pi2 else math.nan

        rows.append(
            {
                "M2_GeV2": float(m2),
                "theta_deg": theta,
                "theta_matrix_median_deg": float(np.median(matrix_angle_samples)),
                "theta_matrix_s0mix_min_deg": float(min(matrix_angle_samples)),
                "theta_matrix_s0mix_max_deg": float(max(matrix_angle_samples)),
                "Pi12_normalized_at_s0mix10": offdiag_norm,
                "s01_fitted_GeV2": s01,
                "s02_fitted_GeV2": s02,
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
        )
    return rows


def pole_fraction(m2: float, s0: float, inp: Inputs) -> float:
    finite = perturbative_matrix(m2, s0, inp)
    # A fixed high endpoint is numerically indistinguishable from infinity in
    # the charm Borel window.
    full = perturbative_matrix(m2, 80.0, inp, n_grid=20000)
    finite_trace = float(np.trace(finite))
    full_trace = float(np.trace(full))
    return finite_trace / full_trace if full_trace else math.nan


def effective_masses(m2: float, s0: float, inp: Inputs, delta_tau: float = 2.0e-4) -> tuple[float, float]:
    tau = 1.0 / m2

    def diagonal_at_tau(value: float) -> np.ndarray:
        local_m2 = 1.0 / value
        matrix, _ = ope_matrix(local_m2, s0, inp)
        _, diag = rotation_from_matrix(matrix)
        return np.diag(diag)

    plus = diagonal_at_tau(tau + delta_tau)
    minus = diagonal_at_tau(tau - delta_tau)
    if np.any(plus <= 0.0) or np.any(minus <= 0.0):
        return math.nan, math.nan
    mass_sq = -(np.log(plus) - np.log(minus)) / (2.0 * delta_tau)
    return tuple(float(math.sqrt(max(x, 0.0))) for x in mass_sq)


def grid_scan(inp: Inputs) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for m2 in np.arange(3.0, 6.0001, 0.1):
        for s0 in np.arange(7.8, 11.0001, 0.1):
            matrix, pieces = ope_matrix(float(m2), float(s0), inp)
            theta, diagonal = rotation_from_matrix(matrix)
            f1, f2 = physical_residues(float(m2), diagonal, inp)
            meff1, meff2 = effective_masses(float(m2), float(s0), inp)
            d5_norm = float(np.linalg.norm(pieces["d5_mixed"]))
            total_norm = float(np.linalg.norm(matrix))
            rows.append(
                {
                    "M2_GeV2": float(m2),
                    "s0_GeV2": float(s0),
                    "theta_deg": theta,
                    "Pi_AA_GeV4": float(matrix[0, 0]),
                    "Pi_AB_GeV4": float(matrix[0, 1]),
                    "Pi_BB_GeV4": float(matrix[1, 1]),
                    "Pi12_rotated_GeV4": float(diagonal[0, 1]),
                    "lambda1_GeV4": float(diagonal[0, 0]),
                    "lambda2_GeV4": float(diagonal[1, 1]),
                    "f1_GeV": f1,
                    "f2_GeV": f2,
                    "m_eff1_GeV": meff1,
                    "m_eff2_GeV": meff2,
                    "pole_fraction_trace": pole_fraction(float(m2), float(s0), inp),
                    "d5_fraction_frobenius": d5_norm / total_norm if total_norm else math.nan,
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, float]], inp: Inputs) -> str:
    accepted = [
        row
        for row in rows
        if row["pole_fraction_trace"] >= 0.50
        and row["d5_fraction_frobenius"] <= 0.15
        and math.isfinite(row["f1_GeV"])
        and math.isfinite(row["f2_GeV"])
        and abs(row["m_eff1_GeV"] - inp.mass_low) <= 0.25
        and abs(row["m_eff2_GeV"] - inp.mass_high) <= 0.25
    ]

    lines = [
        "Ds1 normalized-current AA/AB/BB two-point sum rule",
        "==================================================",
        "OPE: exact-mass LO perturbative + local d=3 through ms^2 + local d=5 mixed condensate.",
        "Ordinary local condensates are included because this is a two-point SVZ sum rule.",
        "No external f1, f2, fA, or fT input is used.",
        "",
        f"central inputs: mc={inp.mc:.4f} GeV, ms={inp.ms:.4f} GeV, "
        f"<ss>={inp.ss:.7f} GeV^3, m0^2={inp.m0_sq:.3f} GeV^2",
        f"scanned points: {len(rows)}; accepted by stated diagnostics: {len(accepted)}",
    ]
    if not accepted:
        lines.append("No point passed every mass/pole/OPE criterion; inspect the full grid before quoting residues.")
        return "\n".join(lines) + "\n"

    for key, label, unit in (
        ("theta_deg", "theta", "deg"),
        ("f1_GeV", "f1", "GeV"),
        ("f2_GeV", "f2", "GeV"),
        ("m_eff1_GeV", "m_eff1", "GeV"),
        ("m_eff2_GeV", "m_eff2", "GeV"),
        ("pole_fraction_trace", "pole fraction", ""),
        ("d5_fraction_frobenius", "d5 fraction", ""),
    ):
        values = np.asarray([row[key] for row in accepted], dtype=float)
        q16, q50, q84 = np.percentile(values, [16.0, 50.0, 84.0])
        lines.append(f"{label}: {q50:.6g} [{q16:.6g}, {q84:.6g}] {unit}".rstrip())

    return "\n".join(lines) + "\n"


def summarize_physical_scan(rows: list[dict[str, float]], inp: Inputs) -> str:
    accepted = [
        row
        for row in rows
        if row["pole_fraction1"] >= 0.50
        and row["pole_fraction2"] >= 0.50
        and row["d5_fraction1"] <= 0.15
        and row["d5_fraction2"] <= 0.15
        and abs(row["m_eff1_GeV"] - inp.mass_low) <= 0.03
        and abs(row["m_eff2_GeV"] - inp.mass_high) <= 0.03
        and math.isfinite(row["f1_GeV"])
        and math.isfinite(row["f2_GeV"])
    ]

    lines = [
        "Ds1 physical-pole residue extraction from the AA/AB/BB matrix",
        "============================================================",
        f"The projection angle is theta_Ds={THETA_DS_DEG:.1f}+-{THETA_DS_SIGMA_DEG:.1f} deg.",
        f"Angle source: {ANGLE_SOURCE}.",
        "The angle from diagonalizing the truncated matrix is diagnostic only.",
        "Separate projected thresholds s01 and s02 are then fitted to the two measured masses.",
        "No external decay constant or overlap parameter enters.",
        "",
        f"Borel points: {len(rows)}; accepted by pole/OPE/mass diagnostics: {len(accepted)}",
    ]
    if not accepted:
        lines.append("No Borel point passed every stated diagnostic.")
        return "\n".join(lines) + "\n"

    for key, label, unit in (
        ("theta_deg", "theta", "deg"),
        ("theta_matrix_median_deg", "truncated-matrix diagnostic theta", "deg"),
        ("Pi12_normalized_at_s0mix10", "normalized off-diagonal residual", ""),
        ("s01_fitted_GeV2", "s01", "GeV^2"),
        ("s02_fitted_GeV2", "s02", "GeV^2"),
        ("f1_GeV", "f1", "GeV"),
        ("f2_GeV", "f2", "GeV"),
        ("pole_fraction1", "pole fraction 1", ""),
        ("pole_fraction2", "pole fraction 2", ""),
        ("d5_fraction1", "d5 fraction 1", ""),
        ("d5_fraction2", "d5 fraction 2", ""),
    ):
        values = np.asarray([row[key] for row in accepted], dtype=float)
        q16, q50, q84 = np.percentile(values, [16.0, 50.0, 84.0])
        lines.append(f"{label}: {q50:.6g} [{q16:.6g}, {q84:.6g}] {unit}".rstrip())

    return "\n".join(lines) + "\n"


def mc_sample(rng: np.random.Generator) -> Inputs:
    mc = float(np.clip(rng.normal(1.27, 0.02), 1.20, 1.34))
    ms = float(np.clip(rng.normal(0.093, 0.011), 0.060, 0.130))
    condensate_scale = float(np.clip(rng.normal(0.240, 0.010), 0.210, 0.270))
    kappa = float(np.clip(rng.normal(0.80, 0.10), 0.50, 1.10))
    m0_sq = float(np.clip(rng.normal(0.80, 0.20), 0.30, 1.30))
    return Inputs(
        mc=mc,
        ms=ms,
        qq=-(condensate_scale**3),
        kappa_s=kappa,
        m0_sq=m0_sq,
    )


def monte_carlo_scan(n_samples: int = 2000) -> list[dict[str, float]]:
    rng = np.random.default_rng(20260722)
    rows: list[dict[str, float]] = []
    delta_tau = 2.0e-4
    threshold_candidates = np.linspace(6.55, 12.00, 219)

    for _ in range(n_samples):
        inp = mc_sample(rng)
        m2 = float(rng.uniform(2.0, 2.7))
        tau = 1.0 / m2
        m2_plus = 1.0 / (tau + delta_tau)
        m2_minus = 1.0 / (tau - delta_tau)
        evaluator = CumulativePerturbative(inp, (m2, m2_plus, m2_minus), n_grid=5000)

        s0_mix = float(rng.uniform(9.0, 11.0))
        mix_matrix = evaluator.matrix(m2, s0_mix) + sum(
            local_components(m2, inp).values(), np.zeros((2, 2), dtype=float)
        )
        theta_matrix, _ = rotation_from_matrix(mix_matrix)
        theta = sample_gaussian_angle(rng, THETA_DS_DEG, THETA_DS_SIGMA_DEG)
        rotation = rotation_matrix(theta)
        rotated_mix = rotation @ mix_matrix @ rotation.T
        offdiag_norm = abs(float(rotated_mix[0, 1])) / math.sqrt(
            abs(float(rotated_mix[0, 0] * rotated_mix[1, 1]))
        )

        fitted: list[tuple[float, float, float, float, float, float]] = []
        failed = False
        for channel, target_mass in enumerate((inp.mass_low, inp.mass_high)):
            best = None
            for s0 in threshold_candidates:
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
                mass_sq = -(math.log(plus) - math.log(minus)) / (2.0 * delta_tau)
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
            pert = evaluator.matrix(m2, s0)
            locals_now = local_components(m2, inp)
            total = pert + sum(locals_now.values(), np.zeros((2, 2), dtype=float))
            pi = projected_value(total, theta, channel)
            pole = projected_value(pert, theta, channel) / projected_value(
                evaluator.full(m2), theta, channel
            )
            d5 = abs(projected_value(locals_now["d5_mixed"], theta, channel)) / abs(pi)
            residue = math.sqrt(
                pi * math.exp(target_mass**2 / m2) / target_mass**2
            ) if pi > 0.0 else math.nan
            fitted.append((s0, meff, pi, pole, d5, residue))

        if failed or len(fitted) != 2:
            continue
        one, two = fitted
        accepted = (
            one[3] >= 0.50
            and two[3] >= 0.50
            and one[4] <= 0.15
            and two[4] <= 0.15
            and abs(one[1] - inp.mass_low) <= 0.03
            and abs(two[1] - inp.mass_high) <= 0.03
            and math.isfinite(one[5])
            and math.isfinite(two[5])
        )
        rows.append(
            {
                "accepted": int(accepted),
                "M2_GeV2": m2,
                "s0_mix_GeV2": s0_mix,
                "mc_GeV": inp.mc,
                "ms_GeV": inp.ms,
                "ss_GeV3": inp.ss,
                "m0_sq_GeV2": inp.m0_sq,
                "theta_deg": theta,
                "theta_matrix_deg": theta_matrix,
                "Pi12_normalized": offdiag_norm,
                "s01_GeV2": one[0],
                "s02_GeV2": two[0],
                "m_eff1_GeV": one[1],
                "m_eff2_GeV": two[1],
                "f1_GeV": one[5],
                "f2_GeV": two[5],
                "pole_fraction1": one[3],
                "pole_fraction2": two[3],
                "d5_fraction1": one[4],
                "d5_fraction2": two[4],
            }
        )
    return rows


def summarize_mc(rows: list[dict[str, float]]) -> str:
    accepted = [row for row in rows if row["accepted"]]
    lines = [
        "Ds1 fixed-angle AA/AB/BB two-point Monte Carlo",
        "================================================",
        "Inputs varied: mc, ms, <qq>, <ss>/<qq>, m0^2, M^2, and s0_mix.",
        f"theta_Ds={THETA_DS_DEG:.1f}+-{THETA_DS_SIGMA_DEG:.1f} deg is sampled from {ANGLE_SOURCE}.",
        "The angle from diagonalizing the truncated matrix is diagnostic only.",
        "The two physical thresholds are refitted to the measured masses in every sample.",
        f"generated samples with valid fits: {len(rows)}; accepted: {len(accepted)}",
    ]
    if not accepted:
        lines.append("No sample passed all diagnostics.")
        return "\n".join(lines) + "\n"
    for key, label, unit in (
        ("theta_deg", "theta input", "deg"),
        ("theta_matrix_deg", "truncated-matrix diagnostic theta", "deg"),
        ("Pi12_normalized", "normalized off-diagonal residual", ""),
        ("s01_GeV2", "s01", "GeV^2"),
        ("s02_GeV2", "s02", "GeV^2"),
        ("f1_GeV", "f1", "GeV"),
        ("f2_GeV", "f2", "GeV"),
    ):
        values = np.asarray([row[key] for row in accepted], dtype=float)
        q16, q50, q84 = np.percentile(values, [16.0, 50.0, 84.0])
        lines.append(f"{label}: {q50:.6g} [{q16:.6g}, {q84:.6g}] {unit}".rstrip())
    return "\n".join(lines) + "\n"


def main() -> None:
    inp = Inputs()
    rows = grid_scan(inp)
    grid_path = OUT / "twopoint_ds1_matrix_grid.csv"
    summary_path = OUT / "twopoint_ds1_matrix_summary.txt"
    write_csv(grid_path, rows)
    summary = summarize(rows, inp)
    summary_path.write_text(summary)
    print(summary)

    physical_rows = physical_threshold_scan(inp)
    physical_grid_path = OUT / "twopoint_ds1_physical_residue_grid.csv"
    physical_summary_path = OUT / "twopoint_ds1_physical_residue_summary.txt"
    write_csv(physical_grid_path, physical_rows)
    physical_summary = summarize_physical_scan(physical_rows, inp)
    physical_summary_path.write_text(physical_summary)
    print(physical_summary)
    print(f"Wrote {grid_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {physical_grid_path}")
    print(f"Wrote {physical_summary_path}")

    mc_rows = monte_carlo_scan()
    mc_path = OUT / "twopoint_ds1_matrix_mc.csv"
    mc_summary_path = OUT / "twopoint_ds1_matrix_mc_summary.txt"
    write_csv(mc_path, mc_rows)
    mc_summary = summarize_mc(mc_rows)
    mc_summary_path.write_text(mc_summary)
    print(mc_summary)
    print(f"Wrote {mc_path}")
    print(f"Wrote {mc_summary_path}")


if __name__ == "__main__":
    main()
