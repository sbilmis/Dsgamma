#!/usr/bin/env python3
"""Direct-Borel evaluator for the all-positive vector radiative G2 sector.

This is intentionally separate from the perturbative vector driver.  It turns
the denominator-reduced all-positive rows in
step3_vec_g2_denominator_reduction.csv into diagonal double-Borel moments in
the convention M1^2=M2^2=2 M^2.

Important convention:
- c_photon: D1,D2 are the active charm-line denominators and D3 is the
  anti-bottom spectator.
- bbar_photon: D3,D2 are the active anti-bottom denominators and D1 is the
  charm spectator.
- Explicit projector factors p2 and pq are evaluated at the physical final
  vector mass and photon kinematics for each channel.
- Single-line and open-pair G2 prefactors follow the Bc-mixing workbench:
  Nc G2/12 and G2*((Nc^2-1)/2)/96, respectively.

The output is a controlled all-positive G2 moment in the same diagonal Borel
calibration used by the Bc radiative pilot.  Contact/crossed rows are not
included here; their support audit is handled separately.
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from bc_vec_complete_analysis import (
    Inputs,
    GEV_TO_KEV,
    gauss_legendre_integral,
    photon_energy,
    summary_stats,
    twopoint_residues,
    bcstar_decay_constant,
    vec_couplings,
    width_vec_keV,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"

NC = 3.0
G2_DEFAULT = 4.0 * math.pi**2 * 0.012
E_C = 2.0 / 3.0
E_BBAR = 1.0 / 3.0

@dataclass(frozen=True)
class ChannelPoint:
    label: str
    mass_initial: float
    mass_final: float
    pq: float


def parse_coeff(expr: str):
    text = expr.replace("^", "**").replace("Pi", "pi")
    text = re.sub(r"\bI\b", "1j", text)
    code = compile(text, "<g2-coefficient>", "eval")
    allowed_names = {"mb", "mc", "p2", "pq", "pi"}
    unknown = set(code.co_names) - allowed_names
    if unknown:
        raise ValueError(f"Unsupported symbols in coefficient {expr}: {sorted(unknown)}")

    def evaluate(mb: float, mc: float, p2: float, pq: float):
        return eval(code, {"__builtins__": {}}, {"mb": mb, "mc": mc, "p2": p2, "pq": pq, "pi": math.pi})

    return evaluate


def kallen(s: float, m1: float, m2: float) -> float:
    return (s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)


def z_limits(s0: float, active_mass: float, spectator_mass: float) -> tuple[float, float] | None:
    lam = kallen(s0, active_mass, spectator_mass)
    if lam <= 0:
        return None
    root = math.sqrt(lam)
    lo = (s0 + active_mass**2 - spectator_mass**2 - root) / (2.0 * s0)
    hi = (s0 + active_mass**2 - spectator_mass**2 + root) / (2.0 * s0)
    lo = max(1.0e-8, lo)
    hi = min(1.0 - 1.0e-8, hi)
    if hi <= lo:
        return None
    return lo, hi


def row_powers(row: dict[str, str]) -> tuple[int, int, int, float, float, float]:
    n1, n2, n3 = (int(row[f"eff_n{i}"]) for i in (1, 2, 3))
    if row["topology"] == "c_photon":
        # active after/before photon, spectator
        return n1, n2, n3, E_C, 1.0, -1.0
    if row["topology"] == "bbar_photon":
        # active after/before photon, spectator.  Charge sign follows the
        # perturbative vector driver: +e_bbar*bbar - e_c*c.
        return n3, n2, n1, E_BBAR, 1.0, 1.0
    raise ValueError(row["topology"])


def class_prefactor(row: dict[str, str], g2_value: float) -> float:
    if row["class"] == "single_line_G2":
        return NC * g2_value / 12.0
    if row["class"] == "open_pair_GG":
        return g2_value * ((NC**2 - 1.0) / 2.0) / 96.0
    raise ValueError(row["class"])


def row_integral(
    row: dict[str, str],
    coeff_func,
    M2: float,
    s0: float,
    inp: Inputs,
    channel: ChannelPoint,
    g2_value: float,
) -> float:
    na, nb, ns, _charge, _unused, sign = row_powers(row)
    if row["topology"] == "c_photon":
        active_mass, spectator_mass = inp.mc, inp.mb
    else:
        active_mass, spectator_mass = inp.mb, inp.mc

    limits = z_limits(s0, active_mass, spectator_mass)
    if limits is None:
        return 0.0
    zlo, zhi = limits
    nsum = na + nb + ns
    gamma_den = math.gamma(na) * math.gamma(nb) * math.gamma(ns)
    pref = class_prefactor(row, g2_value)
    charge = -E_C if row["topology"] == "c_photon" else E_BBAR
    # The -I phase follows the direct-Borel phase in BcMixingMomentum.wl;
    # the loop factor is the usual four-dimensional one-loop 1/(16 pi^2).
    coeff_val = complex(coeff_func(inp.mb, inp.mc, channel.mass_final**2, channel.pq))
    coeff_real = (-1j * coeff_val).real
    overall = sign * charge * pref * coeff_real / (16.0 * math.pi**2 * gamma_den)

    def integrand(z):
        z = np.asarray(z, dtype=float)
        x = 0.5 * (1.0 - z)
        y = x
        tau = 1.0 / (M2 * z * (1.0 - z))
        mass_bar = (1.0 - z) * active_mass**2 + z * spectator_mass**2
        jac = M2 / z
        return (
            overall
            * (x ** (na - 1))
            * (y ** (nb - 1))
            * (z ** (ns - 1))
            * (tau ** (nsum - 3))
            * jac
            * np.exp(-tau * mass_bar)
        )

    return gauss_legendre_integral(integrand, zlo, zhi, n=160)


def load_all_positive_rows():
    rows = list(csv.DictReader((OUT / "step3_vec_g2_denominator_reduction.csv").open()))
    rows = [r for r in rows if all(int(r[f"eff_n{i}"]) > 0 for i in (1, 2, 3))]
    parsed = [(r, parse_coeff(r["coefficient"])) for r in rows]
    return parsed


def g2_moments(M2: float, s0: float, inp: Inputs, channel: ChannelPoint, parsed_rows, g2_value: float):
    totals = {
        "A_single": 0.0,
        "A_open": 0.0,
        "B_single": 0.0,
        "B_open": 0.0,
    }
    for row, coeff in parsed_rows:
        val = row_integral(row, coeff, M2, s0, inp, channel, g2_value)
        key = f"{row['current']}_{'single' if row['class'] == 'single_line_G2' else 'open'}"
        totals[key] += val
    totals["A_total"] = totals["A_single"] + totals["A_open"]
    totals["B_total"] = totals["B_single"] + totals["B_open"]
    return totals


def run_grid():
    inp = Inputs()
    parsed_rows = load_all_positive_rows()
    rows = []
    for M2 in [7.0, 8.0, 9.0]:
        for s0 in [53.0, 54.0, 55.0]:
            res = twopoint_residues(M2, s0, inp)
            fstar = bcstar_decay_constant(6.0, 42.0, inp, 1.0 / 3.0)
            pert = vec_couplings(M2, s0, inp, inp.theta_deg, res["f1_GeV"], res["f2_GeV"], fstar)
            for label, mi, mf, f_i, gkey, wkey in [
                ("Bc1(6743)->Bc* gamma", inp.m_low, inp.m_bcstar, res["f1_GeV"], "g1_GeV_inv", "Gamma1_keV"),
                ("Bc1(6750)->Bc* gamma", inp.m_high, inp.m_bcstar, res["f2_GeV"], "g2_GeV_inv", "Gamma2_keV"),
            ]:
                ch = ChannelPoint(label, mi, mf, 0.5 * (mi**2 - mf**2))
                moments = g2_moments(M2, s0, inp, ch, parsed_rows, G2_DEFAULT)
                th = math.radians(inp.theta_deg)
                st, ct = math.sin(th), math.cos(th)
                if "6743" in label:
                    pi_g2 = st * moments["A_total"] + ct * moments["B_total"]
                    pref = math.exp((mi**2 + mf**2) / (2.0 * M2)) / (mf * fstar * mi * f_i)
                else:
                    pi_g2 = ct * moments["A_total"] - st * moments["B_total"]
                    pref = math.exp((mi**2 + mf**2) / (2.0 * M2)) / (mf * fstar * mi * f_i)
                g_g2 = pref * pi_g2
                g_total = pert[gkey] + g_g2
                rows.append(
                    {
                        "M2": M2,
                        "s0": s0,
                        "channel": label,
                        "f_Bcstar_standard_GeV": fstar,
                        "g_pert_GeV_inv": pert[gkey],
                        "Gamma_pert_keV": pert[wkey],
                        **moments,
                        "Pi_G2_rotated": pi_g2,
                        "g_G2_all_positive_GeV_inv": g_g2,
                        "g_pert_plus_G2ap_GeV_inv": g_total,
                        "Gamma_pert_plus_G2ap_keV": width_vec_keV(g_total, mi, mf),
                        "delta_g_over_g": g_g2 / pert[gkey] if pert[gkey] else 0.0,
                    }
                )
    return rows


def main() -> None:
    rows = run_grid()
    write_csv(OUT / "bc_vec_g2_all_positive_direct_borel_grid.csv", rows)
    summary_rows = []
    for channel in sorted({r["channel"] for r in rows}):
        sub = [r for r in rows if r["channel"] == channel]
        for obs in [
            "g_G2_all_positive_GeV_inv",
            "g_pert_plus_G2ap_GeV_inv",
            "Gamma_pert_plus_G2ap_keV",
            "delta_g_over_g",
        ]:
            summary_rows.append({"channel": channel, "observable": obs, **summary_stats([r[obs] for r in sub])})
    write_csv(OUT / "bc_vec_g2_all_positive_direct_borel_summary.csv", summary_rows)

    lines = [
        "Bc1 -> Bc* gamma all-positive radiative G2 direct-Borel result",
        "================================================================",
        "",
        "Convention:",
        "- diagonal double Borel M1^2=M2^2=2 M^2;",
        "- physical p^2=m_Bc*^2 and p.q=(m_i^2-m_Bc*^2)/2 in projector prefactors;",
        "- single-line prefactor Nc G2/12;",
        "- open-pair prefactor G2*((Nc^2-1)/2)/96;",
        "- contact/crossed rows are not included.",
        "",
    ]
    for row in summary_rows:
        if row["observable"] == "Gamma_pert_plus_G2ap_keV":
            lines.append(
                "{channel}: Gamma(pert + all-positive G2) = {median:.4g} "
                "[{q16:.4g}, {q84:.4g}] keV".format(**row)
            )
        if row["observable"] == "delta_g_over_g":
            lines.append(
                "{channel}: delta g/g from all-positive G2 = {median:.4g} "
                "[{q16:.4g}, {q84:.4g}]".format(**row)
            )
    (OUT / "bc_vec_g2_all_positive_direct_borel_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
