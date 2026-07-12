#!/usr/bin/env python3
"""Direct-Borel evaluator for all-positive Bc1 -> Bc gamma G2 terms."""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from bc_ps_complete_analysis import (
    Inputs,
    gauss_legendre_integral,
    summary_stats,
    twopoint_residues,
    ps_couplings,
    width_ps_keV,
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
    code = compile(text, "<ps-g2-coefficient>", "eval")
    allowed = {"mb", "mc", "p2", "pq", "pi"}
    unknown = set(code.co_names) - allowed
    if unknown:
        raise ValueError(f"Unsupported symbols in coefficient {expr}: {sorted(unknown)}")

    def evaluate(mb: float, mc: float, p2: float, pq: float):
        return eval(code, {"__builtins__": {}}, {"mb": mb, "mc": mc, "p2": p2, "pq": pq, "pi": math.pi})

    return evaluate


def kallen(s: float, m1: float, m2: float) -> float:
    return (s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)


def z_limits(s0: float, active_mass: float, spectator_mass: float):
    lam = kallen(s0, active_mass, spectator_mass)
    if lam <= 0:
        return None
    root = math.sqrt(lam)
    lo = (s0 + active_mass**2 - spectator_mass**2 - root) / (2.0 * s0)
    hi = (s0 + active_mass**2 - spectator_mass**2 + root) / (2.0 * s0)
    lo = max(1.0e-8, lo)
    hi = min(1.0 - 1.0e-8, hi)
    return (lo, hi) if hi > lo else None


def class_prefactor(row: dict[str, str], g2_value: float) -> float:
    if row["class"] == "single_line_G2":
        return NC * g2_value / 12.0
    if row["class"] == "open_pair_GG":
        return g2_value * ((NC**2 - 1.0) / 2.0) / 96.0
    raise ValueError(row["class"])


def mapped_powers(row: dict[str, str]):
    n1, n2, n3 = (int(row[f"eff_n{i}"]) for i in (1, 2, 3))
    if row["topology"] == "c_photon":
        return n1, n2, n3, E_C
    if row["topology"] == "bbar_photon":
        return n3, n2, n1, E_BBAR
    raise ValueError(row["topology"])


def row_integral(row, coeff_func, M2, s0, inp: Inputs, channel: ChannelPoint, g2_value: float):
    na, nb, ns, _charge_abs = mapped_powers(row)
    active_mass, spectator_mass = (inp.mc, inp.mb) if row["topology"] == "c_photon" else (inp.mb, inp.mc)
    limits = z_limits(s0, active_mass, spectator_mass)
    if limits is None:
        return 0.0
    zlo, zhi = limits
    nsum = na + nb + ns
    gamma_den = math.gamma(na) * math.gamma(nb) * math.gamma(ns)
    charge = -E_C if row["topology"] == "c_photon" else E_BBAR
    coeff_val = complex(coeff_func(inp.mb, inp.mc, channel.mass_final**2, channel.pq))
    coeff_real = (-1j * coeff_val).real
    overall = charge * class_prefactor(row, g2_value) * coeff_real / (16.0 * math.pi**2 * gamma_den)

    def integrand(z):
        z = np.asarray(z, dtype=float)
        x = 0.5 * (1.0 - z)
        y = x
        tau = 1.0 / (M2 * z * (1.0 - z))
        mass_bar = (1.0 - z) * active_mass**2 + z * spectator_mass**2
        return (
            overall
            * x ** (na - 1)
            * y ** (nb - 1)
            * z ** (ns - 1)
            * tau ** (nsum - 3)
            * (M2 / z)
            * np.exp(-tau * mass_bar)
        )

    return gauss_legendre_integral(integrand, zlo, zhi, n=160)


def load_rows():
    rows = list(csv.DictReader((OUT / "step3_ps_g2_denominator_reduction.csv").open()))
    rows = [r for r in rows if all(int(r[f"eff_n{i}"]) > 0 for i in (1, 2, 3))]
    return [(r, parse_coeff(r["coefficient"])) for r in rows]


def moments(M2, s0, inp, channel, parsed_rows):
    totals = {"A_single": 0.0, "A_open": 0.0, "B_single": 0.0, "B_open": 0.0}
    for row, coeff in parsed_rows:
        val = row_integral(row, coeff, M2, s0, inp, channel, G2_DEFAULT)
        key = f"{row['current']}_{'single' if row['class'] == 'single_line_G2' else 'open'}"
        totals[key] += val
    totals["A_total"] = totals["A_single"] + totals["A_open"]
    totals["B_total"] = totals["B_single"] + totals["B_open"]
    return totals


def run_grid():
    inp = Inputs()
    parsed = load_rows()
    rows = []
    for M2 in [7.0, 8.0, 9.0]:
        for s0 in [53.0, 54.0, 55.0]:
            res = twopoint_residues(M2, s0, inp)
            pert = ps_couplings(M2, s0, inp, inp.theta_deg, res["f1_GeV"], res["f2_GeV"])
            for label, mi, mf, f_i, gkey, wkey in [
                ("Bc1(6743)->Bc gamma", inp.m_low, inp.m_bc, res["f1_GeV"], "g1_GeV_inv", "Gamma1_keV"),
                ("Bc1(6750)->Bc gamma", inp.m_high, inp.m_bc, res["f2_GeV"], "g2_GeV_inv", "Gamma2_keV"),
            ]:
                ch = ChannelPoint(label, mi, mf, 0.5 * (mi**2 - mf**2))
                mom = moments(M2, s0, inp, ch, parsed)
                th = math.radians(inp.theta_deg)
                st, ct = math.sin(th), math.cos(th)
                if "6743" in label:
                    pi_g2 = st * mom["A_total"] + ct * mom["B_total"]
                else:
                    pi_g2 = ct * mom["A_total"] - st * mom["B_total"]
                pref = math.exp((mi**2 + mf**2) / (2.0 * M2)) * (inp.mb + inp.mc) / (mf**2 * inp.f_bc * mi * f_i)
                g_g2 = pref * pi_g2
                g_total = pert[gkey] + g_g2
                rows.append(
                    {
                        "M2": M2,
                        "s0": s0,
                        "channel": label,
                        "g_pert_GeV_inv": pert[gkey],
                        "Gamma_pert_keV": pert[wkey],
                        **mom,
                        "Pi_G2_rotated": pi_g2,
                        "g_G2_all_positive_GeV_inv": g_g2,
                        "g_pert_plus_G2ap_GeV_inv": g_total,
                        "Gamma_pert_plus_G2ap_keV": width_ps_keV(g_total, mi, mf),
                        "delta_g_over_g": g_g2 / pert[gkey] if pert[gkey] else 0.0,
                    }
                )
    return rows


def main():
    rows = run_grid()
    write_csv(OUT / "bc_ps_g2_all_positive_direct_borel_grid.csv", rows)
    summary_rows = []
    for channel in sorted({r["channel"] for r in rows}):
        sub = [r for r in rows if r["channel"] == channel]
        for obs in ["g_G2_all_positive_GeV_inv", "g_pert_plus_G2ap_GeV_inv", "Gamma_pert_plus_G2ap_keV", "delta_g_over_g"]:
            summary_rows.append({"channel": channel, "observable": obs, **summary_stats([r[obs] for r in sub])})
    write_csv(OUT / "bc_ps_g2_all_positive_direct_borel_summary.csv", summary_rows)
    lines = [
        "Bc1 -> Bc gamma all-positive radiative G2 direct-Borel result",
        "==============================================================",
        "",
    ]
    for row in summary_rows:
        if row["observable"] == "Gamma_pert_plus_G2ap_keV":
            lines.append(
                "{channel}: Gamma(pert + all-positive G2) = {median:.4g} [{q16:.4g}, {q84:.4g}] keV".format(**row)
            )
        if row["observable"] == "delta_g_over_g":
            lines.append("{channel}: delta g/g = {median:.4g} [{q16:.4g}, {q84:.4g}]".format(**row))
    (OUT / "bc_ps_g2_all_positive_direct_borel_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
