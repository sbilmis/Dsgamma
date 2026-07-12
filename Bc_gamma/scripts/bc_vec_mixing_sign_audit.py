#!/usr/bin/env python3
"""Audit sign, mixing and state-assignment choices for Bc1 -> Bc* gamma.

The baseline vector result has a larger Bc1(6750)->Bc* gamma width than most
model calculations.  This script isolates whether the hierarchy comes from
phase space, the mixed-current rotation, a possible relative sign of the
tensor-current radiative amplitude, or an assignment swap of the two mixed
eigenstates.
"""

from __future__ import annotations

import csv
import importlib.util
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


def load_vector_module():
    path = ROOT / "scripts" / "bc_vec_complete_analysis.py"
    spec = importlib.util.spec_from_file_location("bc_vec_complete_analysis_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def standard_fbcstar(mod, inp):
    rows = []
    for M2v in [5.0, 6.0, 7.0, 8.0]:
        for s0v in [40.0, 42.0, 44.0]:
            rows.append(mod.bcstar_decay_constant(M2v, s0v, inp, 1.0 / 3.0))
    rows.sort()
    return rows[len(rows) // 2]


def row_ingredients(mod, inp, M2, s0, fstar):
    res = mod.twopoint_residues(M2, s0, inp)
    theta = math.radians(res["theta_sumrule_deg"])
    st, ct = math.sin(theta), math.cos(theta)
    pq_low = 0.5 * (inp.m_low**2 - inp.m_bcstar**2)
    pq_high = 0.5 * (inp.m_high**2 - inp.m_bcstar**2)
    amp_a = mod.hard_integral(M2, s0, lambda x: mod.rho_vec_A_diag(x, inp), inp)
    amp_b_low = mod.hard_integral(M2, s0, lambda x: mod.rho_vec_B_diag(x, inp, pq_low), inp)
    amp_b_high = mod.hard_integral(M2, s0, lambda x: mod.rho_vec_B_diag(x, inp, pq_high), inp)
    pref_low_f1 = math.exp((inp.m_low**2 + inp.m_bcstar**2) / (2.0 * M2)) / (
        inp.m_bcstar * fstar * inp.m_low * res["f1_GeV"]
    )
    pref_high_f2 = math.exp((inp.m_high**2 + inp.m_bcstar**2) / (2.0 * M2)) / (
        inp.m_bcstar * fstar * inp.m_high * res["f2_GeV"]
    )
    pref_low_f2 = math.exp((inp.m_low**2 + inp.m_bcstar**2) / (2.0 * M2)) / (
        inp.m_bcstar * fstar * inp.m_low * res["f2_GeV"]
    )
    pref_high_f1 = math.exp((inp.m_high**2 + inp.m_bcstar**2) / (2.0 * M2)) / (
        inp.m_bcstar * fstar * inp.m_high * res["f1_GeV"]
    )
    return {
        "M2": M2,
        "s0": s0,
        "theta_deg": res["theta_sumrule_deg"],
        "sin_theta": st,
        "cos_theta": ct,
        "f1_GeV": res["f1_GeV"],
        "f2_GeV": res["f2_GeV"],
        "amp_A": amp_a,
        "amp_B_low": amp_b_low,
        "amp_B_high": amp_b_high,
        "pref_low_f1": pref_low_f1,
        "pref_high_f2": pref_high_f2,
        "pref_low_f2": pref_low_f2,
        "pref_high_f1": pref_high_f1,
    }


def variant_widths(mod, inp, ing, variant):
    st, ct = ing["sin_theta"], ing["cos_theta"]
    a = ing["amp_A"]
    b_low = ing["amp_B_low"]
    b_high = ing["amp_B_high"]

    if variant == "baseline":
        c_low = st * a + ct * b_low
        c_high = ct * a - st * b_high
        g_low = ing["pref_low_f1"] * c_low
        g_high = ing["pref_high_f2"] * c_high
    elif variant == "radiative_JB_sign_flip":
        c_low = st * a - ct * b_low
        c_high = ct * a + st * b_high
        g_low = ing["pref_low_f1"] * c_low
        g_high = ing["pref_high_f2"] * c_high
    elif variant == "alternate_rotation_c_s":
        c_low = ct * a + st * b_low
        c_high = -st * a + ct * b_high
        g_low = ing["pref_low_f1"] * c_low
        g_high = ing["pref_high_f2"] * c_high
    elif variant == "alternate_rotation_c_minus_s":
        c_low = ct * a - st * b_low
        c_high = st * a + ct * b_high
        g_low = ing["pref_low_f1"] * c_low
        g_high = ing["pref_high_f2"] * c_high
    elif variant == "swap_eigenstate_assignment":
        c_low = ct * a - st * b_low
        c_high = st * a + ct * b_high
        g_low = ing["pref_low_f2"] * c_low
        g_high = ing["pref_high_f1"] * c_high
    elif variant == "swap_assignment_and_JB_flip":
        c_low = ct * a + st * b_low
        c_high = st * a - ct * b_high
        g_low = ing["pref_low_f2"] * c_low
        g_high = ing["pref_high_f1"] * c_high
    else:
        raise ValueError(variant)

    w_low = mod.width_vec_keV(g_low, inp.m_low, inp.m_bcstar)
    w_high = mod.width_vec_keV(g_high, inp.m_high, inp.m_bcstar)
    return {
        "g_low_GeV_m2": g_low,
        "g_high_GeV_m2": g_high,
        "Gamma_low_keV": w_low,
        "Gamma_high_keV": w_high,
        "ratio_high_over_low": w_high / w_low if w_low else float("inf"),
    }


def summarize(values):
    vals = sorted(values)
    n = len(vals)
    return {
        "median": vals[n // 2],
        "min": vals[0],
        "max": vals[-1],
    }


def main():
    mod = load_vector_module()
    inp = mod.Inputs()
    fstar = standard_fbcstar(mod, inp)
    variants = [
        "baseline",
        "radiative_JB_sign_flip",
        "alternate_rotation_c_s",
        "alternate_rotation_c_minus_s",
        "swap_eigenstate_assignment",
        "swap_assignment_and_JB_flip",
    ]

    detail_rows = []
    for M2 in [7.0, 8.0, 9.0]:
        for s0 in [53.0, 54.0, 55.0]:
            ing = row_ingredients(mod, inp, M2, s0, fstar)
            for variant in variants:
                detail_rows.append(
                    {
                        "variant": variant,
                        **{k: ing[k] for k in ["M2", "s0", "theta_deg", "f1_GeV", "f2_GeV", "amp_A", "amp_B_low", "amp_B_high"]},
                        **variant_widths(mod, inp, ing, variant),
                    }
                )

    summary_rows = []
    for variant in variants:
        rows = [r for r in detail_rows if r["variant"] == variant]
        for key in ["Gamma_low_keV", "Gamma_high_keV", "ratio_high_over_low"]:
            s = summarize([r[key] for r in rows])
            summary_rows.append({"variant": variant, "observable": key, **s})

    literature = [
        ("Godfrey", 60.0, 11.0),
        ("Ebert-Faustov-Galkin", 78.9, 13.6),
        ("Fulcher", 75.6, 26.2),
        ("Gershtein et al.", 77.8, 8.1),
        ("T. Y. Li et al.", 49.0, 10.5),
        ("Q. Li et al.", 70.0, 40.0),
        ("X. J. Li et al.", 47.8, 25.6),
        ("Eichten-Quigg", 62.5, 7.5),
        ("Bondar-Milstein", 75.0, 2.0),
    ]
    lit_rows = [
        {
            "publication": name,
            "Gamma_low_keV": low,
            "Gamma_high_keV": high,
            "ratio_high_over_low": high / low,
        }
        for name, low, high in literature
    ]

    central = row_ingredients(mod, inp, 8.0, 54.0, fstar)
    st = central["sin_theta"]
    ct = central["cos_theta"]
    central_components = {
        "baseline_low_A_piece": st * central["amp_A"],
        "baseline_low_B_piece": ct * central["amp_B_low"],
        "baseline_high_A_piece": ct * central["amp_A"],
        "baseline_high_B_piece": -st * central["amp_B_high"],
        "baseline_low_sum": st * central["amp_A"] + ct * central["amp_B_low"],
        "baseline_high_sum": ct * central["amp_A"] - st * central["amp_B_high"],
    }
    central_rows = [
        {"variant": variant, **variant_widths(mod, inp, central, variant)}
        for variant in variants
    ]

    for path, rows in [
        (OUT / "bc_vec_mixing_sign_audit_grid.csv", detail_rows),
        (OUT / "bc_vec_mixing_sign_audit_summary.csv", summary_rows),
        (OUT / "bc_vec_mixing_sign_audit_literature.csv", lit_rows),
        (OUT / "bc_vec_mixing_sign_audit_central.csv", central_rows),
    ]:
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    with (OUT / "bc_vec_mixing_sign_audit_components.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(central_components.keys()))
        writer.writeheader()
        writer.writerow(central_components)

    lines = [
        "Bc1 -> Bc* gamma mixing/sign audit",
        "====================================",
        "",
        "Central point: M2=8 GeV^2, s0=54 GeV^2, standard f_Bc* normalization.",
        f"theta = {central['theta_deg']:.6g} deg, f1 = {central['f1_GeV']:.6g} GeV, f2 = {central['f2_GeV']:.6g} GeV",
        f"A = {central['amp_A']:.6g}, B_low = {central['amp_B_low']:.6g}, B_high = {central['amp_B_high']:.6g}",
        "",
        "Baseline mixed-amplitude decomposition:",
        f"low:  sin(theta) A = {central_components['baseline_low_A_piece']:.6g}, "
        f"cos(theta) B = {central_components['baseline_low_B_piece']:.6g}, "
        f"sum = {central_components['baseline_low_sum']:.6g}",
        f"high: cos(theta) A = {central_components['baseline_high_A_piece']:.6g}, "
        f"-sin(theta) B = {central_components['baseline_high_B_piece']:.6g}, "
        f"sum = {central_components['baseline_high_sum']:.6g}",
        "The lower-state amplitude is cancellation-dominated; the higher-state amplitude is constructive.",
        "",
        "Central variant widths:",
    ]
    for row in central_rows:
        lines.append(
            "{variant}: Gamma_low={Gamma_low_keV:.4g} keV, "
            "Gamma_high={Gamma_high_keV:.4g} keV, high/low={ratio_high_over_low:.4g}".format(**row)
        )
    lines.extend(["", "Nine-point grid medians:"])
    for variant in variants:
        rows = [r for r in summary_rows if r["variant"] == variant]
        lookup = {r["observable"]: r for r in rows}
        lines.append(
            f"{variant}: low={lookup['Gamma_low_keV']['median']:.4g} "
            f"[{lookup['Gamma_low_keV']['min']:.4g},{lookup['Gamma_low_keV']['max']:.4g}] keV, "
            f"high={lookup['Gamma_high_keV']['median']:.4g} "
            f"[{lookup['Gamma_high_keV']['min']:.4g},{lookup['Gamma_high_keV']['max']:.4g}] keV, "
            f"ratio={lookup['ratio_high_over_low']['median']:.4g}"
        )
    lines.extend(["", "Literature high/low ratios for Bc* gamma:"])
    for row in lit_rows:
        lines.append(
            f"{row['publication']}: low={row['Gamma_low_keV']:.4g}, "
            f"high={row['Gamma_high_keV']:.4g}, ratio={row['ratio_high_over_low']:.4g}"
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "- Baseline reproduces the current QCDSR result and has high/low > 1.",
            "- All representative model ratios are < 1.",
            "- A radiative J_B sign flip gives the literature-like hierarchy but is not supported by the tensor-current convention shared with the Bc-mixing code.",
            "- Therefore the vector result should remain flagged until the independent vector tensor-current trace/sign convention is rederived.",
        ]
    )
    (OUT / "bc_vec_mixing_sign_audit.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
