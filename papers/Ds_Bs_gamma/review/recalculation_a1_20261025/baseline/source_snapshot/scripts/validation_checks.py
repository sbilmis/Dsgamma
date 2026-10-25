"""Validation and limiting-case checks for the Ds1 -> Ds gamma LCSR setup."""

from __future__ import annotations

from contextlib import contextmanager
import csv
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT.parents[1] / "shared"))
sys.path.insert(0, str(ROOT / "scripts"))

import photon_da as pda
from stage1_axial_g1_baseline import (
    central_inputs,
    colangelo_inputs,
    gauss_legendre_integral,
    g1_stage1,
)
from stage1_tensor_gb_hard_candidate import hard_tensor_gb
from stage1_tensor_gb_soft2p_estimate import gb_soft_two_particle, width_keV
from stage2_axial_g1_three_particle import F1_integral, g1_stage2
from stage2_tensor_gb_three_particle_estimate import gb_em_da, gb_three_particle_from_integral
from final_stage2_uncertainty_scan import matched_f1_integral, precompute_f1_basis
from rohrwild_em_da import P_T4_BRACE, self_check as em_da_self_check


@contextmanager
def photon_params(**updates):
    old = {k: pda.PARAMS[k] for k in updates}
    try:
        for key, value in updates.items():
            err = pda.PARAMS[key][1]
            pda.PARAMS[key] = (value, err)
        yield
    finally:
        for key, value in old.items():
            pda.PARAMS[key] = value


def summarize_range(values):
    arr = np.array(values, dtype=float)
    return arr.min(), arr.max()


def colangelo_axial_benchmark():
    """Reproduce the Colangelo axial-current LCSR band with their inputs."""
    rows = []
    with photon_params(kappa=0.2, kappa_p=0.0, zeta1=0.4, zeta1_p=0.0, zeta2=0.3, zeta2_p=0.0):
        inputs = colangelo_inputs()
        f1_total, _, _ = F1_integral(u0=0.5)
        for s0_root in (2.50, 2.55, 2.60):
            s0 = s0_root * s0_root
            for M2 in np.linspace(3.0, 6.0, 13):
                row = g1_stage2(
                    float(M2),
                    s0,
                    inputs,
                    f1_total,
                    transition_scheme="colangelo_local_benchmark",
                )
                row["s0_root"] = s0_root
                rows.append(row)
    return rows


def da_normalization_checks():
    phi_int = gauss_legendre_integral(pda.phi_gamma, 0.0, 1.0, n=300)
    psi_a_sym = max(abs(float(pda.psi_a(u) - pda.psi_a(1.0 - u))) for u in np.linspace(0, 1, 51))
    psi_v_anti = max(abs(float(pda.psi_v(u) + pda.psi_v(1.0 - u))) for u in np.linspace(0, 1, 51))
    f1_total, f1_d1, f1_d2 = F1_integral(u0=0.5)
    return {
        "int_phi_gamma": phi_int,
        "psi_a_symmetry_max_abs": psi_a_sym,
        "psi_v_antisymmetry_max_abs": psi_v_anti,
        "F1_total": f1_total,
        "F1_domain1": f1_d1,
        "F1_domain2": f1_d2,
        "F1_closure": f1_total - f1_d1 - f1_d2,
    }


def final_stage2_point(inputs, M2=4.5, s0_root=2.55, theta_deg=35.3):
    f1_axial, _, _ = F1_integral(u0=0.5)
    f1_tensor = matched_f1_integral(inputs, precompute_f1_basis())
    fT = 0.256
    fB = fT * inputs["m_ds1"] / (inputs["mc"] + inputs["ms"])
    s0 = s0_root * s0_root
    theta = math.radians(theta_deg)
    axial = g1_stage2(M2, s0, inputs, f1_axial)
    gb_hard, _ = hard_tensor_gb(M2, s0, inputs, fB)
    gb_soft = gb_soft_two_particle(axial, inputs, fB)
    gb_3p = gb_three_particle_from_integral(axial, inputs, fB, f1_tensor)
    gb_em, _ = gb_em_da(axial, inputs, fB)
    GA = axial["g1_stage2_GeV_inv"]
    GB = gb_hard + gb_soft + gb_3p + gb_em
    G2460 = math.sin(theta) * GA + math.cos(theta) * GB
    G2536 = math.cos(theta) * GA - math.sin(theta) * GB
    return {
        "GA": GA,
        "GB": GB,
        "GB_hard": gb_hard,
        "GB_soft2p": gb_soft,
        "GB_3p": gb_3p,
        "GB_em": gb_em,
        "G2460": G2460,
        "G2536": G2536,
        "Gamma2460_keV": width_keV(inputs["m_ds1"], inputs["m_ds"], G2460),
        "Gamma2536_keV": width_keV(2.53511, inputs["m_ds"], G2536),
    }


def limiting_case_checks():
    base = central_inputs()
    mq_eps = 1.0e-4
    cases = []

    def add_case(label, edits):
        inputs = dict(base)
        inputs.update(edits)
        # Keep thresholds above partonic threshold in unphysical toy limits.
        cases.append({"case": label, **final_stage2_point(inputs)})

    add_case("central", {})
    add_case("m_s_to_0_eps", {"ms": mq_eps})
    add_case("equal_parton_masses_ms_eq_mc_toy", {"ms": base["mc"]})
    add_case("heavy_charge_off_toy_ec_0", {"ec": 0.0})
    add_case("strange_charge_off_toy_es_0", {"es": 0.0})

    # Light-quark counterpart diagnostics: OPE-level replacements only.  These
    # are not physical D1 predictions because the nonstrange hadron inputs are
    # not yet in inputs_table.csv.  The tiny mass avoids a 0/0 representation
    # in the current hard strange-line spectral density; the exact massless
    # limit should be taken analytically if these rows become paper results.
    add_case("light_u_toy_mq_eps_eu", {"ms": mq_eps, "es": 2.0 / 3.0, "ss": -(0.240) ** 3})
    add_case("light_d_toy_mq_eps_ed", {"ms": mq_eps, "es": -1.0 / 3.0, "ss": -(0.240) ** 3})
    return cases


def hard_charge_decomposition():
    base = central_inputs()
    full = final_stage2_point(base)
    no_heavy = final_stage2_point({**base, "ec": 0.0})
    no_strange = final_stage2_point({**base, "es": 0.0})
    return {
        "GB_hard_full": full["GB_hard"],
        "GB_hard_no_heavy_charge": no_heavy["GB_hard"],
        "GB_hard_no_strange_charge": no_strange["GB_hard"],
        "GB_hard_heavy_line_estimate": full["GB_hard"] - no_heavy["GB_hard"],
        "GB_hard_strange_line_estimate": full["GB_hard"] - no_strange["GB_hard"],
    }


def transition_scheme_checks():
    inputs = central_inputs()
    f1, _, _ = F1_integral(u0=0.5)
    args = (4.5, 2.55**2, inputs, f1)
    rw = g1_stage2(*args, transition_scheme="rohrwild_nonlocal")
    col = g1_stage2(*args, transition_scheme="colangelo_local_benchmark")
    bare = g1_stage2(*args, transition_scheme="rohrwild_local_free_diagnostic")
    em_da_self_check()
    assert rw["heavy_local"] == 0.0 and rw["em_delta_qcd"] != 0.0
    assert col["heavy_local"] != 0.0 and col["em_delta_qcd"] == 0.0
    assert bare["heavy_local"] == 0.0 and bare["em_delta_qcd"] == 0.0
    return {
        "local_qcd": col["heavy_local"],
        "nonlocal_em_qcd": rw["em_delta_qcd"],
        "em_over_local": rw["em_delta_qcd"] / col["heavy_local"],
        "delta_gA_rw_minus_col": rw["g1_stage2_GeV_inv"] - col["g1_stage2_GeV_inv"],
        "P_T4gamma_brace": P_T4_BRACE,
    }


def main():
    col_rows = colangelo_axial_benchmark()
    col_g = [r["g1_stage2_GeV_inv"] for r in col_rows]
    col_w = [r["width_stage2_keV"] for r in col_rows]
    da = da_normalization_checks()
    limits = limiting_case_checks()
    charge = hard_charge_decomposition()
    schemes = transition_scheme_checks()

    with (OUT / "validation_limiting_cases.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(limits[0].keys()))
        writer.writeheader()
        writer.writerows(limits)

    lines = [
        "Validation checks",
        "=================",
        "Colangelo axial benchmark with paper-like inputs and BBK kappa=0.2:",
        f"  g1 range = {min(col_g):+.4f} to {max(col_g):+.4f} GeV^-1",
        f"  width range = {min(col_w):.2f} to {max(col_w):.2f} keV",
        "  target from hep-ph/0505195: g1 = -0.37 to -0.29 GeV^-1, width = 19 to 29 keV",
        "",
        "Photon DA checks:",
        f"  int phi_gamma du = {da['int_phi_gamma']:.10f}",
        f"  max psi_a symmetry residual = {da['psi_a_symmetry_max_abs']:.3e}",
        f"  max psi_v antisymmetry residual = {da['psi_v_antisymmetry_max_abs']:.3e}",
        f"  F1 total = {da['F1_total']:+.10e}",
        f"  F1 closure residual = {da['F1_closure']:+.3e}",
        "",
        "Mutually exclusive transition-scheme check at the central point:",
        f"  Colangelo local QCD term = {schemes['local_qcd']:+.6e} GeV^3",
        f"  Rohrwild nonlocal EM term = {schemes['nonlocal_em_qcd']:+.6e} GeV^3",
        f"  nonlocal/local ratio = {schemes['em_over_local']:+.6f}",
        f"  gA(RW)-gA(local) = {schemes['delta_gA_rw_minus_col']:+.6e} GeV^-1",
        f"  P[T4gamma] finite brace = {schemes['P_T4gamma_brace']:+.12f}",
        "",
        "Hard charge decomposition at M2=4.5 GeV^2, sqrt(s0)=2.55 GeV:",
    ]
    for key, value in charge.items():
        lines.append(f"  {key} = {value:+.6e}")
    lines.extend(
        [
            "",
            "Limiting-case diagnostics at M2=4.5 GeV^2, sqrt(s0)=2.55 GeV are in:",
            f"  {OUT / 'validation_limiting_cases.csv'}",
            "Light-u/light-d rows are OPE diagnostics only, not physical D1 predictions.",
        ]
    )
    path = OUT / "validation_checks_summary.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
