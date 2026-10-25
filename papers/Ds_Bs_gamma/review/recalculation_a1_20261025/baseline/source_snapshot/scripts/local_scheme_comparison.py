"""Numerical and conceptual comparison of transition-condensate schemes.

The comparison is deliberately made at the axial-current OPE level, where
both prescriptions are defined.  Colangelo et al. did not provide the mixed
tensor-current J_B local term, so no such number is invented here.  The full
mixed-current paper prediction uses the Rohrwild nonlocal scheme only.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from bs1_stage2_baseline import bs1_inputs
from rohrwild_em_da import P_T4_BRACE, SUPPORT_NOTE, tensor_em_qcd
from stage1_axial_g1_baseline import central_inputs
from stage2_axial_g1_three_particle import F1_integral, g1_stage2


CONCEPTUAL_ROWS = [
    {
        "item": "transition correlator ordinary local condensate",
        "rohrwild_nonlocal": "absent",
        "colangelo_local_benchmark": "heavy-charge term retained",
        "paper_choice": "Rohrwild: absent",
    },
    {
        "item": "transition correlator S_gamma and T4_gamma",
        "rohrwild_nonlocal": "included as nonlocal photon DAs",
        "colangelo_local_benchmark": "not included separately",
        "paper_choice": "Rohrwild: included",
    },
    {
        "item": "background-gluon S,S-tilde,T1...T4",
        "rohrwild_nonlocal": "included as nonlocal photon DAs",
        "colangelo_local_benchmark": "included through F1",
        "paper_choice": "included",
    },
    {
        "item": "light-quark mass",
        "rohrwild_nonlocal": "retained in our extension",
        "colangelo_local_benchmark": "retained in benchmark formula",
        "paper_choice": "retained",
    },
    {
        "item": "two-point decay-constant sum rules",
        "rohrwild_nonlocal": "local condensates allowed and required",
        "colangelo_local_benchmark": "local condensates allowed and required",
        "paper_choice": "separate two-point OPE",
    },
    {
        "item": "local plus S_gamma/T4_gamma in one transition result",
        "rohrwild_nonlocal": "forbidden double counting",
        "colangelo_local_benchmark": "forbidden double counting",
        "paper_choice": "never added",
    },
]


def comparison_row(label, scenario, inputs, M2, s0_root, f1_total):
    s0 = s0_root * s0_root
    bare = g1_stage2(
        M2,
        s0,
        inputs,
        f1_total,
        transition_scheme="rohrwild_local_free_diagnostic",
    )
    col = g1_stage2(
        M2,
        s0,
        inputs,
        f1_total,
        transition_scheme="colangelo_local_benchmark",
    )
    roh = g1_stage2(
        M2,
        s0,
        inputs,
        f1_total,
        transition_scheme="rohrwild_nonlocal",
    )
    tem = tensor_em_qcd(M2, s0, inputs)
    local = float(col["heavy_local"])
    em = float(roh["em_delta_qcd"])
    return {
        "sector": label,
        "scenario": scenario,
        "M2_GeV2": M2,
        "sqrt_s0_GeV": s0_root,
        "colangelo_local_qcd_GeV3": local,
        "rohrwild_em_qcd_GeV3": em,
        "em_over_local": em / local if local else math.nan,
        "relative_qcd_difference_percent": 100.0 * (em - local) / local if local else math.nan,
        "gA_without_either_GeV_inv": bare["g1_stage2_GeV_inv"],
        "gA_colangelo_local_GeV_inv": col["g1_stage2_GeV_inv"],
        "gA_rohrwild_nonlocal_GeV_inv": roh["g1_stage2_GeV_inv"],
        "GammaA_colangelo_local_keV": col["width_stage2_keV"],
        "GammaA_rohrwild_nonlocal_keV": roh["width_stage2_keV"],
        "JB_rohrwild_em_qcd_GeV3": tem["em_qcd"],
        "JB_colangelo_local_status": "not_given_for_tensor_current",
        "P_T4gamma_brace": P_T4_BRACE,
    }


def scan_summary(label, scenario, inputs, M2_values, s0_roots, f1_total):
    ratios = []
    dg = []
    dw = []
    for M2 in M2_values:
        for s0_root in s0_roots:
            row = comparison_row(label, scenario, inputs, float(M2), float(s0_root), f1_total)
            ratios.append(row["em_over_local"])
            dg.append(
                row["gA_rohrwild_nonlocal_GeV_inv"]
                - row["gA_colangelo_local_GeV_inv"]
            )
            dw.append(
                row["GammaA_rohrwild_nonlocal_keV"]
                - row["GammaA_colangelo_local_keV"]
            )
    return {
        "sector": label,
        "scenario": scenario,
        "em_over_local_min": min(ratios),
        "em_over_local_max": max(ratios),
        "delta_gA_min_GeV_inv": min(dg),
        "delta_gA_max_GeV_inv": max(dg),
        "delta_GammaA_min_keV": min(dw),
        "delta_GammaA_max_keV": max(dw),
    }


def write_latex(rows):
    lines = [
        r"\begin{table}",
        r"\centering",
        r"\caption{Axial-current comparison of the alternative transition-OPE prescriptions. The entries in the two OPE columns are alternatives and are not added.}",
        r"\begin{tabular}{llrrrr}",
        r"\hline",
        r"Sector & input & $\Pi_{\rm loc}$ & $\Pi_{S_\gamma,T_4^\gamma}$ & $g_A^{\rm loc}$ & $g_A^{\rm RW}$ \\",
        r" & & \multicolumn{2}{c}{[GeV$^3$]} & \multicolumn{2}{c}{[GeV$^{-1}$]} \\",
        r"\hline",
    ]
    for r in rows:
        scenario_tex = r["scenario"].replace("_", r"\_")
        lines.append(
            f"{r['sector']} & {scenario_tex} & "
            f"{r['colangelo_local_qcd_GeV3']:+.5e} & "
            f"{r['rohrwild_em_qcd_GeV3']:+.5e} & "
            f"{r['gA_colangelo_local_GeV_inv']:+.4f} & "
            f"{r['gA_rohrwild_nonlocal_GeV_inv']:+.4f} \\\\"
        )
    lines.extend(
        [
            r"\hline",
            r"\end{tabular}",
            r"\label{tab:local-vs-rohrwild}",
            r"\end{table}",
        ]
    )
    path = OUT / "local_scheme_comparison_table.tex"
    path.write_text("\n".join(lines) + "\n")
    return path


def main():
    f1_total, _, _ = F1_integral(u0=0.5)
    cases = [
        ("Ds", "reference_point", central_inputs(), 4.5, 2.55),
        ("Ds", "selected_window_midpoint", central_inputs(), 3.75, math.sqrt(8.0)),
        (
            "Bs",
            "legacy_chi_selected_midpoint",
            {k: v for k, v in bs1_inputs(5.82870, "legacy_chi_condensate").items() if k not in {"fT", "fperp_s_used"}},
            12.0,
            math.sqrt(39.5),
        ),
        (
            "Bs",
            "lattice_fperp_selected_midpoint",
            {k: v for k, v in bs1_inputs(5.82870, "lattice_fperp_s").items() if k not in {"fT", "fperp_s_used"}},
            12.0,
            math.sqrt(39.5),
        ),
    ]
    rows = [comparison_row(*case, f1_total) for case in cases]

    conceptual_path = OUT / "local_scheme_conceptual_comparison.csv"
    with conceptual_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(CONCEPTUAL_ROWS[0]))
        writer.writeheader()
        writer.writerows(CONCEPTUAL_ROWS)

    numerical_path = OUT / "local_scheme_numerical_comparison.csv"
    with numerical_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summaries = [
        scan_summary(
            "Ds",
            "reference_window",
            central_inputs(),
            np.linspace(3.0, 6.0, 13),
            (2.50, 2.55, 2.60),
            f1_total,
        ),
        scan_summary(
            "Ds",
            "selected_window",
            central_inputs(),
            np.linspace(3.0, 4.5, 13),
            tuple(math.sqrt(x) for x in (7.5, 8.0, 8.5)),
            f1_total,
        ),
        scan_summary(
            "Bs",
            "legacy_chi_selected_window",
            cases[2][2],
            np.linspace(10.0, 14.0, 13),
            tuple(math.sqrt(x) for x in (38.0, 39.5, 41.0)),
            f1_total,
        ),
        scan_summary(
            "Bs",
            "lattice_fperp_selected_window",
            cases[3][2],
            np.linspace(10.0, 14.0, 13),
            tuple(math.sqrt(x) for x in (38.0, 39.5, 41.0)),
            f1_total,
        ),
    ]
    summary_csv = OUT / "local_scheme_scan_summary.csv"
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    tex_path = write_latex(rows)
    lines = [
        "Local contribution comparison: Rohrwild vs Colangelo",
        "=====================================================",
        "The two transition prescriptions are alternatives; they are never added.",
        f"P[T4gamma] brace = {P_T4_BRACE:+.12f}",
        f"Support qualification: {SUPPORT_NOTE}",
        "",
    ]
    for r in rows:
        lines.append(
            f"{r['sector']} ({r['scenario']}): local={r['colangelo_local_qcd_GeV3']:+.6e}, "
            f"nonlocal EM={r['rohrwild_em_qcd_GeV3']:+.6e} GeV^3, "
            f"ratio={r['em_over_local']:+.4f}; "
            f"gA {r['gA_colangelo_local_GeV_inv']:+.5f} -> "
            f"{r['gA_rohrwild_nonlocal_GeV_inv']:+.5f} GeV^-1"
        )
    lines.extend(
        [
            "",
            "No J_B Colangelo-local value is quoted: that tensor-current local OPE was not",
            "given in the reference and is not needed in the paper's Rohrwild scheme.",
            f"Wrote {conceptual_path}",
            f"Wrote {numerical_path}",
            f"Wrote {summary_csv}",
            f"Wrote {tex_path}",
        ]
    )
    summary_path = OUT / "local_scheme_comparison_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
