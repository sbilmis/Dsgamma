"""Pure-current/no-mixing comparison tables.

This is not an angle fit.  It compares the basis-current LCSR result obtained
with Colangelo-like inputs to the corresponding published LCSR benchmark.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
PSEUDOSCALAR_OUT = ROOT.parent / "Ds_Bs_gamma" / "outputs"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def value_range(rows: list[dict[str, str]], key: str) -> tuple[float, float]:
    vals = np.asarray([float(r[key]) for r in rows], dtype=float)
    return float(np.min(vals)), float(np.max(vals))


def fmt_interval(lo: float, hi: float, ndigits: int = 3) -> str:
    return f"[{lo:.{ndigits}f},{hi:.{ndigits}f}]"


def main() -> None:
    ds_rows = [
        r for r in read_csv(PSEUDOSCALAR_OUT / "stage1_axial_g1_baseline.csv")
        if r["input_set"] == "colangelo_inputs_chi_positive_stage1"
    ]
    dsstar_rows = [
        r for r in read_csv(ROOT / "Dsstar_gamma" / "outputs" / "stage1_axial_colangelo_gstar.csv")
        if r["input_set"] == "colangelo_like_original_window_positive_chi"
    ]

    ds_g = value_range(ds_rows, "g1_GeV_inv")
    ds_w = value_range(ds_rows, "width_keV")
    dsstar_g = value_range(dsstar_rows, "gA_star")
    dsstar_w = value_range(dsstar_rows, "width_A_keV")

    rows = [
        {
            "channel": "D_s1(2460) -> D_s gamma",
            "current_limit": "pure axial J_A, no J_B mixing",
            "benchmark": "Colangelo et al. LCSR",
            "benchmark_coupling": "[-0.37,-0.29] GeV^-1",
            "our_coupling": fmt_interval(ds_g[0], ds_g[1]) + " GeV^-1",
            "benchmark_width": "[19,29] keV",
            "our_width": fmt_interval(ds_w[0], ds_w[1], 2) + " keV",
            "comment": "same-method scale reproduced; residual shift from inputs/window and omitted stage-2 term in this baseline",
        },
        {
            "channel": "D_s1 -> D_s* gamma",
            "current_limit": "pure axial J_A, no J_B mixing",
            "benchmark": "Colangelo et al. LCSR",
            "benchmark_coupling": "[-0.18,-0.13]",
            "our_coupling": fmt_interval(dsstar_g[0], dsstar_g[1]) ,
            "benchmark_width": "not emphasized in comparison table",
            "our_width": fmt_interval(dsstar_w[0], dsstar_w[1], 3) + " keV",
            "comment": "direct normalization/sign benchmark; excellent agreement for g_A*",
        },
    ]

    csv_path = OUT / "no_mixing_method_comparison.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    tex_path = OUT / "no_mixing_method_comparison_table.tex"
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\caption{No-mixing/pure-current benchmark against same-method LCSR results.}",
        r"\label{tab:no-mixing-method-check}",
        r"\begin{tabular}{lllll}",
        r"\toprule",
        r"Channel & Current limit & Literature coupling & This work coupling & Width check\\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            f"{row['channel']} & {row['current_limit']} & {row['benchmark_coupling']} & "
            f"{row['our_coupling']} & {row['benchmark_width']} vs. {row['our_width']}\\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    tex_path.write_text("\n".join(lines))

    note_path = ROOT / "notes" / "no_mixing_method_comparison.md"
    note_path.write_text(
        "# No-Mixing Same-Method Comparison\n\n"
        "This is the intended benchmark check: switch off physical mixing and compare the pure basis-current LCSR result with literature calculations using the same current.\n\n"
        "With our convention\n\n"
        "\\[\nJ_{1\\mu}=\\sin\\theta\\,J^A_\\mu+\\cos\\theta\\,J^B_\\mu,\n\\]\n\n"
        "the pure axial-current comparison corresponds to the basis quantity `g_A`, equivalently the lower-state current at `theta=90 deg` before adding tensor-current mixing.\n\n"
        "## Results\n\n"
        f"- `D_s1(2460) -> D_s gamma`: Colangelo-like pure axial baseline gives "
        f"`g_1={fmt_interval(ds_g[0], ds_g[1])} GeV^-1` and "
        f"`Gamma={fmt_interval(ds_w[0], ds_w[1], 2)} keV`, compared with the published "
        "`g_1=[-0.37,-0.29] GeV^-1` and `Gamma=[19,29] keV`.\n"
        f"- `D_s1 -> D_s* gamma`: Colangelo-like pure axial baseline gives "
        f"`g_A*={fmt_interval(dsstar_g[0], dsstar_g[1])}`, compared with the published "
        "`g*=[-0.18,-0.13]`.  This is a direct sign/normalization agreement.\n\n"
        "The comparison table is stored in:\n\n"
        f"- `{csv_path}`\n"
        f"- `{tex_path}`\n"
    )

    print(f"Wrote {csv_path}")
    print(f"Wrote {tex_path}")
    print(f"Wrote {note_path}")


if __name__ == "__main__":
    main()
