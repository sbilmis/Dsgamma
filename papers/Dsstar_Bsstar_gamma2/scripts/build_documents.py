#!/usr/bin/env python3
"""Generate audited LaTeX tables and the final numerical manifest.

The calculation scripts own every number.  This script only formats their CSV
outputs, so the manuscripts cannot silently drift away from the machine-
readable results.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "outputs" / "csv"
GEN = ROOT / "tex" / "generated"
MANIFEST = ROOT / "output" / "manifests"

STATE_LABELS = {
    "Ds1_2460": r"$D_{s1}(2460)$",
    "Ds1_2536": r"$D_{s1}(2536)$",
    "Bs1_lower": r"$B_{s1}^{L}$",
    "Bs1_5830": r"$B_{s1}(5830)$",
}

CHANNEL_LABELS = {
    "Ds1_2460": r"$D_{s1}(2460)^+\to D_s^{\ast+}\gamma$",
    "Ds1_2536": r"$D_{s1}(2536)^+\to D_s^{\ast+}\gamma$",
    "Bs1_lower": r"$B_{s1}^{L\,0}\to B_s^{\ast0}\gamma$ (predicted lower state)",
    "Bs1_5830": r"$B_{s1}(5830)^0\to B_s^{\ast0}\gamma$ (observed higher state)",
}

SOURCE_LABELS = {
    "Colangelo2005": r"Colangelo et al. (2005)~\cite{Colangelo2005}",
    "Godfrey2003": r"Godfrey (2003)~\cite{Godfrey2003}",
    "Bardeen2003": r"Bardeen et al. (2003)~\cite{Bardeen2003}",
    "Chen2020": r"Chen et al. (2020)~\cite{Chen2020}",
    "Radford2009": r"Radford et al. (2009)~\cite{Radford2009}",
    "Green2017": r"Green et al. (2017)~\cite{Green2017}",
    "Godfrey2005": r"Godfrey (2005)~\cite{Godfrey2005}",
    "CloseSwanson2005": r"Close and Swanson (2005)~\cite{CloseSwanson2005}",
    "GoityRoberts2001": r"Goity and Roberts (2001)~\cite{GoityRoberts2001}",
    "Korner1993": r"K{\"o}rner et al. (1993)~\cite{Korner1993}",
    "Wang2008": r"Wang (2008)~\cite{Wang2008}",
    "Fu2022": r"Fu et al. (2022)~\cite{Fu2022}",
    "Li2021": r"Li et al. (2021)~\cite{Li2021}",
    "BondarMilstein2025": (
        r"Bondar and Milstein (2025)~\cite{BondarMilstein2025}"
    ),
    "PDG2026": r"PDG (2026)~\cite{PDG2026}",
}


def fmt(x: float, digits: int = 3) -> str:
    if pd.isna(x):
        return "--"
    if x == 0:
        return "0"
    if abs(x) < 1.0e-3 or abs(x) >= 1.0e4:
        return rf"\num{{{x:.{digits}e}}}"
    return rf"\num{{{x:.{digits}f}}}"


def interval(center: float, lo: float, hi: float, digits: int = 3) -> str:
    return (
        rf"${center:.{digits}f}_{{-{center-lo:.{digits}f}}}"
        rf"^{{+{hi-center:.{digits}f}}}$"
    )


def write(name: str, body: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(body.rstrip() + "\n", encoding="utf-8")


def result_table(mc: pd.DataFrame) -> None:
    lines = [
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"State & $g$ & $\Gamma_\gamma$ [keV] & $f_i$ [GeV] & $f_V$ [GeV]\\",
        r"\midrule",
    ]
    for _, r in mc.iterrows():
        lines.append(
            "{} & {} & {} & {} & {} \\\\".format(
                STATE_LABELS[r.state],
                interval(r.g_median, r.g_p16, r.g_p84, 3),
                interval(r.width_median_keV, r.width_p16_keV, r.width_p84_keV, 3),
                interval(
                    r.f_initial_median_GeV,
                    r.f_initial_p16_GeV,
                    r.f_initial_p84_GeV,
                    3,
                ),
                interval(
                    r.f_vector_median_GeV,
                    r.f_vector_p16_GeV,
                    r.f_vector_p84_GeV,
                    3,
                ),
            )
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("results_table.tex", "\n".join(lines))


def branching_table(mc: pd.DataFrame) -> None:
    rows = mc[mc.branching_central.notna()]
    lines = [
        r"\begin{tabular}{lcc}",
        r"\toprule",
        r"State & central radiative branching fraction & experimental status\\",
        r"\midrule",
    ]
    status = {
        "Ds1_2536": r"possibly seen; no branching fraction",
        "Bs1_5830": r"mode not observed",
    }
    for _, r in rows.iterrows():
        lines.append(
            rf"{STATE_LABELS[r.state]} & "
            rf"\num{{{r.branching_central:.3e}}} & "
            rf"{status[r.state]} \\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("branching_table.tex", "\n".join(lines))


def windows_table(df: pd.DataFrame) -> None:
    lines = [
        r"\begin{tabular}{lcccccc}",
        r"\toprule",
        r"State & $M^2$ [GeV$^2$] & $s_0$ [GeV$^2$] & accepted & PC & "
        r"$|d=5|/|\Pi|$ & $m_{\rm SR}$ [GeV]\\",
        r"\midrule",
    ]
    for _, r in df.iterrows():
        lines.append(
            rf"{STATE_LABELS[r.state]} & {r.M2_min_GeV2:.1f}--{r.M2_max_GeV2:.1f}"
            rf" & {r.s0_min_GeV2:.2f}--{r.s0_max_GeV2:.2f}"
            rf" & {int(r.accepted_points)}/{int(r.M2_points*r.s0_points)}"
            rf" & {r.pc_min:.2f}--{r.pc_max:.2f}"
            rf" & {r.d5frac_max:.3f}"
            rf" & {r.mass_min_GeV:.3f}--{r.mass_max_GeV:.3f} \\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("windows_table.tex", "\n".join(lines))


def stability_table(df: pd.DataFrame) -> None:
    piv = df.pivot(index="state", columns=["observable", "scan"], values="relative_variation")
    lines = [
        r"\begin{tabular}{lrrrr}",
        r"\toprule",
        r"State & $\Delta_{M^2}g$ & $\Delta_{s_0}g$ & "
        r"$\Delta_{M^2}\Gamma$ & $\Delta_{s_0}\Gamma$\\",
        r"\midrule",
    ]
    for state in STATE_LABELS:
        lines.append(
            rf"{STATE_LABELS[state]} & "
            rf"{100*piv.loc[state, ('coupling','M2')]:.1f}\% & "
            rf"{100*piv.loc[state, ('coupling','s0')]:.1f}\% & "
            rf"{100*piv.loc[state, ('width','M2')]:.1f}\% & "
            rf"{100*piv.loc[state, ('width','s0')]:.1f}\% \\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("stability_table.tex", "\n".join(lines))


def transition_table(df: pd.DataFrame) -> None:
    piv = df.pivot_table(
        index="state", columns="term", values="value_GeV4", aggfunc="first"
    )
    order = [
        "hard",
        "tw2",
        "tw3_2p",
        "tw3_3p",
        "tw4_2p",
        "tw4_3p",
        "local_heavy_diagnostic",
    ]
    lines = [
        r"\begin{tabular}{lrrrrrrr}",
        r"\toprule",
        r"State & hard & tw2 & tw3(2p) & tw3(3p) & tw4(2p) & tw4(3p)"
        r" & local diag.\\",
        r"\midrule",
    ]
    for state in STATE_LABELS:
        lines.append(
            STATE_LABELS[state]
            + " & "
            + " & ".join(fmt(piv.loc[state, term], 3) for term in order)
            + r" \\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("transition_table.tex", "\n".join(lines))


def validation_table(df: pd.DataFrame) -> None:
    gate_tex = {
        "rho_AB_equals_rho_BA": r"$\rho_{AB}=\rho_{BA}$",
        "spectral_matrix_positive": "spectral positivity",
        "mass_dimensions": "mass dimensions",
        "physical_threshold": "physical threshold",
        "Ward_identity": "Ward identity",
        "no_local_transition_condensate": "no local transition condensate",
        "two_point_local_d3_d5": r"two-point $d=3,5$",
        "pure_A_limit": "pure A limit",
        "pure_B_limit": "pure B limit",
        "equal_Borel_limit": "equal Borel limit",
        "zero_ms_threshold_limit": r"$m_s\to0$ threshold",
        "charge_switch": "charge switch",
        "python_mathematica_regression": "Python--Mathematica regression",
        "subleading_tensor_kernels": "subleading tensor kernels",
        "alpha_s_transition": r"transition $\alpha_s$",
    }
    lines = [
        r"\begin{longtable}{p{0.29\linewidth}p{0.13\linewidth}p{0.50\linewidth}}",
        r"\toprule",
        r"Gate & status & evidence\\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Gate & status & evidence\\",
        r"\midrule",
        r"\endhead",
    ]
    for _, r in df.iterrows():
        evidence = (
            str(r.evidence)
            .replace("&", r"\&")
            .replace("%", r"\%")
            .replace("_", r"\_")
            .replace("^", r"\textasciicircum{}")
            .replace("#", r"\#")
            .replace("<", r"$<$")
            .replace(">", r"$>$")
        )
        lines.append(
            rf"{gate_tex.get(r.gate, r.gate.replace('_', r'\_'))} & "
            rf"{r.status} & {evidence}\\"
        )
    lines += [r"\bottomrule", r"\end{longtable}"]
    write("validation_table.tex", "\n".join(lines))


def regression_table(df: pd.DataFrame) -> None:
    selected = df[
        df.quantity.str.contains(
            r"rho_|T_quoted|_g$|width_keV$", regex=True
        )
    ].copy()
    lines = [
        r"\begin{longtable}{llrrr}",
        r"\toprule",
        r"sector & quantity & Python & Mathematica & relative difference\\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"sector & quantity & Python & Mathematica & relative difference\\",
        r"\midrule",
        r"\endhead",
    ]
    for _, r in selected.iterrows():
        quantity = str(r.quantity).replace("_", r"\_")
        lines.append(
            rf"{r.sector} & {quantity} & {fmt(r.value_python, 6)} & "
            rf"{fmt(r.value_mathematica, 6)} & {fmt(r.relative_difference, 2)}\\"
        )
    lines += [r"\bottomrule", r"\end{longtable}"]
    write("regression_table.tex", "\n".join(lines))


def literature_width(value: object) -> str:
    text = str(value)
    if "+/-" in text:
        left, right = text.split("+/-", maxsplit=1)
        return rf"${left}\pm {right}$"
    if "-" in text:
        left, right = text.split("-", maxsplit=1)
        return rf"${left}\,\text{{--}}\,{right}$"
    return f"${text}$"


def comparison_table(
    mc: pd.DataFrame, literature: pd.DataFrame, experiment: pd.DataFrame
) -> None:
    mc_by_state = mc.set_index("state")
    experiment_by_state = experiment.set_index("state")
    experimental_result = {
        "Ds1_2460": (
            r"$\mathcal B(D_s^\ast\gamma)<8\%$ (90\% CL); "
            r"$\Gamma_{\rm tot}<3.5$ MeV"
        ),
        "Ds1_2536": (
            r"mode possibly seen; no $\mathcal B$ or partial width; "
            r"$\Gamma_{\rm tot}=0.92\pm0.05$ MeV"
        ),
        "Bs1_lower": r"state unobserved; no measurement or limit",
        "Bs1_5830": r"radiative mode unobserved; no measurement or limit",
    }
    lines = [
        r"{\footnotesize",
        r"\setlength{\tabcolsep}{4pt}",
        r"\renewcommand{\arraystretch}{1.10}",
        r"\begin{longtable}{>{\raggedright\arraybackslash}p{0.13\linewidth}"
        r">{\raggedright\arraybackslash}p{0.23\linewidth}"
        r">{\raggedright\arraybackslash}p{0.27\linewidth}"
        r">{\raggedright\arraybackslash}p{0.28\linewidth}}",
        (
            r"\caption{Radiative-width predictions and current experimental "
            r"status.  Width entries are in keV.  ``This work'' gives the "
            r"median and 16th--84th percentile interval; experimental rows "
            r"report direct branching-fraction, total-width, or observation "
            r"information and must not be read as radiative widths.}"
            r"\label{tab:literature-full}\\"
        ),
        r"\toprule",
        r"Entry & width [keV] or experimental status & method/status & source and note\\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{4}{l}{\footnotesize\itshape Table "
        r"\thetable\ continued from the previous page}\\",
        r"\toprule",
        r"Entry & width [keV] or experimental status & method/status & source and note\\",
        r"\midrule",
        r"\endhead",
        r"\midrule",
        r"\multicolumn{4}{r}{\footnotesize Continued on the next page}\\",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
    ]
    for state in STATE_LABELS:
        current = mc_by_state.loc[state]
        experiment_row = experiment_by_state.loc[state]
        if state == "Bs1_lower":
            lines.append(r"\pagebreak")
        lines.append(r"\addlinespace[3pt]")
        lines.append(
            rf"\multicolumn{{4}}{{l}}{{\textbf{{{CHANNEL_LABELS[state]}}}}}\\"
        )
        lines.append(
            rf"\textbf{{This work}} & "
            rf"\textbf{{{interval(current.width_median_keV, current.width_p16_keV, current.width_p84_keV, 3)}}} "
            rf"& mixed-current external-photon LCSR "
            rf"& this calculation\\"
        )
        for _, row in literature[literature.state.eq(state)].iterrows():
            note = str(row.comment).replace("&", r"\&").replace("%", r"\%")
            note = (
                note.replace("Bs1(5748)", r"$B_{s1}(5748)$")
                .replace("Bs1(5829)", r"$B_{s1}(5829)$")
            )
            if note.startswith("comparison value collected"):
                note = ""
            method = str(row.method).replace("&", r"\&")
            source = SOURCE_LABELS[str(row.reference)]
            source_note = source + (f"; {note}" if note else "")
            lines.append(
                rf"Published & {literature_width(row.literature_value)} "
                rf"& {method} & {source_note}\\"
            )
        source = SOURCE_LABELS[str(experiment_row.source)]
        lines.append(
            rf"\textbf{{Experiment}} & \textbf{{{experimental_result[state]}}} "
            rf"& {str(experiment_row.status).replace('-', '--')} "
            rf"& {source}\\"
        )
    lines += [r"\end{longtable}", r"}"]
    write("comparison_table.tex", "\n".join(lines))


def numerical_macros(mc: pd.DataFrame, regression: pd.DataFrame) -> None:
    by = mc.set_index("state")
    material = regression[
        regression[["value_python", "value_mathematica"]]
        .abs()
        .max(axis=1)
        .gt(2.0e-8)
    ]
    macros = [
        rf"\newcommand{{\RegressionPass}}{{{int(regression['pass'].sum())}/{len(regression)}}}",
        rf"\newcommand{{\RegressionMaxRel}}{{{material.relative_difference.max():.2e}}}",
        rf"\newcommand{{\RegressionMaxAbs}}{{{regression.absolute_difference.max():.2e}}}",
        r"\newcommand{\NMonteCarlo}{1200}",
    ]
    short = {
        "Ds1_2460": "DsLow",
        "Ds1_2536": "DsHigh",
        "Bs1_lower": "BsLow",
        "Bs1_5830": "BsHigh",
    }
    for state, key in short.items():
        r = by.loc[state]
        macros += [
            rf"\newcommand{{\{key}G}}{{{r.g_median:.3f}}}",
            rf"\newcommand{{\{key}Width}}{{{r.width_median_keV:.3f}}}",
        ]
    write("numerical_macros.tex", "\n".join(macros))


def final_summary(mc: pd.DataFrame) -> None:
    out = mc.copy()
    out.insert(1, "uncertainty_definition", "16th-84th Monte Carlo percentiles")
    out.to_csv(CSV / "final_summary.csv", index=False)


def file_manifest() -> None:
    rows = []
    manifest_path = MANIFEST / "file_manifest.csv"
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        if any(part in {".git", "__pycache__", "tmp"} for part in path.parts):
            continue
        if path == manifest_path:
            continue
        if path.parent == ROOT / "output" / "rendered":
            continue
        if path.suffix in {
            ".aux",
            ".bbl",
            ".blg",
            ".fdb_latexmk",
            ".fls",
            ".log",
            ".out",
            ".toc",
        }:
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": digest,
            }
        )
    MANIFEST.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(manifest_path, index=False)


def main() -> None:
    mc = pd.read_csv(CSV / "monte_carlo_summary.csv")
    regression = pd.read_csv(CSV / "python_mathematica_regression.csv")
    result_table(mc)
    branching_table(mc)
    windows_table(pd.read_csv(CSV / "two_point_windows.csv"))
    stability_table(pd.read_csv(CSV / "stability_summary.csv"))
    transition_table(pd.read_csv(CSV / "transition_term_breakdown.csv"))
    validation_table(pd.read_csv(CSV / "validation_status.csv"))
    regression_table(regression)
    comparison_table(
        mc,
        pd.read_csv(CSV / "literature_comparison.csv"),
        pd.read_csv(CSV / "experimental_comparison.csv"),
    )
    numerical_macros(mc, regression)
    final_summary(mc)
    file_manifest()
    print(f"Wrote generated LaTeX to {GEN}")


if __name__ == "__main__":
    main()
