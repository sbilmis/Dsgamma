"""Create publication tables for vector-final-state radiative widths.

The literature values are the comparison entries collected by Bondar and
Milstein for Ds1 -> Ds* gamma and Bs1 -> Bs* gamma.  The "This work" values are
the preferred lattice-photon, fixed-theta calibrated LCSR results from the
completed Dsstar_gamma and Bsstar_gamma analyses.
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


ROWS = [
    {
        "sector": "Ds",
        "channel_plain": "Ds1(2536)->Ds* gamma",
        "channel_tex": r"$D_{s1}(2536)\to D_s^\ast\gamma$",
        "Godfrey2005": r"$5.6$",
        "GoityRoberts2001": r"$14.6$--$22.8$",
        "GreenRepkoRadford2017": r"$9.21$",
        "RadfordRepkoSaelim2009": r"$8.9$",
        "Chen2020": r"$2.96$--$3.02$",
        "KornerPirjolSchilcher1993": r"$10.4\pm1.0$",
        "LiNiZhong2021": "",
        "GodfreyMoatsSwanson2016": "",
        "LuPanWang2016": "",
        "BondarMilstein2025": r"$29$",
        "This_work_lattice_theta_fixed": r"$23.0\,[16.5,35.1]$",
        "This_work_lattice_theta_scan": r"$22.3\,[10.6,41.8]$",
        "This_work_legacy_theta_fixed": r"$9.74\,[5.78,16.5]$",
        "Experimental_status": r"no direct radiative measurement; $\Gamma_{\rm tot}(D_{s1}(2536))=0.92\pm0.05$ MeV",
    },
    {
        "sector": "Ds",
        "channel_plain": "Ds1(2460)->Ds* gamma",
        "channel_tex": r"$D_{s1}(2460)\to D_s^\ast\gamma$",
        "Godfrey2005": r"$5.5$",
        "GoityRoberts2001": r"$14.0$--$25.1$",
        "GreenRepkoRadford2017": r"$17.4$",
        "RadfordRepkoSaelim2009": r"$15.5$",
        "Chen2020": r"$4.74$--$4.79$",
        "KornerPirjolSchilcher1993": "",
        "LiNiZhong2021": "",
        "GodfreyMoatsSwanson2016": "",
        "LuPanWang2016": "",
        "BondarMilstein2025": r"$104$",
        "This_work_lattice_theta_fixed": r"$92.3\,[59.7,138.7]$",
        "This_work_lattice_theta_scan": r"$93.1\,[57.1,145.1]$",
        "This_work_legacy_theta_fixed": r"$24.2\,[8.63,53.7]$",
        "Experimental_status": r"no direct $D_s^\ast\gamma$ width; BaBar measured $D_s\gamma$ and $D_s^\ast\pi^0$ modes",
    },
    {
        "sector": "Bs",
        "channel_plain": "Bs1(5748/5750)->Bs* gamma",
        "channel_tex": r"$B_{s1}(5750)\to B_s^\ast\gamma$",
        "Godfrey2005": "",
        "GoityRoberts2001": "",
        "GreenRepkoRadford2017": "",
        "RadfordRepkoSaelim2009": "",
        "Chen2020": "",
        "KornerPirjolSchilcher1993": "",
        "LiNiZhong2021": r"$56$",
        "GodfreyMoatsSwanson2016": r"$36.9$",
        "LuPanWang2016": r"$39.5$",
        "BondarMilstein2025": r"$20$",
        "This_work_lattice_theta_fixed": r"$10.2\,[6.63,16.24]$",
        "This_work_lattice_theta_scan": r"$11.5\,[6.62,17.79]$",
        "This_work_legacy_theta_fixed": r"$1.89\,[0.471,5.11]$",
        "Experimental_status": r"no established lower $B_{s1}$ state or radiative measurement",
    },
    {
        "sector": "Bs",
        "channel_plain": "Bs1(5829/5830)->Bs* gamma",
        "channel_tex": r"$B_{s1}(5830)\to B_s^\ast\gamma$",
        "Godfrey2005": "",
        "GoityRoberts2001": "",
        "GreenRepkoRadford2017": "",
        "RadfordRepkoSaelim2009": "",
        "Chen2020": "",
        "KornerPirjolSchilcher1993": "",
        "LiNiZhong2021": r"$53$",
        "GodfreyMoatsSwanson2016": r"$57.3$",
        "LuPanWang2016": r"$98.8$",
        "BondarMilstein2025": r"$41$",
        "This_work_lattice_theta_fixed": r"$7.46\,[3.30,24.0]$",
        "This_work_lattice_theta_scan": r"$2.55\,[0.805,9.75]$",
        "This_work_legacy_theta_fixed": r"$4.33\,[2.05,12.05]$",
        "Experimental_status": r"no direct radiative measurement; observed in $B^{(*)}K$ spectroscopy",
    },
]


def write_csv():
    path = OUT / "vector_final_state_literature_comparison.csv"
    fields = list(ROWS[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(ROWS)
    return path


def make_latex():
    path = OUT / "vector_final_state_literature_comparison_table.tex"
    lines = [
        r"\begin{widetext}",
        r"\begin{table}[t]",
        r"\caption{Representative comparison for vector-final-state radiative widths. Widths are in keV. The last column gives the present calibrated LCSR result with the preferred lattice photon normalization and fixed $\theta=35.3^\circ$; percentile intervals are shown in square brackets.}",
        r"\label{tab:vector-final-literature}",
        r"\begin{ruledtabular}",
        r"\begin{tabular}{lcccccccccp{0.25\textwidth}}",
        r"Channel & \cite{Godfrey:2005ww} & \cite{Goity:2000dk} & \cite{Green:2016occ} & \cite{Radford:2009bs} & \cite{Chen:2020ejk} & \cite{Korner:1992pz} & \cite{Li:2021qgz,Godfrey:2016nwn,Lu:2016bbk} & \cite{Bondar:2025smw} & This work & Experiment\\",
        r"\hline",
    ]
    for r in ROWS:
        bottom_group = ", ".join(x for x in (r["LiNiZhong2021"], r["GodfreyMoatsSwanson2016"], r["LuPanWang2016"]) if x)
        if not bottom_group:
            bottom_group = r"$\cdots$"
        values = [
            r["channel_tex"],
            r["Godfrey2005"] or r"$\cdots$",
            r["GoityRoberts2001"] or r"$\cdots$",
            r["GreenRepkoRadford2017"] or r"$\cdots$",
            r["RadfordRepkoSaelim2009"] or r"$\cdots$",
            r["Chen2020"] or r"$\cdots$",
            r["KornerPirjolSchilcher1993"] or r"$\cdots$",
            bottom_group,
            r["BondarMilstein2025"],
            r["This_work_lattice_theta_fixed"],
            r["Experimental_status"],
        ]
        lines.append(" & ".join(values) + r"\\")
    lines.extend(
        [
            r"\end{tabular}",
            r"\end{ruledtabular}",
            r"\end{table}",
            r"\end{widetext}",
            "",
            r"% For the bottom rows, the combined column lists, in order,",
            r"% \cite{Li:2021qgz}, \cite{Godfrey:2016nwn}, and \cite{Lu:2016bbk}.",
            r"% The machine-readable CSV stores the theta-scan and legacy-photon alternatives.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")
    return path


def make_markdown_note():
    path = OUT / "vector_final_state_literature_comparison_summary.md"
    lines = [
        "# Vector-final-state literature comparison",
        "",
        "Widths are in keV. `This work` uses the preferred lattice photon normalization.",
        "",
        "| Channel | Literature pattern | Bondar-Milstein | This work | Experimental status |",
        "|---|---:|---:|---:|---|",
    ]
    for r in ROWS:
        if r["sector"] == "Ds":
            lit = ", ".join(
                x
                for x in (
                    f"Godfrey {r['Godfrey2005']}",
                    f"Goity-Roberts {r['GoityRoberts2001']}",
                    f"Green-Repko-Radford {r['GreenRepkoRadford2017']}",
                    f"Radford-Repko-Saelim {r['RadfordRepkoSaelim2009']}",
                    f"Chen et al. {r['Chen2020']}",
                    f"Korner-Pirjol-Schilcher {r['KornerPirjolSchilcher1993']}",
                )
                if not x.endswith(" ")
            )
        else:
            lit = (
                f"Li-Ni-Zhong {r['LiNiZhong2021']}, "
                f"Godfrey-Moats-Swanson {r['GodfreyMoatsSwanson2016']}, "
                f"Lu-Pan-Wang {r['LuPanWang2016']}"
            )
        lines.append(
            f"| {r['channel_plain']} | {lit} | {r['BondarMilstein2025']} | "
            f"{r['This_work_lattice_theta_fixed']} | {r['Experimental_status']} |"
        )
    lines.extend(
        [
            "",
            "The CSV also stores the theta-scan preferred result and the legacy-photon result.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")
    return path


def main():
    paths = [write_csv(), make_latex(), make_markdown_note()]
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
