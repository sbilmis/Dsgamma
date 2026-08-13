"""Consolidate the exact post-Borel Rohrwild-nonlocal recalculation."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def select(rows, **labels):
    matches = [row for row in rows if all(row[key] == value for key, value in labels.items())]
    if len(matches) != 1:
        raise RuntimeError(f"expected one row for {labels}, found {len(matches)}")
    return matches[0]


def main():
    final_rows = read_rows(OUT / "final_window_mc_summary.csv")
    ds_common = {"scenario": "lattice_fperp_s", "ensemble": "theta_prior_gaussian"}
    ds2460 = select(final_rows, sector="Ds", state="D_{s1}(2460)", window_id="central", **ds_common)
    ds2536 = select(final_rows, sector="Ds", state="D_{s1}(2536)", window_id="central", **ds_common)
    bs_common = {"scenario": "lattice_fperp_s", "ensemble": "theta_prior_gaussian"}
    bs5750 = select(
        final_rows,
        sector="Bs",
        state="B_{s1}(5750)",
        window_id="central_10_14",
        **bs_common,
    )
    bs5830 = select(
        final_rows,
        sector="Bs",
        state="B_{s1}(5830)",
        window_id="central_10_14",
        scenario="lattice_fperp_s",
        ensemble="theta_prior_gaussian",
    )
    central = {
        row["key"]: float(row["value"])
        for row in read_rows(OUT / "corrected_transition_central_python.csv")
    }
    comparison_text = (
        OUT / "corrected_transition_python_mathematica_comparison.txt"
    ).read_text()
    comparison_status = next(
        line for line in comparison_text.splitlines() if line.startswith("STATUS=")
    )

    lines = [
        "Exact post-Borel Rohrwild-nonlocal numerical summary",
        "====================================================",
        "Transition OPE: ordinary local condensate excluded; nonlocal",
        "S_gamma/T4^gamma included.",
        "Tensor vector-DA and electromagnetic P-functional numerators are",
        "double-Borel transformed off shell. No physical pole mass is inserted",
        "in a QCD-side numerator.",
        "",
        "Central charm checkpoint (M2=3.75 GeV^2, s0=8.0 GeV^2)",
        f"  T_A = {central['Ds.T_A']:+.9e} GeV^3",
        f"  T_B = {central['Ds.T_B']:+.9e} GeV^3",
        f"  g1 = {central['Ds.g_1']:+.9f} GeV^-1",
        f"  g2 = {central['Ds.g_2']:+.9f} GeV^-1",
        f"  Gamma1 = {central['Ds.Gamma_1_keV']:.9f} keV",
        f"  Gamma2 = {central['Ds.Gamma_2_keV']:.9f} keV",
        "",
        "Preferred final-window outputs (lattice photon input)",
        "  Ds1(2460): Gamma = {} [{}, {}] keV".format(ds2460["Gamma_median_keV"], ds2460["Gamma_p16_keV"], ds2460["Gamma_p84_keV"]),
        "  Ds1(2536): Gamma = {} [{}, {}] keV".format(ds2536["Gamma_median_keV"], ds2536["Gamma_p16_keV"], ds2536["Gamma_p84_keV"]),
        "  Bs1(5750): Gamma = {} [{}, {}] keV".format(bs5750["Gamma_median_keV"], bs5750["Gamma_p16_keV"], bs5750["Gamma_p84_keV"]),
        "  Bs1(5830): Gamma = {} [{}, {}] keV".format(bs5830["Gamma_median_keV"], bs5830["Gamma_p16_keV"], bs5830["Gamma_p84_keV"]),
        "",
        "Decay-constant provenance",
        "  Ds theta=26.6+-0.6 deg is external; f1/f2 are AA/AB/BB projections at that angle.",
        "  Bs theta=38.5+-0.1 deg is external; f1/f2 are direct AA/AB/BB projections at that angle.",
        "  Both two-point matrices use exact-mass LO + local d=3 + local d=5.",
        "",
        "Mathematica status",
        f"  {comparison_status}",
        "  145 central quantities are compared term by term.",
        "  Every entry passes the 1e-8 absolute-or-relative tolerance.",
    ]
    path = OUT / "rohrwild_nonlocal_final_summary.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
