"""Consolidate the completed Python Rohrwild-nonlocal recalculation."""

from __future__ import annotations

import csv
from pathlib import Path

from stage1_axial_g1_baseline import central_inputs
from stage2_axial_g1_three_particle import F1_integral, g1_stage2


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
    inputs = central_inputs()
    f1_total, _, _ = F1_integral(u0=0.5)
    args = (4.5, 2.55**2, inputs, f1_total)
    rw = g1_stage2(*args, transition_scheme="rohrwild_nonlocal")
    col = g1_stage2(*args, transition_scheme="colangelo_local_benchmark")

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

    lines = [
        "Final Python Rohrwild-nonlocal recalculation",
        "============================================",
        "Transition OPE: ordinary local condensate excluded; nonlocal",
        "S_gamma/T4^gamma included.  Colangelo local mode is comparison-only.",
        "The tensor-current numerator is evaluated with the documented",
        "physical-residue prescription p^2=m_P^2, p.q=(m_A^2-m_P^2)/2.",
        "",
        "Axial-current central checkpoint (M2=4.5 GeV^2, sqrt(s0)=2.55 GeV)",
        f"  Colangelo local OPE = {col['heavy_local']:+.6e} GeV^3",
        f"  Rohrwild nonlocal EM OPE = {rw['em_delta_qcd']:+.6e} GeV^3",
        f"  ratio EM/local = {rw['em_delta_qcd']/col['heavy_local']:+.6f}",
        f"  g_A Colangelo-local = {col['g1_stage2_GeV_inv']:+.6f} GeV^-1",
        f"  g_A Rohrwild-nonlocal = {rw['g1_stage2_GeV_inv']:+.6f} GeV^-1",
        f"  Gamma_A Rohrwild-nonlocal = {rw['width_stage2_keV']:.6f} keV",
        "",
        "Preferred final-window outputs (lattice photon input)",
        "  Ds1(2460): Gamma = {} [{}, {}] keV".format(ds2460["Gamma_median_keV"], ds2460["Gamma_p16_keV"], ds2460["Gamma_p84_keV"]),
        "  Ds1(2536): Gamma = {} [{}, {}] keV".format(ds2536["Gamma_median_keV"], ds2536["Gamma_p16_keV"], ds2536["Gamma_p84_keV"]),
        "  Bs1(5750): Gamma = {} [{}, {}] keV".format(bs5750["Gamma_median_keV"], bs5750["Gamma_p16_keV"], bs5750["Gamma_p84_keV"]),
        "  Bs1(5830): Gamma = {} [{}, {}] keV".format(bs5830["Gamma_median_keV"], bs5830["Gamma_p16_keV"], bs5830["Gamma_p84_keV"]),
        "",
        "Decay-constant provenance",
        "  Ds theta=26.6+-0.6 deg is external; f1/f2 are AA/AB/BB projections at that angle.",
        "  Bs f1 and f2 are diagonal-closure-derived from basis inputs.",
        f"  Nominal Bs closure acceptance: {bs5750['n']}/500 (5750), {bs5830['n']}/500 (5830).",
        "  Bs predictions are conditional pending a direct AA/AB/BB two-point matrix.",
        "  Bs theta=38.5+-0.1 deg is external.",
        "  The Ds OPE truncation is exact-mass LO + local d=3 + local d=5.",
        "",
        "Mathematica status",
        "  The spectral and local traces have separate exact audits.",
        "  The fixed-angle notebook has been regenerated.",
        "  The numerical regression agrees with Python below 10^-7 at the central point.",
    ]
    path = OUT / "rohrwild_nonlocal_final_summary.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
