"""Consolidate the provisional Rohrwild local-free recalculation.

The report intentionally separates completed numerical propagation from the
remaining current-specific S_gamma/T4^gamma and tensor double-Borel work.
"""

from __future__ import annotations

import csv
from pathlib import Path

from stage1_axial_g1_baseline import central_inputs, g1_stage1
from stage2_axial_g1_three_particle import F1_integral, g1_stage2, width_keV


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def select(rows: list[dict[str, str]], **labels: str) -> dict[str, str]:
    matches = [row for row in rows if all(row[key] == value for key, value in labels.items())]
    if len(matches) != 1:
        raise RuntimeError(f"expected one row for {labels}, found {len(matches)}")
    return matches[0]


def main() -> None:
    inputs = central_inputs()
    m2 = 4.5
    s0 = 2.55**2
    f1_total, _, _ = F1_integral(u0=0.5)
    local_free = g1_stage2(m2, s0, inputs, f1_total)
    legacy_stage1 = g1_stage1(
        m2,
        s0,
        **inputs,
        transition_scheme="colangelo_local_benchmark",
    )
    legacy_g = legacy_stage1["g1_GeV_inv"] + local_free["F1_delta_g1"]
    legacy_width = width_keV(inputs["m_ds1"], inputs["m_ds"], legacy_g)

    ds_summary = read_rows(OUT / "stage3_ds1_physical_residue_summary.csv")
    bs_summary = read_rows(OUT / "stage3_bs1_physical_decay_constant_summary.csv")
    ds_lattice = select(ds_summary, scenario="lattice_fperp_s", ensemble="theta_fixed")
    bs_lattice = select(
        bs_summary,
        state="Bs1_5830_observed_high_combo",
        scenario="lattice_fperp_s",
        ensemble="theta_fixed",
    )

    lines = [
        "Provisional Rohrwild local-free recalculation",
        "==============================================",
        "",
        "What has been propagated",
        "  - the standalone Colangelo heavy-charge local condensate is zero in",
        "    every default transition call;",
        "  - hard, two-particle, and background-gluon three-particle terms were",
        "    regenerated through the Ds and Bs mixed-current scans;",
        "  - the legacy local term remains available only as a benchmark switch.",
        "",
        "What is not complete",
        "  - current-specific S_gamma and T4^gamma electromagnetic DA kernels;",
        "  - exact tensor-current double-Borel reduction (the existing GB sector",
        "    remains residue-matched and is therefore an estimate).",
        "",
        "Axial-current central checkpoint (M2=4.5 GeV^2, sqrt(s0)=2.55 GeV)",
        f"  Colangelo-local benchmark g_A = {legacy_g:+.6f} GeV^-1",
        f"  Colangelo-local benchmark width = {legacy_width:.6f} keV",
        f"  local-free provisional g_A = {local_free['g1_stage2_GeV_inv']:+.6f} GeV^-1",
        f"  local-free provisional width = {local_free['width_stage2_keV']:.6f} keV",
        f"  removed local OPE term = {local_free['legacy_heavy_local_diagnostic']:+.6e} GeV^3",
        "",
        "Mixed-current provisional outputs (lattice normalization, fixed angle)",
        "  Ds1(2460): Gamma = {} [{}, {}] keV".format(
            ds_lattice["Gamma2460_median"],
            ds_lattice["Gamma2460_p16"],
            ds_lattice["Gamma2460_p84"],
        ),
        "  Bs1(5830): Gamma = {} [{}, {}] keV".format(
            bs_lattice["Gamma_quoted_phys_keV_median"],
            bs_lattice["Gamma_quoted_phys_keV_p16"],
            bs_lattice["Gamma_quoted_phys_keV_p84"],
        ),
        "",
        "Decay-constant status",
        "  Ds f1 is an external anchor and f2 is overlap-derived.",
        "  Bs f1 and f2 are diagonal-closure-derived from basis inputs.",
        "  Neither pair is an independent complete two-point-QCDSR prediction.",
        "",
        "Do not use this report as the final paper table until the two missing",
        "transition-kernel items above have been completed and re-propagated.",
    ]
    path = OUT / "rohrwild_local_free_recalculation_summary.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
