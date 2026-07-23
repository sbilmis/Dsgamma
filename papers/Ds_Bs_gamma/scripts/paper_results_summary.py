"""Build paper-facing tables from the final-window Monte Carlo summary."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def photon_label(scenario):
    if scenario == "lattice_fperp_s":
        return "f_gamma_s^perp(1 GeV)=-55.1+-1.9 MeV"
    return "chi_s(1 GeV)=3.15+-0.30 GeV^-2"


def provenance(sector):
    if sector == "Ds":
        return "prior_theta_plus_AA_AB_BB_projection_LO_d3_d5"
    return "basis_inputs_plus_diagonal_closure"


def status(row):
    if row["sector"] == "Bs" and "5750" in row["state"]:
        return "unobserved_lower_state_diagnostic"
    if (
        row["scenario"] == "lattice_fperp_s"
        and row["ensemble"] == "theta_prior_gaussian"
    ):
        return "preferred_photon_and_mixing_scenario"
    return "comparison_scenario"


def selected_final_rows(rows):
    selected = []
    for row in rows:
        if row["sector"] == "Bs" and row["window_id"] != "central_10_14":
            continue
        if row["sector"] == "Ds" and row["window_id"] != "central":
            continue
        selected.append(row)
    return selected


def main():
    source = selected_final_rows(read(OUT / "final_window_mc_summary.csv"))
    compact = []
    for row in source:
        compact.append(
            {
                "sector": row["sector"],
                "state": row["state"],
                "photon_input": row["scenario"],
                "theta_treatment": row["ensemble"],
                "window_id": row["window_id"],
                "transition_scheme": "rohrwild_nonlocal",
                "decay_constant_provenance": provenance(row["sector"]),
                "f1_median_GeV": row["f1_median_GeV"],
                "f1_p16_GeV": row["f1_p16_GeV"],
                "f1_p84_GeV": row["f1_p84_GeV"],
                "f2_median_GeV": row["f2_median_GeV"],
                "f2_p16_GeV": row["f2_p16_GeV"],
                "f2_p84_GeV": row["f2_p84_GeV"],
                "G_median_GeV_inv": row["g_median_GeV_inv"],
                "G_p16_GeV_inv": row["g_p16_GeV_inv"],
                "G_p84_GeV_inv": row["g_p84_GeV_inv"],
                "Gamma_median_keV": row["Gamma_median_keV"],
                "Gamma_p16_keV": row["Gamma_p16_keV"],
                "Gamma_p84_keV": row["Gamma_p84_keV"],
                "status": status(row),
            }
        )
    compact.sort(
        key=lambda r: (
            r["sector"], r["state"], r["photon_input"], r["theta_treatment"]
        )
    )

    compact_path = OUT / "paper_final_results_rohrwild_nonlocal.csv"
    with compact_path.open("w", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=list(compact[0]))
        writer.writeheader()
        writer.writerows(compact)

    fields = [
        "sector", "state", "photon_input", "mixing_input", "status",
        "normalization_input", "residue_1_GeV", "residue_1_p16_GeV",
        "residue_1_p84_GeV", "residue_2_GeV", "residue_2_p16_GeV",
        "residue_2_p84_GeV", "effective_mixed_residue_GeV",
        "effective_mixed_residue_p16_GeV", "effective_mixed_residue_p84_GeV",
        "G_GeV_minus1_median", "G_GeV_minus1_p16", "G_GeV_minus1_p84",
        "Gamma_keV_median", "Gamma_keV_minus", "Gamma_keV_plus",
        "Gamma_keV_p16", "Gamma_keV_p84", "comment",
    ]
    combined = []
    for row in compact:
        median = float(row["Gamma_median_keV"])
        p16 = float(row["Gamma_p16_keV"])
        p84 = float(row["Gamma_p84_keV"])
        combined.append(
            {
                "sector": row["sector"],
                "state": row["state"],
                "photon_input": photon_label(row["photon_input"]),
                "mixing_input": row["theta_treatment"],
                "status": row["status"],
                "normalization_input": row["decay_constant_provenance"],
                "residue_1_GeV": row["f1_median_GeV"],
                "residue_1_p16_GeV": row["f1_p16_GeV"],
                "residue_1_p84_GeV": row["f1_p84_GeV"],
                "residue_2_GeV": row["f2_median_GeV"],
                "residue_2_p16_GeV": row["f2_p16_GeV"],
                "residue_2_p84_GeV": row["f2_p84_GeV"],
                "effective_mixed_residue_GeV": "",
                "effective_mixed_residue_p16_GeV": "",
                "effective_mixed_residue_p84_GeV": "",
                "G_GeV_minus1_median": row["G_median_GeV_inv"],
                "G_GeV_minus1_p16": row["G_p16_GeV_inv"],
                "G_GeV_minus1_p84": row["G_p84_GeV_inv"],
                "Gamma_keV_median": median,
                "Gamma_keV_minus": median - p16,
                "Gamma_keV_plus": p84 - median,
                "Gamma_keV_p16": p16,
                "Gamma_keV_p84": p84,
                "comment": (
                    "final selected window; Rohrwild nonlocal transition OPE; "
                    "ordinary local term excluded; " + row["decay_constant_provenance"]
                ),
            }
        )
    combined_path = OUT / "combined_recommended_results_table.csv"
    with combined_path.open("w", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=fields)
        writer.writeheader()
        writer.writerows(combined)

    preferred = [r for r in compact if r["status"] == "preferred_photon_and_mixing_scenario"]
    print(f"Wrote {compact_path}")
    print(f"Updated {combined_path}")
    print("Preferred lattice-photon rows:")
    for row in preferred:
        print(
            f"  {row['state']}: G={float(row['G_median_GeV_inv']):+.4g} "
            f"[{float(row['G_p16_GeV_inv']):+.4g},{float(row['G_p84_GeV_inv']):+.4g}], "
            f"Gamma={float(row['Gamma_median_keV']):.4g} "
            f"[{float(row['Gamma_p16_keV']):.4g},{float(row['Gamma_p84_keV']):.4g}] keV"
        )


if __name__ == "__main__":
    main()
