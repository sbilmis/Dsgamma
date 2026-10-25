"""Paired audit of charm transition/residue sampling; no production writes.

Replays the existing seed and compares to the saved 500-sample lattice run.
Then changes only mc, ms, ss on the transition side to the corresponding
residue sample, keeping the sampled lattice photon normalization fixed.
This is a diagnostic variant, not an adopted replacement for paper tables.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import final_window_monte_carlo as mc


def main():
    rng = np.random.default_rng(mc.SEED)
    residues = mc.read_accepted_twopoint_samples()
    conv = mc.precompute_convolutions()
    baseline, paired, inputs = [], [], []
    for scenario in ("legacy_chi_condensate", "lattice_fperp_s"):
        for i in range(mc.N_POINTS):
            vals = mc.sample_ds_inputs(rng, scenario)
            index = int(rng.integers(0, len(residues)))
            residue = residues[index]
            m2 = float(rng.uniform(mc.DS_WINDOW["M2_min"], mc.DS_WINDOW["M2_max"]))
            s01 = float(rng.uniform(mc.DS_WINDOW["s0_2460_min"], mc.DS_WINDOW["s0_2460_max"]))
            s02 = float(rng.uniform(mc.DS_WINDOW["s0_2536_min"], mc.DS_WINDOW["s0_2536_max"]))
            if scenario != "lattice_fperp_s":
                continue
            baseline.extend(mc.evaluate_ds_physical(vals, m2, s01, s02, residue,
                            scenario, "theta_prior_gaussian", conv))
            shared = dict(vals)
            shared.update(mc=float(residue["mc_GeV"]), ms=float(residue["ms_GeV"]),
                          ss=float(residue["ss_GeV3"]))
            shared["chi"] = shared["fperp_s_used"] / shared["ss"]
            paired.extend(mc.evaluate_ds_physical(shared, m2, s01, s02, residue,
                          scenario, "theta_prior_gaussian", conv))
            inputs.append(dict(sample=i, residue_index=index, M2=m2, s01=s01, s02=s02,
                               transition_mc=vals["mc"], residue_mc=shared["mc"],
                               transition_ms=vals["ms"], residue_ms=shared["ms"],
                               transition_ss=vals["ss"], residue_ss=shared["ss"]))
            if (i+1) % 100 == 0:
                print(f"Finished {i+1}/{mc.N_POINTS} paired samples", flush=True)
    with (ROOT / "outputs/final_window_mc_scan.csv").open() as f:
        saved = [r for r in csv.DictReader(f) if r["sector"] == "Ds"
                 and r["scenario"] == "lattice_fperp_s"
                 and r["ensemble"] == "theta_prior_gaussian"]
    assert len(saved) == len(baseline) == len(paired) == 1000
    max_errors = {}
    for key in ["M2", "s0", "theta_deg", "f1_GeV", "f2_GeV", "T_A", "T_B", "g_GeV_inv", "Gamma_keV"]:
        max_errors[key] = max(abs(float(a[key])-float(b[key])) for a,b in zip(saved,baseline))
    reproduced = max(max_errors.values()) < 1e-8
    summaries = []
    for state in ["Ds1_2460", "Ds1_2536"]:
        for observable in ["Gamma_keV", "g_GeV_inv"]:
            for mode, rows in [("saved_sampling_replay",baseline), ("shared_QCD_inputs",paired)]:
                values = [float(r[observable]) for r in rows if r["state_key"] == state]
                if observable == "g_GeV_inv":
                    values = np.abs(values).tolist()
                summaries.append(dict(state=state, observable=observable, mode=mode,
                                      **mc.summarize(values)))
    for name, rows in [("ds_shared_input_summary.csv", summaries),
                       ("ds_shared_input_draws.csv",inputs),
                       ("ds_saved_sampling_replay.csv",baseline),
                       ("ds_shared_input_variant.csv",paired)]:
        mc.write_csv(HERE / name, rows)
    sources = [ROOT/"scripts/final_window_monte_carlo.py",
               ROOT/"scripts/lattice_photon_normalization_comparison.py",
               ROOT/"outputs/final_window_mc_scan.csv",
               ROOT/"outputs/twopoint_ds1_matrix_mc.csv"]
    result = dict(baseline_reproduced=reproduced, baseline_max_absolute_errors=max_errors,
                  samples=mc.N_POINTS, seed=mc.SEED,
                  changed="Transition mc, ms, ss use selected residue sample; fperp fixed; chi=fperp/ss",
                  unchanged="Same angle, residues, photon-shape and other draws, transition M2/s0",
                  status="Audit variant only; original scripts and outputs unchanged",
                  sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  summaries=summaries)
    (HERE/"ds_shared_input_check.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    assert reproduced, "Replay mismatch: interpret no paired effect until provenance is reconciled"


if __name__ == "__main__":
    main()
