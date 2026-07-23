# Outputs Folder

Generated numerical tables, symbolic reductions, Monte Carlo samples, and plots
live here.  The manuscript and notes refer to these files with relative paths, so
avoid moving individual output files unless the TeX sources are updated too.

## Mathematica Audit Outputs

- `mathematica_current_wick_trace_audit.txt` - currents through E1 cores and
  two-particle soft traces; 19 exact checks.
- `mathematica_correlator_routing_audit.txt` - standard versus legacy routing
  and fermion-line Ward identities.
- `mathematica_hard_loop_parameterization.txt` - shifted heavy- and
  strange-emission triangle denominators; 4 exact checks.

## Headline Tables

- `paper_final_results_rohrwild_nonlocal.csv` - authoritative compact table;
  explicitly labels the transition scheme and decay-constant provenance.
- `combined_recommended_results_table.csv` - the same regenerated Ds1/Bs1
  values in the established manuscript schema.
- `local_scheme_numerical_comparison.csv`, `local_scheme_scan_summary.csv`, and
  `local_scheme_comparison_table.tex` - Colangelo-local versus
  Rohrwild-nonlocal comparison.  The two contributions are alternatives.
- `rohrwild_nonlocal_final_summary.txt` - central check and final qualifications.
- `bondar_literature_comparison_table.csv` - comparison with representative
  literature values.
- `paper_input_summary_table.csv` - compact input table for the manuscript.
- `input_citation_map.csv` - citation bookkeeping for inputs.
- `experimental_comparison_table.csv` - interpretation against experiment.

## Headline Plots

- `stage2_stability_central.pdf` - Ds1 Borel/threshold stability.
- `bs1_5830_stability.pdf` - Bs1 Borel/threshold stability.
- `lattice_fperp_mc_width_histograms_theta_fixed.pdf` - Ds1 preferred MC widths.
- `stage3_bs1_physical_width_histograms.pdf` - Bs1 physical-current MC widths.

## Current Bs1 Normalization Outputs

- `stage3_bs1_physical_decay_constant_summary.csv`
- `stage3_bs1_physical_decay_constant_summary.txt`
- `stage3_bs1_physical_decay_constant_scan.csv`

The older `stage3_bs1_pz_*` files are retained as a basis-normalization
diagnostic, not as the final physical-current result.
