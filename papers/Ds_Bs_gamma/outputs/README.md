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
- `corrected_transition_python_mathematica_comparison.txt` - term-by-term
  comparison of all 145 corrected central transition quantities.
- `step14_explicit_double_borel_forms.txt` - raw \(E_Q\) double-Borel
  identities, explicitly separated from continuum subtraction.
- `step15_complete_three_particle_borel.txt` - final three-particle
  invariants, with \(E_Q\) for the channel-specific gluonic term and
  \(E_Q-E_0\) for the gauge-completion electromagnetic term.

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

- `final_angle_stability.pdf` - preferred-angle coupling and width stability
  versus \(M^2\) and \(s_0\).
- `final_window_mc_width_histograms.pdf` - final Ds1/Bs1 Monte Carlo widths.
- `final_window_bs_window_crosscheck.pdf` - selected versus cross-check bottom
  windows.

## Current Bs1 Normalization Outputs

- `twopoint_bs1_matrix_mc_summary.txt`
- `twopoint_bs1_matrix_mc.csv`
- `twopoint_bs1_physical_residue_grid.csv`

The older `stage3_bs1_pz_*` and
`stage3_bs1_physical_decay_constant_*` files are retained as
basis-normalization diagnostics, not as the final physical-current result.
