# Outputs Folder

Generated numerical tables, symbolic reductions, Monte Carlo samples, and plots
live here.  The manuscript and notes refer to these files with relative paths, so
avoid moving individual output files unless the TeX sources are updated too.

## Headline Tables

- `combined_recommended_results_table.csv` - final combined Ds1/Bs1 values.
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
