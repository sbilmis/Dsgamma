# Scripts Folder

This folder contains the symbolic and numerical workflow.

## Main Numerical Scripts

- `stage3_ds1_physical_residue_normalization.py` - Ds1 physical-current
  normalization.
- `stage3_bs1_physical_decay_constants.py` - Bs1 physical-current normalization
  from Pullin-Zwicky basis-current inputs and the mixed-current diagonalization
  closure.
- `lattice_photon_normalization_comparison.py` - legacy susceptibility versus
  lattice photon normalization.
- `stage2_stability_plots.py` - Borel and threshold stability plots.
- `redo_stability_windows.py` - redo Borel and continuum-threshold stability
  plots with \(s_0\) written as the squared threshold in GeV^2.
- `window_diagnostics.py` - compares candidate \(M^2,s_0\) windows using
  perturbative pole-fraction proxies, photon-DA continuum retention, OPE
  hierarchy, and numerical stability.
- `publication_stability_selected_windows.py` - publication-style selected
  window plots, varying \(M^2\) at fixed \(s_0\) and varying \(s_0\) at fixed
  \(M^2\).
- `publication_plots.py` and `uncertainty_publication_plots.py` - publication
  plotting helpers.

## Symbolic/FeynCalc Scripts

The `.wl` files store the Mathematica/FeynCalc trace and kernel reductions.  On
this machine, run them with the direct Wolfram kernel path:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/step2_current_building_blocks.wl
```

## Generated Files

Scripts write tables, summaries, and figures into `../outputs/`.
