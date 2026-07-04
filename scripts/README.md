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
