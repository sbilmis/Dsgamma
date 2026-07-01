# Dsgamma

Light-cone QCD sum-rule analysis for radiative transitions

- `Ds1(2460) -> Ds gamma`
- `Ds1(2536) -> Ds gamma`

The calculation is organized around a staged OPE:

- Stage 1: hard photon emission plus two-particle photon distribution amplitudes.
- Stage 2: add three-particle photon DA corrections from the background-gluon insertion.

## Key Files

- `conventions.md` - metric, gamma-matrix, current, and amplitude conventions.
- `inputs_table.csv` - numerical inputs used in the analysis.
- `photon_da.py` - photon DA functions through twist 4.
- `roadmap_Ds1_radiative_LCSR.md` - project roadmap and calculation plan.
- `notes/` - LaTeX notes and reports compiled to PDF.
- `scripts/` - Mathematica/FeynCalc scripts for symbolic traces.
- `outputs/` - generated symbolic outputs.

## Mathematica/FeynCalc

Scripts are run with the direct Wolfram kernel path:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/step2_current_building_blocks.wl
```

Plain `wolframscript` does not reliably find the kernel on this machine.
