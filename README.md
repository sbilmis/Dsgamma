# Dsgamma

Light-cone QCD sum-rule analysis for radiative transitions

- `Ds1(2460) -> Ds gamma`
- `Ds1(2536) -> Ds gamma`
- `Bs1(5830) -> Bs gamma`
- a lower, not-yet-established `Bs1 -> Bs gamma` diagnostic channel

The calculation is organized around a staged OPE:

- Stage 1: hard photon emission plus two-particle photon distribution amplitudes.
- Stage 2: add three-particle photon DA corrections from the background-gluon insertion.
- Stage 3: convert basis-current amplitudes to physical-current normalizations using
  two-point residue inputs or calibrated residue scenarios.

The preferred numerical tables use the modern lattice photon-normalization input
for `f_gamma,s^perp`, while the legacy Rohrwild/BBK susceptibility scenario is kept
as a comparison against the older photon-DA literature.

## Key Files

- `conventions.md` - metric, gamma-matrix, current, and amplitude conventions.
- `inputs_table.csv` - numerical inputs used in the analysis.
- `photon_da.py` - photon DA functions through twist 4.
- `roadmap_Ds1_radiative_LCSR.md` - project roadmap and calculation plan.
- `notes/` - active calculation notes, redo-analysis plan, and reports compiled to PDF.
- `draft_prd/` - preserved PRD-style draft with the fuller theory-framework version.
- `manuscript/` - archived local manuscript drafts; not the active editing target.
- `scripts/` - Mathematica/FeynCalc and Python scripts for symbolic and numerical work.
- `outputs/` - generated symbolic outputs, Monte Carlo scans, plots, citation maps,
  and publication-ready tables.

## Active Working Entry Points

- Use `notes/Ds1_to_Ds_gamma_LCSR_notes.tex` as the active calculation and redo-analysis record.
- Use `notes/Ds1_to_Ds_gamma_LCSR_notes.pdf` to read the compiled notes.
- Keep `draft_prd/main.tex` only as a reference copy for the detailed PRD-style derivation.
- The paper itself will be edited elsewhere; `manuscript/` is retained only as an archived snapshot.

## Current Headline Outputs

- `outputs/combined_recommended_results_table.csv` - combined Ds1 and Bs1 form-factor
  and width table.
- `outputs/ds1_recommended_results_table.csv` - Ds1 Stage-3 physical-current table.
- `outputs/stage3_bs1_physical_decay_constant_summary.csv` - Bs1 physical-current
  normalization using Pullin-Zwicky basis-current inputs and the mixed-current
  diagonalization closure.
- `outputs/experimental_comparison_table.csv` - interpretation against available
  experimental information.
- `outputs/input_citation_map.csv` - citation map for the numerical inputs.
- `notes/Ds1_to_Ds_gamma_LCSR_notes.pdf` - full working notes.

## Mathematica/FeynCalc

Scripts are run with the direct Wolfram kernel path:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/step2_current_building_blocks.wl
```

Plain `wolframscript` does not reliably find the kernel on this machine.
