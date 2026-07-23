# \(D_s\gamma\) and \(B_s\gamma\) paper

Light-cone QCD sum-rule analysis for radiative transitions

- `Ds1(2460) -> Ds gamma`
- `Ds1(2536) -> Ds gamma`
- `Bs1(5830) -> Bs gamma`
- a lower, not-yet-established `Bs1 -> Bs gamma` diagnostic channel

The calculation is organized around a staged OPE:

- Stage 1: hard photon emission plus two-particle photon distribution amplitudes.
- Stage 2: add three-particle photon DA corrections from the background-gluon insertion.
- Stage 3: convert the charm basis-current amplitudes using the complete
  normalized-current (AA/AB/BB) two-point QCD sum rule.  The older external-
  (f_1)/overlap closure is retained only as a legacy comparison.  The bottom
  sector still uses its separately documented basis-input closure.

The preferred numerical tables use the modern lattice photon-normalization input
for `f_gamma,s^perp`, while the legacy Rohrwild/BBK susceptibility scenario is kept
as a comparison against the older photon-DA literature.

## Key Files

- `../../shared/conventions.md` - shared metric, gamma-matrix, current, and amplitude conventions.
- `inputs_table.csv` - numerical inputs used in the analysis.
- `../../shared/photon_da.py` - shared photon DA functions through twist 4.
- `roadmap_Ds1_radiative_LCSR.md` - project roadmap and calculation plan.
- `notes/` - active calculation notes, redo-analysis plan, and reports compiled to PDF.
- `notes/DsBs_gamma_Mathematica_derivation.pdf` - clean standalone record of
  the new Mathematica/FeynCalc derivation.
- `notes/Ds1_AA_AB_BB_two_point_sumrule.pdf` - clean standalone record of the
  charm two-point matrix normalization and its Python/Mathematica comparison.
- `draft_prd/` - preserved PRD-style draft with the fuller theory-framework version.
- `manuscript/` - archived local manuscript drafts; not the active editing target.
- `scripts/` - Mathematica/FeynCalc and Python scripts for symbolic and numerical work.
- `notebooks/DsBs_gamma_symbolic_derivation.nb` - canonical cell-by-cell
  Mathematica restart from currents through hard-loop Feynman parameters.
- `notebooks/Ds1_AA_AB_BB_two_point_sumrule.nb` - self-contained cell-by-cell
  Mathematica calculation of the charm mixing angle and both physical residues.
- `outputs/` - generated symbolic outputs, Monte Carlo scans, plots, citation maps,
  and publication-ready tables.

## Active Working Entry Points

- Use `notes/Ds1_to_Ds_gamma_LCSR_notes.tex` as the active calculation and redo-analysis record.
- Use `notes/Ds1_to_Ds_gamma_LCSR_notes.pdf` to read the compiled notes.
- Use `notes/DsBs_gamma_Mathematica_derivation.tex` only for the independent
  symbolic restart and its continuation.
- Keep `draft_prd/main.tex` only as a reference copy for the detailed PRD-style derivation.
- The paper itself will be edited elsewhere; `manuscript/` is retained only as an archived snapshot.

## Current Headline Outputs

- `outputs/paper_final_results_rohrwild_nonlocal.csv` - compact authoritative
  Ds/Bs table with explicit transition-scheme and decay-constant provenance.
- `outputs/twopoint_ds1_matrix_mc_summary.txt` - charm (AA/AB/BB) two-point
  Monte Carlo result for θ, (f_1), and (f_2).
- `outputs/mathematica_twopoint_ds1_matrix_check.txt` - independent central-point
  Mathematica/Python regression.
- `outputs/combined_recommended_results_table.csv` - the same regenerated Ds1
  and Bs1 values in the established manuscript schema.
- `outputs/local_scheme_numerical_comparison.csv` and
  `outputs/local_scheme_comparison_table.tex` - Rohrwild-nonlocal versus
  Colangelo-local transition-OPE comparison (the terms are alternatives).
- `outputs/rohrwild_nonlocal_final_summary.txt` - final Python-side checkpoint
  and decay-constant provenance statement.
- `outputs/ds1_recommended_results_table.csv` - older pre-completion table;
  superseded by the two regenerated tables above.
- `outputs/stage3_bs1_physical_decay_constant_summary.csv` - Bs1 physical-current
  normalization using Pullin-Zwicky basis-current inputs and the mixed-current
  diagonalization closure.
- `outputs/experimental_comparison_table.csv` - interpretation against available
  experimental information.
- `outputs/input_citation_map.csv` - citation map for the numerical inputs.
- `notes/Ds1_to_Ds_gamma_LCSR_notes.pdf` - full working notes.

## Mathematica/FeynCalc

For interactive work, open
`notebooks/DsBs_gamma_symbolic_derivation.nb` and evaluate its 35 input cells
from top to bottom.  The corresponding batch audits are documented in
`scripts/README.md`.

Scripts are run with the direct Wolfram kernel path:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/step2_current_building_blocks.wl
```

Plain `wolframscript` does not reliably find the kernel on this machine.
Run these commands from `papers/Ds_Bs_gamma/`.
