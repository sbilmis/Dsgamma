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
  sector is now obtained from the direct normalized-current AA/AB/BB
  continuation as well.

The continuum prescription is stated term by term.  The leading twist-2
contribution carries
\(\Delta E_Q=e^{-m_Q^2/M^2}-e^{-s_0/M^2}\), while twist-3, twist-4, and the
gluonic three-particle term retain \(E_Q=e^{-m_Q^2/M^2}\), as in the
channel-specific \(D_{s1}\to D_s\gamma\) sum rule.  The separate
gauge-completion electromagnetic three-particle term uses Rohrwild's
\(I_F\)-type \(\Delta E_Q\) subtraction.

The preferred numerical tables use the modern lattice photon-normalization input
for `f_gamma,s^perp`, while the legacy Rohrwild/BBK susceptibility scenario is kept
as a comparison against the older photon-DA literature.

## Key Files

- `arxiv/main-3.pdf` - final Prism PDF selected for arXiv submission.
- `arxiv/figures/` - final figure PDFs used with the Prism/arXiv version.
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
- `manuscript/dsbs_radiative_lcsr_polished.pdf` - current paper PDF.
- `scripts/` - Mathematica/FeynCalc and Python scripts for symbolic and numerical work.
- `notebooks/DsBs_gamma_symbolic_derivation.nb` - canonical cell-by-cell
  Mathematica restart from currents through hard-loop Feynman parameters.
- `notebooks/Ds1_AA_AB_BB_two_point_sumrule.nb` - self-contained cell-by-cell
  Mathematica calculation of the charm mixing angle and both physical residues.
- `outputs/` - generated symbolic outputs, Monte Carlo scans, plots, citation maps,
  and publication-ready tables.
- `archive/legacy_latex/` - preserved older LaTeX paper drafts; these are not
  the current Prism source of truth.

## Active Working Entry Points

- Use `arxiv/main-3.pdf` as the final paper snapshot selected for arXiv.
- Use `notes/Ds1_to_Ds_gamma_LCSR_notes.tex` as the active calculation and redo-analysis record.
- Use `notes/Ds1_to_Ds_gamma_LCSR_notes.pdf` to read the compiled notes.
- Use `notes/DsBs_gamma_Mathematica_derivation.tex` only for the independent
  symbolic restart and its continuation.
- Keep `archive/legacy_latex/draft_prd/main.tex` only as a reference copy for
  the older detailed PRD-style derivation.
- Keep `archive/legacy_latex/manuscript/` as a preserved pre-Prism manuscript
  snapshot.

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
- `outputs/twopoint_bs1_matrix_mc_summary.txt` - direct Bs1 AA/AB/BB
  two-point normalization and accepted-window Monte Carlo.
- `outputs/experimental_comparison_table.csv` - interpretation against available
  experimental information.
- `outputs/input_citation_map.csv` - citation map for the numerical inputs.
- `notes/Ds1_to_Ds_gamma_LCSR_notes.pdf` - full working notes.

## Mathematica/FeynCalc

For interactive work, open
`notebooks/DsBs_gamma_symbolic_derivation.nb` and evaluate its 53 input cells
from top to bottom.  The corresponding batch audits are documented in
`scripts/README.md`.

Scripts are run with the direct Wolfram kernel path:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/step2_current_building_blocks.wl
```

Plain `wolframscript` does not reliably find the kernel on this machine.
Run these commands from `papers/Ds_Bs_gamma/`.
