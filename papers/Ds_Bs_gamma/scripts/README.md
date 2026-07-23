# Scripts Folder

This folder contains the symbolic and numerical workflow.

## Main Numerical Scripts

- `twopoint_ds1_matrix_sumrule.py` - preferred Ds1 (AA/AB/BB) two-point
  normalization, including the deterministic scan and input Monte Carlo.
- `stage3_ds1_physical_residue_normalization.py` - legacy external-(f_1)/overlap
  normalization retained only for reproducibility and comparison.
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
- `paper_results_summary.py` - regenerates the authoritative compact paper
  table and the established combined-table schema from the Stage-3 summaries.

## Symbolic/FeynCalc Scripts

The `.wl` files store the Mathematica/FeynCalc trace and kernel reductions.  On
this machine, run them with the direct Wolfram kernel path:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/step2_current_building_blocks.wl
```

For the consolidated independent check from currents through Wick contraction,
propagators, hard traces, E1 projection, and two-particle Fierz traces, run:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/mathematica_current_wick_trace_audit.wl
```

The audit writes `outputs/mathematica_current_wick_trace_audit.txt` and exits
with a failure status if any exact symbolic assertion fails.

The standard-correlator momentum routing and the fermion-line Ward identities
are checked independently by:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/mathematica_correlator_routing_audit.wl
```

The next hard-loop stage, including the exact shifted Feynman-parameter
denominators, is checked by:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/mathematica_hard_loop_parameterization.wl
```

The physical-state rotation, residue-weighted coupling, Colangelo-style OPE
component assembly, mixed-current decay constants, radiative width, and
pure-axial Colangelo limit are checked by:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/mathematica_state_mixing_projection.wl
```

The complete normalized-current charm two-point matrix has three independent
audits:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/twopoint_dirac_spectral_audit.wl
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/twopoint_local_condensate_audit.wl
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/mathematica_twopoint_ds1_matrix_check.wl
```

For interactive evaluation open
`../notebooks/Ds1_AA_AB_BB_two_point_sumrule.nb`. Regenerate it with
`scripts/build_twopoint_ds1_matrix_notebook.wl`.

For interactive, cell-by-cell evaluation, open
`../notebooks/DsBs_gamma_symbolic_derivation.nb`.  It contains 35 input cells
covering the standard correlator, Wick contraction, propagator routing, Ward
identities, the explicit vacuum/background-field propagator decomposition,
axial/tensor traces, E1 cores, Feynman parameterization, and soft two-particle
Fierz traces, followed by the state-mixing, explicit axial OPE assembly, and
Colangelo-limit audit.
Regenerate it after editing the source cells with:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt -script scripts/build_dsbs_symbolic_notebook.wl
```

The older `step2_*` and `step3_*` hard-trace scripts are retained as legacy
comparison artifacts.  Their strange-line `k-p+q` routing is not the routing
used in the restarted notebook or the three audits above.

The completed Python Rohrwild-scheme implementation is tracked by:

- `step12_em_da_current_kernels.wl`, which derives the current-specific
  FeynCalc kernels for the electromagnetic (S_\gamma,T_4^\gamma) sector and
  compares them with Rohrwild's vector-current reference;
- `rohrwild_em_da.py`, which implements the analytic symmetric-double-Borel
  S_gamma/T4_gamma kernels, including the explicitly documented support
  correction in the printed P-functional;
- `local_scheme_comparison.py`, which generates conceptual, numerical, scan,
  and LaTeX tables comparing the alternative Rohrwild-nonlocal and
  Colangelo-local prescriptions;
- `rohrwild_nonlocal_final_summary.py`, which consolidates the regenerated Ds
  and Bs results;
- `decay_constant_provenance_audit.py`, which distinguishes the anchored or
  closure-derived (f_1,f_2) values from a genuine two-point-QCDSR result.

`step12_em_da_current_kernels.wl` is intentionally kept as a standalone batch
cell source.  The next stage is to duplicate the finalized Python formulas
and numerical tables independently in the canonical `.nb`.

## Generated Files

Scripts write tables, summaries, and figures into `../outputs/`.
