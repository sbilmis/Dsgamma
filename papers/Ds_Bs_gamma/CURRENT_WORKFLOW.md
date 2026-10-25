# Current calculation map — updated 25 October 2026

This map separates the saved publication baseline, diagnostics, proposed revisions, and historical experiments. It does not certify the physical assumptions of the calculation.

The canonical forward task queue is now [project.org](project.org). A fresh
Codex session must ask for the author's explicit approval before starting
that isolated recalculation. The source map below describes the preserved
production baseline; it is not an instruction to run those entry points.

## Numerical pipeline

| Role | Source | Data or result |
| --- | --- | --- |
| Published angle inputs | `scripts/mixing_angle_inputs.py` | External 26.6 ± 0.6 degrees and 38.5 ± 0.1 degrees |
| Ds two-point matrix and residues | `scripts/twopoint_ds1_matrix_sumrule.py` | `outputs/twopoint_ds1_matrix_mc.csv` |
| Bs two-point matrix and residues | `scripts/twopoint_bs1_matrix_sumrule.py` | `outputs/twopoint_bs1_matrix_mc.csv` |
| Transition invariants and state projection | `scripts/rohrwild_transition_exact.py` | Called by central checks, scans and figure builders |
| Publication uncertainty sampling | `scripts/final_window_monte_carlo.py` | `outputs/final_window_mc_scan.csv`, `outputs/final_window_mc_summary.csv` |
| Paper-table formatting | `scripts/paper_results_summary.py` | `outputs/paper_final_results_rohrwild_nonlocal.csv`, `outputs/combined_recommended_results_table.csv` |
| Central Python/Mathematica comparison | `scripts/corrected_transition_central_python.py`, `scripts/mathematica_corrected_transition_numerics.wl`, `scripts/compare_corrected_transition_numerics.py` | `outputs/corrected_transition_*` |
| Consolidated checkpoint | `scripts/rohrwild_nonlocal_final_summary.py` | `outputs/rohrwild_nonlocal_final_summary.txt` |

Most entry points overwrite their established outputs when run normally. Audit reruns belong in a separate directory; do not replace the baseline with an unreviewed diagnostic variant.

The preferred table/MC selection is:

- `scenario = lattice_fperp_s` (called `photon_input` in the compact table);
- `ensemble = theta_prior_gaussian` (called `theta_treatment` in that table);
- `window_id = central` for Ds and `central_10_14` for Bs.

`paper_final_results_rohrwild_nonlocal.csv` includes comparison scenarios. Its first row is not necessarily the preferred result.

## Figures

| Result | Builder | Inputs |
| --- | --- | --- |
| Current manuscript Figure 1: individual f1/f2 | `review/referee1_revision/build_decay_constant_figures.py` | Representative inputs shared with the angular-width scan |
| Current manuscript Figure 2: Borel dependence | `figures/figure1_borel_form_factors/build_figure1_data.py` | Accepted residue samples plus transition code |
| Current manuscript Figure 3: angular widths | `figures/figure2_mixing_angle_widths/build_figure2_data.py` | Accepted residue samples; reprojects f1/f2 at every angle |
| Current manuscript Figure 4: distributions | `figures/figure5_monte_carlo_widths/build_figure5_data.py` | Preferred rows of the final-window MC scan; historical folder number |
| Response-only f2/f1 | `review/referee1_revision/build_decay_constant_figures.py` | Same representative inputs as current Figure 3 |
| Additional signed Pi_A/Pi_B/Pi_2 panel | `review/referee_a_major_ii/plot_signed_invariants.py` | Saved old-prescription Ds2536 diagnostic, to be regenerated after validation |

The decay-constant curves are fixed-input illustrations, not marginal
distributions or uncertainty bands. Figure numbering above follows the
current r1-2.pdf, not the historical builder folder names.
`figures/ope_topologies/` remains in place for ongoing referee work.
`arxiv/` preserves submission snapshots; Prism is the live manuscript.

## Audit status

Read [the audit report](review/correctness_audit_20261019/REPORT.md) before using the referee reply as final.

- The original numerical implementation of PRD 111, 014005 has not been located. The current `.nb` file is not that implementation.
- The printed 2025 appendix has a leading-term label/sign mismatch relative to the independently traced current matrix. A marked reconstruction is saved in the audit; it has not established the original code or published uncertainty extraction.
- The leading mixed-condensate matrix has now been derived independently. The
  current code omits the heavy-line background-gluon contribution to AB and BB.
  [The resolution package](review/mixed_condensate_resolution/README.md) contains
  an isolated implementation, fixed-point impact checks, and an unapplied
  consistent patch. Production sources and accepted ensembles still contain
  the original light-only d5 contribution; do not mistake them for the completed
  implementation or silently mix their outputs with the proposed correction.
- The current charm MC independently samples some QCD inputs on the transition and residue sides. A paired shared-input calculation is saved in the audit and leaves the baseline untouched.
- Algebraic and numerical checks do not establish exclusive coupling of each rotated current to one physical state. That hadronic assumption remains a separate theoretical question.

## Dependencies and historical files

The static local Python import inventory is in `review/correctness_audit_20261019/python_dependency_inventory.json`. It is supplemented by source inspection; it does not automatically detect all dynamic imports, Wolfram calls or data dependencies.

Old-looking modules still used transitively include `final_stage2_uncertainty_scan.py`, `stage1_axial_g1_baseline.py`, `stage1_tensor_gb_hard_candidate.py`, `stage1_tensor_gb_soft2p_estimate.py`, `stage2_axial_g1_three_particle.py`, `stage2_tensor_gb_three_particle_estimate.py`, `step10_stage2_f1_split.py`, `step11_tensor_xgamma_ratio_integral.py` and `bs1_stage2_baseline.py`. Retaining a helper module does not make every standalone calculation in that module a preferred prediction.

Three inactive figure bundles are preserved under `archive/inactive_figures_20261019/`. Core scripts and `outputs/` were not moved. Ambiguous historical outputs are retained until their complete file dependencies are established.
