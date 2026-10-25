# Correctness and workspace audit — 19 October 2026

**Later same-day resolution of the leading mixed-condensate coefficient:**
the author requested a direct derivation. Both tree-level d5 contributions
have now been evaluated independently. The current implementation omits the
heavy-line background-gluon contribution to AB and BB (zero in AA). The
complete leading coefficient also differs from the earlier printed expression.
See [the derivation, checks, isolated impact calculation, and proposed patch](../mixed_condensate_resolution/README.md).
This supersedes Section 3's unresolved-coefficient wording at leading d5.
The historical comparison below is retained as provenance; the original full
numerical extraction and its uncertainty calculation remain unidentified.

**Status: the present numerical implementation has passed several checks, but the original 2025 mixing-angle calculation has not yet been reproduced. Do not finalize Referee 1 item 1 or claim negligible input sensitivity from this audit.**

The earlier 30.05-degree and 39.43-degree values came from the current LO+d3+d5 residue matrix with earlier inputs. They were not a rerun of the original angle calculation. This audit found specific differences explaining why that comparison cannot be treated as a direct test of the published angles.

All numerical variants and rerun outputs are in this directory. The production scripts, baseline output tables, proposed reply, and `final_manuscript/` were not modified by these calculations. Documentation and an inactive-figure archive were updated separately.

## 1. Provenance of the two calculations

The current Python two-point scripts and Mathematica notebook take the published angles as external inputs. They implement an exact-mass LO perturbative matrix plus the local quark-condensate and mixed-condensate terms described in the present manuscript. Their computed matrix angle is explicitly diagnostic.

The original 2025 numerical notebook has not been positively identified in the project scripts/notebooks, the surrounding Dsgamma tree, or the additional research locations searched. Searches elsewhere were incomplete inventories of accessible local files, not proof that the original code no longer exists. The user has been asked for another filename or location.

For the printed formula comparison, both the [published paper](https://doi.org/10.1103/PhysRevD.111.014005) and [arXiv v2 source](https://arxiv.org/abs/2408.08014v2) were obtained. Copies, extracted text, and rendered appendix pages are under `references/`. The final PDF confirms the same expressions and labels as the arXiv source.

The author-confirmed normalization `1/(m_Q+m_s)` in the original calculations remains an established fact. No angle rescaling has been applied on account of its omission from the printed current definition.

**Subsequent source-location follow-up:** the user supplied
`/Users/sbilmis/developer/physics/Bc_mixing`. Its Bc perturbative matrix was
checked algebraically and agrees with the current one up to a common factor
after matching flavor roles. It does not contain the required heavy-light
quark-condensate implementation. A nearby original D/B manuscript project was
also found, with plotting routines and older raw curves, but the final angle
calculation is still unidentified. See [the detailed follow-up](bc_folder_followup/README.md).

## 2. Independent leading-order matrix check: passed

`check_perturbative_matrix.py` reconstructs the cut-loop numerator with explicit 4×4 Dirac matrices, the Dirac adjoint at the second vertex, and the two-body phase-space factor. It does not reuse the implemented spectral polynomials when constructing its reference result.

At 24 physical points covering charm/bottom, zero/nonzero light mass, and different distances from threshold, all AA, AB and BB entries agree with the current code. The largest relative matrix discrepancy is 2.01×10^-14. Positivity and AB symmetry hold at those points. The existing FeynCalc spectral audit was also rerun in isolation and agrees algebraically.

Evidence: `perturbative_trace_check.csv`, `perturbative_trace_check.json`, and `isolated_checks/outputs/twopoint_dirac_spectral_audit.txt`.

This validates the tested leading-order implementation, not the complete OPE, the physical current-state ansatz, or the original 2025 uncertainty calculation.

## 3. The printed 2025 appendix needs reconciliation

Let m=m_Q and take m_s=0, with condensates set to zero. In the positive spin-one spectral convention used by the current code, the independently checked expressions are

\[
\rho_{AA}=\frac{(s-m^2)^2(2s+m^2)}{8\pi^2s^2},\qquad
\rho_{AB}=\frac{3(s-m^2)^2}{8\pi^2s},\qquad
\rho_{BB}=\frac{(s-m^2)^2(s+2m^2)}{8\pi^2m^2s}.
\]

Denote the printed right-hand sides of Eqs. (A1) and (A2) by F_A1 and F_A2, including their displayed integral signs. Their perturbative integrands satisfy

\[
F_{A1}^{\rm pert}\leftrightarrow-\rho_{AB},\qquad
F_{A2}^{\rm pert}\leftrightarrow\rho_{AA}-\rho_{BB}.
\]

The left-hand labels in the paper are instead AA−BB and AB, respectively. A common overall transverse-structure sign cannot exchange the off-diagonal combination and the diagonal difference. This is an apparent printed label/sign mismatch. It does **not** identify what the original numerical code evaluated.

To explore this specific lead, `check_printed_appendix.py` transcribes all displayed terms, including gluon condensates, and evaluates them with the earlier central inputs and window midpoints. It records both the literal printed-label angle and a separate hypothesis using `tan(2 theta) = 2 F_A1/F_A2`, motivated by the checked leading-term assignment. No fit or tuning to the quoted angles was performed.

| Sector | M², s0 [GeV²] | Published angle | Current truncated-matrix diagnostic with old inputs | Printed full expressions with the leading-term assignment hypothesis |
| --- | --- | --- | --- | --- |
| Ds | 4.5, 10 | 26.6 ± 0.6 degrees | 30.04965 degrees | 26.73607 degrees |
| Bs | 11, 43 | 38.5 ± 0.1 degrees | 39.43370 degrees | 38.32237 degrees |

These latter values are closer to the published results, but are **not a validated reproduction**. In particular, the Bs midpoint differs from the quoted central value by about 0.178 degrees; a midpoint is not the fitted ensemble central value. The original scale treatment, finite-mass expansion, implementation, and uncertainty extraction still need checking. Numerical integration of the transcribed expressions was stable when increasing Gauss–Legendre nodes from 256 to 512.

There is a further explicit mismatch beyond the labels. Setting m_s and gluon condensates to zero, the mixed-condensate part of printed F_A2 is

\[
-\frac{m\,m_0^2\langle\bar ss\rangle}{M^2}\,e^{-m^2/M^2},
\]

whereas the current code's mixed-condensate contribution to AA−BB is half this expression. Printed F_A1's mixed term agrees with minus the current AB mixed term in that limit. This was checked by subtracting the z=0 result from the z=m0² result and independently evaluating the current local components. See `printed_appendix_massless_d5_check.json`.

This difference is not resolved by relabeling A1 and A2. At this stage of the
audit the coordinate-space derivative alone did not prove completeness of the
dimension-five coefficient. The subsequent independent derivation linked above
identified the second contribution and computed the complete leading result.
Its addition leaves AA unchanged and modifies AB and BB. The original program's
implementation remains a distinct provenance question.

## 4. Current transition numerical implementation: passed

The central Python and Mathematica implementations were rerun into `isolated_checks/outputs/`. All 187 matched quantities pass the established absolute-or-relative 10^-8 criterion; there are no missing keys. The largest absolute difference is 6.58×10^-8 in a bottom width, with relative difference below 5.26×10^-9. The Python rerun exactly reproduces all saved central Python values.

This checks numerical consistency between implementations. It is not an independent proof of the physical transition OPE or the assumed separation of the two hadronic states.

## 5. Confirmed charm sampling inconsistency

`evaluate_ds_physical` in `scripts/final_window_monte_carlo.py` selects f1, f2 and theta from an accepted two-point sample, but uses independently drawn mc, ms and the strange condensate for the transition side. The bottom evaluator explicitly shares those QCD inputs with the selected residue sample.

`check_ds_shared_inputs.py` reproduces the existing preferred 500-sample charm run exactly, including all checked kinematics, residues, transition invariants, couplings and widths. The paired variant then changes only the transition mc, ms and strange condensate to the selected residue sample. The sampled lattice f_gamma,s^perp is held fixed and chi is recomputed as fperp/ss. All other draws, windows, angles and residues are held fixed.

| Width [keV] | Saved baseline median [p16,p84] | Shared-QCD-input variant median [p16,p84] |
| --- | --- | --- |
| Ds1(2460) → Ds γ | 45.9157 [37.4092,56.0037] | 45.5342 [36.9479,56.2219] |
| Ds1(2536) → Ds γ | 0.55034 [0.07087,2.28056] | 0.65247 [0.04959,2.38470] |

The median shifts lie inside the existing intervals. This does not justify retaining inconsistent sampling. The more suppressed channel is visibly sensitive in relative terms. These are paired 500-sample comparisons, not a converged new uncertainty recommendation. Sharing accepted two-point draws also changes the transition-side marginal distribution to the accepted ensemble; the effect should not be described as a pure covariance-only change.

Evidence: `ds_shared_input_check.json`, `ds_shared_input_summary.csv`, the raw draw CSV, and both full paired-output CSVs. No published table or manuscript number has been replaced.

Before a production fix is adopted, preserve the primitive photon input for each scenario: fperp in the lattice scenario, chi in the legacy susceptibility scenario. The current bottom code's unconditional fperp-preserving update also warrants review for the legacy comparison; the preferred lattice calculation already uses the intended choice.

## 6. Workspace organization completed

- Created `CURRENT_WORKFLOW.md` and corrected stale paths and entry-point descriptions in the root and scripts READMEs.
- Recorded static local Python imports and the transitive dependencies of the selected current entry points in `python_dependency_inventory.json`. Some stage-named modules remain active helpers.
- Moved 27 files in three inactive figure bundles to `archive/inactive_figures_20261019/`, preserving their bytes. The archive README and manifest give reasons, original paths and restoration guidance.
- Retained all active scripts, baseline outputs, current figure bundles, ongoing OPE-topology work, submission snapshots and manuscript files. No files were deleted.
- Stored a SHA-256 manifest of the protected manuscript and a broader preorganization manifest for final verification.

Final verification passed: all 11 protected manuscript files are unchanged,
all 27 archived files match their original hashes, and all 21 identified
active Python source paths remain present. Against the 398-file preorganization
manifest, the only content edit is `scripts/README.md`; scientific sources and
baseline outputs are unchanged. See `final_integrity_check.json`.

## Required next steps

1. Obtain the original 2025 numerical implementation, or independently reconstruct and validate its full OPE if that code cannot be recovered. The leading mixed-condensate coefficient has now been derived; verify the remaining labels, finite-mass/higher-order terms, and actual uncertainty procedure without reopening that completed calculation.
2. Only after establishing that baseline, perform Referee 1's old-versus-new input comparison at matched scales and auxiliary parameters. No negligible-sensitivity conclusion is currently justified.
3. Adopt a consistent shared-input sampler in a separate reviewed calculation branch, regenerate affected summaries/figures there, and then decide manuscript changes with the author.
4. Address the hadronic current-state isolation assumption separately. Agreement of matrix arithmetic or Python/Mathematica numerics does not establish it.
5. Keep the Referee 1 response as a draft until these findings are discussed. Keep `final_manuscript/` read-only.

## Reproduction commands

From the project root, the new audit drivers write only into this audit directory:

```bash
python3 review/correctness_audit_20261019/check_perturbative_matrix.py
python3 review/correctness_audit_20261019/check_printed_appendix.py
python3 review/correctness_audit_20261019/check_ds_shared_inputs.py
```

Existing Wolfram trace scripts were copied unchanged into `isolated_checks/scripts/` to redirect their path-derived output directories. The central Mathematica transition script was run by absolute filename with `isolated_checks/` as its working directory. The Python central and comparison modules were imported with only their output-path globals redirected to `isolated_checks/outputs/`.
