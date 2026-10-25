# Referee 1, item 1: condensate-only comparison

19 October 2026. **Completed conditional numerical check; original-program
validation remains open.** The author requested this comparison for the referee
response only. No manuscript or published input value has been replaced.

**Later same-day update:** the leading mixed-condensate issue has been derived
independently in [mixed_condensate_resolution/](../mixed_condensate_resolution/README.md).
The current code omits the heavy-line contribution to AB and BB; the complete
leading coefficient also differs from the printed appendix used here. The
tables below preserve the original conditional reconstruction, not a result
validated by the new derivation. Keep them distinct from the completed
current-OPE sensitivity calculation in the new directory.

## Result and its scope

Keeping all other inputs fixed, changing the central light-quark condensate
from -0.0165 to -(0.240)^3 = -0.013824 GeV cubed produces the following shifts
in the previously documented reconstruction of the earlier appendix:

| Sector | M², s0 [GeV²] | Earlier-input angle | Updated-condensate angle | Shift | Quoted angle error | Shift/error |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Ds | 4.5, 10 | 26.736070° | 26.921234° | +0.185164° | 0.6° | 0.309 |
| Bs | 11, 43 | 38.322373° | 38.379469° | +0.057096° | 0.1° | 0.571 |

The old and new entries are reconstructed midpoint values, **not** replacements
for the published Monte Carlo central values 26.6° and 38.5°. Do not add these
shifts to the published values and call the sums a new extraction.

A 21 × 21 grid for each sector, including endpoints and midpoints, gives:

| Sector | M² range [GeV²] | s0 range [GeV²] | Sampled shift range | Range/quoted error |
| --- | --- | --- | --- | --- |
| Ds | 3–6 | 9–11 | +0.118680° to +0.368437° | 0.198–0.614 |
| Bs | 9–13 | 42–44 | +0.044416° to +0.079634° | 0.444–0.796 |

The shifts decrease with both auxiliary parameters on this grid. The largest
occurs at the lower/lower corner and the smallest at the upper/upper corner.
These are sampled extrema, not a proof of continuous bounds. The window
spread is not a statistical uncertainty.

The supported description is **smaller than the quoted angle uncertainties
throughout the sampled windows**. Calling the effect simply negligible would
conceal its size relative to the narrow bottom uncertainty. This is a central
input sensitivity check, not a fresh propagation of all uncertainties or a
validation of the published uncertainty estimates.

## Inputs held fixed

Both condensate values refer to 1 GeV. The magnitude decreases by 16.2182%.
Use kappa_s = 0.8 and m0² = 0.8 GeV² for both calculations, so that:

| Condensate | Earlier | Updated |
| --- | ---: | ---: |
| <qbar q> [GeV³] | -0.0165 | -0.013824 |
| <sbar s> = 0.8 <qbar q> [GeV³] | -0.0132 | -0.0110592 |
| <sbar g sigma.G s> = m0² <sbar s> [GeV⁵] | -0.01056 | -0.00884736 |

All terms containing the strange condensate, including mixed and products
with gluon condensates, change consistently. Fixed values are mc = 1.27 GeV,
mb = 4.18 GeV, ms(1 GeV) = 0.126 GeV, <g²G²> = 4 pi² × 0.012 GeV⁴,
and <g³G³> = 0.57 GeV⁶. In particular, ms is not independently changed to its
2 GeV value. Neither a new kappa_s error nor new heavy-quark errors are
propagated in this comparison.

The source of the earlier inputs, formulas, and windows is
[PRD 111, 014005 / arXiv:2408.08014v2](https://arxiv.org/html/2408.08014v2).
Local source/PDF hashes are in `comparison_summary.json`. The current
condensate value is from the author-designated `/Users/sbilmis/Downloads/r1.pdf`,
as recorded in the revision notes.

## Reconstruction limitation

Let F_A1 and F_A2 denote the printed right-hand sides of the earlier appendix.
This check uses tan(2 theta) = 2 F_A1/F_A2, with theta on the 0–90° branch,
following the assignment F_A1 = -AB and F_A2 = AA-BB identified from the
independently checked leading perturbative terms. This is the same explicit
hypothesis used in the [correctness audit](../correctness_audit_20261019/REPORT.md).
No coefficients were fitted or tuned to the published angles.

The original Ds/Bs calculation code has not been recovered. A subsequent
independent derivation has resolved the leading mixed-condensate coefficient
and found an omitted heavy-line contribution in the current implementation;
the complete result also differs from this printed-appendix reconstruction.
Numerical agreement between the transcriptions below establishes arithmetic
consistency, not completeness of the printed OPE or equivalence to the original
program. The response draft remains provisional until the prescription for the
reported comparison is selected and validated. The author-confirmed tensor
normalization is not in question; no normalization conversion is applied here.

As a separate historical diagnostic, the uncompleted current LO+d3+d5 matrix gave
midpoint shifts +0.177103° (Ds) and +0.044706° (Bs). Across nine endpoint/midpoint
locations its shifts span 0.117107–0.335123° and 0.035968–0.059038°, respectively.
After the leading-d5 completion, midpoint shifts in that current OPE are
+0.174981 degrees and +0.045126 degrees; the new directory records the full
21 by 21 grid. These checks do not replace the missing original baseline or
define a systematic error band.
The current diagnostic angles are not included in the referee draft.

## Numerical verification and mechanism

- A separate Wolfram transcription groups the expressions by powers of the
  condensate and evaluates the integrals with 50-digit working precision.
  All 126 compared values at 18 points agree: maximum angle difference
  3.20 × 10^-14 degrees, maximum shift difference 1.57 × 10^-14 degrees.
- Increasing Python quadrature from 256 to 512 nodes changes the checked
  angles/shifts by at most 6.40 × 10^-14 degrees.
- The affine condensate dependence and the branch-consistent direct angle
  difference were checked. All 11 protected manuscript files and input source
  hashes are unchanged.

`signed_midpoint_contributions.csv` separates perturbative, q-independent
gluon, and all strange-condensate-dependent terms. At the midpoints the
perturbative terms dominate both reconstructed combinations. The condensate
update changes the off-diagonal magnitude more than the diagonal difference;
the ratio therefore changes slightly rather than cancelling exactly.

`qq_only_endpoint_envelopes.csv` additionally varies only the quoted light
condensate interval at fixed kappa_s. These are endpoint envelopes, not total
one-sigma angle intervals. They must not be combined with the quoted errors
or interpreted as having already propagated the changed kappa_s uncertainty.

## Files and reproduction

- `condensate_shift_grid.csv`: all 882 paired points.
- `midpoint_comparison.csv`, `comparison_summary.json`: central results and provenance.
- `current_truncation_crosscheck.csv`: separate current-OPE diagnostic.
- `verify_condensate_inputs.wl`, `wolfram_high_precision.json`: independent evaluation.
- `verify_numerical_results.py`, `verification.json`: cross-language checks.
- `referee_response_condensate_check.tex`: provisional response-only wording.

From the project root:

```bash
python3 -B review/referee1_item1/compare_condensate_inputs.py
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noinit -noprompt -script review/referee1_item1/verify_condensate_inputs.wl
python3 -B review/referee1_item1/verify_numerical_results.py
```

These commands write only into this comparison directory. The existing full
referee-response PDF has not been rebuilt or silently finalized with these
conditional results.
