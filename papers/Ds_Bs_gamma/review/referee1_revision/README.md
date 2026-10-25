# Referee 1 revision package

Prepared 19 October 2026 using `/Users/sbilmis/Downloads/r1.pdf` as the manuscript reference. These are proposed changes for the author to insert into Prism. No manuscript source or existing figure has been edited.

The two figures and the reply sections for comments 2 and 3 are prepared. The reply to comment 1 includes the method and normalization clarification, but its quantitative input-sensitivity conclusion is **pending**. The reply PDF explicitly marks that unfinished part and is an author-review draft, not a submission-ready response.

## What to upload and where to insert it

Page, equation and figure numbers below refer to the original `r1.pdf`; they will change after insertion. Use the quoted text anchors to find the right locations in Prism.

| Proposed change | Location in `r1.pdf` | File to use |
| --- | --- | --- |
| Mixing-angle method and normalization footnote | Page 11, immediately after the paragraph following Eqs. (37)-(38) ending “their expressions are presented in Appendix C.” Insert before “Substituting the residues obtained from the two-point sum rules...” | `manuscript_item1_method.tex` |
| New individual-decay-constant figure, caption and discussion | Page 12, after the paragraph giving the four decay constants and the theta = 90-degree check, ending “compatible with those reported in [34].” Insert before “Because the initial and final mesons in each channel...” | `manuscript_item3_figure.tex` and `figures/decay_constants_vs_theta.pdf` |
| Clarify the existing width-versus-angle figure description | Pages 13-14: replace the sentence beginning “The remaining inputs are fixed at their central values in this figure...” and the final sentence of the original Figure 2 caption | `manuscript_figure2_clarification.tex` contains the two separate replacements |
| Vector-final-state outlook | Page 17, end of Conclusion, after “Searches for the lower bottom-strange 1+ state would provide an additional test of this picture.” Before Appendix A | `manuscript_item2_outlook.tex` |

Upload `decay_constants_vs_theta.pdf` to Prism under the name used in the LaTeX snippet, or adjust the image path. The new figure uses label `fig:decay-constants-theta`. Let LaTeX renumber the figures automatically; the old Figure 2 retains its existing label and image. Do not add the ratio figure to the manuscript.

The key `Aliev:2025mixing` is present in the local bibliography, but the live Prism bibliography has not been inspected. Retain the existing Prism citation key for Ref. [29] if it differs. In the Figure 2 replacement text, preserve the existing label-based reference to Table III where available. The new equations have unique `eq:r1-...` labels.

If the earlier method paragraph has already been inserted in Prism, replace or reconcile it with the supplied method snippet rather than duplicating it.

## Figures and numerical meaning

- **Manuscript:** `figures/decay_constants_vs_theta.pdf`, with PNG preview. Two panels, each showing f1 and f2, for charm-strange and bottom-strange mesons.
- **Response only:** `figures/decay_constant_ratio_vs_theta.pdf`, with PNG preview. It shows **f2/f1 only**, following the user's final correction. Both sectors use the same representative samples as the individual-constant plot.
- **Reply:** `referee1_reply.tex` and its compiled PDF. The ratio figure is included. The wording about manuscript changes assumes the author has applied the accompanying snippets before submission.

These are illustrative angle scans from 25 to 45 degrees. QCD inputs, masses, two-point Borel parameters and thresholds are fixed; the residues are recomputed at every angle. No uncertainty bands or new statistical ratio estimates are claimed. A fixed representative sample need not reproduce the marginal median of every quantity.

At the preferred angles, the plotted representative values are:

| Sector | Angle (degrees) | f1 (MeV) | f2 (MeV) | f2/f1 |
| --- | ---: | ---: | ---: | ---: |
| Ds | 26.6 | 404.997 | 169.094 | 0.417520 |
| Bs | 38.5 | 547.834 | 91.094 | 0.166280 |

They do not replace the manuscript's propagated decay-constant estimates. In particular, the representative Bs f1 differs from the quoted median 536 MeV while lying within its quoted interval. The chosen inputs reproduce the existing Figure 2 exactly; they were not rescaled to make the curves pass through the quoted medians.

The fixed two-point auxiliary values are:

| Sector | M² (GeV²) | s01 (GeV²) | s02 (GeV²) |
| --- | ---: | ---: | ---: |
| Ds | 2.38020685125 | 7.950 | 10.025 |
| Bs | 6.95283998335 | 44.750 | 43.225 |

These are **two-point residue inputs**, not the transition windows in Table II. Full selected input rows, source hashes and validation results are stored in `provenance.json`.

## Why the referee asks for the ratio

The referee's motivation is an application to the relative production of the two charm-strange axial states in factorized B decays. The weak amplitude involves the axial-current decay constants. The manuscript's f1 and f2 are defined using rotated interpolating currents; they must not be identified unconditionally with weak-current constants.

Under the adopted exclusive-current-coupling ansatz, the inverse rotation gives fA,1 = sin(theta) f1 and fA,2 = cos(theta) f2. Thus the weak-current ratio is cot(theta) f2/f1. Phase-space and B-to-D form-factor factors are also needed for a branching-fraction ratio. The response makes this distinction but does not perform a nonleptonic branching-fraction calculation. The conditional relation does not establish the state-isolation assumption raised by Referee 2.

The experimental reference is the [LHCb amplitude analysis, arXiv:2411.03399](https://arxiv.org/html/2411.03399v3). Its invariant-mass fits contain both axial-state components. We do not independently assert that this particular article measures both absolute secondary branching fractions; that statement appears only in the referee's reproduced comment.

## What remains for comment 1

See `input_check_status.md`. The current scripts take the published mixing angles as external inputs. Their optional truncated-matrix diagnostic has not been validated as the original 2025 calculation. Its numerical output must not be reported as a recomputation of the published angles or as evidence that the earlier paper is incorrect.

The established normalization correction is already settled: the earlier calculations included 1/(mQ+ms), and only the printed definition omitted it. This is separate from input sensitivity.

## Reproduction and checks

Run from the project root:

```bash
python3 review/referee1_revision/build_decay_constant_figures.py
```

The driver selects the same accepted representative samples as the existing Figure 2 builder, recomputes both residues for all 161 angles per sector, and compares them with both state rows in the saved Figure 2 CSV. All 1,288 comparisons agree to better than 1e-12 GeV; the observed maximum difference was 1.11e-16 GeV. All plotted residues are positive and finite. The preferred-angle benchmark values are directly evaluated, not interpolated from the 0.125-degree grid.

Compile the reply from this directory with:

```bash
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build referee1_reply.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build referee1_reply.tex
```

The `build` directory contains temporary compiler and visual-check files. The delivered PDF is copied to this directory after compilation. Both plots are vector PDFs with embedded fonts. `decay_constants_theta.csv` and `benchmark_values.csv` provide the numbers for independent reuse.

Checksums confirm that `final_manuscript/` remained unchanged during figure generation; the complete manifest is recorded in `provenance.json`.
