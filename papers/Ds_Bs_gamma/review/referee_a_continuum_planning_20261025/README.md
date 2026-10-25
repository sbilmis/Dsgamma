# Continuum discussion and recalculation roadmap — 25 October 2026

The author asked whether the unequal-Borel diagnostic requires manuscript
changes, requested perturbative/total ratios at the working points, and
asked to discuss the continuum-subtraction comment before further width
recalculations. This package contains working-point diagnostics and a
proposed roadmap. No new continuum variant, width ensemble, manuscript
change, or production patch has been adopted or run here.

## Decision about the manuscript

Do not replace a manuscript prediction by the provisional 5.18362-keV
fixed-input diagnostic. The general-Borel convention and continuum treatment
remain unvalidated. Conversely, if a consistently derived calculation in
admissible windows confirms a material shift or larger systematic uncertainty,
the affected manuscript numbers, errors and claims must reflect it. Such a
change cannot be confined to the reply letter. This need not mean changing
every central value; the outcome determines which values and uncertainties
change. The earlier preference for response-only robustness checks does not
justify retaining an unsupported precision claim.

## Two different ratios

For each invariant, `perturbative_to_total.csv` records the SIGNED ratio
`Pi_pert/Pi_total`, and separately an absolute-component measure. The physical
invariant is the appropriate mixed combination, using exactly the saved
midpoint decompositions from the previous diagnostic.

Each cell below is symmetric -> mass-based, with M^2 and s0 in GeV^2:

| Channel | M^2 | s0 | pert/total A | pert/total B | pert/total physical |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ds2460 | 3.75 | 9.0 | 2.486 -> 1.970 | 1.563 -> 1.478 | 1.750 -> 1.593 |
| Ds2536 | 3.75 | 9.5 | 2.396 -> 1.872 | 1.549 -> 1.469 | 19.379 -> 3.344 |
| Bs5750 | 12.0 | 40.0 | -0.118 -> -0.118 | 0.310 -> 0.312 | 0.166 -> 0.167 |
| Bs5830 | 12.0 | 41.0 | -0.108 -> -0.109 | 0.312 -> 0.314 | 2.056 -> 2.087 |

Ratios above one or negative are possible because these terms have signs.
They are not pole fractions or probabilities, and a perturbative fraction
alone does not establish light-cone OPE convergence. The leading photon DA
is part of the leading light-cone approximation even though it is nonperturbative.

For continuum diagnostics the signed ratio of hard integrals below s0 and
to infinity is also saved. In the bottom sector it becomes negative or
exceeds one because the hard spectral density itself changes sign. It is
not a useful positive pole probability. The bounded absolute-weight diagnostic is

    R_abs = Integral_sth^s0 exp(-s/M^2) abs(rho(s)) ds
            / Integral_sth^infinity exp(-s/M^2) abs(rho(s)) ds.

Its midpoint values are:

| Channel | R_abs,A | R_abs,B |
| --- | ---: | ---: |
| Ds2460 | 0.69353 | 0.77054 |
| Ds2536 | 0.72584 | 0.79927 |
| Bs5750 | 0.51309 | 0.78085 |
| Bs5830 | 0.53135 | 0.79760 |

Over the Table II windows the corresponding A/B ranges are:

- Ds2460: 0.57502--0.81656 / 0.67130--0.86767.
- Ds2536: 0.61129--0.84126 / 0.70595--0.88803.
- Bs5750: 0.39639--0.65644 / 0.67369--0.87384.
- Bs5830: 0.41320--0.67408 / 0.69289--0.88674.

For a nonnegative absolute density the extremal ratios lie at the stated
opposite corners of M^2/s0. They were recomputed by integration split at the
spectral zeros; controls are in `numerical_check.json`. These quantities describe
only the existing diagonal perturbative spectral model. They are independent
of r at fixed effective M^2 in that model. They do NOT validate the full
two-variable continuum domain, the soft continuum prescription, or the
individual axial and pseudoscalar Borel windows.

## What the continuum comment requires

The referee's specific comment 1 asks for the duality domain in both dispersion
variables, a derivation including boundary terms, comparison with Ref. [25],
and a quantitative effect on the suppressed channel. Merely naming duality
or citing a prescription does not supply this derivation.

Verified primary sources:

- Colangelo, De Fazio, Ozpineci, hep-ph/0505195, Eq. (4.4): the leading-twist
  term has E_Q-E_0; its higher-twist terms carry E_Q without that subtraction.
  https://arxiv.org/pdf/hep-ph/0505195
- Rohrwild, arXiv:0708.1405, Eq. (20): its displayed two-particle soft block
  carries E_Q-E_0, including higher-twist terms. The surrounding discussion
  explicitly ties the subtraction prescription to the accompanying sum rules.
  This is precedent, not a derivation for the present axial/tensor correlator.
  https://arxiv.org/pdf/0708.1405
- Kim, Lee, Oka, hep-ph/0011072, Sec. 2: subtraction terms and a truncated
  double dispersion can generate spurious continuum pieces. It supports
  checking contact/surface terms; it does not by itself choose the prescription
  for this correlator. https://arxiv.org/pdf/hep-ph/0011072

Illustration only (NOT the requested full double-dispersion derivation):
the one-variable Borel image M^2 exp(-m^2/M^2) can come from a constant
density above m^2, whose cutoff gives M^2(E_m-E_0). An exp(-m^2/M^2)
image can instead come from delta(s-m^2); if the pole is inside the cutoff,
it is unchanged. Therefore an exponential appearing in two terms does not
establish an identical subtraction rule. The actual decision must be made
from each term's two-variable distribution and duality domain.

The earlier local `main3_savci_222_finalization_review.md` contains categorical
statements that the uniform soft subtraction is correct. Agreement between
numerical implementations verifies the chosen rule, not its double-dispersion
derivation. Treat that older methodological conclusion as unestablished until
the derivation below is completed; no replacement prescription is adopted yet.

## Proposed dependency order for recalculation

1. **Freeze one consistent correlator and DA convention.** Align current
   positions, Fourier phases, p/pPrime, tensor momentum, two-particle fractions,
   three-particle support, antisymmetric DA orientation, integrated vector DA,
   the B(u) coefficient, and the B4 auxiliary integral. Derive which apparently
   omitted DA projections actually vanish. Supply the nonlocal definitions.
   This joins specific comments 3, 5 and 6 with the asymmetry prerequisite.
   Deliver a dictionary and tests that work away from u0=1/2.

2. **Derive continuum subtraction at general Borel parameters.** Obtain the
   hard and soft double-spectral distributions, including delta derivatives
   and single-variable/contact pieces. State the hadronic duality domain
   explicitly, e.g. a rectangle with channel thresholds or a justified averaged
   domain. Determine its relation to the Table II single effective s0.
   Retain surface terms when integrating by parts. Check the uncut limit and
   symmetric limit. Do not assume that the existing tensor diagonal hard mapping
   exhausts the general density. This is specific comment 1 and part of major ii.
   Start with I_n[F]=Integral_0^1 da F(a)/[mQ^2-a*pPrime^2-(1-a)*p^2]^n
   for n=1,2,3, including the numerator p.q in the vector/tensor kernels.
   Recover the raw Schwinger/Borel localization first; then impose the selected
   spectral cutoff, keeping the DA polynomial components and their surface terms.
   Extend the verified operation to the three-particle line functionals.

3. **Resolve inputs and physical-state extraction before the production run.**
   Check the finite-ms photon normalization/subtraction scheme and scale
   evolution (specific comment 2) together with the hard/soft separation.
   Include the already derived leading-d5 two-point correction in the proposed
   regenerated variant, with its stated approximation, before adopting it.
   Check the off-diagonal two-point correlator/pole separation (major i); retain
   both poles if the one-pole identification is not justified. The tensor-current
   normalization typo is settled and must not be reopened.

4. **Run controlled fixed-input comparisons and reassess the existing windows.**
   Compare the old uniform-soft subtraction, a leading-twist-only-subtracted
   sensitivity variant, and the derived prescription, keeping other terms and
   inputs identical. This isolates continuum sensitivity; do not simultaneously
   add the alternative local-condensate treatment of another paper. Separately
   benchmark the pure axial transition kernel against Eq. (4.4) with common
   inputs and explicitly matched conventions (specific comment 7). Set M^2,
   thresholds and allowed Borel-ratio ranges from continuum/OPE and pole tests,
   before favoring any width. Check both basis and physical combinations.

5. **Propagate uncertainties coherently.** Use the same sampled QCD inputs and
   angle for residues and transitions, preserving their correlations. The current
   charm sampling issue is already documented; its repair is part of this run.
   Existing angle plots already recompute residues with theta. Re-evaluate them
   with the final kernels and inputs. Vary a normalized leading-twist DA shape,
   independent twist-four parameters with their constraints, and the allowed
   Borel asymmetry. Resolve the rounded twist-four EOM-constraint issue first.
   Assess omitted radiative/higher-twist contributions on an absolute amplitude
   scale and compare them with the residual; an estimate is not a calculation
   of an uncomputed order. Separate input, auxiliary-parameter and prescription
   uncertainties and do not assign unjustified Gaussian errors to cancellations.

6. **Update the response and affected manuscript results together.** Present
   signed components, window criteria, a compact uncertainty budget, and the
   difference between hierarchy robustness and precision of the small width.
   Regenerate all affected tables/figures consistently in both sectors. Proposed
   edits stay outside final_manuscript for author approval and Prism insertion.

The next substantive task is the combined convention/continuum derivation,
not another large numerical ensemble. The roadmap is a recommendation under
discussion; this turn does not launch those recalculations or adopt patches.

## Pullin--Zwicky prescription review, later on 25 October 2026

The author decided to keep the two Borel parameters independent for the
forthcoming analysis and requested discussion/agreement on the continuum
formula BEFORE recalculation. The published source is accessible:
https://doi.org/10.1007/JHEP09(2021)023, arXiv:2106.13617v3.
The published PDF and a text extraction are saved in `references/`.
No named copy was found in the searched project/Downloads paths before
downloading the published PDF from INSPIRE. No author upload is needed.

Verified source locations: Sec. 3.4.1, Eqs. (3.17)--(3.19); Appendix C,
Eqs. (C.3)--(C.4); Appendix D, Eqs. (D.1)--(D.18). The general expressions
retain independent Borel parameters; Eq. (3.12) is the authors' numerical
equal-Borel choice, not a restriction on the generic construction.

The duality region is s,t >= mQ^2 with (s/S)^a+(t/T)^a <= 1, and its
intersection with s=t is sigma=S*T/(S^a+T^a)^(1/a). For a=1 the upper
boundary is t=T*(1-s/S) and sigma=S*T/(S+T). Boundary terms from the
derivatives of delta(t-s) must be retained. D.11 alone is not the subtracted
answer; D.14--D.15 supply its boundary contribution. In the symmetric
linear-domain limit, D.18 gives

    Fhat_n,0 = K^(2-n)/Gamma(n) f_n(1/2)
               * [exp(-mQ^2/K)-delta_n1*exp(-sigma/K)].

This statement concerns the stated class of F_n,0 kernels, not all terms
of a particular twist. Terms with a numerator virtuality, F_n,1, have a
different rule. The polynomial representation also has endpoint assumptions;
logarithmic terms require a justified treatment (see source footnote 16).

Our raw n=2 double-Borel maps match D.16 after exchanging the hadronic
virtuality labels. Our production assembly does NOT yet implement the full
PZ boundary prescription: `scripts/rohrwild_transition_exact.py` uses the
common delta_exp for every soft component. This is a source/code comparison,
not a new width determination or a claim that all previous results are wrong.

Notation: paper p_B is the pseudoscalar momentum and its q is the axial/vector
hadronic momentum, NOT our real-photon momentum. Thus paper M1^2=M_P^2,
M2^2=M_A^2. For our kernel mQ^2-x*P_A^2-(1-x)*p_P^2,
x0=M_P^2/(M_A^2+M_P^2), while paper u0=1-x0. All DA and support relabelling
must be performed consistently, especially for antisymmetric functions.

### Independent boundary-formula controls

`check_pz_boundary_formulas.wl` and `pz_boundary_checks.json` are analytic
test-kernel checks only. For a linear boundary and independent inverse
Borel parameters, exact distribution differentiation agrees with printed
C.4 for n=1,2,3,4; all four symbolic residuals vanish. The test kernel
(s-1)^4 exp(-x*s-y*t) removes lower-endpoint complications to those orders.

There is an apparent defect in using printed C.4 literally for curved
domains at higher delta-derivative order. An explicit control uses n=3,
rho=(s-1)^2, g=1, s,t>=1, s^2+t^2<=8. Set v=t-s; the upper endpoint is
S(v)=(-v+sqrt(16-v^2))/2. Its diagonal value is 2 and S''(0)=-1/8.
The direct distribution integral is G''(0)/Gamma(3)=3/16, with
G(v)=[(S(v)-1)^3-(max(1,1-v)-1)^3]/3. Printed C.4 gives 1/4; the missing
curvature term is -1/16. Exact symbolic algebra and finite-difference
evaluation agree. This check does NOT invalidate the linear a=1 formula.
Do not claim that the printed general-a boundary formula has been verified
for arbitrary curved domains or use it for an a-scan without repairing this.

### Recommendation for discussion, not yet adopted

Use the linear a=1 duality domain, keeping both Borel parameters general.
Start from the EXISTING accepted/adopted Table II effective-M^2 and s0
windows; identify sigma with their effective s0 for a controlled comparison.
The ratio T/S is an additional domain choice not fixed by that single s0.
A concrete candidate is T/S=m_A^2/m_P^2, as a modelling prescription, not
an identity. At the mass-proportional Borel choice this also makes the Borel
exponent constant along the linear boundary. For Borel-ratio variations,
keep the selected domain fixed so that Borel and duality variations can be
distinguished. Check any chosen domain against the pole/daughter and OPE
diagnostics before promoting it to the final numerical prescription.

These recommendations remain for author discussion. No continuum variant,
production source, residue fit, uncertainty ensemble, or manuscript has
been changed by this source review.
