# Unequal-Borel sensitivity calculation — 25 October 2026

The author requested this calculation for the referee response only. All
scripts, data, checks, and proposed text belong in this folder. The production
baseline and `final_manuscript/` remain unchanged.

Status: the model-preserving fixed-input scan and numerical checks are complete.
This is an author-review diagnostic. It is not yet a validated general
unequal-Borel extraction or submission-ready referee evidence; the scope and
the manuscript/derivation discrepancy below must be retained with the results.

## Main result

Using the stored off-shell fraction convention and the existing hard spectral
and soft continuum prescriptions, the squared-mass Borel ratio does **not**
produce a negligible shift in the suppressed charm channel.

| Channel | r_mass | g(r=1), GeV^-1 | g(r_mass), GeV^-1 | width(r=1), keV | width(r_mass), keV |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ds1(2460) | 1.561310 | -0.474360 | -0.489020 | 47.2989 | 50.2675 |
| Ds1(2536) | 1.658794 | -0.024265 | -0.129243 | 0.182716 | 5.183620 |
| Bs1(5750) | 1.147851 | +0.312485 | +0.307496 | 12.0623 | 11.6803 |
| Bs1(5830) | 1.179499 | -0.202114 | -0.195630 | 8.66890 | 8.12156 |

These are paired values at the Table II M^2/s0 midpoints, with fixed QCD
inputs and residues from the existing Borel-figure reference. They are NOT
the manuscript's Monte Carlo medians (in particular, 0.182716 is not 0.55).
Both columns use the same Table I meson masses. The baseline-reproduction
check separately uses the original plotting masses, which differ in the
last few digits for the higher states.

At the mass-based choice the Ds1(2536) width spans 1.69265--10.10340 keV
on the 41 x 21 effective-M^2/s0 grid. The symmetric scan spans zero to
2.54841 keV (the signed coupling changes sign). Thus, in this retained
model the suppression relative to Ds1(2460) survives, but the sub-keV
midpoint does not. This is a sensitivity result, not proof of an incorrect
published prediction and not a newly estimated uncertainty interval.

The extended scan covers r in [0.85, 1.15*r_mass], with additional points
between 1 and r_mass and exact benchmark nodes. Each channel has 264 ratio
points and 41 x 21 M^2/s0 points, 227304 points per channel. The modest
r=1.15 benchmark gives Ds1(2536) width 0.657801 keV at the midpoint.
All sampled extrema, not statistical intervals, are in the accompanying CSVs.

## Files

- `benchmark_midpoints.csv`: paired signed couplings, widths, and every OPE component.
- `benchmark_window_ranges.csv`: M^2/s0 envelopes for each of the two choices.
- `ratio_scan.csv`: midpoint curves and M^2/s0 envelopes as functions of r.
- `scan_envelopes.csv`: full sampled ranges and extremum locations.
- `borel_ratio_signed_g2.png` / `.svg`: higher-state coupling plots for both sectors.
- `borel_comparison_author_review.tex`: numerical table, explicitly for author review.
- `fraction_convention_diagnostic.csv`: alternative three-particle argument, not an error band.
- `Ds2460_scan.npz`, `Ds2536_scan.npz`, `Bs5750_scan.npz`, `Bs5830_scan.npz`:
  full component arrays. Axis ordering is documented in the Python builder.
- `inputs.json`, `numerical_checks.json`, `symbolic_checks.json` and the check CSVs:
  reproducibility information and independent controls.

Run `python3 -B review/referee_a_borel_asymmetry_20261025/calculate_borel_sensitivity.py`
and then the adjacent `present_results.py` from the project root. The optional
Wolfram script derives the electromagnetic line-density derivative and checks
the Borel saddle independently. No production entry point is executed.

The test compares `r=M1^2/M2^2=1` with `r=m_A^2/m_P^2`, keeping the effective
`M^2=M1^2*M2^2/(M1^2+M2^2)` and all other inputs paired. The separate pole
exponents must be retained. The completed scan includes both prescriptions,
the Table II effective-Borel and threshold regions, and moderate departures
on either side of the symmetric point.

## Fraction convention and what remains to establish

The two-particle stored kernel uses D=mQ^2-a*pPrime^2-(1-a)*p^2. Its Schwinger
constraints are tau*a=1/M1^2 and tau*(1-a)=1/M2^2, hence a=1/(1+r).
The stored three-particle derivation defines its heavy-line fraction by
a=alpha_qbar+v*alpha_g and uses the same D. Its line projection therefore
has z=a. This is the primary branch in the new calculation.

The older derivation note calls this fraction a0, whereas Eq. (24) of
r1-2.pdf calls it u0. The older note's u0 is its complement. This naming
change must not silently reflect only one sector of the numerical code.

Two qualifications were identified during inspection:

1. The manuscript r1-2.pdf uses `u0=1/(1+r)` in its two-particle terms,
   but `z=1-u0` in its three-particle terms. The derivation scripts instead
   use a single heavy-line fraction `a` in both sectors. Equal Borel
   parameters conceal this difference. The argument conventions must be
   reconciled before quoting a general-Borel formula in the reply. The new
   primary branch follows the common fraction in the stored off-shell
   derivation. Reflecting only the three-particle argument gives 0.053050 keV
   instead of 5.183620 keV for the higher charm channel at r_mass. Both
   reproduce the equal-Borel number, so that check alone cannot choose them.
   The alternative is saved only to expose the issue; it is not a second
   physical prediction or an uncertainty estimate.
2. The production hard tensor density uses a diagonal double-dispersion
   model. Holding it fixed at fixed effective M^2 is mathematically consistent
   with that model, but is not a new derivation of the full unequal-Borel
   hard triangle. The soft-sector continuum subtraction is also retained
   as the production prescription, not independently validated by this scan.

There is also a DA-definition detail to reconcile with the printed appendix:
the shared numerical T_i and T4gamma functions are written with first minus
second argument, while the manuscript displays alpha_sbar-alpha_s. A relabelling
requires consistent changes to the support constraint and integration variables.
No sign or DA parameter has been changed in this calculation, and the
alternative-argument diagnostic is not a literal transcription of all printed
Appendix A/B formulas. Do not infer a physical sign correction from this
observation alone. The manuscript's B(u) coefficient also differs from the
production BBK dictionary (3 versus 2 multiplying kappa-plus); this scan keeps
the established production dictionary, as recorded in the existing source.

At fixed effective M^2 the hard term is independent of r within its assumed
diagonal density rho(s1)*delta(s1-s2), because the weight is exp(-s/M^2).
This establishes the implementation of that model, not the absence of other
double-spectral terms in a complete hard calculation. No separate final-channel
continuum threshold or general two-variable duality domain has been fitted.
The individual-channel OPE/continuum admissibility of the mass-based Borel
pair has not been established by merely retaining the Table II effective window.

## Why the suppressed amplitude moves

At the Ds1(2536) midpoint, the physical numerator Pi_2 changes from
-0.001844197 to -0.010687041 GeV^3. The perturbative contribution stays
-0.035738745 GeV^3 under the diagonal model. The leading-twist contribution
decreases from +0.028127806 to +0.026400913, the twist-three contribution
from +0.002054102 to +0.000007156, and the electromagnetic three-particle
term changes from +0.001098517 to -0.004068427 GeV^3. The separate pole
normalization decreases by about 8.1%; it does not cause the enhancement.
The effect comes from changing the small residual of several contributions.

## Verification

All 372 saved symmetric Borel points agree: maximum coupling difference
8.44e-15 GeV^-1 and width difference 1.65e-12 keV. Two independent
three-particle integration parameterizations agree to 2.17e-9; increasing
the original quadrature order changes results by at most 4.27e-14. The
analytic electromagnetic line-density derivative agrees with direct numerical
differentiation to 5.12e-13. Raising the hard quadrature from 1600 to 2400
nodes changes the checked invariants by at most 9.44e-12 GeV^3.
Four independent symbolic Borel-kinematic identities vanish exactly.
These checks verify the numerical evaluation of the selected expressions;
they do not remove the physical/prescription qualifications above.

All protected manuscript, production source and selected baseline hashes
remain unchanged. No full Monte Carlo ensemble, decay-constant fit, DA-shape
scan, twist-four parameter scan, or manuscript edit was performed.

Do not promote a model-preserving diagnostic into a complete validation of
unequal-Borel sum rules or into a new statistical uncertainty interval.
