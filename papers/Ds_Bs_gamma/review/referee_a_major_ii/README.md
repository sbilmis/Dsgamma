# Referee A, major comment (ii): signed contribution figure

## Mass-proportional Borel prescription — latest discussion

The author accepted the revised Pi_A/Pi_B/Pi_2 plot and asked to verify
M1^2/M2^2=m_A^2/m_P^2 for the unequal-Borel test, with its numerical results
intended for the referee response rather than a new manuscript table/figure.
This is a motivated auxiliary-parameter prescription, not an exact QCD
constraint. J. Rohrwild, Phys. Rev. D 75, 074025 (2007),
https://arxiv.org/pdf/hep-ph/0701085 (v3), Eq. (49), explicitly uses
M1^2/M2^2=m_Delta^2/m_N^2; the subsequent text and Figures 6--7 examine
variations of the ratio. Do not cite Eq. (50) from the older v1 numbering.

With the manuscript's assignment M1 to the axial channel and M2 to the
pseudoscalar channel, the prescription equates m_A^2/M1^2 and m_P^2/M2^2.
Keeping the effective M^2 fixed gives M1^2=(1+r)M^2,
M2^2=(1+r)M^2/r, and u0=1/(1+r). Using Table I central masses:

| Channel | r_mass=m_A^2/m_P^2 | u0 |
| --- | ---: | ---: |
| Ds1(2460) -> Ds gamma | 1.561309509 | 0.390425287 |
| Ds1(2536) -> Ds gamma | 1.658793763 | 0.376110405 |
| Bs1(5750) -> Bs gamma | 1.147850813 | 0.465581685 |
| Bs1(5830) -> Bs gamma | 1.179499153 | 0.458821009 |

These are kinematic ratios only; no unequal-Borel OPE or widths have been
calculated. The Table I masses differ slightly in their last quoted digits
from the preserved plotting reference; use paired identical inputs in the
eventual comparison. A response-only table comparing r=1, r=r_mass, and
an explicitly defined ratio-scan envelope is proposed. Evaluate both
prescriptions over the working M^2 and threshold regions and check the
validity of the individual Borel channels. A second fixed ratio alone is
a benchmark, not a demonstration of stability over a region. The DA-shape
and twist-four-input tests remain separate. Do not predetermine that all
effects are small or that the manuscript conclusions remain unchanged.

## Latest figure choice, later on 25 October

The author proposed replacing g2 in the third panel by
`Pi_2 = cos(theta) Pi_A - sin(theta) Pi_B` and matching the existing
manuscript plot style. The revised figure is
[ds2536_signed_invariants.png](ds2536_signed_invariants.png), with a
[vector SVG](ds2536_signed_invariants.svg),
[proposed LaTeX caption](ds2536_signed_invariants_caption.tex), and
[midpoint table in GeV^3](ds2536_signed_invariants_midpoint.csv).
Reproduce it with `python3 -B review/referee_a_major_ii/plot_signed_invariants.py`.
This reads the unchanged saved data and does not repeat the OPE evaluation.
The earlier g2 figure below remains available for comparison.

All three revised panels are in GeV^3. The physical numerator is plotted
directly, not obtained by relabelling the g2 axis. Since the normalization
between Pi_2 and g2 depends on M^2, their Borel shapes differ, while their
relative component sizes at any fixed point are unchanged. The midpoint
Pi_2 is -0.0018441970136229627 GeV^3. Algebraic projection, component-sum,
and g2-normalization checks pass to roundoff. Data and manuscript hashes are
unchanged. The PNG was visually inspected against manuscript Figure 2:
sans-serif math, matching blue/red/black colours, boxed axes, inward ticks,
no grid, and compact labels; explanatory details are in the caption.

The proposed presentation for the remaining stability checks is a compact
table of variation ranges and the resulting widths in both sectors, plus
a focused plot of signed g2 versus the Borel ratio if needed to show motion
of its zero. The DA shape, twist-four, and asymmetric-Borel scans have NOT
been performed. Use a normalized Gegenbauer deformation for leading twist;
choose justified ranges and keep normalization fixed for the shape-only
test. Vary independent twist-four parameters coherently across all affected
DAs. BBK (hep-ph/0207307), Eq. (4.32), constrains
zeta1+11*zeta2-2*zeta2_plus=7/2. Before selecting an EOM-consistent scan,
reconcile the rounded baseline values (0.4, 0.3, 0), which give 3.7,
with that exact constraint; no parameter was changed in these plots.

For unequal Borel parameters, retain the effective M^2 and the manuscript's
actual convention in Eq. (24): u0=M2^2/(M1^2+M2^2), a0=1-u0. With
r=M1^2/M2^2, M1^2=(1+r)*M^2, M2^2=(1+r)*M^2/r, and u0=1/(1+r).
The full kernels, continuum treatment, DA arguments, and separate hadronic
exponents must be checked consistently before varying r. Changing only the
current code's a0 argument is not a validated asymmetric-Borel calculation.
No proposed sensitivity envelope should be labelled a new 68% interval.

25 October 2026. The author requested a three-panel figure after discussing
the distinction between convergence of the basis amplitudes and the size of
the suppressed physical combination. This package covers
`Ds1(2536) -> Ds gamma` only. Other channels and the requested DA, asymmetric
Borel, continuum, and omitted-term sensitivity analyses remain pending.

## Figure and reproduction

- [Three-panel PNG](ds2536_signed_contributions.png)
- [Vector SVG](ds2536_signed_contributions.svg)
- [All sampled contributions](ds2536_signed_contributions_data.csv)
- [Midpoint decomposition](ds2536_signed_contributions_midpoint.csv)
- [Inputs, source hashes, and numerical checks](ds2536_signed_contributions_provenance.json)
- [Builder](build_signed_contribution_figure.py)

Run `python3 -B review/referee_a_major_ii/build_signed_contribution_figure.py`
from the project root. Outputs remain in this review directory.

## Prescription

The reference manuscript is `/Users/sbilmis/Downloads/r1-2.pdf`, Table II.
The three panels show signed contributions to Pi_A, Pi_B, and g2. Each has
perturbative, twist-two, twist-three, combined twist-four, and total curves.
Combined twist four includes the retained two-particle, gluonic
three-particle, and electromagnetic three-particle terms; these are also
saved separately in the CSV files. The midpoint CSV contains both the three
separate twist-four entries and their subtotal, which must not be counted
twice when summing.

The scan uses 61 M^2 points from 3 to 4.5 GeV^2 and 21 s0 points from 9 to
10 GeV^2. Lines use s0=9.5 GeV^2; shading is the pointwise sampled s0
envelope, not a confidence interval. Theta=26.6 degrees, u0=1/2, and
M1^2=M2^2=2M^2. All other inputs are the fixed reference inputs of the
existing manuscript Borel builder. In particular f2=0.16667036765413243
GeV is fixed, as in that builder. The transition threshold being scanned
is separate from the two-point thresholds used to obtain f2.

For every retained contribution k, the third panel uses

    g2[k] = N2 * (cos(theta) * Pi_A[k] - sin(theta) * Pi_B[k]).

It therefore decomposes the amplitude, not a sum of separate widths.
The axial twist-three term is zero at this symmetric point in the current
implementation; the tensor derivative term survives.

## Midpoint result and interpretation

At M^2=3.75 GeV^2 and s0=9.5 GeV^2:

| Contribution | Pi_A (GeV^3) | Pi_B (GeV^3) | g2 (GeV^-1) |
| --- | ---: | ---: | ---: |
| Perturbative | -0.10398970 | -0.12784573 | -0.47022825 |
| Twist 2 | +0.05897475 | +0.05495080 | +0.37008824 |
| Twist 3 | 0 | -0.00458752 | +0.02702660 |
| Twist 4, combined | +0.00161247 | -0.00507157 | +0.04884861 |
| Total | -0.04340247 | -0.08255401 | -0.02426480 |

The total g2 varies from -0.09061996 to +0.05990042 GeV^-1 over the sampled
window. At the midpoint both projected twist-three and combined twist-four
contributions exceed the magnitude of the residual. This explains why
small terms in the basis amplitudes can matter for the suppressed channel.
It does not by itself establish failure of the OPE or quantify missing
terms. The point value is not the Monte Carlo median in the manuscript.

## Checks and limitations

All 93 Ds1(2536) points in the saved manuscript Borel scan and all 42 saved
central numerical entries are reproduced exactly. The largest discrepancy
between the sum of plotted components and the total is 3.68e-16. The PNG
was visually inspected and its title/footer spacing adjusted.

This is a decomposition of the preserved production baseline. It adopts
neither the pending two-point mixed-condensate patch nor changes to the
continuum prescription or the shared-input Monte Carlo treatment. It is
not a new uncertainty estimate or a resolution of the other referee items.
Source/reference hashes and a before/after `final_manuscript/` manifest
confirm that the protected manuscript and baseline references were not
modified.

Major comment (i) remains separate: the discussion of arXiv:1007.0814
explains the diagonal-current construction but has not supplied the
requested present-OPE off-diagonal check or a quantitative pole-contamination
bound. Do not mark that item complete when using this figure.
