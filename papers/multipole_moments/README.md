<!-- LATTICE UPDATE -->
**3 October 2026:** [Lattice-normalized results with both twist-2 signs](black_box_comparison/recalculation/lattice_results.md). Includes nominal and refitted 90-degree windows, window intervals and separate contribution medians. Frozen outputs and their SHA-256 manifest remain unchanged.

# Electromagnetic multipole moments of the physical \(D_{s1}\) and \(B_{s1}\) mesons

**Current numerical analysis (25 September 2026):** [independent recalculation](black_box_comparison/recalculation/RECALCULATION_REPORT.md) and [updated PDF](output/pdf/multipole_moments_recalculated_2026_09_25.pdf). This replaces the earlier numerical recommendations and uncertainty estimates below. It includes the repaired OPE and residues, clean multipole projections, new window scans, angle curves and independent consistency checks. The step-numbered records remain historical diagnostics.

The current report also contains the [susceptibility-sign and charge derivation](black_box_comparison/recalculation/susceptibility_charge_audit/EXPLANATION.md) and [all documented structure extractions at nominal mixing and 90 degrees](black_box_comparison/recalculation/structure_tables/STRUCTURE_TABLES.md), with frozen and refitted windows distinguished. Use its corrected conventions and extraction formulas where they supersede the historical text below. Alternative structure results still differ, including in the pure tensor-current channel, and are identified as diagnostics rather than equivalent predictions.

## Aim

This study will calculate the static electromagnetic properties of the four
physical strange heavy-light axial-vector states

\[
D_{s1}(2460)^+,\qquad D_{s1}(2536)^+,\qquad
B_{s1}(5750),\qquad B_{s1}(5830).
\]

For each state the targets are the magnetic dipole form factor \(F_M(0)\)
and the electric quadrupole form factor \(F_D(0)\).  The charge form factor
\(F_C(0)\) is retained as a normalization and Ward-identity diagnostic.

This file fixes the conventions and records the current OPE status. Detailed
derivation notes, symbolic scripts, and numerical outputs are added as the
calculation progresses. Any table labelled `hard-only` or `hard + 2p partial`
is a
diagnostic, not a final multipole prediction.

The shareable computational-notes PDF numbers every displayed equation by
note, for example (4.7), (11.14), and (12.3). Note 12 records the recommended
paper appendix for the explicit \(AA\), \(BB\), and physical \(AB+BA\)
invariant amplitudes.

## Literature decision: which form factor is the charge?

The standard spin-one electromagnetic vertex is written in terms of three
covariant form factors \(G_1,G_2,G_3\):

\[
\begin{aligned}
\langle H(p',\eta')|j_\alpha^{\rm em}|H(p,\eta)\rangle
=-\eta^\rho\eta'^{*\beta}\bigg[
&G_1(Q^2)g_{\rho\beta}(p+p')_\alpha\\
&+G_2(Q^2)(q_\rho g_{\alpha\beta}-q_\beta g_{\alpha\rho})\\
&+\frac{G_3(Q^2)}{2m_H^2}q_\rho q_\beta(p+p')_\alpha
\bigg],
\end{aligned}
\]

where \(q=p'-p\), \(Q^2=-q^2\), and
\(\eta=Q^2/(4m_H^2)\).  The Sachs-type multipole form factors are

\[
F_C=G_1+\frac{2}{3}\eta F_D,\qquad
F_M=G_2,\qquad
F_D=G_1-G_2+(1+\eta)G_3.
\]

At \(Q^2=0\), the static quantities in the convention used here are

\[
F_C(0)=Q_H/e,\qquad
\mu_H=\frac{e}{2m_H}F_M(0),\qquad
\mathcal D_H=\frac{e}{m_H^2}F_D(0).
\]

Therefore the handwritten note should read

\[
F_C^{D_{s1}(2460)^+}(0)=F_C^{D_{s1}(2536)^+}(0)=+1,
\qquad
F_C^{B_{s1}(5750)}(0)=F_C^{B_{s1}(5830)}(0)=0,
\]

not \(F_D(0)=1\) or \(0\).  The quadrupole form factor is not fixed by the
electric charge.  In the hadronic tensor it naturally occurs through
\(F_C(0)+F_D(0)\), so \(F_D(0)\) is obtained by subtracting the independently
projected \(F_C(0)\).

## Why physical-state mixing is required

The 2020 LCSR calculation used the single current
\(\bar s\gamma_\mu\gamma_5 Q\) and quoted one \(D_{s1}\) and one \(B_{s1}\)
result.  A single current does not separate the two physical \(1^+\) states.
The present study explicitly asks for both partners in each flavor sector;
therefore the two-current description is necessary.

Use the normalized basis currents

\[
J_{A\mu}=\bar s\gamma_\mu\gamma_5 Q,
\qquad
J_{B\mu}=\frac{i}{m_Q+m_s}\bar s\sigma_{\mu\alpha}P^\alpha\gamma_5Q,
\]

where \(Q=c,b\) and \(P\) is the momentum entering the current.  The physical
currents are

\[
J_{1\mu}=\sin\theta\,J_{A\mu}+\cos\theta\,J_{B\mu},\qquad
J_{2\mu}=\cos\theta\,J_{A\mu}-\sin\theta\,J_{B\mu}.
\]

State 1 denotes the lower member and state 2 the higher member:

| Sector | State 1 | State 2 | Nominal angle |
|---|---|---|---:|
| charm-strange | \(D_{s1}(2460)^+\) | \(D_{s1}(2536)^+\) | \(26.6\pm0.6{}^\circ\) |
| bottom-strange | \(B_{s1}(5750)\) | \(B_{s1}(5830)\) | \(38.5\pm0.1{}^\circ\) |

These are the values quoted in the 2025 mixing-angle paper,
\(\theta_{D_s}=26.6\pm0.6^\circ\) and
\(\theta_{B_s}=38.5\pm0.1^\circ\). Its displayed tensor current is

\[
J_{B\mu}^{\rm paper}=i\bar s\sigma_{\mu\alpha}P^\alpha\gamma_5Q,
\]

without the \(1/(m_Q+m_s)\) factor. The authors have clarified that this
factor was omitted only in the written current definition and was included
in the actual calculation. The published angles therefore belong to the
normalized basis used here and are the nominal inputs.

## Important difference from a transition correlator

This is an elastic current-current correlator.  Both external meson currents
are mixed, so the physical OPE is quadratic in the rotation coefficients.
For \(s_\theta=\sin\theta\), \(c_\theta=\cos\theta\),

\[
\begin{aligned}
\Pi_{11}&=s_\theta^2\Pi_{AA}+c_\theta^2\Pi_{BB}
          +s_\theta c_\theta(\Pi_{AB}+\Pi_{BA}),\\
\Pi_{22}&=c_\theta^2\Pi_{AA}+s_\theta^2\Pi_{BB}
          -s_\theta c_\theta(\Pi_{AB}+\Pi_{BA}).
\end{aligned}
\]

Keeping only a linear combination of an axial and a tensor contribution would
be appropriate when only one end of a transition correlator is mixed, but it
would be incorrect here.  The \(AA,AB,BA,BB\) terms must be derived and checked
separately before rotation.

## Hadronic projections

After summing the initial and final spin-one polarizations and imposing
\(q^2=0\) and \(q\cdot\epsilon_\gamma=0\), the double-pole numerator contains

\[
\begin{aligned}
&2F_C(0)(p\cdot\epsilon_\gamma)
  \left(g_{\mu\nu}-\frac{p_\mu p_\nu}{m_H^2}
  -\frac{p_\mu q_\nu}{m_H^2}\right),\\
&F_M(0)\left(q_\mu\epsilon_{\gamma\nu}
  -q_\nu\epsilon_{\gamma\mu}
  -\frac{p\cdot\epsilon_\gamma}{m_H^2}
   (p_\mu q_\nu-p_\nu q_\mu)\right),\\
&-[F_C(0)+F_D(0)]
  \frac{p\cdot\epsilon_\gamma}{m_H^2}q_\mu q_\nu.
\end{aligned}
\]

This printed tensor is being treated as a result to test, not an identity to
assume.  The first FeynCalc audit found that a literal polarization sum of the
displayed covariant vertex gives

\[
\begin{aligned}
H_{\mu\nu}^{\rm direct}={}&-2F_C(p\!\cdot\!\epsilon)
 \left[g_{\mu\nu}-\frac{p_\mu p_\nu}{m_H^2}
 -\frac{q_\mu p_\nu}{m_H^2}
 -\frac{3q_\mu q_\nu}{2m_H^2}\right]\\
&+F_M\left[\epsilon_\mu q_\nu-q_\mu\epsilon_\nu
 -\frac{p\!\cdot\!\epsilon}{m_H^2}
  (p_\mu q_\nu-q_\mu p_\nu)\right]
 -F_D\frac{p\!\cdot\!\epsilon}{m_H^2}q_\mu q_\nu,
\end{aligned}
\]

for the momentum and overall-vertex conventions stated above.  This direct
tensor obeys all three required identities

\[
H_{\mu\nu}q^\alpha\big|_{\epsilon_\alpha\to q_\alpha}=0,
\qquad p^\mu H_{\mu\nu}=0,
\qquad (p+q)^\nu H_{\mu\nu}=0.
\]

The printed 2020 expression is electromagnetic-gauge invariant but, under the
same routing, gives nonzero meson-transversality residuals proportional to
\(F_C+F_M\).  Double-projecting it does not remove the mismatch.  This is now
an explicit audit item: the OPE projectors will not be fixed until momentum
routing, index ordering, and the published reduction have been independently
reconciled.  Alternative projectors will also be checked because the 2009
light-meson analysis emphasizes contact-term-free structures, while the 2020
heavy-meson paper favors maximal-momentum structures for OPE convergence.

The subsequent coefficient audit resolves how the direct tensor must be
projected. In the maximal-momentum basis, define \(h_{pp},h_{qp},h_{qq}\)
as the coefficients of
\((p\cdot\epsilon)p_\mu p_\nu/m_H^2\),
\((p\cdot\epsilon)q_\mu p_\nu/m_H^2\), and
\((p\cdot\epsilon)q_\mu q_\nu/m_H^2\). The exact result is

\[
h_{pp}=2F_C,\qquad
h_{qp}=2F_C+F_M,\qquad
h_{qq}=3F_C-F_D,
\]

so that

\[
F_C=\frac{h_{pp}}2,\qquad
F_M=h_{qp}-h_{pp},\qquad
F_D=\frac32h_{pp}-h_{qq}.
\]

The determinant of this map is \(-2\), and substitution recovers all three
form factors with zero symbolic residual. This demonstrates why the OPE
invariant amplitudes cannot be assigned to \(F_C,F_M,F_D\) one-to-one before
the full matching matrix is established.

For the magnetic moment, the final projector choice is instead the
contact-term-free structure identified in the 2009 spin-one analysis.  Its
correlator orders the two current indices oppositely to ours, so its
\(q_\mu\epsilon_\nu\) structure maps to
\(\epsilon_\mu q_\nu\) here.  The direct tensor gives

\[
H_{\mu\nu}\supset F_M
(\epsilon_\mu q_\nu-q_\mu\epsilon_\nu),
\]

and the charge-calibrated invariant convention therefore gives

\[
F_M(0)=-\frac{e^{m_H^2/M_B^2}}
{f_H^2m_H^2}\Pi_{\epsilon_\mu q_\nu}^{\rm OPE}.
\]

The contact-free route is retained as the selected working prescription.
Historical comparisons are isolated in `verification_prd101_054009` and are
not used to define the present result. The main theory/numerical draft now
also compares the independent \(\epsilon_\mu q_\nu\) and
\(q_\mu\epsilon_\nu\) contact-free projections with the rematched maximal
basis on identical correlated input rows. This spread is reported as a
structure-choice diagnostic, not an uncertainty estimate.

## Reproducibility checks already completed

On 2026-08-13 and 2026-08-14:

- `/Applications/Wolfram.app/Contents/MacOS/WolframKernel` evaluated `2+2` as
  `4`.
- FeynCalc 10.2.1 loaded from
  `/Users/sbilmis/Library/Wolfram/Applications/FeynCalc/`.
- The smoke test
  \(\operatorname{Tr}(\gamma^\mu\gamma^\nu)=4g^{\mu\nu}\) passed exactly.
- The seven-element static tensor decomposition and both three-structure
  matching matrices passed with zero residuals.
- The perturbative hard-photon numerators have been generated separately for
  heavy- and strange-line emission in all \(AA,AB,BA,BB\) channels. All eight
  fermion-line Ward residuals and the \(q=0\) transposed \(AB-BA\) residual
  vanish.
- The alternative full-tensor rematching of the published
  \(\Pi_1,\Pi_2,\Pi_3\) amplitudes recovers the three static form factors with
  zero symbolic residual. Its disagreement with the printed one-to-one
  assignment remains an explicit verification question.
- The exact-mass perturbative triangles have been reduced to Gram-safe
  one-dimensional double-Borel kernels for all eight channel/topology pairs.
  Every rank-four numerator residual vanishes. The Minkowski numerator
  moments use the Wick-corrected rule
  \(\langle\ell_M^2\rangle=-2/T\); this affects only the \(BB\) channel.
- The pure-\(AA\) perturbative charge coefficient exactly reproduces the 2020
  Appendix-B normalization, fixing the common factor \(N_c/(16\pi^2)\) and
  the relative heavy/strange emission sign.
- All 20 two-particle Fierz traces pass an independent cyclic-trace check.
- The full two-particle photon-DA sector through twist 4 has been Fourier
  transformed and double-Borel transformed in all \(AA,AB,BA,BB\) channels.
  Its leading-twist entries reproduce the independent Step-9 derivation with
  zero residual, and its pure-\(AA\) twist-3/twist-4 entries reproduce the
  corresponding 2020 Appendix-B terms with zero residual.
- The eight three-particle photon-DA operator tensors pass gauge-invariance,
  gluon-index antisymmetry, tensor-bilinear antisymmetry, and cyclic-trace
  checks.  The source-paper names of \(\mathcal T_1\) and
  \(\mathcal T_2\) are used consistently.
- Both terms of the background-gluon heavy propagator have been projected in
  all \(AA,AB,BA,BB\) channels.  The pure-\(AA\)
  \(\mathcal A,\mathcal V\) and \(\mathcal T_{1\ldots4}\) results reproduce
  Appendix B with exact zero residuals.  The pure-\(AA\)
  \(\mathcal S,\widetilde{\mathcal S}\) maximal-momentum coefficients vanish
  exactly, as required.
- The complete tensor current matrix contains all 16 channel/operator
  entries, no free Lorentz index, no weighted-theta remainder, and no
  \(h_{pp}\) contribution.
- All native three-particle simplex moments at \(u_0=1/2\) have been
  integrated analytically, including the delta-functional Jacobian.
- The hard numerical scan contains 16 state/Borel rows. Its 192-versus-384
  point quadrature comparison has a maximum relative difference
  \(1.62\times10^{-9}\).
- The perturbative charge audit proves exact pointwise equal-charge
  cancellation in \(AA,AB,BA\). The \(BB\) sum is exactly a total
  \(r\)-derivative. Removing its symmetric finite-continuum contact boundary
  restores the neutral hard charge to \(1.4\times10^{-17}\) over all 16 rows.
- The correctly index-matched contact-free
  \(\epsilon_\mu q_\nu\) coefficient has been derived independently in the
  hard, two-particle, and three-particle photon-DA sectors. Ordinary local
  vacuum-condensate insertions are excluded from the external-photon LCSR. The
  selected 16-row partial magnetic scan is stored in
  `outputs/step42_epsq_contact_free_full_magnetic_scan.csv`.
- The multipole scans use the reported accepted-scan median
  \(s_0=9.925\,\mathrm{GeV}^2\) for \(D_{s1}(2536)\).  The nearby
  \(9.975\,\mathrm{GeV}^2\) number is only a single-point Mathematica
  regression value and is not substituted for the published result.
- The selected magnetic and quadrupole-diagnostic chains have also been run
  on the two-point-accepted correlated rows.  These contain 14 charm points
  over \(2.05\leq M_B^2\leq2.70\,\mathrm{GeV}^2\) and 17 bottom points over
  \(6.50\leq M_B^2\leq7.30\,\mathrm{GeV}^2\).  At every point the fitted
  \((s_0,f_H)\) pair from the same two-point row is used.  Step 53 verifies
  identical provenance for all 62 magnetic and quadrupole rows.  A later
  state-indexed audit found that all 17 legacy bottom pairs obeyed
  \(s_{0,2}\leq s_{0,1}\), opposite to the expected physical-pole ordering.
  The bottom rows in these legacy tables are therefore provisional and must
  be rerun with the ordered thresholds before they are quoted as final
  physical-state results.  The independent ordered refit is documented in
  `verification_prd101_054009/` and does not overwrite the legacy artifacts.
- The resulting selected no-\(G^2\) median magnetic form factors are
  \(-1.995\), \(-2.280\), \(+0.125\), and \(-1.761\) for
  \(D_{s1}(2460)^+\), \(D_{s1}(2536)^+\), \(B_{s1}(5750)\), and
  \(B_{s1}(5830)\), respectively.  These medians and their accepted-window
  ranges are deterministic Borel diagnostics, not full statistical errors.
  In particular, the \(B_{s1}(5750)\) result is a strong hard--photon-DA
  cancellation and is not yet a precision prediction.
- The state-ordered rerun is stored separately in Steps 80--88.  Its median
  magnetic form factors are \(-1.995\), \(-2.280\), \(+0.122\), and
  \(-1.724\), and its selected no-\(G^2\) quadrupole form factors are
  \(+4.587\), \(+2.252\), \(+2.360\), and \(-8.804\), in the same state
  order.  These deterministic Borel-window diagnostics supersede the legacy
  bottom medians for state-indexed reporting; they are not an uncertainty
  analysis and remain incomplete without the gauge-complete vacuum-\(G^2\)
  sector.
- Mixing-angle endpoint propagation and Monte Carlo input variations are
  deferred until the central ordered-input calculation and the
  gauge/contact completion are settled.
- The quadrupole tensor
  \((p\!\cdot\!\epsilon)q_\mu q_\nu/m_H^2\) is symmetric under current-index
  exchange and is unchanged by the \(pp\)-only charge-contact subtraction.
  Step 59 proves exact hadronic recovery and reconstructs all 62 values with
  maximum residual \(1.8\times10^{-15}\).  The charge-normalized Step-52
  \(F_D\) column is therefore the selected no-\(G^2\) quadrupole result.
  No quadrupole uncertainty interval is assigned in this version.
- The candidate vacuum-gluon insertion inventory contains all 48
  single-line/open-pair objects, reproduces all 24 maximal and both sets of
  eight contact-free unit-power hard kernels, and has zero rank-six
  residuals.  It is not used numerically: the pure-\(AA\) spectral audit of
  both selected \(\epsilon_\mu q_\nu\) and \(h_{qq}\) coefficients retain
  nonzero lower-endpoint poles.  This proves that mixed photon--gluon
  gauge/contact terms are still missing.  The failure persists after the
  Wick-sign and current-index corrections.
- A deterministic no-mixing counterfactual has been evaluated with
  \(\theta=0\), so state 1 is pure \(J_B\) and state 2 is pure \(J_A\).
  Thresholds and residues are refitted at every reference Borel point. A
  common charm window survives, but the pure-\(J_A\) assignment to
  \(B_{s1}(5830)\) has no accepted two-point window and leaves no common
  bottom window. The four numerical entries are therefore retained as a
  diagnostic table in Note 10, with the bottom results explicitly excluded
  from the class of accepted no-mixing predictions. No uncertainty analysis
  is attached to this counterfactual.

The detailed derivations are in
[`notes/01_scope_and_hadronic_audit.md`](notes/01_scope_and_hadronic_audit.md),
[`notes/02_static_projector_matching.md`](notes/02_static_projector_matching.md),
[`notes/03_elastic_ope_routing.md`](notes/03_elastic_ope_routing.md),
[`notes/04_ope_building_blocks_and_borel.md`](notes/04_ope_building_blocks_and_borel.md),
[`notes/05_two_particle_photon_das.md`](notes/05_two_particle_photon_das.md),
[`notes/06_three_particle_photon_das.md`](notes/06_three_particle_photon_das.md),
[`notes/08_charge_normalization_and_continuum.md`](notes/08_charge_normalization_and_continuum.md),
[`notes/09_contact_free_magnetic_projector.md`](notes/09_contact_free_magnetic_projector.md),
[`notes/10_no_mixing_comparison.md`](notes/10_no_mixing_comparison.md),
and
[`notes/11_explicit_spin1_sum_rules.md`](notes/11_explicit_spin1_sum_rules.md),
with the current-basis appendix in
[`notes/12_current_basis_invariant_amplitudes.md`](notes/12_current_basis_invariant_amplitudes.md),
with explicit two-point spectral densities in
[`notes/14_explicit_twopoint_spectral_densities.md`](notes/14_explicit_twopoint_spectral_densities.md).
The paper-style hadronic representation and explicit spin-one form-factor
sum rules are collected in Note 11.

The complete shareable notes, including the explicit spin-one sum rules, the
zero-mixing comparison, the independent PRD-2020 verification appendix, and
the two original calculation sketches, are compiled as
[`output/pdf/multipole_moments_computational_notes.pdf`](output/pdf/multipole_moments_computational_notes.pdf).
The shorter theory-and-numerical checking draft is
[`output/pdf/multipole_moments_theory_numerical_check.pdf`](output/pdf/multipole_moments_theory_numerical_check.pdf).
Regenerate the document with
`scripts/build_computational_notes_pdf.sh`.

The independent convention-matched comparison with the 2020 single-current
calculation, including the ordered-threshold and \(f_1/f_2\) audits, is kept
separately in [`verification_prd101_054009/`](verification_prd101_054009/).
This comparison deliberately does not assume that either projector
prescription is authoritative.

## Primary references

1. T. M. Aliev, S. Bilmis, and M. Savci, *Multipole Moments of Heavy Vector
   and Axial-Vector Mesons in QCD*, Phys. Rev. D **101**, 054009 (2020),
   [arXiv:1912.05988](https://arxiv.org/abs/1912.05988).
2. T. M. Aliev, A. Ozpineci, and M. Savci, *Magnetic and quadrupole moments
   of light spin-1 mesons in light cone QCD sum rules*,
   [arXiv:0902.4627](https://arxiv.org/abs/0902.4627).
3. T. M. Aliev, S. Bilmis, and M. Savci, *Mixing angle between
   \(^{3}P_1\) and \(^{1}P_1\) states in heavy axial vector mesons within the
   QCD sum rules framework*, Phys. Rev. D **111**, 014005 (2025),
   [arXiv:2408.08014](https://arxiv.org/abs/2408.08014).

## Complete structure calculation (Steps 98--100, 2026-09-16)

All seven nonzero tensors of the static pole numerator are now covered. The
previous statement that the raw `(p.epsilon) p_mu q_nu` hard term vanishes
was incorrect: all eight current-channel/emission kernels are nonzero.
Steps 98--99 independently project this tensor and the metric tensor, with
full two-particle and gluonic three-particle DA contributions in the existing
OPE inventory. Step 100 compares their results with the existing projections
on identical 62 ordered correlated input rows, including charge and
quadrupole diagnostics. These remain partial OPE/continuum diagnostics;
vacuum G2 and gauge/contact completion are not supplied by this extension.

- `outputs/step100_all_structure_tables.md`: readable results and formulas.
- `outputs/step100_all_structure_rows.csv`: row-level invariants and moments.
- `outputs/step100_all_structure_summary.csv`: medians and full row ranges.
- `outputs/step98_{pq,g}_*_kernels.wl`: independently extracted kernels.
- `outputs/step99_hard_complete_basis_audit.txt`: exact symbolic regressions.

Reproduce with `python3 scripts/step98_all_structure_projections.py`, then
WolframKernel `-script scripts/step99_hard_complete_basis_audit.wl`, and
`python3 scripts/step100_all_structure_summary.py`. Step 98 uses isolated
copies under `outputs/step98_work`; its manifest records the source hashes
and the internal legacy-slot mapping. It never overwrites the original
projection artifacts. The PDFs contain the updated comparisons.

## All form factors at nominal angles and 90 degrees (Steps 101--102)

The complete comparison is `outputs/step102_all_form_factors_angles_tables.md`,
with row-level data in `step102_all_form_factors_angles_rows.csv` and medians
and accepted-row extrema in `step102_all_form_factors_angles_summary.csv`
(both in `outputs`). It contains FC, FM, FD, mu/muN, and D/(e fm^2).
Each extraction is displayed at both the nominal angle and theta=90 degrees.
The latter uses pure AA for state 1 and pure BB for state 2, with consistently
rotated two-point residues at fixed nominal M^2, s0, and mass. It is not an
independently refitted pure-current window.

Step 90 now accepts `MULTIPOLE_EXTRA_STRUCTURE=pq` or `g`, producing the
additional current-basis files `step101_pq_angle_basis.csv` and
`step101_g_angle_basis.csv` without changing the default Step-90 result.
Then run `python3 scripts/step102_all_form_factors_angles.py`. The new
calculation checks the nominal results row by row and reproduces existing
90-degree form factors and residues. `step102_validation.json` records the
residuals and source hashes. These remain partial OPE/continuum diagnostics.
