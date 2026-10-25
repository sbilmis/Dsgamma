# Mixed-condensate coefficient: independent derivation

19 October 2026. The author asked to resolve the coefficient issue rather
than leave it as a condition in the referee-response draft.

**Finding:** the current two-point implementation contains the light-bilocal
contribution to the mixed condensate but omits the contribution with a
background gluon on the heavy-quark line. The omitted contribution is zero
in AA and nonzero in AB and BB. Both contributions have now been derived
explicitly at tree level and leading order in the light mass. This resolves
the leading dimension-five coefficient used in the current approximation.
It does not reconstruct the original 2025 numerical program or its larger OPE.

An isolated implementation and an **unapplied** four-file source patch are
provided here. Production sources, accepted samples, figures, and
`final_manuscript/` have not been replaced. The author-confirmed tensor-current
normalization remains unchanged.

## Conventions and result

Use the present basis

\[
J_A^\mu=\bar s\gamma^\mu\gamma_5 Q,\qquad
J_B^\mu=\frac{i}{d}\bar s\sigma^{\mu\nu}p_\nu\gamma_5 Q,
\qquad d=m_Q+m_s.
\]

Write \(m=m_Q\), \(K=M^2\),
\(C_5=\langle\bar s g_s\sigma\!\cdot\!G s\rangle\), and
\(E=\exp(-m^2/K)\). The invariant is the coefficient of
\(-g^{\mu\nu}+p^\mu p^\nu/p^2\), with the same sign as the positive
spectral matrix in the current code. The Wilson coefficients below are
evaluated at \(m_s=0\), while retaining the specified current normalization
\(d=m_Q+m_s\). Terms proportional to \(m_s C_5\), radiative corrections,
and gluon-only condensates are outside this calculation.

The missing Borel-transformed matrix is

\[
\boxed{\displaystyle
\delta\widehat{\boldsymbol\Pi}_5
=C_5 E\left(1-\frac{m^2}{K}\right)
\begin{pmatrix}
0&\dfrac{1}{12d}\\[3pt]
\dfrac{1}{12d}&\dfrac{m}{6d^2}
\end{pmatrix}.}
\]

After adding it, the complete leading mixed-condensate coefficients are

\[
\widehat\Pi_{AA}^{(5)}=-C_5E\frac{m^3}{4K^2},\qquad
\widehat\Pi_{AB}^{(5)}=\widehat\Pi_{BA}^{(5)}
=C_5E\frac{4K^2+2Km^2-3m^4}{12dK^2},
\]

\[
\widehat\Pi_{BB}^{(5)}
=C_5E\frac{m(2K^2+4Km^2-3m^4)}{12d^2K^2}.
\]

Consequently, the existing pure axial-current check is unchanged. Agreement
in that limit cannot detect this omitted contribution to AB and BB.

## Where the additional term comes from

There are two contributions to the same local operator at this order:

1. The light-bilocal term \(S_s(x)|_{C_5}=-x^2C_5/192\), already present
   in the code. Fourier transformation gives \(x^2\to-\Box_k\); its
   coefficient after the color trace is
   \(C_5\operatorname{Tr}[\Gamma_L\Box_k H(k)\Gamma_R]/64\), with
   \(H(k)=(\slashed k+m)/(k^2-m^2)\).
2. The heavy propagator with one background gluon, contracted with the
   open local quark-gluon operator. In the consistently stripped propagator
   convention used in the trace scripts,

\[
S_Q^{(G)}(p)=-\frac{g_s t^a G^a_{\alpha\beta}}4
\frac{\{\sigma^{\alpha\beta},\slashed p+m\}}{(p^2-m^2)^2},
\qquad
\langle g_s s_\alpha^i\bar s_\beta^jG_{\mu\nu}^a\rangle
=-\frac{C_5}{192}(t^a)_{ij}(\sigma_{\mu\nu})_{\alpha\beta}.
\]

The projector normalization is checked using
\(\sum_a\operatorname{Tr}_c(t^at^a)=4\) and
\(\operatorname{Tr}_D(\sigma^{\mu\nu}\sigma_{\mu\nu})=48\).
Their product is 192. The two factors above give
\(C_5\operatorname{Tr}[\Gamma_L\{\sigma^{\alpha\beta},\slashed p+m\}
\Gamma_R\sigma_{\alpha\beta}]/[192(p^2-m^2)^2]\).
The right vertex is the Dirac adjoint of the conjugate current; the explicit
rest-frame vertices used in the symbolic calculation obey this relation.

Before the Borel transform, the omitted transverse entries, with \(C_5\)
factored out and \(s=p^2\), are

\[
\delta\Pi_{AA}=0,\qquad
\delta\Pi_{AB}=\delta\Pi_{BA}=-\frac{s}{12d(m^2-s)^2},\qquad
\delta\Pi_{BB}=-\frac{m\,s}{6d^2(m^2-s)^2}.
\]

The Borel convention is

\[
\mathcal B\big[(m^2-p^2)^{-n}\big]
=\frac{E}{(n-1)!K^{n-1}},\qquad
\frac{s}{(m^2-s)^2}=\frac{m^2}{(m^2-s)^2}-\frac1{m^2-s}.
\]

This gives the boxed matrix directly; no fit to any published angle is used.

The two mixed-condensate topologies and their background-field construction
are described in Section 6 of [Jamin and Münz, hep-ph/9208201](https://arxiv.org/abs/hep-ph/9208201).
The axial, vector, and pseudoscalar control expressions are available in
Eqs. (14), (18), and (23) of
[Wang, Eur. Phys. J. C 75, 427 (2015)](https://link.springer.com/article/10.1140/epjc/s10052-015-3653-9).
Local copies and the inspected pages are in `references/`. The tensor-current
matrix above is the independent calculation in this directory, not a formula
quoted from those papers.

## What this says about the earlier factor-of-two observation

For \(m_s=0\), so \(d=m\), the three expressions for the mixed-condensate
part of \(\widehat\Pi_{AA}-\widehat\Pi_{BB}\) are:

| Expression | Coefficient after factoring out \(C_5 E\) |
| --- | --- |
| Current implementation, light-bilocal term only | \(-m/(2K)\) |
| Complete leading coefficient derived here | \(-m/(3K)-1/(6m)\) |
| Printed earlier appendix, with the leading-term label assignment used in the audit | \(-m/K\) |

Thus choosing one of the original two coefficients is insufficient. The
complete result has an additional Borel dependence. The earlier printed
off-diagonal expression also matches only the light-bilocal contribution in
this limit. The current missing topology is now identified; whether the
original numerical code used the printed expressions remains unknown.
This distinction must not be collapsed into a claim that the published
angles themselves have been disproved.

The earlier response draft's +0.185-degree and +0.057-degree condensate shifts
were obtained from the printed-appendix reconstruction. They remain results
of that reconstruction, not values certified by this new derivation. Do not
silently replace them with the different-OPE check below or call either a
verified rerun of the original angle extraction.

## Verification

- `derive_dimension5.wl` uses explicit 4 by 4 gamma matrices. All seven
  symbolic checks pass: AB=BA, reproduction of the three existing light-only
  entries, and the standard axial, vector, and pseudoscalar controls.
- In units of \(C_5E\), the controls are axial \(-m^3/(4K^2)\), vector
  \(+m^3/(4K^2)\), and pseudoscalar \(-m/(2K)+m^3/(4K^2)\).
  The heavy-line term is essential for the first term of the pseudoscalar
  result; this tests its relative normalization and sign.
- `check_covariant_traces.py` independently contracts explicit color and
  gamma matrices with the transverse projector, in rest and boosted frames.
  All 24 charm/bottom and normalization checks pass, with maximum absolute
  difference \(1.11\times10^{-15}\).
- The proposed Python patch is evaluated in memory at four points and agrees
  with the isolated analytic addition to \(2.17\times10^{-19}\). The other
  local contributions are unchanged. The patch incorporates the correction
  into the existing `d5_mixed` component, so OPE diagnostics will include both
  contributions rather than testing the light-only term.
- All three proposed Wolfram sources parse successfully without running their
  output-writing routines. Evaluating only each `localPieces` definition gives
  exact symbolic agreement with the independently derived complete d5 matrix;
  see `proposed_wolfram_patch_check.json`.
- Hash checks confirm that 11 manuscript files, the two production Python
  sources, the two accepted-residue CSVs, and the preferred transition scan
  were unchanged during the impact calculation.

## Fixed-point numerical impact

Use the same representative accepted input samples as the angle-plot builder,
and the external angles 26.6 degrees and 38.5 degrees. Keep the input masses,
Borel parameters, and thresholds fixed. This gives:

| Sector/state | Existing f [MeV] | With complete leading d5 [MeV] | Change in f | Width change at fixed transition numerator |
| --- | ---: | ---: | ---: | ---: |
| Ds, 1 | 404.9967 | 404.4145 | -0.1438% | +0.2881% |
| Ds, 2 | 169.0940 | 169.3804 | +0.1693% | -0.3378% |
| Bs, 1 | 547.8345 | 547.9954 | +0.0294% | -0.0587% |
| Bs, 2 | 91.0942 | 90.9877 | -0.1169% | +0.2341% |

These are representative values, not marginalized predictions. The width
comparison uses only the two-point denominator change,
\(\Gamma_{\rm completed}/\Gamma_{\rm existing}=(f_{\rm existing}/f_{\rm completed})^2\).
No transition OPE change is inferred from this two-point calculation.
Thresholds have not been refitted, accepted samples regenerated, or total
uncertainties recomputed. In particular, the effective-mass acceptance
criterion must be reevaluated when adopting the correction.

A separate paired condensate check in the **completed current LO+d3+d5 OPE**
changes \(\langle\bar qq\rangle=-0.0165\) to \(-0.013824\,\mathrm{GeV}^3\),
holding \(m_s=0.126\,\mathrm{GeV}\), \(\kappa_s=0.8\), and
\(m_0^2=0.8\,\mathrm{GeV}^2\) fixed:

| Sector | Midpoint M², s0 [GeV²] | Condensate-induced angle shift | Range on 21 by 21 grid |
| --- | --- | ---: | --- |
| Ds | 4.5, 10 | +0.174981 degrees | +0.115984 to +0.331162 degrees |
| Bs | 11, 43 | +0.045126 degrees | +0.036148 to +0.060047 degrees |

The grids are M²=3--6, s0=9--11 GeV² for Ds and M²=9--13, s0=42--44 GeV²
for Bs. Relative to the uncompleted current code, the midpoint sensitivities
change by -0.002122 degrees and +0.000420 degrees. This quantifies the influence
of the missing topology on that check. The full original OPE and its
uncertainty procedure have not thereby been validated.

## Implementation and next steps

- `candidate_dimension5.py` evaluates the completed matrix in isolation.
- `proposed_implementation.patch` adds the same expression to the Python
  two-point source, the Ds/Bs Wolfram checks, and the notebook builder. The Bs
  Python code imports the common local components, so it inherits the fix.
- `dimension5_correction.tex` is a derivation note for author review; it is not
  a claim that the manuscript has been changed.
- Before replacing numerical results, adopt the proposed change explicitly
  and regenerate the accepted ensembles, affected plots, and summaries in a
  separate output directory. Treat the previously identified charm shared-input
  sampling issue separately and record the effects of each change.
- For Referee 1 item 1, select and validate the OPE prescription used for the
  reported input comparison. The leading d5 derivation is complete; historical
  printed labels, higher-order terms, and the original uncertainty procedure
  remain separate provenance questions.

Reproduce from the project root:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noinit -noprompt -script review/mixed_condensate_resolution/derive_dimension5.wl
python3 -B review/mixed_condensate_resolution/check_covariant_traces.py
python3 -B review/mixed_condensate_resolution/evaluate_impact.py
python3 -B review/mixed_condensate_resolution/prepare_implementation_patch.py
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noinit -noprompt -script review/mixed_condensate_resolution/verify_proposed_wolfram_patch.wl
```

These commands write audit products in this review directory and do not apply
the proposed patch. JSON files retain exact numerical results and check status;
the CSV files retain all sampled points.
