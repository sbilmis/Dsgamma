# Master prompt for mixed-current radiative QCD sum-rule studies

Copy the prompt below into a new task and replace every item in angle brackets.
Do not remove the validation gates.  The prompt is designed for studies such as
\(D_{s1}\to D_s^{(*)}\gamma\), \(B_{s1}\to B_s^{(*)}\gamma\), and radiative
\(B_c\) transitions, but it requires an explicit applicability decision before
the OPE is constructed.

---

## Copy-ready prompt

We are studying the radiative transition

\[
  \langle\text{INITIAL STATE, }J_i^{P_i},\text{ mass}\rangle
  \longrightarrow
  \langle\text{FINAL STATE, }J_f^{P_f},\text{ mass}\rangle+\gamma
\]

with valence flavours `<INITIAL FLAVOURS>` and `<FINAL FLAVOURS>`.  Use the
same disciplined framework as our \(D_s/B_s\gamma\) calculation, generalized
to the spins and flavours of this channel.  The final result must be
independently reproducible in Python and Mathematica.

### Inputs that must be fixed before calculation

Use or determine the following, and create an input-provenance table:

- transition and charge-conjugate convention: `<TRANSITION>`;
- quark flavours, masses, electric charges, and renormalization scale:
  `<QUARK INPUTS>`;
- interpolating currents, including all normalization factors:
  `<BASIS CURRENTS>`;
- physical-state rotation and sign convention: `<ROTATION>`;
- mixing angle and uncertainty: `<ANGLE AND SOURCE>`;
- whether the angle is an external result from an earlier study or is to be
  determined here;
- two-point OPE truncation: `<PERTURBATIVE ORDER AND LOCAL OPERATORS>`;
- transition OPE truncation: `<HARD, TWO-PARTICLE DA, THREE-PARTICLE DA>`;
- photon-DA set and scale: `<PHOTON DA INPUTS>`;
- preliminary Borel and continuum-threshold search ranges: `<SEARCH RANGES>`;
- experimental masses, total widths, branching fractions, limits, and their
  primary sources: `<EXPERIMENTAL INPUTS>`;
- literature calculations to be compared: `<REFERENCE LIST>`.

Never silently import a decay constant, mixing angle, overlap, or continuum
threshold.  Label every quantity as calculated here, externally fixed, fitted,
or used only as a diagnostic.

### 1. Applicability audit

Before deriving formulas, decide which framework applies.

1. For a heavy--light channel with a light valence quark, use the
   Rohrwild-style external-photon LCSR.  Include the light-quark mass unless a
   controlled limit is being tested.
2. Do not insert ordinary local vacuum-condensate propagator pieces as
   separate transition-LCSR terms.  In this scheme the relevant long-distance
   photon coupling is represented by nonlocal vacuum-to-photon matrix elements
   and photon distribution amplitudes.  Local \(d=3,d=5,\ldots\) condensates
   remain allowed and normally required in the separate two-point SVZ sum
   rules used to calculate residues.
3. For a double-heavy channel such as a valence \(B_c=\bar b c\) system, do
   not copy the light-quark photon-DA terms.  First determine whether the
   appropriate object is a background-field/two-point sum rule or a
   three-point sum rule dominated by hard photon emission from the two heavy
   lines.  State which nonperturbative operators survive and derive that OPE
   independently.
4. If the requested channel cannot consistently use the same LCSR, stop the
   algebraic reuse, explain the necessary modification, and continue with the
   correct channel-specific correlator.

Include this applicability decision explicitly in both final documents.

### 2. Conventions and theoretical-framework opening

Start the paper's theoretical framework from first principles:

- define \(p'=p+q\), with \(p'\) the initial momentum, \(p\) the final momentum,
  and \(q^2=0\);
- give the metric, \(\gamma_5\), \(\sigma_{\mu\nu}\), Levi--Civita, Fourier,
  electric-charge, and Borel-transform conventions;
- write every basis current \(J_A,J_B,\ldots\) explicitly;
- state the current dimensions and explain any mass normalization;
- write the physical-current rotation matrix \(R(\theta)\) and its inverse;
- define all decay constants/residues;
- parameterize the physical radiative matrix element using the complete
  gauge-invariant Lorentz basis appropriate to \(J_i^{P_i}\to J_f^{P_f}\gamma\);
- identify the multipole/form factors being calculated;
- derive the decay-width formula in the same convention.

Use \(p'\), not an ambiguous capital \(P\), for the initial momentum.  Keep a
one-page sign-and-normalization sheet and use it in every script.

### 3. Complete \(AA/AB/BB\) two-point calculation

For all basis currents \(i,j\in\{A,B\}\), start from

\[
 \Pi_{ij}^{\mu\nu}(p')=
 i\int d^4x\,e^{ip'\cdot x}
 \langle0|T\{J_i^\mu(x)J_j^{\nu\dagger}(0)\}|0\rangle .
\]

Perform and display the following steps.

1. Give the Wick contractions, momentum routing, colour factor, full Dirac
   traces, and projection onto the invariant used for the \(1^+\) pole.
2. Reduce the traces analytically with exact quark masses.
3. Derive the dispersion representation and the physical threshold
   \(s_{\rm th}\).
4. Derive and print fully explicit perturbative densities

   \[
     \rho_{AA}(s),\qquad \rho_{BB}(s),\qquad
     \rho_{AB}(s)=\rho_{BA}(s).
   \]

   Do not leave them as `DiracTrace`, `Im Pi`, an unevaluated Feynman-parameter
   integral, or an undefined helper function.  If a compact kernel is useful
   in the paper, also provide its completely expanded form in the detailed
   calculation document and as a machine-readable Mathematica expression.
5. Give the pre-Borel invariant amplitudes as

   \[
     \Pi_{ij}(p'^2)=
     \int_{s_{\rm th}}^\infty
       \frac{\rho_{ij}(s)}{s-p'^2}\,ds+
     \Pi_{ij}^{\rm local}(p'^2)
   \]

   with every retained local \(d=3,d=5,\ldots\) contribution explicit.
6. Apply the Borel transform term by term and show the continuum-subtracted
   result

   \[
     \widehat\Pi_{ij}(M^2,s_0)=
     \int_{s_{\rm th}}^{s_0}ds\,e^{-s/M^2}\rho_{ij}(s)
     +\widehat\Pi_{ij}^{\rm local}(M^2).
   \]

   Strictly speaking, the \(\rho_{ij}\) are the pre-Borel spectral densities;
   the requested “after-Borel expressions” are the explicit
   \(\widehat\Pi_{ij}\), including the exponential weights and local matrices.
7. Print the complete \(2\times2\) Borel matrix before rotation and the
   rotated matrix

   \[
     \widehat\Pi^{\rm phys}=R(\theta)
     \widehat\Pi^{\rm basis}R^T(\theta).
   \]
8. Write \(\widehat\Pi_{11}\), \(\widehat\Pi_{22}\), and
   \(\widehat\Pi_{12}\) explicitly, including every \(AA,BB,AB\) term.
9. Derive, rather than merely quote,

   \[
     f_i^2m_i^2e^{-m_i^2/M^2}
       =\widehat\Pi_{ii}(M^2,s_0^{(i)};\theta),
     \qquad
     f_i=
     \frac{e^{m_i^2/(2M^2)}}{m_i}
     \sqrt{\widehat\Pi_{ii}} .
   \]

   Adapt the pole normalization if the chosen current definition requires it.
10. Give explicit numerical substitutions for \(f_1\) and \(f_2\), not only
    their final values.

Use the externally established mixing angle as the nominal input when one
exists.  Matrix diagonalization may be reported as an OPE-truncation
diagnostic, but it must not silently replace that angle.

Report separate thresholds when the two physical poles require them.  Provide
mass-reproduction, pole-contribution, OPE-convergence, positivity,
off-diagonal-residual, Borel-window, threshold, and angle-stability tables.

### 4. Transition correlator and Wick contraction

Write the Rohrwild-style external-photon correlator for the selected channel,
with its operator ordering explicit:

\[
 \Pi_{\cdots}(p,q)=
 i\int d^4x\,e^{ip\cdot x}
 \langle\gamma(q,\varepsilon)|
 T\{J_{\rm initial}(x)J_{\rm final}^\dagger(0)\}
 |0\rangle .
\]

Then:

1. perform the Wick contractions with all global fermionic signs;
2. show photon emission from every electrically charged quark line;
3. give the momentum routing for each diagram;
4. write the light propagator in Rohrwild notation, including the light-quark
   mass, but without ordinary local vacuum-condensate pieces;
5. write the heavy propagator, electromagnetic background insertion, and
   background-gluon insertion explicitly;
6. list every retained two- and three-particle photon-DA matrix element in the
   same notation as the chosen reference;
7. warn prominently if any ordinary local transition condensate appears, and
   do not use it unless a channel-specific derivation proves that it belongs;
8. calculate and display the Dirac traces for every basis current;
9. decompose the result into a complete gauge-invariant Lorentz basis and
   identify the invariant amplitudes used for the physical form factors;
10. verify the Ward identity before numerical evaluation.

### 5. Explicit transition expressions before and after double Borel

For each basis current, give the complete pre-Borel invariant,

\[
 T_A(p^2,p'^2),\qquad T_B(p^2,p'^2),
\]

split into clearly named contributions, for example

\[
 T_i=T_i^{\rm pert}+T_i^{\rm tw2}+T_i^{\rm tw3}
     +T_i^{\rm tw4}+T_i^{\rm 3p},
\]

with notation adapted to the actual channel and DA basis.  Every term must be
explicit: denominators, Feynman parameters, DA arguments, masses, charges,
integration limits, and numerator factors.

Keep \(p^2\) and \(p'^2\) independent on the QCD side.  Do not replace them by
physical pole masses, and do not set \(p\cdot q\) to its on-shell hadronic
value before the double Borel transformation.

State and derive the double-Borel identities.  Define

\[
 u_0=\frac{M_1^2}{M_1^2+M_2^2},\qquad
 \mathcal M^2=\frac{M_1^2M_2^2}{M_1^2+M_2^2},
\]

and only set \(u_0=1/2\) after the general result has been obtained if equal
Borel parameters are chosen.

Print the fully reduced post-Borel expressions

\[
 \widehat T_A(M_1^2,M_2^2,s_0),\qquad
 \widehat T_B(M_1^2,M_2^2,s_0)
\]

term by term.  Do not leave a statement such as
`\(\mathcal B_{12}[T]\)` without evaluating it.  If expressions are too long
for the main text, use named kernels in the paper, but print the expanded
kernels in the appendix, detailed Mathematica PDF, notebook, and a plain-text
or `.wl` export.

Continuum subtraction must be stated for each class of term rather than
assumed universally.  In a Rohrwild-style external-photon prescription
define

\[
 E_Q=e^{-m_Q^2/\mathcal M^2},\qquad
 E_0=e^{-s_0/\mathcal M^2},\qquad
 \Delta E_Q=E_Q-E_0 .
\]

Distinguish the two stages explicitly:

- the raw double-Borel identities contain \(E_Q\);
- the continuum factor in the final sum rule is channel- and term-specific.

For the \(D_{s1}\to D_s\gamma\) framework used here, implement

\[
\text{twist-2}:\Delta E_Q,\qquad
\text{twist-3, twist-4, 3p gluonic}:E_Q,\qquad
\text{gauge-completion 3p electromagnetic}:\Delta E_Q .
\]

Do not infer a common factor merely because the raw kernels have the same
Borel exponential.  Derive or cite the subtraction prescription for every
class of term, and flag alternative duality prescriptions as systematic
variants.  Distinguish pole-side physical-mass substitutions from QCD-side
variables.

### 6. Physical projection, couplings, and widths

Project the basis amplitudes with exactly the same rotation used in the
two-point problem.  Give the analytic physical couplings explicitly in terms
of:

- \(\widehat T_A,\widehat T_B\);
- \(f_1,f_2\) obtained from the direct \(AA/AB/BB\) sum rule;
- final-state decay constant/residue;
- masses, Borel exponentials, charges, and mixing angle.

Display separately the individual \(A\)- and \(B\)-current contributions to
each physical amplitude so constructive or destructive interference is
visible before squaring.  Locate any cancellation angle analytically or
numerically, and show the width as a function of the mixing angle.

Calculate the form factors/couplings, partial widths, branching fractions when
the total width is known, and all propagated uncertainties.  Use a correlated
Monte Carlo where the same QCD input enters the two-point and transition sum
rules.

### 7. Stability-window selection and mandatory plots

Choose the working \(M^2\) and \(s_0\) windows from standard QCD sum-rule
criteria, not by selecting the narrowest interval that makes the result look
flat.

#### Window-selection rules

1. Determine the lower edge of the Borel window from OPE convergence.  Show
   the fractional contribution of every retained perturbative, condensate,
   twist, and three-particle term.  The highest-dimensional or highest-twist
   retained terms must remain subleading.  State the numerical convergence
   criterion used and justify it for the channel.
2. Determine the upper edge from pole dominance.  Define the pole contribution
   explicitly, for example

   \[
     {\rm PC}(M^2,s_0)=
     \frac{\displaystyle\int_{s_{\rm th}}^{s_0}
       ds\,e^{-s/M^2}\rho(s)+\widehat\Pi^{\rm local}}
     {\displaystyle\int_{s_{\rm th}}^\infty
       ds\,e^{-s/M^2}\rho(s)+\widehat\Pi^{\rm local}},
   \]

   with the expression adapted when local terms or a double dispersion
   relation require a different treatment.  Quote the adopted minimum pole
   fraction.  A typical \(40\%-50\%\) condition may be used as guidance, but it
   must be checked rather than imposed blindly.
3. Select \(s_0\) using the physical threshold, the expected first-excitation
   gap, mass reproduction, pole dominance, and stability.  Do not tune \(s_0\)
   solely to reproduce the desired coupling.
4. Use broad, rounded, physically readable working windows.  For example, if
   the diagnostics support approximately
   \(2.73\lesssim M^2\lesssim2.93~{\rm GeV}^2\), normally report and scan a
   transparent interval such as
   \(2.5\leq M^2\leq3.0~{\rm GeV}^2\), provided the endpoint diagnostics still
   pass.  Do not quote artificially precise endpoints unless they correspond
   to a genuine physical boundary.
5. Sample each scan on an evenly spaced grid.  Choose a rounded step appropriate
   to the scale, such as \(0.1\), \(0.25\), or \(0.5~{\rm GeV}^2\), and include
   at least six \(M^2\) points and at least five \(s_0\) points.  State the
   endpoints, step, and number of points in the caption and output table.
6. Use the same declared grid in Python and Mathematica.  Do not use a fine
   hidden optimization grid in one implementation and a different publication
   grid in the other without documenting both.
7. If rounded endpoints fail the convergence or pole-dominance requirements,
   move to the nearest rounded admissible endpoints and explain the choice.
   Physics criteria take precedence over cosmetic round numbers.

#### Mandatory one-dimensional scans

For every physical state and every reported coupling, produce:

1. \(G\) versus \(M^2\) at fixed central \(s_0\);
2. \(G\) versus \(s_0\) at fixed central \(M^2\);
3. \(\Gamma\) versus \(M^2\) at fixed central \(s_0\);
4. \(\Gamma\) versus \(s_0\) at fixed central \(M^2\).

The fixed value must be written in the axis caption or legend.  Show the full
accepted working interval, not only a small central portion.  Mark the central
point and, where useful, display curves for the lower, central, and upper
representative values of the fixed parameter.

Every plot must have:

- physical units on both axes;
- the state and transition in the title or caption;
- the fixed \(s_0\) or \(M^2\) value;
- the accepted window indicated by shading or vertical boundaries;
- an uncertainty band when input uncertainties are propagated;
- a companion CSV containing every plotted point;
- a quantitative stability summary giving the minimum, maximum, relative
  variation, and central slope over the accepted window.

Also produce the following diagnostic plots:

- pole contribution versus \(M^2\) for the tested \(s_0\) values;
- relative OPE/twist contribution versus \(M^2\);
- reproduced mass versus \(M^2\) for the two-point sum rule;
- \(f_1\) and \(f_2\) versus \(M^2\) and \(s_0\);
- coupling and width versus the mixing angle when mixed currents are used;
- separate \(A\)- and \(B\)-current amplitude components versus the mixing
  angle for interference-sensitive states.

Do not infer width stability from coupling stability alone.  Generate and
inspect the width scans independently because phase space, normalization, and
the squaring of a cancellation-sensitive amplitude can change the apparent
stability.

Include the principal four stability plots in the paper PDF, with additional
diagnostics in an appendix or supplement.  Include all plots, numerical grids,
selection criteria, and pass/fail diagnostics in the detailed Mathematica PDF.

### 8. Independent Python and Mathematica implementations

Maintain two independent calculation paths.

**Python**

- numerical integrations, threshold fits, window scans, Monte Carlo,
  uncertainty propagation, figures, and publication tables;
- save accepted samples and rejection reasons;
- use a fixed documented random seed for reproducibility.

**Mathematica**

- one readable `.nb` notebook that evaluates section by section;
- companion `.wl` scripts suitable for noninteractive regression;
- explicit currents, contractions, Dirac traces, spectral densities,
  pre-Borel expressions, Borel rules, post-Borel expressions, physical
  projection, \(f_1,f_2\), couplings, and widths;
- use FeynCalc where useful, but export fully reduced scalar expressions;
- do not import Python's final numbers as Mathematica “results.”

Create a term-by-term comparison table.  Compare spectral densities at test
points, every Borel-matrix entry, \(f_1,f_2\), every contribution to
\(\widehat T_A,\widehat T_B\), physical couplings, and widths.  State absolute
and relative tolerances.  A comparison is a pass only if every required entry
meets the declared tolerance.

### 9. Analytic and numerical validation gates

The work cannot be called complete until all applicable checks pass:

- \(\rho_{AB}=\rho_{BA}\);
- correct mass dimensions and physical thresholds;
- gauge invariance/Ward identities;
- cancellation of spurious Lorentz structures;
- pure-\(A\), pure-\(B\), zero-light-mass, equal-Borel, heavy-quark, and
  charge-switch limits where meaningful;
- reproduction of the corresponding Rohrwild or Colangelo formula in the
  precisely matched limit, without claiming that a partial comparison is a
  full reproduction;
- no physical pole mass inserted into a QCD-side pre-Borel numerator;
- no double counting of photon DAs and background-field propagator terms;
- no ordinary local condensate in the heavy--light transition LCSR;
- local condensates retained and identified in the two-point SVZ sum rule;
- Borel stability, threshold stability, pole dominance, OPE convergence, and
  mixing-angle sensitivity;
- rounded, evenly spaced \(M^2\) and \(s_0\) scans with the mandatory coupling
  and width stability plots;
- explicit interference table for every physical state;
- independent Python--Mathematica regression.

If a check fails, retain the failed result as a diagnostic, explain it, and do
not label the calculation final.

### 10. Literature and experimental comparison

Search current primary sources and verify bibliographic details.  Create a
consolidated table with one row per state and literature prediction.  Include
at least:

| State/transition | Observable | This work | Experimental status/value/limit | Literature result | Method | Reference/year | Compatibility/comment |
|---|---|---:|---|---:|---|---|---|

Distinguish direct measurements, upper limits, inferred limits, unobserved
channels, and model-state predictions.  Do not write “agrees with experiment”
when only an upper limit exists; write “compatible with the current limit.”
Explain large deviations through normalization, mixing, interference,
phase-space, DA inputs, or OPE content rather than only listing numbers.

### 11. Required deliverables

Produce two polished PDFs and all sources.

#### A. Paper PDF

Create `<CHANNEL>_paper.tex` and `<CHANNEL>_paper.pdf` with:

1. abstract;
2. introduction and experimental motivation;
3. theoretical framework beginning with conventions, currents, physical
   rotation, matrix element, and width;
4. transition correlator and OPE;
5. direct \(AA/AB/BB\) decay-constant calculation;
6. numerical inputs and provenance;
7. results, the four mandatory \(M^2/s_0\) coupling and width stability plots,
   pole/OPE diagnostics, interference discussion, and uncertainties;
8. literature/experimental comparison table;
9. conclusions;
10. appendices containing the explicit analytic
    \(\rho_{AA},\rho_{BB},\rho_{AB}\), pre-/post-Borel two-point formulas,
    explicit pre-/post-double-Borel transition formulas, and projection
    formulas.

The main text may use compact kernels, but each kernel must be expanded in an
appendix or an explicitly cited supplemental source.

#### B. Detailed Mathematica-calculation PDF

Create `<CHANNEL>_Mathematica_derivation.tex` and
`<CHANNEL>_Mathematica_derivation.pdf`, generated from and cross-referenced to
`<CHANNEL>_derivation.nb`.  It must show, step by step:

1. conventions and assumptions;
2. currents and correlators;
3. Wick contractions and momentum routing;
4. propagators and photon DAs;
5. every Dirac trace;
6. invariant projections and Ward checks;
7. explicit \(\rho_{AA},\rho_{BB},\rho_{AB}\);
8. pre-Borel and post-Borel two-point matrices;
9. explicit \(f_1,f_2\) formulas and numerical substitutions;
10. pre-Borel \(T_A,T_B\);
11. every double-Borel rule and fully reduced post-Borel
    \(\widehat T_A,\widehat T_B\);
12. physical rotation and interference;
13. numerical couplings and widths;
14. evenly spaced \(M^2\) and \(s_0\) grids, pole/OPE window diagnostics, and
    all mandatory coupling, width, residue, and mixing-angle stability plots;
15. Python--Mathematica comparison tables;
16. a final status table marking each analytic and numerical component as
    derived, checked, not applicable, or still open.

Also deliver:

- Python scripts and Mathematica `.nb`/`.wl` files;
- CSV tables for inputs, term breakdowns, every plotted \(M^2/s_0\) grid,
  accepted windows, Monte Carlo, literature, and experiment;
- PDF/PNG stability and interference plots;
- a README containing exact reproduction commands;
- a manifest mapping every number and table in the paper to the script and
  output file that generated it.

### 12. Final reporting rules

At the end, report:

- the physical results with uncertainties;
- whether \(f_1,f_2\) were fully derived or externally anchored;
- the selected angle and its provenance;
- which local and nonlocal contributions were included in each correlator;
- the dominant uncertainty and any cancellation sensitivity;
- the rounded \(M^2\) and \(s_0\) windows, grid increments, pole fractions,
  OPE-convergence criteria, and quantitative coupling/width stability;
- experimental compatibility;
- Python--Mathematica agreement;
- all remaining limitations, such as missing \(\alpha_s\) or higher-twist
  corrections.

Do not call the study complete merely because numerical values exist.  It is
complete only when the explicit analytic formulas, independent numerical
implementations, validation gates, two PDFs, notebook, tables, and
reproduction instructions are all present.

---

## Recommended use for the next projects

For \(D_s^{(*)}/B_s^{(*)}\) channels, retain the heavy--light branch and change
the final-state current and Lorentz/multipole decomposition.  Re-derive the
transition traces; do not infer them from the pseudoscalar channel by replacing
a polarization vector or mass.

For \(B_c\) channels, begin with the double-heavy applicability audit.  Reuse
the conventions, two-point-matrix logic, rotation, validation, software, and
document structure, but do not assume that the light-quark photon-DA OPE or its
condensate prescription survives.
