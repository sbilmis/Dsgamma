# Roadmap for \(D_{s1}\to D_s^\ast\gamma\)

## 1. Physics Scope

The channels missing from the present \(D_s\gamma\) paper are

\[
\Gamma[D_{s1}(2460)\to D_s^\ast\gamma],
\qquad
\Gamma[D_{s1}(2536)\to D_s^\ast\gamma].
\]

These are not just a trivial replacement \(D_s\to D_s^\ast\).  The final meson
has spin one, so the correlator has two external vector polarizations: the
final \(D_s^\ast\) polarization and the photon polarization.  Consequently the
hadronic tensor has more independent gauge-invariant structures than the
\(1^+\to0^-\gamma\) case.

## 2. Whether To Integrate Into The Current Paper

My recommendation is:

- Do not add this immediately to the nearly finished \(D_s\gamma\) paper.
- Develop it in this folder as a companion calculation.
- After we have stable \(D_s^\ast\gamma\) results, decide whether to merge.

Reasons:

- The current paper already has a clean story: suppression of
  \(D_{s1}(2536)\to D_s\gamma\) and bottom-strange analogues.
- \(D_s^\ast\gamma\) requires a new hadronic decomposition and new traces.
- The literature comparison table from Bondar--Milstein includes these
  channels, so a complete companion paper can be well motivated.
- If both \(D_s\gamma\) and \(D_s^\ast\gamma\) are included in one paper, the
  theoretical section becomes much longer and less focused.

## 3. Correlator

Use the vacuum-to-photon correlator

\[
\Pi_{\mu\nu}^{(i)}(p,q)
= i\int d^4x\,e^{ip\cdot x}
\langle \gamma(q,\varepsilon)|
T\{J_\nu^{D_s^\ast}(x) J_{i\mu}^\dagger(0)\}
|0\rangle ,
\]

where

\[
J_\nu^{D_s^\ast}=\bar s\gamma_\nu Q,
\]

and

\[
J_{1\mu}=\sin\theta\,J_\mu^A+\cos\theta\,J_\mu^B,
\qquad
J_{2\mu}=\cos\theta\,J_\mu^A-\sin\theta\,J_\mu^B.
\]

The basis currents are the same as in the \(D_s\gamma\) paper:

\[
J_\mu^A=\bar s\gamma_\mu\gamma_5 Q,\qquad
J_\mu^B=\frac{i}{m_Q+m_s}\bar s\sigma_{\mu\alpha}P^\alpha\gamma_5 Q.
\]

## 4. Hadronic Matrix Elements

The vector decay constant is defined by

\[
\langle 0|J_\nu^{D_s^\ast}|D_s^\ast(p,\xi)\rangle
=m_{D_s^\ast}f_{D_s^\ast}\xi_\nu .
\]

For the initial axial state,

\[
\langle 0|J_{i\mu}|D_{s1}^{(i)}(P,\eta)\rangle
=m_i f_i \eta_\mu .
\]

The radiative matrix element

\[
\langle D_s^\ast(p,\xi)\gamma(q,\varepsilon)|D_{s1}(P,\eta)\rangle
\]

must be decomposed into gauge-invariant tensors.  This is the first real
derivation step.  We should not assume the same single form factor as in
\(1^+\to0^-\gamma\).

## 5. Width Formula

For \(1^+\to1^-\gamma\), the width is obtained by summing/averaging over the
initial axial, final vector, and photon polarizations:

\[
\Gamma =
\frac{|\vec q\,|}{8\pi m_A^2}\,
\frac{1}{3}
\sum_{\eta,\xi,\varepsilon}
|\mathcal M|^2,
\qquad
|\vec q\,|=
\frac{m_A^2-m_V^2}{2m_A}.
\]

The explicit expression depends on the chosen form-factor basis.  This must be
derived before quoting widths.

## 6. OPE Strategy

Reuse from the \(D_s\gamma\) analysis:

- external electromagnetic-field quark propagators;
- hard photon emission from heavy and strange lines;
- two-particle photon DAs;
- three-particle photon DAs;
- lattice \(f_{\gamma,s}^{\perp}\) normalization as preferred input;
- same mixing-angle scan \(25^\circ\le\theta\le45^\circ\).

Redo from scratch:

- Dirac traces with \(J_\nu^{D_s^\ast}=\bar s\gamma_\nu Q\);
- projection onto the selected \(1^+\to1^-\gamma\) tensor structures;
- hadronic normalization and decay-width formula;
- stability windows and uncertainty scan.

## 7. First Controlled Baseline

Stage 1 should include:

- basis-current amplitudes for \(J_\mu^A\) and \(J_\mu^B\);
- hard photon emission;
- leading two-particle photon DA contribution;
- strange-quark mass terms retained;
- no premature use of full three-particle terms until the tensor basis is
  tested.

Stage 2 then adds:

- all two-particle twist corrections;
- three-particle photon DAs;
- final Monte Carlo uncertainty analysis;
- comparison with Bondar--Milstein and other \(D_s^\ast\gamma\) predictions.

## 8. Inputs Needed

At minimum:

- \(m_{D_s^\ast}\);
- \(f_{D_s^\ast}\);
- \(m_{D_{s1}(2460)}\), \(m_{D_{s1}(2536)}\);
- physical-current residues \(f_1,f_2\) already used in the \(D_s\gamma\)
  paper;
- photon DA inputs;
- \(m_c,m_s,\langle\bar s s\rangle\), and \(f_{\gamma,s}^{\perp}\).

## 9. Paper Strategy

Possible title if separate:

> Radiative \(D_{s1}\to D_s^\ast\gamma\) transitions in light-cone QCD sum rules

Possible combined-paper strategy:

- Add one section after the \(D_s\gamma\) results.
- Present \(D_s^\ast\gamma\) only if the derivation is compact and the same
  mixing/cancellation story remains clear.

At this stage, a separate companion paper is cleaner.

