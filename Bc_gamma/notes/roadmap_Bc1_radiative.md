# Roadmap for `Bc1 -> Bc^(*) gamma`

## 1. Scope

We treat the four radiative transitions
\[
B_{c1}(6743)\to B_c\gamma,\qquad
B_{c1}(6743)\to B_c^\ast\gamma,\qquad
B_{c1}(6750)\to B_c\gamma,\qquad
B_{c1}(6750)\to B_c^\ast\gamma .
\]

This should be developed as a separate paper from the strange-meson analysis.
The reason is structural: the \(D_s\) and \(B_s\) calculations contain a light
strange spectator and therefore have soft photon-distribution-amplitude (DA)
terms such as \(\chi_s\langle \bar s s\rangle\). In the \(B_c\) system both
constituents are heavy, so the leading QCD-sum-rule baseline is a heavy-heavy
hard-photon calculation. The photon couples perturbatively to the \(c\) and
\(\bar b\) lines.

## 2. Currents and states

For \(B_c^+\sim c\bar b\), a convenient current convention is
\[
j_5(x)=\bar b(x)i\gamma_5 c(x),\qquad
j_\nu^V(x)=\bar b(x)\gamma_\nu c(x).
\]
The two axial basis currents may be chosen in analogy with the strange-sector
analysis,
\[
J_\mu^A(x)=\bar b(x)\gamma_\mu\gamma_5 c(x),
\qquad
J_\mu^B(x)=\frac{i}{m_b+m_c}\bar b(x)\sigma_{\mu\alpha}P^\alpha\gamma_5c(x),
\]
where \(P\) is the axial-meson momentum. The physical currents are then written
as
\[
J_{1\mu}=\sin\theta\,J_\mu^A+\cos\theta\,J_\mu^B,\qquad
J_{2\mu}=\cos\theta\,J_\mu^A-\sin\theta\,J_\mu^B .
\]

The precise relation between these local currents and the nonrelativistic
\({}^1P_1\)-\({}^3P_1\) basis must be checked carefully. For this reason the
Stage-0 script below is only a quark-model sanity check, not the final sum-rule
result.

## 3. Correlation functions

For the pseudoscalar final state we use
\[
\Pi_\mu^{(5,i)}(p,q)=i\int d^4x\, e^{ip\cdot x}
\langle\gamma(q,\epsilon)|T\{j_5(x)J_{i\mu}^\dagger(0)\}|0\rangle ,
\]
and for the vector final state
\[
\Pi_{\mu\nu}^{(V,i)}(p,q)=i\int d^4x\, e^{ip\cdot x}
\langle\gamma(q,\epsilon)|T\{j_\nu^V(x)J_{i\mu}^\dagger(0)\}|0\rangle .
\]

At leading order the photon is attached to each heavy propagator. Schematically,
\[
S_Q^{(\gamma)}(x,0)=
-ie_Q\int d^4z\,S_Q^{(0)}(x,z)\gamma_\rho A^\rho(z)S_Q^{(0)}(z,0),
\]
or equivalently in a constant external electromagnetic field. There is no
replacement of a light-quark line by photon DAs in this heavy-heavy baseline.

## 4. Hadronic normalization

The hadronic sides require the decay constants
\[
\langle 0|j_5|B_c(p)\rangle =
\frac{m_{B_c}^2 f_{B_c}}{m_b+m_c},
\qquad
\langle 0|j_\nu^V|B_c^\ast(p,\eta)\rangle =
m_{B_c^\ast}f_{B_c^\ast}\eta_\nu,
\]
and
\[
\langle0|J_{1\mu}|B_{c1}(6743)(P,\eta)\rangle =
m_1 f_1\eta_\mu,\qquad
\langle0|J_{2\mu}|B_{c1}(6750)(P,\eta)\rangle =
m_2 f_2\eta_\mu .
\]

These constants should be obtained consistently from heavy-heavy two-point sum
rules, rather than imported from the strange-sector analysis.

## 5. Literature benchmarks

Bondar and Milstein quote the following comparison table for
\(B_{c1}\to B_c^{(*)}\gamma\) widths in keV:

| Channel | Godfrey [22] | Ebert et al. [23] | Fulcher [24] | Gershtein et al. [25] | Li et al. [26] | Li et al. [27] | Li et al. [28] | Eichten-Quigg [29] | Bondar-Milstein |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| \(B_{c1}(6743)\to B_c\gamma\) | 13 | 18.4 | 32.4 | 11.6 | 16.6 | 35 | 30.1 | 9.9 | 22 |
| \(B_{c1}(6743)\to B_c^\ast\gamma\) | 60 | 78.9 | 75.6 | 77.8 | 49 | 70 | 47.8 | 62.5 | 75 |
| \(B_{c1}(6750)\to B_c\gamma\) | 80 | 132 | 127.8 | 131.1 | 66.6 | 74 | 64 | 92.3 | 47 |
| \(B_{c1}(6750)\to B_c^\ast\gamma\) | 11 | 13.6 | 26.2 | 8.1 | 10.5 | 40 | 25.6 | 7.5 | 2 |

Bondar-Milstein use
\[
M(B_c^\ast)-M(B_c)=61~{\rm MeV},\quad
M(B_{c1}(6750))-M(B_{c1}(6743))=8~{\rm MeV},\quad
M(B_{c1}(6743))-M(B_c)=497~{\rm MeV}.
\]
They also state that there are currently no experimental data for these
radiative \(B_{c1}\) decays.

## 6. Immediate tasks

1. Derive the two heavy-heavy perturbative spectral densities for
   \(J_A\) and \(J_B\) into \(B_c\gamma\).
2. Derive the corresponding tensor decomposition and spectral densities for
   \(B_c^\ast\gamma\).
3. Build two-point sum rules for \(f_{B_c}\), \(f_{B_c^\ast}\), \(f_1\), and
   \(f_2\), or document which external inputs are used if we decide not to
   recalculate a given decay constant.
4. Choose \(M^2\) and \(s_0\) windows by pole dominance, OPE convergence, and
   form-factor stability.
5. Run uncertainty analysis over masses, thresholds, Borel parameters, decay
   constants, and the mixing angle.
6. Compare with the literature table above.

## 7. Current artifacts

The first FeynCalc trace pass is
`Bc_gamma/scripts/step1_hard_photon_numerators.wl`. It generates
`Bc_gamma/outputs/step1_hard_photon_numerators.txt`.

This pass keeps only free heavy-quark propagator numerators and inserts the
photon vertex on either the \(c\) line or the \(\bar b\) line. The term counts
from the successful run are

| Current/final channel | \(c\)-line terms | \(\bar b\)-line terms |
|---|---:|---:|
| \(J_A\to B_c\gamma\) | 14 | 19 |
| \(J_B\to B_c\gamma\) | 26 | 28 |
| \(J_A\to B_c^\ast\gamma\) | 16 | 17 |
| \(J_B\to B_c^\ast\gamma\) | 49 | 52 |

The next symbolic step is to project these numerator traces onto a minimal set
of gauge-invariant Lorentz structures and then perform the double-dispersion
reduction. For \(B_c^\ast\gamma\), the tensor basis is larger than for the
pseudoscalar final state, so this should be handled as a separate script rather
than folded into the \(B_c\gamma\) projector.

The auxiliary file `Bc_gamma/scripts/stage0_e1_benchmark.py` produces a
nonrelativistic E1 sanity check in `Bc_gamma/outputs/stage0_e1_benchmark.md`.
It is useful only as a hierarchy check. It should not be quoted as the QCD
sum-rule prediction.
