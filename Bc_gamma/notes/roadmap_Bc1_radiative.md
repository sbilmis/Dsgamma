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

## 8. Stage-1 symbolic projection results

### 8.1 \(B_{c1}\to B_c\gamma\)

The script `Bc_gamma/scripts/step2_ps_e1_triangle_reduction.wl` projects the
pseudoscalar final-state correlator onto
\[
S_{\mu\nu}=p_\nu q_\mu-(p\cdot q)g_{\mu\nu},\qquad
P=p+q ,
\]
with projector \(S_{\mu\nu}/[2(p\cdot q)^2]\).  The FeynCalc check gives
\[
S_{\mu\nu}{\cal P}^{\mu\nu}=1.
\]
We denote \(p^2=p_2\) and \(p\cdot q=pq\).  The projected triangle-core
numerators, after setting the inverse denominators to zero, are
\[
\begin{aligned}
C_A^{(c)} &=
-4i\,m_c+\frac{4i\,m_bm_c^2}{pq}-\frac{4i\,m_c^3}{pq},\\
C_A^{(\bar b)} &=
-4i\,m_b+\frac{4i\,m_b^3}{pq}-\frac{4i\,m_b^2m_c}{pq},
\end{aligned}
\]
and therefore
\[
C_A=e_c C_A^{(c)}+e_{\bar b}C_A^{(\bar b)} .
\]
For the tensor current,
\[
\begin{aligned}
C_B^{(c)} &=
\frac{4i\,m_bm_c-8i\,m_c^2}{m_b+m_c}
-\frac{4i\,m_c^2p_2}{(m_b+m_c)pq},\\
C_B^{(\bar b)} &=
-\frac{4i\,m_bm_c}{m_b+m_c}
-\frac{4i\,m_b^2p_2}{(m_b+m_c)pq},
\end{aligned}
\]
and
\[
C_B=e_c C_B^{(c)}+e_{\bar b}C_B^{(\bar b)} .
\]
The same output file also stores the denominator-cancellation terms. These
terms should not be discarded automatically; they correspond to reduced
two-denominator/contact contributions and must be either included in the
dispersion analysis or shown to vanish under the chosen Borel projection.

### 8.2 \(B_{c1}\to B_c^\ast\gamma\)

For the vector final state the trace contains Levi-Civita structures.  The
working gauge-invariant tensor is
\[
V_{\mu\nu\rho}
=p_\nu\,\epsilon_{\mu\rho p q}
-(p\cdot q)\,\epsilon_{\mu\nu\rho p}.
\]
It is transverse in the photon index for \(q^2=0\).  FeynCalc gives
\[
V_{\mu\nu\rho}V^{\mu\nu\rho}=-4p_2(pq)^2,\qquad
V_{\mu\nu\rho}{\cal P}^{\mu\nu\rho}=1 .
\]

The split scripts
`Bc_gamma/scripts/step2_vec_epsilon_A_reduction.wl` and
`Bc_gamma/scripts/step2_vec_epsilon_B_reduction.wl` should be used instead of
the all-in-one vector script because the tensor-current trace is much larger.
The projected \(J_A\) triangle cores are
\[
\begin{aligned}
C_{A,V}^{(c)} &=
-\frac{4i\,m_bm_c}{p_2}
+\frac{2i\,m_c^2}{p_2}
+\frac{4i\,m_c^2}{pq},\\
C_{A,V}^{(\bar b)} &=
-\frac{2i\,m_b^2}{p_2}
+\frac{4i\,m_bm_c}{p_2}
+\frac{4i\,m_b^2}{pq}.
\end{aligned}
\]
The tensor-current vector cores are
\[
\begin{aligned}
C_{B,V}^{(c)} &=
\frac{2i\,m_c}{m_b+m_c}
+\frac{2i\,m_b^2m_c-4i\,m_bm_c^2+2i\,m_c^3}
       {(m_b+m_c)p_2}\\
&\quad
-\frac{4i\,m_bm_c^2-4i\,m_c^3}{(m_b+m_c)pq}
+\frac{2i\,m_c\,pq}{(m_b+m_c)p_2},\\
C_{B,V}^{(\bar b)} &=
\frac{2i\,m_b}{m_b+m_c}
+\frac{-6i\,m_b^3+4i\,m_b^2m_c+2i\,m_bm_c^2}
       {(m_b+m_c)p_2}\\
&\quad
-\frac{4i\,m_b^3-4i\,m_b^2m_c}{(m_b+m_c)pq}
+\frac{2i\,m_b\,pq}{(m_b+m_c)p_2}.
\end{aligned}
\]
Again \(C_{A,V}=e_c C_{A,V}^{(c)}+e_{\bar b}C_{A,V}^{(\bar b)}\) and
\(C_{B,V}=e_c C_{B,V}^{(c)}+e_{\bar b}C_{B,V}^{(\bar b)}\).

In the double-dispersion variables one should use
\[
P^2=p^2+2p\cdot q,\qquad pq=\frac{P^2-p^2}{2},
\]
before Borel transformation.  The diagonal single-variable reduction used in
the strange-sector scripts is a useful calibration tool, but for \(B_c^\ast\)
the additional \(1/p_2\) and \(pq/p_2\) structures mean that the double
dispersion should be written explicitly before numerical work.

## 9. First numerical hard-photon pass

The script `Bc_gamma/scripts/bc_radiative_pilot_sumrule.py` turns the compact
triangle cores into a first numerical estimate.  The input window is inherited
from the submitted \(B_c\)-mixing analysis,
\[
M^2=\{7,8,9\}\ {\rm GeV}^2,\qquad s_0=\{53,54,55\}\ {\rm GeV}^2,
\]
with \(m_b=4.18~{\rm GeV}\), \(m_c=1.27~{\rm GeV}\), and
\(\theta=43.3^\circ\).  The decay constants used in this first pass are
\[
f_{B_c}=0.371~{\rm GeV},\qquad
f_{B_c^\ast}=0.384~{\rm GeV},\qquad
f_{B_{c1}}=0.373~{\rm GeV}.
\]
The first two are standard external inputs/benchmarks, while the last is a
common axial-vector benchmark value.  We have not yet extracted separate
physical \(f_1\) and \(f_2\) residues from a two-point \(B_c\) mixed-current
sum rule in this folder.

The resulting grid summary is

| Channel | Status | Coupling/projection coefficient | Width |
|---|---|---:|---:|
| \(B_{c1}(6743)\to B_c\gamma\) | hard-QCDSR baseline | \(+0.0219\,[0.0209,0.0232]~{\rm GeV}^{-1}\) | \(0.108\,[0.098,0.121]~{\rm keV}\) |
| \(B_{c1}(6750)\to B_c\gamma\) | hard-QCDSR baseline | \(-0.354\,[ -0.387,-0.331]~{\rm GeV}^{-1}\) | \(29.6\,[25.9,35.3]~{\rm keV}\) |
| \(B_{c1}(6743)\to B_c^\ast\gamma\) | diagonal vector-pilot | \(-0.0620\,[-0.0669,-0.0583]\) | \(46.3\,[41.0,53.8]~{\rm keV}\) |
| \(B_{c1}(6750)\to B_c^\ast\gamma\) | diagonal vector-pilot | \(+0.362\,[0.335,0.401]\) | \(1.67\,[1.43,2.04]\times10^3~{\rm keV}\) |

Only the \(B_c\gamma\) entries should be regarded as a controlled hard-photon
baseline at this stage.  The \(B_c^\ast\gamma\) entries are diagnostic because
they use the \(V_{\mu\nu\rho}\) gauge-invariant projection and a diagonal
single-variable reduction.  Before they are promoted to paper results we must
match the \(V_{\mu\nu\rho}\) coefficient to the physical \(1^+\to1^-\gamma\)
tensor basis and perform the full double-Borel audit, including the
denominator-cancellation/contact terms.
