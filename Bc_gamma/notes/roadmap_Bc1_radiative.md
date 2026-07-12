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

This mixing is not optional in the physical amplitudes.  We derive the OPE for
the basis currents \(J_A\) and \(J_B\) only because the traces are simpler and
can be checked independently.  The physical invariant amplitudes are obtained
afterward by the same rotation,
\[
\widehat\Pi_1=\sin\theta\,\widehat\Pi_A+\cos\theta\,\widehat\Pi_B,\qquad
\widehat\Pi_2=\cos\theta\,\widehat\Pi_A-\sin\theta\,\widehat\Pi_B .
\]
For \(B_{c1}\to B_c\gamma\),
\[
g_i=
\frac{m_b+m_c}{m_{B_c}^2 f_{B_c}\,m_i f_i}
\exp\!\left(\frac{m_i^2}{M_1^2}+\frac{m_{B_c}^2}{M_2^2}\right)
\widehat\Pi_i^{(P)},
\]
and for \(B_{c1}\to B_c^\ast\gamma\),
\[
g_i^\ast=
\frac{1}{m_{B_c^\ast} f_{B_c^\ast}\,m_i f_i}
\exp\!\left(\frac{m_i^2}{M_1^2}+\frac{m_{B_c^\ast}^2}{M_2^2}\right)
\widehat\Pi_i^{(V)}.
\]
The pilot script applies this \(\sin\theta,\cos\theta\) rotation.  Its current
simplification is not the mixing; it is the normalization, where a common
benchmark \(f_{B_{c1}}=0.373~{\rm GeV}\) is used instead of separately
extracted physical \(f_1\) and \(f_2\).

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

## 3.1 Step-by-step derivation checklist

This subsection records the calculation in the order a student should repeat
it.

**Step 1: choose the physical currents.**  Work with the two basis currents
\(J_A,J_B\), but remember that the physical currents are the rotated
combinations \(J_1,J_2\).  The basis-current OPEs are intermediate objects; the
physical form factors use \(\widehat\Pi_1\) and \(\widehat\Pi_2\).

**Step 2: write the vacuum-to-photon correlators.**  For the \(B_c\gamma\)
channel the object is
\[
\Pi_\mu^{(5,K)}(p,q)=
i\int d^4x\,e^{ipx}
\langle\gamma(q,\epsilon)|
T\{\bar b(x)i\gamma_5 c(x)\,J_{K\mu}^\dagger(0)\}|0\rangle ,
\qquad K=A,B .
\]
For the \(B_c^\ast\gamma\) channel replace the pseudoscalar current by the
vector current:
\[
\Pi_{\mu\nu}^{(V,K)}(p,q)=
i\int d^4x\,e^{ipx}
\langle\gamma(q,\epsilon)|
T\{\bar b(x)\gamma_\nu c(x)\,J_{K\mu}^\dagger(0)\}|0\rangle .
\]

**Step 3: perform the Wick contractions.**  Since both valence constituents are
heavy, no light-cone photon DA is inserted.  Instead the photon is emitted from
the charm line or from the anti-bottom line:
\[
\Pi^{(K)}=\Pi^{(K,c{\rm -emission})}+\Pi^{(K,\bar b{\rm -emission})}.
\]
In momentum space the free propagator is
\[
S_Q^{(0)}(k)=\frac{i(\slashed{k}+m_Q)}{k^2-m_Q^2+i0},
\]
and the photon insertion is
\[
S_Q^{(\gamma)}(k,k+q)=
S_Q^{(0)}(k+q)\,i e_Q\slashed{\epsilon}\,S_Q^{(0)}(k).
\]
For the anti-bottom line the code uses the effective antiquark charge
\(e_{\bar b}=+1/3\), while \(e_c=+2/3\).

**Step 4: calculate the Dirac traces.**  For each final state and each basis
current calculate two traces:
\[
{\rm Tr}[J_K\,S_c^{(\gamma)}\,j_{\rm final}\,S_b^{(0)}],
\qquad
{\rm Tr}[J_K\,S_c^{(0)}\,j_{\rm final}\,S_b^{(\gamma)}].
\]
The FeynCalc files in `Bc_gamma/scripts` do precisely this and store the raw
and reduced projections in `Bc_gamma/outputs`.

**Step 5: project onto gauge-invariant structures.**  For \(B_c\gamma\) use
\[
S_{\mu\nu}=p_\nu q_\mu-(p\cdot q)g_{\mu\nu},
\qquad
{\cal P}^{\mu\nu}=\frac{S^{\mu\nu}}{2(p\cdot q)^2},
\]
so that \(S_{\mu\nu}{\cal P}^{\mu\nu}=1\).  For \(B_c^\ast\gamma\) the current
working tensor is
\[
V_{\mu\nu\rho}=p_\nu\epsilon_{\mu\rho p q}
-(p\cdot q)\epsilon_{\mu\nu\rho p},
\]
with \(V_{\mu\nu\rho}V^{\mu\nu\rho}=-4p^2(p\cdot q)^2\).  This vector-channel
projection is still a diagnostic basis; it must be matched to the physical
\(1^+\to1^-\gamma\) tensor basis before final publication.

**Step 6: reduce loop scalar products to denominators.**  The charm-emission
routing uses
\[
D_1=(k+q)^2-m_c^2,\quad D_2=k^2-m_c^2,\quad
D_3=(k-p)^2-m_b^2,
\]
while the anti-bottom-emission routing uses
\[
D_1=k^2-m_c^2,\quad D_2=(k-p)^2-m_b^2,\quad
D_3=(k-p+q)^2-m_b^2.
\]
Terms proportional to the full product \(D_1D_2D_3\) give the triangle
contribution.  Terms where a numerator cancels one denominator give
two-denominator/contact contributions.  These must not be silently dropped;
the present vector-channel pilot is marked unfinished partly for this reason.

**Step 7: write the double dispersion and Borel match.**  The true
three-point sum rule depends on \(p^2\) and \(P^2=(p+q)^2\).  The relation
\[
p\cdot q=\frac{P^2-p^2}{2}
\]
must be used before the double Borel transformation.  The pilot numerical
script uses a diagonal reduction \(M_1^2=M_2^2=2M^2\), which gives the
exponential factor used in the code,
\[
\exp\!\left[\frac{m_i^2+m_f^2}{2M^2}\right].
\]
This is acceptable for the first \(B_c\gamma\) hard baseline, but the
\(B_c^\ast\gamma\) channel needs a full double-Borel audit.

**Step 8: rotate to the physical states.**  After obtaining \(\widehat\Pi_A\)
and \(\widehat\Pi_B\), form
\[
\widehat\Pi_{6743}=\sin\theta\,\widehat\Pi_A+\cos\theta\,\widehat\Pi_B,
\qquad
\widehat\Pi_{6750}=\cos\theta\,\widehat\Pi_A-\sin\theta\,\widehat\Pi_B .
\]
The numerical pilot uses \(\theta=43.3^\circ\).

**Step 9: convert the form factor to a width.**  For the pseudoscalar final
state we use
\[
\Gamma(1^+\to0^-\gamma)=
\frac{\alpha_{\rm em}}{3}\,g^2 E_\gamma^3,\qquad
E_\gamma=\frac{m_i^2-m_f^2}{2m_i}.
\]
For the vector final state the current numbers are obtained from the diagnostic
polarization-summed \(V_{\mu\nu\rho}\) basis, so they should be treated as
checks, not final quoted widths.

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

For comparison with other predictions we should cite each calculation
separately.  Bondar and Milstein collect several model results for
\(B_{c1}\to B_c^{(*)}\gamma\); the original sources are listed explicitly
below.  All entries are widths in keV.

| Publication | \(6743\to B_c\gamma\) | \(6743\to B_c^\ast\gamma\) | \(6750\to B_c\gamma\) | \(6750\to B_c^\ast\gamma\) |
|---|---:|---:|---:|---:|
| Godfrey, Phys. Rev. D 70, 054017 (2004) | 13 | 60 | 80 | 11 |
| Ebert, Faustov, and Galkin, Phys. Rev. D 67, 014027 (2003) | 18.4 | 78.9 | 132 | 13.6 |
| Fulcher, Phys. Rev. D 60, 074006 (1999) | 32.4 | 75.6 | 127.8 | 26.2 |
| Gershtein, Kiselev, Likhoded, and Tkabladze, Phys. Rev. D 51, 3613 (1995) | 11.6 | 77.8 | 131.1 | 8.1 |
| T. Y. Li et al., Phys. Rev. D 108, 034019 (2023) | 16.6 | 49 | 66.6 | 10.5 |
| Q. Li et al., Phys. Rev. D 99, 096020 (2019) | 35 | 70 | 74 | 40 |
| X. J. Li et al., Eur. Phys. J. C 83, 1080 (2023) | 30.1 | 47.8 | 64 | 25.6 |
| Eichten and Quigg, Phys. Rev. D 99, 054025 (2018) | 9.9 | 62.5 | 92.3 | 7.5 |
| Bondar and Milstein, Phys. Rev. D 111, 114019 (2025) | 22 | 75 | 47 | 2 |

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

For a student reproducing the controlled \(B_c\gamma\) baseline, the minimal
workflow is:

1. Run `Bc_gamma/scripts/step1_hard_photon_numerators.wl` to generate the raw
   hard-photon traces.
2. Run `Bc_gamma/scripts/step2_ps_e1_triangle_reduction.wl` for the
   \(B_c\gamma\) E1 projector.
3. Run `Bc_gamma/scripts/step2_vec_epsilon_A_reduction.wl` and
   `Bc_gamma/scripts/step2_vec_epsilon_B_reduction.wl` for the diagnostic
   \(B_c^\ast\gamma\) tensor projection.
4. Check the projector normalizations:
   \(S_{\mu\nu}{\cal P}^{\mu\nu}=1\) and
   \(V_{\mu\nu\rho}{\cal P}^{\mu\nu\rho}=1\).
5. Insert the projected cores into
   `Bc_gamma/scripts/bc_ps_complete_analysis.py`.
6. Use \(M^2=\{7,8,9\}\,{\rm GeV}^2\), \(s_0=\{53,54,55\}\,{\rm GeV}^2\),
   \(m_b=4.18\,{\rm GeV}\), \(m_c=1.27\,{\rm GeV}\), and the
   \(B_c\)-mixing angle distribution to reproduce the current grid and Monte
   Carlo sample.
7. Extract \(f_1\) and \(f_2\) from the diagonalized two-point moments, not
   from a common \(f_{B_{c1}}\) benchmark.
8. Treat the \(B_c\gamma\) and \(B_c^\ast\gamma\) entries as final
   \(G^2_{\rm ap}\)-included hard-QCDSR predictions, with the crossed-contact
   rows retained only in the audit trail.

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

## 9. Controlled \(B_{c1}\to B_c\gamma\) pass

The script `Bc_gamma/scripts/bc_ps_complete_analysis.py` turns the compact
\(B_c\gamma\) triangle cores into the controlled pseudoscalar final-state
analysis.  The input window is inherited from the submitted \(B_c\)-mixing
analysis,
\[
M^2=\{7,8,9\}\ {\rm GeV}^2,\qquad s_0=\{53,54,55\}\ {\rm GeV}^2,
\]
with \(m_b=4.18~{\rm GeV}\), \(m_c=1.27~{\rm GeV}\), and
\(\theta=43.3^\circ\).  More precisely, the submitted \(B_c\)-mixing
Monte Carlo gives
\[
\theta_{\rm mean}=43.295^\circ,\qquad
\sigma_\theta=0.151^\circ,\qquad
\theta_{16\%}=43.147^\circ,\qquad
\theta_{84\%}=43.455^\circ .
\]
The Monte Carlo samples this mixing-angle distribution.  The pseudoscalar
decay-constant benchmark is taken as
\[
f_{B_c}=0.371\pm0.037~{\rm GeV}.
\]

The physical \(B_{c1}\) residues are no longer set equal to a common
benchmark.  They are extracted from the diagonalized two-point heavy-heavy
moments:
\[
\Pi_1=s_\theta^2\Pi_{AA}+2s_\theta c_\theta\Pi_{AB}
      +c_\theta^2\Pi_{BB},\qquad
\Pi_2=c_\theta^2\Pi_{AA}-2s_\theta c_\theta\Pi_{AB}
      +s_\theta^2\Pi_{BB},
\]
\[
f_i(M^2,s_0)=
\left[
\frac{e^{m_i^2/M^2}\Pi_i(M^2,s_0)}{m_i^2}
\right]^{1/2}.
\]
Using the perturbative heavy-heavy two-point spectral densities gives

| Quantity | Median | \(16\)--\(84\%\) | Grid range |
|---|---:|---:|---:|
| \(\theta\) | \(43.29^\circ\) | \(43.10\)--\(43.44^\circ\) | \(42.97\)--\(43.58^\circ\) |
| \(f_1\) | \(0.264~{\rm GeV}\) | \(0.259\)--\(0.269~{\rm GeV}\) | \(0.255\)--\(0.273~{\rm GeV}\) |
| \(f_2\) | \(0.722~{\rm GeV}\) | \(0.709\)--\(0.738~{\rm GeV}\) | \(0.695\)--\(0.750~{\rm GeV}\) |

After folding in the explicitly integrated all-positive \(G^2\) sector, the
Monte Carlo summary for the controlled \(B_c\gamma\) channel is

| Channel | Status | Coupling/projection coefficient | Width |
|---|---|---:|---:|
| \(B_{c1}(6743)\to B_c\gamma\) | final \(G^2_{\rm ap}\)-included result | \(+0.0487\,[0.0414,0.0572]~{\rm GeV}^{-1}\) | \(0.532\,[0.385,0.734]~{\rm keV}\) |
| \(B_{c1}(6750)\to B_c\gamma\) | final \(G^2_{\rm ap}\)-included result | \(-0.206\,[-0.230,-0.188]~{\rm GeV}^{-1}\) | \(10.04\,[8.29,12.45]~{\rm keV}\) |

For the vector final state we then recomputed the normalization instead of
using a fixed \(f_{B_c^\ast}\).  The vector two-point current is
\[
  j_\rho^V=\bar c\gamma_\rho b,\qquad
  \langle0|j_\rho^V|B_c^\ast(p,\eta)\rangle
  =m_{B_c^\ast}f_{B_c^\ast}\eta_\rho .
\]
With the same spin--1 invariant convention used in the mixing calculation,
the leading perturbative vector spectral density can be written as
\[
  \rho_V(s)=
  \frac{3}{16\pi^2}\frac{\sqrt{\lambda(s,m_b^2,m_c^2)}}{s}
  \left[
    -\frac{8(k\cdot p)^2}{s}-4m_c^2+12(k\cdot p)+12m_bm_c
  \right],
\]
where \(k\cdot p=(s+m_c^2-m_b^2)/2\).  In the equal-mass limit this reduces to
the expected \(S\)-wave vector-current form proportional to
\((s+2m_Q^2)\sqrt{1-4m_Q^2/s}\).  The sum rule is
\[
  f_{B_c^\ast}^2m_{B_c^\ast}^2e^{-m_{B_c^\ast}^2/M_V^2}
  =
  \int_{(m_b+m_c)^2}^{s_{0,V}} ds\,e^{-s/M_V^2}\rho_V(s).
\]
Using \(M_V^2=5\)--\(8~{\rm GeV}^2\) and
\(s_{0,V}=40\)--\(44~{\rm GeV}^2\), the direct spin-projector density gives
\[
  f_{B_c^\ast}^{\rm projector}=0.753\,[0.693,0.808]~{\rm GeV}.
\]
This number is a useful diagnostic because it uses exactly the projector
normalization of the \(B_c\)-mixing notes.  It is not the decay constant
normalization to insert into
\(\langle0|\bar c\gamma_\rho b|B_c^\ast\rangle=m_{B_c^\ast}f_{B_c^\ast}\eta_\rho\).
For that standard vector-invariant convention one divides the density by the
factor of three associated with the transverse spin trace.  The physical
two-point result becomes
\[
  f_{B_c^\ast}^{\rm standard}=0.435\,[0.400,0.467]~{\rm GeV}.
\]
This is close to the frequently used phenomenological scale
\(f_{B_c^\ast}\simeq0.384~{\rm GeV}\), and removes the apparent factor-of-two
normalization tension.  We therefore keep the standard result as the
self-contained QCDSR normalization and show the reference value as a
sensitivity check:

| Channel | Normalization | Coupling/projection coefficient | Width |
|---|---|---:|---:|
| \(B_{c1}(6743)\to B_c^\ast\gamma\) | standard vector-invariant \(f_{B_c^\ast}\), \(G^2_{\rm ap}\) included | \(-0.0764\,[-0.0821,-0.0717]~{\rm GeV}^{-1}\) | \(70.3\,[61.9,81.2]~{\rm keV}\) |
| \(B_{c1}(6750)\to B_c^\ast\gamma\) | standard vector-invariant \(f_{B_c^\ast}\), \(G^2_{\rm ap}\) included | \(+0.173\,[0.161,0.185]~{\rm GeV}^{-1}\) | \(382\,[331,437]~{\rm keV}\) |
| \(B_{c1}(6743)\to B_c^\ast\gamma\) | reference \(f_{B_c^\ast}=0.384\pm0.038~{\rm GeV}\) | \(-0.0877\,[-0.0986,-0.0795]~{\rm GeV}^{-1}\) | \(92.6\,[76.1,117]~{\rm keV}\) |
| \(B_{c1}(6750)\to B_c^\ast\gamma\) | reference \(f_{B_c^\ast}=0.384\pm0.038~{\rm GeV}\) | \(+0.187\,[0.169,0.210]~{\rm GeV}^{-1}\) | \(447\,[365,564]~{\rm keV}\) |

These \(B_c^\ast\gamma\) entries now include the explicit all-positive
radiative \(G^2\) sector.  The vector contact audit leaves no ordinary
two-channel contribution to add to the standard physical double-Borel sum.  The
projector-diagnostic
normalization gives \(24.8[21.8,28.5]\) and \(120[104,137]\) keV for the two
channels, respectively, and is kept in the CSV outputs only as a convention
cross-check.

The vector stability plots are produced by
`Bc_gamma/scripts/bc_vec_complete_analysis.py`.  The main output figures are

- `Bc_gamma/outputs/bc_vec_fbcstar_standard_M2_stability.pdf`;
- `Bc_gamma/outputs/bc_vec_standard_g1_M2_stability.pdf`;
- `Bc_gamma/outputs/bc_vec_standard_g2_M2_stability.pdf`;
- `Bc_gamma/outputs/bc_vec_standard_gamma1_M2_stability.pdf`;
- `Bc_gamma/outputs/bc_vec_standard_gamma2_M2_stability.pdf`.

The standard-normalization grid is stable enough to quote as a leading
baseline: the \(6743\) width varies from about \(66.7\) to \(80.9\) keV on the
central \(3\times3\) grid, while the \(6750\) width varies from about \(310\)
to \(406\) keV.  The larger \(6750\) width is therefore not a numerical
instability of the chosen Borel window; it follows from the leading vector
projection and mixing combination.

### 9.2 Superseded vector-channel \(G^2\)/contact screening diagnostic

For \(B_c^\ast\gamma\) the denominator of the sum rule contains both the
physical axial residue \(f_i\) and the final-state vector decay constant
\(f_{B_c^\ast}\).  Therefore the conservative condensate screen is slightly
larger than in \(B_c\gamma\).  Using the same completed two-point
\(G^2/{\rm pert}\) ratios from the \(B_c\)-mixing analysis, we estimate
\[
  \left|\frac{\delta g}{g}\right|
  \lesssim
  \left|\frac{\delta\Pi^{(3)}}{\Pi^{(3)}}\right|
  +
  \left|\frac{\delta f_i}{f_i}\right|
  +
  \left|\frac{\delta f_{B_c^\ast}}{f_{B_c^\ast}}\right|
  \simeq
  2\,\max\left|\frac{\Pi_{G^2}}{\Pi_{\rm pert}}\right|,
\]
and hence
\[
  \left|\frac{\delta\Gamma}{\Gamma}\right|
  \lesssim
  4\,\max\left|\frac{\Pi_{G^2}}{\Pi_{\rm pert}}\right|.
\]
This was only a screening envelope; it is kept here as a diagnostic that
motivated the explicit calculation.

For the preferred standard vector normalization the screen gives

| Channel | perturbative width | median \(G^2\)/contact width envelope | max envelope |
|---|---:|---:|---:|
| \(B_{c1}(6743)\to B_c^\ast\gamma\) | \(74.3[65.3,85.5]~{\rm keV}\) | \(\pm 11.1~{\rm keV}\) | \(\pm 27.7~{\rm keV}\) |
| \(B_{c1}(6750)\to B_c^\ast\gamma\) | \(360[312,412]~{\rm keV}\) | \(\pm 53.6~{\rm keV}\) | \(\pm 139~{\rm keV}\) |

For the reference \(f_{B_c^\ast}=0.384\pm0.038~{\rm GeV}\) normalization the
corresponding median envelopes are \(\pm14.3\) keV and \(\pm68.8\) keV.  These
numbers were large enough to justify the explicit all-positive \(G^2\)
calculation.  The final numbers are the \(G^2_{\rm ap}\)-included results quoted
above and in Sec. 12.5.

The script writes:

- `Bc_gamma/outputs/bc1_twopoint_residue_grid.csv`
- `Bc_gamma/outputs/bc1_twopoint_residue_summary.csv`
- `Bc_gamma/outputs/bc_ps_complete_grid.csv`
- `Bc_gamma/outputs/bc_ps_complete_grid_summary.csv`
- `Bc_gamma/outputs/bc_ps_complete_monte_carlo.csv`
- `Bc_gamma/outputs/bc_ps_complete_monte_carlo_summary.csv`
- `Bc_gamma/outputs/bc_ps_g1_M2_stability.pdf`
- `Bc_gamma/outputs/bc_ps_g2_M2_stability.pdf`
- `Bc_gamma/outputs/bc_ps_gamma1_mc_hist.pdf`
- `Bc_gamma/outputs/bc_ps_gamma2_mc_hist.pdf`
- `Bc_gamma/outputs/bc_vec_fbcstar_twopoint_grid.csv`
- `Bc_gamma/outputs/bc_vec_fbcstar_twopoint_summary.csv`
- `Bc_gamma/outputs/bc_vec_complete_grid.csv`
- `Bc_gamma/outputs/bc_vec_complete_grid_summary.csv`
- `Bc_gamma/outputs/bc_vec_complete_monte_carlo.csv`
- `Bc_gamma/outputs/bc_vec_complete_monte_carlo_summary.csv`
- `Bc_gamma/outputs/bc_vec_complete_summary.txt`

## 10. Comparison with experiment and other theory

### Experimental status

LHCb has observed a broad peaking structure in the \(B_c^+\gamma\) mass
spectrum with a significance above seven standard deviations.  A minimal
two-peak description gives
\[
M_1=6704.8\pm5.5\pm2.8\pm0.3~{\rm MeV},\qquad
M_2=6752.4\pm9.5\pm3.1\pm0.3~{\rm MeV},
\]
where the uncertainties are statistical, systematic, and from the \(B_c^+\)
mass.  LHCb does not measure the individual radiative partial widths.  The
natural widths are expected to be only a few hundred keV and are neglected in
the experimental fit because they are small compared with the photon-energy
resolution.

Therefore the direct comparison is qualitative:

| Quantity | Experiment | Present calculation |
|---|---:|---:|
| observed \(B_c(1P)^+\)-like structure | LHCb PRL 2025: yes, \(>7\sigma\) in \(B_c^+\gamma\) | channels studied explicitly |
| lower peak location | \(6704.8\pm6.2~{\rm MeV}\) total in quadrature, plus correlations | input state \(6743~{\rm MeV}\) |
| higher peak location | \(6752.4\pm10.0~{\rm MeV}\) total in quadrature, plus correlations | input state \(6750~{\rm MeV}\) |
| partial widths | not measured | \(B_c\gamma\): \(0.532\) and \(10.04~{\rm keV}\); \(B_c^\ast\gamma\): \(70.3\) and \(382~{\rm keV}\) |

The lower experimental peak should not be identified one-to-one with the
input \(6743~{\rm MeV}\) mass without care.  In the measured \(B_c\gamma\)
spectrum, decays through \(B_c^\ast\gamma\), followed by an unreconstructed
soft \(B_c^\ast\to B_c\gamma\), can shift visible peaks by the unknown
hyperfine splitting.

### Theory comparison

The comparison should cite each paper separately, not compress different
calculations into a single range.  The following table uses the values collected
in Bondar--Milstein, with the original model papers shown explicitly.

| Publication | \(6743\to B_c\gamma\) | \(6743\to B_c^\ast\gamma\) | \(6750\to B_c\gamma\) | \(6750\to B_c^\ast\gamma\) |
|---|---:|---:|---:|---:|
| Godfrey, Phys. Rev. D 70, 054017 (2004) | 13 | 60 | 80 | 11 |
| Ebert, Faustov, and Galkin, Phys. Rev. D 67, 014027 (2003) | 18.4 | 78.9 | 132 | 13.6 |
| Fulcher, Phys. Rev. D 60, 074006 (1999) | 32.4 | 75.6 | 127.8 | 26.2 |
| Gershtein, Kiselev, Likhoded, and Tkabladze, Phys. Rev. D 51, 3613 (1995) | 11.6 | 77.8 | 131.1 | 8.1 |
| T. Y. Li et al., Phys. Rev. D 108, 034019 (2023) | 16.6 | 49 | 66.6 | 10.5 |
| Q. Li et al., Phys. Rev. D 99, 096020 (2019) | 35 | 70 | 74 | 40 |
| X. J. Li et al., Eur. Phys. J. C 83, 1080 (2023) | 30.1 | 47.8 | 64 | 25.6 |
| Eichten and Quigg, Phys. Rev. D 99, 054025 (2018) | 9.9 | 62.5 | 92.3 | 7.5 |
| Bondar and Milstein, Phys. Rev. D 111, 114019 (2025) | 22 | 75 | 47 | 2 |
| This work | \(0.532[0.385,0.734]\) | \(70.3[61.9,81.2]\) | \(10.04[8.29,12.45]\) | \(382[331,437]\) |

A paper-style LaTeX version is

```latex
\begin{table}[t]
\centering
\caption{Comparison of the present \(B_{c1}\) radiative results with
experiment and representative model predictions.  Experimental peak locations
are from LHCb~\cite{LHCbBc1PObservation2025,LHCbBc1PStudy2025}.  The
comparison values are collected in Bondar--Milstein~\cite{BondarMilstein2025};
the original model papers are cited in the first column.}
\begin{tabular}{lcccc}
\toprule
Publication
& \(6743\to B_c\gamma\)
& \(6743\to B_c^\ast\gamma\)
& \(6750\to B_c\gamma\)
& \(6750\to B_c^\ast\gamma\)\\
\midrule
Godfrey~\cite{Godfrey2004} & 13 & 60 & 80 & 11\\
Ebert--Faustov--Galkin~\cite{Ebert2003} & 18.4 & 78.9 & 132 & 13.6\\
Fulcher~\cite{Fulcher1999} & 32.4 & 75.6 & 127.8 & 26.2\\
Gershtein et al.~\cite{Gershtein1995} & 11.6 & 77.8 & 131.1 & 8.1\\
T. Y. Li et al.~\cite{LiTang2023} & 16.6 & 49 & 66.6 & 10.5\\
Q. Li et al.~\cite{LiZhong2019} & 35 & 70 & 74 & 40\\
X. J. Li et al.~\cite{LiLiu2023} & 30.1 & 47.8 & 64 & 25.6\\
Eichten--Quigg~\cite{EichtenQuigg2018} & 9.9 & 62.5 & 92.3 & 7.5\\
Bondar--Milstein~\cite{BondarMilstein2025} & 22 & 75 & 47 & 2\\
This work & \(0.532[0.385,0.734]\) & \(70.3[61.9,81.2]\)
& \(10.04[8.29,12.45]\) & \(382[331,437]\)\\
\bottomrule
\end{tabular}
\label{tab:bc1-radiative-comparison}
\end{table}
```

The most important conclusion is the hierarchy, not point-by-point numerical
agreement.  The \(B_c\gamma\) and \(B_c^\ast\gamma\) channels are controlled
hard-QCDSR predictions including the explicit all-positive \(G^2\) sector.  The
contact-support audits leave no ordinary two-channel contact contribution to
add in the standard double-Borel projection.

## 11. Perturbative and condensate content

The current controlled \(B_c\gamma\) table contains the leading hard-photon
perturbative three-point baseline, normalized with perturbative two-point
physical residues.  More explicitly:

- The photon is emitted perturbatively from the \(c\) or \(\bar b\) heavy-quark
  line.
- The leading radiative three-point OPE uses free heavy-quark propagators.
- The perturbative spectral density is the hard triangle contribution generated
  by the projected Dirac traces.
- There is no photon DA, \(\chi\langle\bar q q\rangle\), or light-quark
  condensate term because the \(B_c\) system has no light valence quark.
- Heavy-heavy power corrections from the background-gluon expansion of the
  heavy propagators are included through the explicit all-positive \(G^2\)
  sector.
- Reduced two-denominator/contact terms from denominator cancellation have
  been audited in the \(B_c^\ast\gamma\) channel and do not give ordinary
  two-channel double-pole contributions.

So the present \(B_c\gamma\) and \(B_c^\ast\gamma\) results should be described
as controlled hard-QCDSR predictions with the explicit all-positive \(G^2\)
sector included.  The contact audits leave only crossed-contact or
single-virtuality candidates outside the standard two-channel projection.

## 12. Dimension-4 radiative gluon-condensate workbench

The heavy-heavy case has no photon DA or light-quark condensate contribution,
but the background-gluon expansion of the heavy propagators can generate
dimension-4 gluon-condensate terms in the radiative three-point correlator.  I
therefore split the calculation into small reproducible Wolfram/Python chunks
instead of asking FeynCalc to reduce the full expression at once.

The topology inventory is:

- photon emitted from the charm line:
  \(S_c(k+q)\,\gamma_\mathrm{em}\,S_c(k)\,\gamma_5\,S_b(k-p)\);
- photon emitted from the anti-bottom line:
  \(S_c(k)\,\gamma_5\,S_b(k-p)\,\gamma_\mathrm{em}\,S_b(k-p+q)\);
- six single-line insertions \(S_Q^{(G^2)}\), one on each heavy propagator
  segment in the two photon-emission diagrams;
- six open-gluon-pair insertions \(S^{(G)}S^{(G)}\), one for each pair of
  propagator segments inside each photon-emission diagram.

After projecting onto the pseudoscalar-final-state E1 tensor
\[
S_{\mu\nu}=p_\nu q_\mu-(p\cdot q)g_{\mu\nu},
\qquad
{\cal P}^{\mu\nu}={S^{\mu\nu}\over 2(p\cdot q)^2},
\]
the workbench gives 24 current-split projected rows: \(J_A\) and \(J_B\) for
each of the 12 gluon-condensate topologies.  The term-count audit is:

- single-line \(S_Q^{(G^2)}\): 6 topology rows, 12 current-split projections,
  221 projected terms;
- open-pair \(S^{(G)}S^{(G)}\): 6 topology rows, 12 current-split projections,
  164 projected terms;
- total: 24 projected rows and 385 projected terms.

The relevant scripts are:

- `Bc_gamma/scripts/step3_ps_g2_topology_inventory.py`;
- `Bc_gamma/scripts/step3_ps_g2_single_line_cphoton.wl`;
- `Bc_gamma/scripts/step3_ps_g2_single_line_bbarphoton.wl`;
- `Bc_gamma/scripts/step3_ps_g2_open_pair_cphoton.wl`;
- `Bc_gamma/scripts/step3_ps_g2_open_pair_bbarphoton.wl`;
- `Bc_gamma/scripts/step3_ps_g2_projection_summary.py`;
- `Bc_gamma/scripts/step3_ps_g2_denominator_bookkeeping.py`.

The corresponding raw and summary outputs are:

- `Bc_gamma/outputs/step3_ps_g2_topology_inventory.{txt,csv}`;
- `Bc_gamma/outputs/step3_ps_g2_single_line_cphoton.txt`;
- `Bc_gamma/outputs/step3_ps_g2_single_line_bbarphoton.txt`;
- `Bc_gamma/outputs/step3_ps_g2_open_pair_cphoton.txt`;
- `Bc_gamma/outputs/step3_ps_g2_open_pair_bbarphoton.txt`;
- `Bc_gamma/outputs/step3_ps_g2_projection_summary.{txt,csv}`;
- `Bc_gamma/outputs/step3_ps_g2_denominator_bookkeeping.{txt,csv}`.

The denominator bookkeeping uses Schwinger parameters \(a,b,r\).  For charm-line
photon emission,
\[
\Phi_c={r(a+b)\over A}p^2+{2ar\over A}(p\cdot q)
-(a+b)m_c^2-rm_b^2
={rb\over A}p^2+{ar\over A}P^2-(a+b)m_c^2-rm_b^2 ,
\]
where \(A=a+b+r\), \(a\) belongs to \(S_c(k+q)\), \(b\) to \(S_c(k)\), and
\(r\) to \(S_b(k-p)\).  For anti-bottom-line photon emission,
\[
\Phi_{\bar b}={a(b+r)\over A}p^2+{2r(b+r)\over A}(p\cdot q)
-a m_c^2-(b+r)m_b^2 ,
\]
or, after \(P^2=p^2+2p\cdot q\),
\[
\Phi_{\bar b}={(b+r)(a-r)\over A}p^2+{r(b+r)\over A}P^2
-a m_c^2-(b+r)m_b^2 .
\]
For \(S_Q^{(G^2)}\) the selected propagator power is four and the other two
powers are one.  For \(S^{(G)}S^{(G)}\), the selected pair has powers two and
two, while the third propagator has power one.

This completes the numerator-projection stage of the radiative
gluon-condensate workbench and attaches the denominator powers to every row.
It is not yet a numerical correction to the decay width.  The next required
step is double-dispersion/double-Borel reduction and numerical integration over
the same Borel and continuum-threshold windows used in the perturbative
baseline.

### 12.1 Conservative \(G^2\) size screen

Before doing the full derivative-delta reduction, I made a conservative size
screen using the already completed \(B_c\)-mixing OPE convergence table in the
same numerical window,
\[
M^2=\{7,8,9\}\ {\rm GeV}^2,\qquad s_0=\{53,54,55\}\ {\rm GeV}^2.
\]
The script is

- `Bc_gamma/scripts/bc_ps_g2_screening_estimate.py`,

and the outputs are

- `Bc_gamma/outputs/bc_ps_g2_screening_estimate.csv`;
- `Bc_gamma/outputs/bc_ps_g2_screening_estimate.txt`.

The screen uses
\[
\Delta_{\rm 2pt}^{G^2}=
\max\left(
\left|{\Pi_{AA}^{G^2}\over\Pi_{AA}^{\rm pert}}\right|,
\left|{\Pi_{AB}^{G^2}\over\Pi_{AB}^{\rm pert}}\right|,
\left|{\Pi_{BB}^{G^2}\over\Pi_{BB}^{\rm pert}}\right|
\right)
\]
from the \(B_c\)-mixing calculation.  Since the residue scales as
\(f\sim \sqrt{\Pi}\), the residue-normalization shift is estimated as
\(\Delta_f\simeq \Delta_{\rm 2pt}^{G^2}/2\).  The deliberately conservative
coupling envelope is then
\[
\left|{\delta g\over g}\right|_{\rm env}
=\Delta_{\rm 2pt}^{G^2}+{1\over2}\Delta_{\rm 2pt}^{G^2},
\qquad
\left|{\delta\Gamma\over\Gamma}\right|_{\rm env}
\simeq 2\left|{\delta g\over g}\right|_{\rm env}.
\]

The result is:

- median conservative \(|\delta g/g|\) envelope: \(0.058\);
- maximum conservative \(|\delta g/g|\) envelope: \(0.128\);
- median conservative \(|\delta\Gamma/\Gamma|\) envelope: \(0.115\);
- maximum conservative \(|\delta\Gamma/\Gamma|\) envelope: \(0.257\).

This is not small enough to justify silently dropping \(G^2\) at publication
level.  It does not prove that the radiative \(G^2\) correction is this large;
it says that, based on the same heavy-heavy OPE in the same Borel window, a
few-percent amplitude correction is plausible and should be either calculated
explicitly or bounded more sharply.  This diagnostic motivated the explicit
all-positive \(G^2\) integration in Sec. 12.5 and is superseded there.

The channel-by-channel summary is produced by

- `Bc_gamma/scripts/bc_ps_g2_corrected_summary.py`,

with outputs

- `Bc_gamma/outputs/bc_ps_g2_corrected_summary.csv`;
- `Bc_gamma/outputs/bc_ps_g2_corrected_summary.txt`.

For the current perturbative Monte Carlo central values, the conservative
\(G^2\) envelope is

| Channel | perturbative \(g\) \({\rm GeV}^{-1}\) | perturbative \(\Gamma\) keV | median \(G^2\) width envelope | max \(G^2\) width envelope |
|---|---:|---:|---:|---:|
| \(B_{c1}(6743)\to B_c\gamma\) | \(0.0316[0.0269,0.0370]\) | \(0.223[0.162,0.308]\) | \(\pm0.025\) keV | \(\pm0.060\) keV |
| \(B_{c1}(6750)\to B_c\gamma\) | \(-0.184[-0.205,-0.167]\) | \(8.01[6.60,9.90]\) | \(\pm0.914\) keV | \(\pm2.30\) keV |

This table should be read as a conservative systematic bound, not as the final
explicit radiative \(G^2\) correction.  It is enough, however, to decide that
we should not simply drop \(G^2\) without comment.

### 12.2 Denominator-reduced \(G^2\) inventory

I also generated the denominator-reduced form of the complete projected
radiative \(G^2\) numerator.  The scripts and outputs are

- `Bc_gamma/scripts/step3_ps_g2_denominator_reduction.wl`;
- `Bc_gamma/scripts/step3_ps_g2_denominator_reduction_summary.py`;
- `Bc_gamma/outputs/step3_ps_g2_denominator_reduction.csv`;
- `Bc_gamma/outputs/step3_ps_g2_denominator_reduction_summary.{txt,csv}`.

The reduction expands each projected numerator in inverse propagator
denominators \(d_1,d_2,d_3\), subtracting the numerator powers from the base
propagator powers.  The audit gives

- total denominator-reduced terms: \(747\);
- terms with all effective denominator powers positive: \(375\);
- contact/derivative-sector terms with at least one nonpositive effective
  denominator power: \(372\).

The contact/derivative sector is therefore not a small bookkeeping detail; it
is roughly half of the complete reduced expression.  A numerical result based
only on the all-positive sector would not be defensible without a contact
support audit.  The refined audit in Sec. 12.3 shows that these contact rows do
not have ordinary physical two-channel support, and Sec. 12.5 gives the final
all-positive \(G^2\) integration.

### 12.3 Vector \(B_c^\ast\gamma\) denominator-reduced \(G^2\) inventory

For the vector final state I repeated the same background-field numerator
reduction with the gauge-invariant structure
\[
  V_{\mu\nu\rho}
  =
  p_\nu\,\epsilon_{\mu\rho pq}
  -
  (p\cdot q)\,\epsilon_{\mu\nu\rho p},
  \qquad
  V_{\mu\nu\rho}V^{\mu\nu\rho}=-4p^2(p\cdot q)^2 .
\]
The working projector is \(V_{\mu\nu\rho}/[V_{\alpha\beta\sigma}
V^{\alpha\beta\sigma}]\).  This is the same tensor normalization used in the
leading vector hard-photon calculation.

The reproducible scripts are

- `Bc_gamma/scripts/step3_vec_g2_denominator_reduction.wl`;
- `Bc_gamma/scripts/step3_vec_g2_denominator_reduction_summary.py`.

The generated outputs are

- `Bc_gamma/outputs/step3_vec_g2_denominator_reduction.csv`;
- `Bc_gamma/outputs/step3_vec_g2_denominator_reduction_summary.csv`;
- `Bc_gamma/outputs/step3_vec_g2_denominator_reduction_summary.txt`.

The vector reduction gives

| Sector | Terms |
|---|---:|
| total reduced terms | \(1095\) |
| all-positive effective denominator powers | \(538\) |
| contact/derivative terms | \(557\) |

Broken down by current:

| Current | all-positive | contact/derivative |
|---|---:|---:|
| \(J^A_\mu\) | \(200\) | \(148\) |
| \(J^B_\mu\) | \(338\) | \(409\) |

Broken down by background-field class:

| Class | all-positive | contact/derivative |
|---|---:|---:|
| single-line \(G^2\) | \(209\) | \(368\) |
| open-pair \(GG\) | \(329\) | \(189\) |

This is an important qualitative result.  In the vector channel the contact
sector is not a correction to a correction; it is slightly larger than the
ordinary all-positive sector.  The tensor current \(J^B_\mu\) is especially
sensitive, with \(409\) contact/derivative rows.  The refined support audit
below shows that these rows do not enter the ordinary two-channel double-Borel
sum, so the final numerical correction is the all-positive sector of Sec. 12.5.

I then applied a support-level double-Borel audit using

- `Bc_gamma/scripts/step3_vec_g2_contact_support_audit.py`.
- `Bc_gamma/scripts/step3_vec_g2_contact_borel_admissibility.py`.

The outputs are

- `Bc_gamma/outputs/step3_vec_g2_contact_support_audit.csv`;
- `Bc_gamma/outputs/step3_vec_g2_contact_support_summary.csv`;
- `Bc_gamma/outputs/step3_vec_g2_contact_support_audit.txt`;
- `Bc_gamma/outputs/step3_vec_g2_contact_active_candidates.csv`.
- `Bc_gamma/outputs/step3_vec_g2_contact_borel_admissibility.csv`;
- `Bc_gamma/outputs/step3_vec_g2_contact_borel_admissibility_summary.csv`;
- `Bc_gamma/outputs/step3_vec_g2_contact_borel_admissibility.txt`.

This audit asks which contact rows still have denominator support in both
external virtualities after denominator cancellation.  The result is

| Contact support class | Rows |
|---|---:|
| double-Borel-active candidates | \(131\) |
| single-virtuality subtraction/contact rows | \(309\) |
| no-external-virtuality rows | \(117\) |

All \(131\) first-pass double-Borel-active candidates come from the
anti-bottom photon topology with remaining denominator support \(\{d_1,d_3\}\).
For this support the combined denominator depends on
\[
  (p-q)^2=2p^2-P^2,
\]
not on independent initial- and final-channel pole variables.  Therefore these
rows are crossed-contact candidates, not ordinary physical double-pole
contributions.  The refined Borel-admissibility audit gives

| Refined contact support class | Rows |
|---|---:|
| crossed single-combination \((p-q)^2\) rows | \(131\) |
| single-virtuality subtraction/contact rows | \(309\) |
| no-external-virtuality rows | \(117\) |
| ordinary two-channel contact rows | \(0\) |

Thus none of the \(557\) contact/derivative rows is added directly to the
physical double-Borel numerical sum.  The \(131\) crossed rows are kept in the
audit trail because they would require a separate crossed-contact prescription
if one wanted to study them beyond the standard two-channel sum-rule
projection.

For a student redoing this part, the checklist is:

1. Reproduce the perturbative vector projector check
   \(V_{\mu\nu\rho}V^{\mu\nu\rho}=-4p^2(p\cdot q)^2\).
2. Replace one or two propagator numerators by the background-field pieces
   \(S_Q^{(G^2)}\) and \(S_Q^{(G)}S_{Q'}^{(G)}\), keeping the two photon
   topologies separate.
3. Project onto \(V_{\mu\nu\rho}\).
4. Rewrite scalar products using the denominator rules for the charm-photon
   topology,
   \[
   k^2=d_2+m_c^2,\quad
   k\cdot q={d_1-d_2\over2},\quad
   k\cdot p={d_2+m_c^2+p^2-m_b^2-d_3\over2},
   \]
   and the corresponding anti-bottom-photon rules.
5. For each term, subtract numerator powers of \(d_1,d_2,d_3\) from the base
   propagator powers.  Positive effective powers are ordinary Schwinger terms;
   nonpositive powers are contact/derivative terms.
6. Run the support audit.  Rows with no external virtuality or only one
   virtuality vanish under the double Borel transform.
7. Run the refined Borel-admissibility audit.  The apparent \(131\) active
   rows reduce to crossed \((p-q)^2\) support and are stored in
   `step3_vec_g2_contact_active_candidates.csv` only as audit candidates.
8. Do not add those \(131\) crossed rows to the physical double-pole sum unless
   a separate crossed-contact prescription is introduced and justified.

### 12.5 Explicit all-positive \(G^2\) integration

After the contact-support audits, the ordinary physical double-Borel
dimension-four correction is the all-positive sector.  This means that, after
rewriting the projected numerator in powers of the inverse denominators
\(d_1,d_2,d_3\), each surviving term has
\[
  n_1>0,\qquad n_2>0,\qquad n_3>0 .
\]
Such terms have standard Schwinger representations and can be double-Borel
transformed without derivative-delta contact prescriptions.

The scripts used for this final step are

- `Bc_gamma/scripts/bc_ps_g2_all_positive_direct_borel.py`;
- `Bc_gamma/scripts/bc_ps_g2_all_positive_mc_update.py`;
- `Bc_gamma/scripts/bc_vec_g2_all_positive_direct_borel.py`;
- `Bc_gamma/scripts/bc_vec_g2_all_positive_mc_update.py`;
- `Bc_gamma/scripts/step3_ps_g2_contact_borel_admissibility.py`.

The corresponding output files are

- `Bc_gamma/outputs/bc_ps_g2_all_positive_direct_borel_grid.csv`;
- `Bc_gamma/outputs/bc_ps_g2_all_positive_monte_carlo_summary.csv`;
- `Bc_gamma/outputs/bc_vec_g2_all_positive_direct_borel_grid.csv`;
- `Bc_gamma/outputs/bc_vec_g2_all_positive_monte_carlo_summary.csv`;
- `Bc_gamma/outputs/step3_ps_g2_contact_borel_admissibility_summary.csv`.

The convention follows the submitted \(B_c\)-mixing analysis:
\[
  S_Q^{(G^2)}:\quad {G^2\over 12},\qquad
  S_Q^{(G)}S_{Q'}^{(G)}:\quad
  {G^2\over 96}{N_c^2-1\over 2},\qquad N_c=3 ,
\]
where \(G^2\) denotes the same gluon-condensate input used in the two-point
residue analysis.  The direct double-Borel check uses the diagonal choice
\[
  M_1^2=M_2^2=2M^2 ,
\]
so that the exponent agrees with the working three-point sum rule.  For charm
photon emission the active denominators are the charm-line denominators, while
for anti-bottom photon emission the active denominators are the bottom-line
denominators.  The final-state virtuality is fixed to
\(p^2=m_{B_c}^2\) or \(p^2=m_{B_c^\ast}^2\), and
\[
  p\cdot q={m_i^2-m_f^2\over 2}
\]
for each physical initial state \(i=1,2\).

The pseudoscalar contact audit gives

| Refined contact support class | Rows |
|---|---:|
| crossed single-combination rows | \(80\) |
| single-physical-virtuality rows | \(220\) |
| no-external-virtuality rows | \(72\) |
| ordinary two-channel contact rows | \(0\) |

Together with the vector audit, this means that no contact/derivative term is
added to the standard physical double-Borel sum.  Crossed rows are retained in
the audit trail only; they would require a separate crossed-contact prescription
and are not part of the present sum rule.

The all-positive direct-Borel correction shifts the form factor by

| Channel | median \(\delta g/g\) from \(G^2_{\rm ap}\) |
|---|---:|
| \(B_{c1}(6743)\to B_c\gamma\) | \(+0.543\,[+0.537,+0.546]\) |
| \(B_{c1}(6750)\to B_c\gamma\) | \(+0.120\,[+0.112,+0.130]\) |
| \(B_{c1}(6743)\to B_c^\ast\gamma\) | \(-0.0245\,[-0.0328,-0.0190]\) |
| \(B_{c1}(6750)\to B_c^\ast\gamma\) | \(+0.0298\,[+0.0264,+0.0347]\) |

Folding these shifts into the same Monte Carlo ensemble gives the final
preferred results

| Channel | \(g\,[{\rm GeV}^{-1}]\) | \(\Gamma\,[{\rm keV}]\) |
|---|---:|---:|
| \(B_{c1}(6743)\to B_c\gamma\) | \(+0.0487[0.0414,0.0572]\) | \(0.532[0.385,0.734]\) |
| \(B_{c1}(6750)\to B_c\gamma\) | \(-0.206[-0.230,-0.188]\) | \(10.04[8.29,12.45]\) |
| \(B_{c1}(6743)\to B_c^\ast\gamma\) | \(-0.0764[-0.0821,-0.0717]\) | \(70.3[61.9,81.2]\) |
| \(B_{c1}(6750)\to B_c^\ast\gamma\) | \(+0.173[0.161,0.185]\) | \(382[331,437]\) |

Thus the old \(G^2\) size screen should now be read only as a diagnostic that
motivated the explicit calculation.  It is superseded by the all-positive
\(G^2_{\rm ap}\) results above for the final numerical tables.

### 12.6 Vector projector and width-normalization audit

Because \(B_{c1}\to B_c^\ast\gamma\) uses the three-index tensor
\[
  V_{\mu\nu\rho}
  =p_\nu\epsilon_{\mu\rho pq}
  -(p\cdot q)\epsilon_{\mu\nu\rho p},
\]
the extracted coefficient is not dimensionally the same as the
\(B_c\gamma\) E1 coupling.  The audit script is

- `Bc_gamma/scripts/bc_vec_projector_width_audit.py`,

with outputs

- `Bc_gamma/outputs/bc_vec_projector_width_audit.csv`;
- `Bc_gamma/outputs/bc_vec_projector_width_audit.txt`.

The explicit polarization sum gives
\[
  {1\over3}\sum_{\rm pol}
  |\eta_A^\mu \epsilon_\gamma^\nu \eta_V^\rho V_{\mu\nu\rho}|^2
  =
  {2\over3}m_i^2E_\gamma^2(2m_f^2+E_\gamma^2).
\]
Therefore
\[
  \Gamma(1^+\to1^-\gamma)
  =
  {\alpha_{\rm em}\over3}\,
  g_V^2(2m_f^2+E_\gamma^2)E_\gamma^3 .
\]
Equivalently, one may define a one-momentum E1-normalized coupling
\[
  g_{E1}=\sqrt{2m_f^2+E_\gamma^2}\,g_V,
  \qquad
  \Gamma={\alpha_{\rm em}\over3}g_{E1}^2E_\gamma^3.
\]
For the two vector channels this conversion factor is \(8.969\).  Thus the
final vector couplings may also be quoted as

| Channel | \(g_V\,[{\rm GeV}^{-2}]\) | \(g_{E1}\,[{\rm GeV}^{-1}]\) |
|---|---:|---:|
| \(B_{c1}(6743)\to B_c^\ast\gamma\) | \(-0.0764[-0.0821,-0.0717]\) | \(-0.685[-0.736,-0.643]\) |
| \(B_{c1}(6750)\to B_c^\ast\gamma\) | \(+0.173[0.161,0.185]\) | \(+1.55[1.45,1.66]\) |

The projector normalization and width conversion are therefore internally
consistent.  They do not explain the large \(6750/6743\) vector-width ratio;
that ratio is driven by the mixed-current amplitude itself and must be checked
through the mixing/sign/current-assignment audit.

### 12.7 Vector mixing/sign and state-assignment audit

The vector hierarchy is the remaining nontrivial issue.  The audit script is

- `Bc_gamma/scripts/bc_vec_mixing_sign_audit.py`,

with outputs

- `Bc_gamma/outputs/bc_vec_mixing_sign_audit.txt`;
- `Bc_gamma/outputs/bc_vec_mixing_sign_audit_grid.csv`;
- `Bc_gamma/outputs/bc_vec_mixing_sign_audit_summary.csv`;
- `Bc_gamma/outputs/bc_vec_mixing_sign_audit_central.csv`;
- `Bc_gamma/outputs/bc_vec_mixing_sign_audit_components.csv`;
- `Bc_gamma/outputs/bc_vec_mixing_sign_audit_no_mixing.csv`;
- `Bc_gamma/outputs/bc_vec_mixing_sign_audit_literature.csv`.

At the central point \(M^2=8~{\rm GeV}^2\), \(s_0=54~{\rm GeV}^2\), the vector
basis amplitudes are
\[
  \widehat\Pi_A=0.00637,\qquad
  \widehat\Pi_B\simeq -0.00848 .
\]
With \(\theta=43.29^\circ\), the baseline physical combinations are
\[
\begin{aligned}
\widehat\Pi_{6743}^{(V)}
&=\sin\theta\,\widehat\Pi_A+\cos\theta\,\widehat\Pi_B\\
&=0.00437-0.00617=-0.00180,\\
\widehat\Pi_{6750}^{(V)}
&=\cos\theta\,\widehat\Pi_A-\sin\theta\,\widehat\Pi_B\\
&=0.00464+0.00582=0.01046 .
\end{aligned}
\]
Thus the lower vector amplitude is cancellation-dominated, while the higher
one is constructive.  The central baseline vector widths are
\[
  \Gamma_{6743}^{(V)}=70.2~{\rm keV},\qquad
  \Gamma_{6750}^{(V)}=338~{\rm keV},
\]
before the small all-positive \(G^2\) fold-in.

The variant scan gives the following nine-point grid medians:

| Convention test | \(\Gamma_{6743}\) keV | \(\Gamma_{6750}\) keV | high/low |
|---|---:|---:|---:|
| baseline | \(70.2\) | \(338\) | \(4.82\) |
| radiative \(J_B\) sign flip | \(2402\) | \(4.29\) | \(0.0018\) |
| alternate rotation \(c,s\) | \(29.8\) | \(344\) | \(11.5\) |
| alternate rotation \(c,-s\) | \(2361\) | \(10.1\) | \(0.0043\) |
| swap eigenstate assignment | \(316\) | \(75.3\) | \(0.239\) |
| swap assignment and \(J_B\) flip | \(3.99\) | \(2572\) | \(645\) |

The no-mixing check is especially useful.  Here the two basis currents are not
rotated into the physical mixed states; this is only a diagnostic, not a final
prediction.  At the same central point it gives

| No-mixing diagnostic | \(\Gamma_{6743}\) keV | \(\Gamma_{6750}\) keV | high/low |
|---|---:|---:|---:|
| pure \(J_A\) | \(878\) | \(126\) | \(0.143\) |
| pure \(J_B\) | \(1554\) | \(222\) | \(0.143\) |
| \(\theta=0\) baseline assignment: low \(=J_B\), high \(=J_A\) | \(1554\) | \(126\) | \(0.0808\) |
| \(\theta=90^\circ\) baseline assignment: low \(=J_A\), high \(=-J_B\) | \(878\) | \(222\) | \(0.253\) |

Thus the no-mixing limits have the same qualitative lower-state-dominant
pattern as the representative model literature.  The reversal in the baseline
mixed-current QCDSR result is therefore not a simple phase-space or mass
effect.  It is produced by destructive interference in the \(6743\) amplitude
and constructive interference in the \(6750\) amplitude.

Representative model calculations all have high/low \(<1\) for
\(B_c^\ast\gamma\): Godfrey \(0.18\), Ebert--Faustov--Galkin \(0.17\), Fulcher
\(0.35\), Gershtein et al. \(0.10\), T. Y. Li et al. \(0.21\), Q. Li et al.
\(0.57\), X. J. Li et al. \(0.54\), Eichten--Quigg \(0.12\), and
Bondar--Milstein \(0.027\).

The only simple convention changes that restore high/low \(<1\) either make the
lower width unrealistically large or require a relative radiative \(J_B\) sign
flip not supported by the tensor-current convention shared with the
two-point \(B_c\)-mixing code.  Therefore the \(B_c^\ast\gamma\) numbers should
not yet be called final physical predictions.  The correct paper wording is:
controlled \(B_c\gamma\) prediction, and \(B_c^\ast\gamma\) baseline under an
independent tensor-current sign and state-assignment audit.

## 13. Final comparison table artifact

The current final comparison table is stored as

- `Bc_gamma/notes/Bc1_final_comparison_table.tex`;
- `Bc_gamma/notes/Bc1_final_comparison_table.pdf`.

It contains:

1. the controlled \(B_{c1}\to B_c\gamma\) predictions including the explicit
   all-positive \(G^2\) correction;
2. the LHCb experimental status, emphasizing that the \(B_c^+\gamma\) structure
   is observed but the individual radiative partial widths are not measured;
3. separate theoretical values from Godfrey, Ebert--Faustov--Galkin, Fulcher,
   Gershtein et al., T. Y. Li et al., Q. Li et al., X. J. Li et al.,
   Eichten--Quigg, and Bondar--Milstein;
4. a separate status table for \(B_{c1}\to B_c^\ast\gamma\), where our QCDSR
   calculation uses the standard \(f_{B_c^\ast}\) normalization and includes
   the explicit all-positive \(G^2\) sector; the reference \(f_{B_c^\ast}\)
   normalization remains only a sensitivity cross-check.
