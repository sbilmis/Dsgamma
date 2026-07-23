# Stage 1 Axial Benchmark For \(D_{s1}\to D_s^\ast\gamma\)

This note records the first numerical \(D_s^\ast\gamma\) result.  The purpose is
not yet to claim a final prediction, but to establish a controlled benchmark
against a known LCSR calculation before adding the tensor-current
\(J_\mu^B\) contribution.

## Correlator And Structure

For the axial basis current

\[
J_\mu^A=\bar s\gamma_\mu\gamma_5 Q
\]

and the final vector current

\[
J_\nu^{D_s^\ast}=\bar s\gamma_\nu Q,
\]

the vacuum-to-photon correlator is

\[
\Pi_{\mu\nu}^A(p,q)=i\int d^4x\,e^{ip\cdot x}
\langle\gamma(q,\varepsilon)|
T\{\bar s(x)\gamma_\nu Q(x)\,\bar Q(0)\gamma_\mu\gamma_5s(0)\}
|0\rangle .
\]

The Colangelo--De Fazio--Ozpineci calculation projects the invariant amplitude
proportional to

\[
\epsilon_{\mu\nu\rho\sigma}\varepsilon^\rho q^\sigma.
\]

At the hadronic level this corresponds to

\[
\langle \gamma(q,\varepsilon)D_s^\ast(p,\xi)|
D_{s1}(P,\eta)\rangle
=i e\,g^\ast\,
\epsilon_{\alpha\beta\rho\sigma}
\eta^\alpha \xi^{\ast\beta}\varepsilon^{\ast\rho}q^\sigma .
\]

The coupling \(g^\ast\) is dimensionless.

## Sum Rule Implemented

The implemented expression is Eq. (5.4) of
Colangelo, De Fazio and Ozpineci, Phys. Rev. D72, 074004 (2005), including

- the perturbative spectral density,
- the local heavy-quark photon-emission term,
- the leading-twist \(\chi\phi_\gamma\) term,
- the two-particle twist-4 term involving \(A(u_0)\), \(H_\gamma(u_0)\), and
  \(\bar H_\gamma(u_0)\),
- the two-particle \(f_{3\gamma}\) term,
- the three-particle \(F_2=S+\widetilde S+T_1-T_2-T_3+T_4\) term,
- the three-particle \(F_3=A+V\) term.

The script is

```bash
python3 Dsstar_gamma/scripts/stage1_axial_colangelo_gstar.py
```

and it writes

```text
Dsstar_gamma/outputs/stage1_axial_colangelo_gstar.csv
Dsstar_gamma/outputs/stage1_axial_colangelo_gstar_summary.txt
```

## Sign Check For \(\chi\)

The project convention stores

\[
\chi(1\,{\rm GeV})=+3.15~{\rm GeV}^{-2}
\]

as a positive magnitude.  The older printed sum-rule convention effectively
uses the opposite sign for the leading-twist photon matrix element.  Therefore
the script runs two diagnostics:

- `positive_chi`: our project convention;
- `printed_chi_sign_diagnostic`: literal insertion into the printed sign.

Only the project convention produces the expected cancellation.  The printed
sign diagnostic gives \(g_A^\ast\sim +1\), which is far too large and therefore
serves as a useful warning check.

## Benchmark Against Colangelo

Using Colangelo-like inputs and their original window

\[
4.0\le M^2\le 6.0~{\rm GeV}^2,\qquad
s_0=(2.50^2,2.55^2,2.60^2)~{\rm GeV}^2,
\]

the script obtains

\[
g_A^\ast=-0.192\ldots -0.141 .
\]

The quoted result in that paper is approximately

\[
g^\ast=-0.18\ldots -0.13 .
\]

This agreement is good enough for a normalization and sign benchmark.  The
remaining difference is consistent with input choices such as \(m_c\), \(m_s\),
\(f_{D_s^\ast}\), and the exact photon-DA parameter set.

## Project Axial-Only Baseline

For the project inputs and the rounded exploratory window

\[
3.0\le M^2\le 6.0~{\rm GeV}^2,\qquad
s_0=7.5,8.0,8.5~{\rm GeV}^2,
\]

with the calculated two-point input

\[
f_{D_s^\ast}=0.227\pm0.028~{\rm GeV},
\]

the basis-current result is

\[
g_A^\ast=-0.250\ldots -0.129 .
\]

At the central point

\[
M^2=4.5~{\rm GeV}^2,\qquad s_0=8.0~{\rm GeV}^2,
\]

the contribution breakdown is unchanged on the QCD side, while the extracted
coupling changes through the hadronic normalization.  The QCD-side contribution
breakdown is

\[
\begin{array}{c|r}
\hbox{term} & \hbox{QCD-side contribution}\\
\hline
\rho^{\rm pert} & +0.01989\\
e_c m_c\langle\bar ss\rangle & -0.00653\\
\chi\phi_\gamma & -0.05274\\
{\rm twist\mbox{-}4\ two\mbox{-}particle} & -0.00099\\
f_{3\gamma}\ {\rm two\mbox{-}particle} & +0.01380\\
F_2\ {\rm three\mbox{-}particle} & -0.00051\\
F_3\ {\rm three\mbox{-}particle} & +0.01380
\end{array}
\]

Thus the small value comes from the same cancellation emphasized in the
literature: the leading \(\chi\phi_\gamma\) term is partially cancelled by the
perturbative and \(f_{3\gamma}\) terms.

## Mapping To Physical Mixed Currents

At this stage only \(J_\mu^A\) is included.  Therefore the physical-current
mapping is only an axial baseline:

\[
f_1 g_1^\ast=\sin\theta\, f_A g_A^\ast,
\qquad
f_2 g_2^\ast=\cos\theta\, f_A g_A^\ast .
\]

Using

\[
\theta=35.3^\circ,\qquad f_A=0.225~{\rm GeV},\qquad
f_{D_s^\ast}=0.227~{\rm GeV},\qquad f_1=0.345~{\rm GeV},\qquad
f_2=0.379~{\rm GeV},
\]

the project-window axial-only ranges are

\[
g_1^\ast(2460)=-0.0940\ldots -0.0486,
\qquad
\Gamma_{2460}^{A{\rm -only}}=0.075\ldots0.282~{\rm keV},
\]

\[
g_2^\ast(2536)=-0.1209\ldots -0.0625,
\qquad
\Gamma_{2536}^{A{\rm -only}}=0.210\ldots0.786~{\rm keV}.
\]

At \(M^2=4.5~{\rm GeV}^2,\ s_0=8.0~{\rm GeV}^2\),

\[
g_1^\ast(2460)=-0.0606,\qquad
\Gamma_{2460}^{A{\rm -only}}=0.117~{\rm keV},
\]

\[
g_2^\ast(2536)=-0.0780,\qquad
\Gamma_{2536}^{A{\rm -only}}=0.327~{\rm keV}.
\]

## Width Formula

For

\[
\mathcal M=i e\,g^\ast S_\epsilon,\qquad
S_\epsilon=
\epsilon_{\alpha\beta\rho\sigma}
\eta^\alpha\xi^\beta\varepsilon^\rho q^\sigma ,
\]

the polarization-summed result from the FeynCalc/rest-frame check gives

\[
\Gamma=
4\pi\alpha_{\rm em}
\frac{m_A^2-m_V^2}{16\pi m_A^3}
\frac{(m_A^2-m_V^2)^2(m_A^2+m_V^2)}
{6m_A^2m_V^2}
|g^\ast|^2 .
\]

This reproduces the keV scale of the Colangelo benchmark.

## Next Calculation

The final mixed-current prediction requires the tensor-current contribution

\[
f_1 g_1^\ast
=\sin\theta\,f_A g_A^\ast+\cos\theta\,f_B g_B^\ast,
\]

\[
f_2 g_2^\ast
=\cos\theta\,f_A g_A^\ast-\sin\theta\,f_B g_B^\ast .
\]

The next step is therefore deriving \(g_B^\ast\) for
\(J_\mu^B=i\bar s\sigma_{\mu\alpha}P^\alpha\gamma_5Q/(m_Q+m_s)\), using the
FeynCalc trace skeleton already stored in
`Dsstar_gamma/outputs/current_building_blocks_dsstar.txt`.

## Axial-Only Monte Carlo Check

After calculating \(f_{D_s^\ast}\) from the two-point vector-current sum rule,
we performed a Monte Carlo uncertainty scan for the completed axial-current
piece.  The script is

```bash
python3 Dsstar_gamma/scripts/axial_only_monte_carlo.py
```

and it writes

```text
Dsstar_gamma/outputs/axial_only_monte_carlo.csv
Dsstar_gamma/outputs/axial_only_monte_carlo_summary.csv
Dsstar_gamma/outputs/axial_only_monte_carlo_summary.txt
Dsstar_gamma/outputs/axial_only_monte_carlo_histograms.pdf
```

The scan uses \(N=500\) samples per scenario and per mixing-angle treatment.  It
samples

\[
3.0\le M^2\le4.5~{\rm GeV}^2,\qquad
7.5\le s_0\le8.5~{\rm GeV}^2,
\]

as well as \(m_c\), \(m_s\), \(m_{D_{s1}}\), \(m_{D_s^\ast}\), \(f_A\),
\(f_{D_s^\ast}=0.227\pm0.028~{\rm GeV}\), \(f_1\), \(f_2\),
\(\langle\bar ss\rangle\), \(f_{3\gamma}\), \(\omega_\gamma^A\), and
\(\omega_\gamma^V\).  Two photon-normalization scenarios are shown:

- legacy \(\chi_s\langle\bar ss\rangle\);
- lattice \(f_{\gamma,s}^{\perp}\).

For fixed \(\theta=35.3^\circ\), the percentile results are

\[
\begin{array}{c|cc}
\hbox{photon input}
& \Gamma_{2460}^{A{\rm -only}}~({\rm keV})
& \Gamma_{2536}^{A{\rm -only}}~({\rm keV})\\
\hline
\chi_s\langle\bar ss\rangle
& 0.412~[0.049,1.389]
& 1.164~[0.146,3.988]\\
f_{\gamma,s}^{\perp}\ {\rm lattice}
& 1.718~[0.502,3.502]
& 4.990~[1.397,10.08]
\end{array}
\]

When the angle is scanned over \(25^\circ\le\theta\le45^\circ\), the
corresponding results are

\[
\begin{array}{c|cc}
\hbox{photon input}
& \Gamma_{2460}^{A{\rm -only}}~({\rm keV})
& \Gamma_{2536}^{A{\rm -only}}~({\rm keV})\\
\hline
\chi_s\langle\bar ss\rangle
& 0.410~[0.036,1.323]
& 1.228~[0.112,4.033]\\
f_{\gamma,s}^{\perp}\ {\rm lattice}
& 1.588~[0.422,3.745]
& 4.820~[1.366,10.82]
\end{array}
\]

These numbers should not be quoted as the final \(D_s^\ast\gamma\) prediction
until the tensor-current contribution \(g_B^\ast\) is included.  Their role is
to show that the axial-current part is numerically under control and that the
dominant uncertainty again comes from the leading photon normalization.
