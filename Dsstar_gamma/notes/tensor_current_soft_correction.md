# Tensor-Current Soft Contribution For \(D_{s1}\to D_s^\ast\gamma\)

This note records the first \(J_\mu^B\) tensor-current contribution in the
\(D_s^\ast\gamma\) channel.  It should be read together with
`stage1_axial_colangelo_gstar.md`.

## Current

The second basis current is

\[
J_\mu^B=
\frac{i}{m_Q+m_s}\,
\bar s\sigma_{\mu\alpha}P^\alpha\gamma_5 Q .
\]

The physical mixed-current normalization is

\[
f_1g_1^\ast
=\sin\theta\,f_Ag_A^\ast+\cos\theta\,f_Bg_B^\ast ,
\]

\[
f_2g_2^\ast
=\cos\theta\,f_Ag_A^\ast-\sin\theta\,f_Bg_B^\ast .
\]

For the basis-current decay constant we use the same convention as in the
\(D_s\gamma\) analysis,

\[
f_B=f_T\,\frac{m_A}{m_Q+m_s}.
\]

## FeynCalc Projection

The script

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
  -script Dsstar_gamma/scripts/soft_2p_tensor_current_projection.wl
```

projects the soft two-particle photon-DA traces onto

\[
\epsilon_{\mu\nu\beta\sigma}q^\sigma ,
\]

where \(\mu\) is the initial axial-current index, \(\nu\) is the final vector
current index, and \(\beta\) is the photon-polarization index.

Since the direct norm of this structure is proportional to \(q^2=0\), the
projection uses the dual tensor \(\epsilon_{\mu\nu\beta\sigma}p^\sigma\).  This
gives a nonzero normalization proportional to \(p\cdot q\).

The output is

```text
Dsstar_gamma/outputs/soft_2p_tensor_current_projection.txt
```

The relevant projected traces are

\[
\Pi_A^{\rm tensor~DA}\propto -8m_Q ,
\]

and

\[
\frac{\Pi_B^{\rm tensor~DA}}{\Pi_A^{\rm tensor~DA}}
=
\frac{11(k\cdot p)+(k\cdot P)+7(p\cdot q)}
{12m_Q(m_Q+m_s)} .
\]

For the vector DA,

\[
\frac{\Pi_B^{\rm vector~DA}}{\Pi_A^{\rm vector~DA}}
=
\frac{m_Q(p\cdot P)}
{(m_Q+m_s)(k\cdot p)} .
\]

At the two-particle Borel saddle

\[
k=p+\bar u q,\qquad \bar u=1-u_0=\frac12,
\]

these become

\[
R_T=
\frac{12m_V^2+14(p\cdot q)}
{12m_Q(m_Q+m_s)},
\qquad
R_V=
\frac{m_Q(m_V^2+p\cdot q)}
{(m_Q+m_s)(m_V^2+\frac12p\cdot q)} .
\]

For the charm-strange central input,

\[
R_T\simeq3.1,\qquad R_V\simeq1.0 .
\]

## Numerical Implementation

The script

```bash
python3 Dsstar_gamma/scripts/tensor_current_soft_correction.py
```

uses the axial-current QCD-side pieces from the completed \(g_A^\ast\) sum rule
and constructs

\[
\Pi_B^{\rm soft}
=R_T\,\Pi_A^{\rm tensor~family}
+R_V\,\Pi_A^{f_{3\gamma}{\rm ~family}} .
\]

Here

\[
\Pi_A^{\rm tensor~family}
=\Pi_A^{\chi\phi_\gamma}
+\Pi_A^{{\rm twist}\mbox{-}4,2p}
+\Pi_A^{F_2},
\]

and

\[
\Pi_A^{f_{3\gamma}{\rm ~family}}
=\Pi_A^{f_{3\gamma},2p}
+\Pi_A^{F_3}.
\]

Then

\[
g_B^\ast=
\frac{e^{(m_A^2+m_V^2)/(2M^2)}}
{m_A f_B\,m_V f_{D_s^\ast}}
\Pi_B^{\rm soft}.
\]

The scan uses the same uncertainty sampling as the axial-only Monte Carlo,
including \(f_{D_s^\ast}=0.227\pm0.028~{\rm GeV}\).

## Results

For fixed \(\theta=35.3^\circ\), the soft \(A+B\) results are

\[
\begin{array}{c|cc}
\hbox{photon input}
& \Gamma_{2460}^{A+B_{\rm soft}}~({\rm keV})
& \Gamma_{2536}^{A+B_{\rm soft}}~({\rm keV})\\
\hline
\chi_s\langle\bar ss\rangle
& 35.0~[13.0,70.2]
& 15.6~[9.48,25.5]\\
f_{\gamma,s}^{\perp}\ {\rm lattice}
& 113.4~[71.9,161.2]
& 31.0~[21.9,47.6]
\end{array}
\]

When \(25^\circ\le\theta\le45^\circ\) is scanned,

\[
\begin{array}{c|cc}
\hbox{photon input}
& \Gamma_{2460}^{A+B_{\rm soft}}~({\rm keV})
& \Gamma_{2536}^{A+B_{\rm soft}}~({\rm keV})\\
\hline
\chi_s\langle\bar ss\rangle
& 35.8~[15.7,73.7]
& 14.8~[7.82,27.2]\\
f_{\gamma,s}^{\perp}\ {\rm lattice}
& 111.6~[71.2,169.3]
& 28.2~[15.7,51.6]
\end{array}
\]

The median projection ratios in the scan are

\[
R_T\simeq3.12,\qquad R_V\simeq1.01 .
\]

The lattice-input result is encouraging because it lands near the expected
literature scale for the \(D_s^\ast\gamma\) channels, unlike the axial-only
baseline.

## Relation To The Full Calibrated Result

This note was the intermediate soft-only tensor-current checkpoint.  The hard
photon-emission tensor-current contribution has now been added in
`tensor_current_hard_calibrated.md`.  The soft piece remains useful because it
shows where the dominant enhancement enters: the tensor photon-DA family is
multiplied by \(R_T\simeq3.1\).  The hard contribution is subleading and partly
reduces the soft-only widths.
