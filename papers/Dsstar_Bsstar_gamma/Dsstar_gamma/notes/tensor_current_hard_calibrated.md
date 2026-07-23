# Hard Tensor-Current Contribution For \(D_{s1}\to D_s^\ast\gamma\)

This note completes the \(J_\mu^B\) hard-photon contribution in the
\(D_s^\ast\gamma\) channel at the same level of approximation used for the
soft tensor-current correction.

## Why Calibration Is Needed

The FeynCalc script

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
  -script Dsstar_gamma/scripts/hard_tensor_triangle_reduction.wl
```

projects the hard-photon triangle numerators onto the epsilon structure

\[
\epsilon_{\mu\nu\rho\sigma}\varepsilon^\rho q^\sigma .
\]

The output is

```text
Dsstar_gamma/outputs/hard_tensor_triangle_reduction.txt
```

The raw triangle cores contain constant terms, \(p\cdot q\) terms, and
\(1/(p\cdot q)\) terms.  This is not a problem by itself, but the published
axial-current hard spectral density for this channel is purely logarithmic:

\[
\rho_A^{\rm hard}(s)
=\frac{3m_sm_Q}{4\pi^2}
\left[
e_s\ln\frac{s-m_Q^2+m_s^2-\lambda^{1/2}(s,m_Q^2,m_s^2)}
{s-m_Q^2+m_s^2+\lambda^{1/2}(s,m_Q^2,m_s^2)}
+e_Q\ln\frac{s-m_s^2+m_Q^2-\lambda^{1/2}(s,m_s^2,m_Q^2)}
{s-m_s^2+m_Q^2+\lambda^{1/2}(s,m_s^2,m_Q^2)}
\right].
\]

Thus the \(p\cdot q\) and \(1/(p\cdot q)\) artifacts must cancel against
denominator-cancellation terms when the full hard diagram is reduced.  We use
this known axial result as the calibration condition.

## Calibrated Tensor Density

The script

```bash
python3 Dsstar_gamma/scripts/tensor_current_hard_calibrated.py
```

splits the known axial density into strange-emission and heavy-emission pieces,

\[
\rho_A^{\rm hard}=\rho_{A,s}^{\rm hard}+\rho_{A,Q}^{\rm hard}.
\]

Then it multiplies each line by the corresponding tensor/axial hard-core ratio
after removing the \(p\cdot q\)-artifact pieces that are absent in the
calibrated axial density:

\[
\rho_B^{\rm hard}(s)
=R_s^{\rm hard}(s)\rho_{A,s}^{\rm hard}(s)
+R_Q^{\rm hard}(s)\rho_{A,Q}^{\rm hard}(s).
\]

The diagonal double-dispersion identification \(p^2\to s\) is used in the
remaining log-producing terms, as in the earlier \(D_s\gamma\) hard-candidate
cross-check.

The hard contribution to the tensor-current basis coupling is then

\[
g_{B,{\rm hard}}^\ast
=
\frac{e^{(m_A^2+m_V^2)/(2M^2)}}
{m_A f_B\,m_V f_{D_s^\ast}}
\int_{(m_Q+m_s)^2}^{s_0}ds\,e^{-s/M^2}\rho_B^{\rm hard}(s).
\]

Finally,

\[
g_B^\ast=g_{B,{\rm soft}}^\ast+g_{B,{\rm hard}}^\ast.
\]

## Monte Carlo Results

The full calibrated scan uses the same uncertainty sampling as the soft-only
scan.  It writes

```text
Dsstar_gamma/outputs/tensor_full_calibrated_monte_carlo.csv
Dsstar_gamma/outputs/tensor_full_calibrated_monte_carlo_summary.csv
Dsstar_gamma/outputs/tensor_full_calibrated_monte_carlo_summary.txt
Dsstar_gamma/outputs/tensor_full_calibrated_monte_carlo_histograms.pdf
```

For fixed \(\theta=35.3^\circ\), the result is

\[
\begin{array}{c|cc}
\hbox{photon input}
& \Gamma_{2460}^{A+B}~({\rm keV})
& \Gamma_{2536}^{A+B}~({\rm keV})\\
\hline
\chi_s\langle\bar ss\rangle
& 24.2~[8.63,53.7]
& 9.74~[5.78,16.5]\\
f_{\gamma,s}^{\perp}\ {\rm lattice}
& 92.3~[59.7,138.7]
& 23.0~[16.5,35.1]
\end{array}
\]

When \(25^\circ\le\theta\le45^\circ\) is scanned,

\[
\begin{array}{c|cc}
\hbox{photon input}
& \Gamma_{2460}^{A+B}~({\rm keV})
& \Gamma_{2536}^{A+B}~({\rm keV})\\
\hline
\chi_s\langle\bar ss\rangle
& 23.9~[6.55,54.4]
& 9.80~[4.23,18.0]\\
f_{\gamma,s}^{\perp}\ {\rm lattice}
& 93.1~[57.1,145.1]
& 22.3~[10.6,41.8]
\end{array}
\]

The hard contribution is subleading but not negligible:

\[
g_{B,{\rm hard}}^\ast\simeq +0.16 ,
\]

while the soft tensor-current contribution is negative and larger in magnitude.
Therefore the hard piece reduces the soft-only widths, but the final scale
remains close to the expected literature values for the \(D_s^\ast\gamma\)
channels.

## Interpretation

This completes the \(J_\mu^B\) tensor-current contribution in the calibrated
LCSR setup used here:

- the soft photon-DA part is obtained directly from FeynCalc projection ratios;
- the hard photon-emission part is obtained by a line-by-line calibration to
  the known axial hard spectral density;
- both pieces are included in the final Monte Carlo scan.

The most conservative wording in a paper is to call this a **calibrated
hard-density implementation** and to explain that the axial channel is used as
the normalization check for the hard triangle reduction.

## Literature Table

A publication-ready comparison table for the vector-final-state channels is
stored in

```text
outputs/vector_final_state_literature_comparison_table.tex
outputs/vector_final_state_literature_comparison.csv
```

The table compares the present \(D_s^\ast\gamma\) results with the values
collected by Bondar--Milstein from
Godfrey, Goity--Roberts, Green--Repko--Radford, Radford--Repko--Saelim,
Chen et al., and Korner--Pirjol--Schilcher.
