# Full Calibrated \(B_{s1}\to B_s^\ast\gamma\) Results

This note records the bottom-strange vector-final-state continuation of the
completed \(D_{s1}\to D_s^\ast\gamma\) calculation.

## Inputs

The final vector current is

\[
j_\nu^{B_s^\ast}=\bar s\gamma_\nu b,
\]

with

\[
\langle0|j_\nu^{B_s^\ast}|B_s^\ast(p,\xi)\rangle
=m_{B_s^\ast}f_{B_s^\ast}\xi_\nu .
\]

The vector decay constant is calculated from the same LO two-point
vector-current baseline used for \(D_s^\ast\):

\[
f_{B_s^\ast}=0.2615\pm0.0336~{\rm GeV}.
\]

The scan uses

\[
m_b=4.18\pm0.03~{\rm GeV},\qquad
m_s=0.093\pm0.011~{\rm GeV},
\]

\[
f_A(B_{s1})=0.305\pm0.027~{\rm GeV},\qquad
f_T(B_{s1})=0.285\pm0.027~{\rm GeV}.
\]

The effective tensor-current decay constant is

\[
f_B=f_T\,\frac{m_{B_{s1}}}{m_b+m_s}.
\]

The physical-current decay constants \(f_1\) and \(f_2\) are obtained from the
same diagonal two-point closure used in the earlier \(B_{s1}\to B_s\gamma\)
analysis.

## Sum Rule Content

The script

```bash
python3 Bsstar_gamma/scripts/bsstar_full_calibrated_monte_carlo.py
```

uses the completed \(D_s^\ast\gamma\) machinery with

\[
c\to b,\qquad e_c\to e_b=-\frac13 .
\]

It includes

- the axial-current \(g_A^\ast\) LCSR,
- the soft tensor-current contribution \(g_{B,{\rm soft}}^\ast\),
- the calibrated hard tensor-current contribution \(g_{B,{\rm hard}}^\ast\),
- the mixed-current normalization

\[
f_1g_{\rm low}^\ast
=\sin\theta\,f_Ag_A^\ast+\cos\theta\,f_Bg_B^\ast,
\]

\[
f_2g_{\rm high}^\ast
=\cos\theta\,f_Ag_A^\ast-\sin\theta\,f_Bg_B^\ast.
\]

The Borel and continuum windows are

\[
10\le M^2\le14~{\rm GeV}^2,
\]

\[
36.5\le s_0\le40.5~{\rm GeV}^2
\quad\hbox{for }B_{s1}(5750),
\]

\[
38.0\le s_0\le41.0~{\rm GeV}^2
\quad\hbox{for }B_{s1}(5830).
\]

## Results

For fixed \(\theta=35.3^\circ\), the preferred lattice photon-normalization
scenario gives

\[
\Gamma(B_{s1}(5750)\to B_s^\ast\gamma)
=10.2\,[6.63,16.24]~{\rm keV},
\]

\[
\Gamma(B_{s1}(5830)\to B_s^\ast\gamma)
=7.46\,[3.30,24.0]~{\rm keV}.
\]

When \(25^\circ\le\theta\le45^\circ\) is scanned,

\[
\Gamma(B_{s1}(5750)\to B_s^\ast\gamma)
=11.5\,[6.62,17.79]~{\rm keV},
\]

\[
\Gamma(B_{s1}(5830)\to B_s^\ast\gamma)
=2.55\,[0.805,9.75]~{\rm keV}.
\]

The legacy \(\chi_s\langle\bar ss\rangle\) scenario gives smaller central
values:

\[
\Gamma(B_{s1}(5750)\to B_s^\ast\gamma)
=1.89\,[0.471,5.11]~{\rm keV},
\]

\[
\Gamma(B_{s1}(5830)\to B_s^\ast\gamma)
=4.33\,[2.05,12.05]~{\rm keV}
\]

for fixed \(\theta=35.3^\circ\).

## Interpretation

The bottom-strange \(B_s^\ast\gamma\) widths are smaller than the
\(D_s^\ast\gamma\) widths mainly because of the reduced phase space and the
bottom-sector mixed-current normalization.  The \(B_{s1}(5830)\) distribution
has a long tail when \(f_2\) becomes small, so the percentile interval should be
quoted rather than the mean and standard deviation.

The outputs are

```text
Bsstar_gamma/outputs/bsstar_full_calibrated_monte_carlo.csv
Bsstar_gamma/outputs/bsstar_full_calibrated_monte_carlo_summary.csv
Bsstar_gamma/outputs/bsstar_full_calibrated_monte_carlo_summary.txt
Bsstar_gamma/outputs/bsstar_full_calibrated_histograms.pdf
```

The publication comparison table is stored in

```text
outputs/vector_final_state_literature_comparison_table.tex
outputs/vector_final_state_literature_comparison.csv
```

For the bottom rows the literature comparison uses Li--Ni--Zhong,
Godfrey--Moats--Swanson, Lu--Pan--Wang--Wang--Li, and Bondar--Milstein.
