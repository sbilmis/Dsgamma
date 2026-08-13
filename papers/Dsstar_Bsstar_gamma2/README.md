# Mixed-current \(D_{s1},B_{s1}\to D_s^\ast,B_s^\ast\gamma\)

This directory is a clean-room recalculation.  Nothing in the sibling
`Dsstar_Bsstar_gamma` directory is imported by the scripts or documents.

The four reported channels are

1. \(D_{s1}(2460)^+\to D_s^{\ast+}\gamma\);
2. \(D_{s1}(2536)^+\to D_s^{\ast+}\gamma\);
3. the predicted lower \(B_{s1}^{(L)0}\to B_s^{\ast0}\gamma\);
4. \(B_{s1}(5830)^0\to B_s^{\ast0}\gamma\).

The two basis currents are

\[
J_A^\mu=\bar s\gamma^\mu\gamma_5 Q,\qquad
J_B^\mu={i\over m_Q+m_s}\bar s\sigma^{\mu\nu}p'_\nu\gamma_5Q ,
\]

and the physical rotation is

\[
\binom{J_L}{J_H}=
\begin{pmatrix}\sin\theta&\cos\theta\\
\cos\theta&-\sin\theta\end{pmatrix}
\binom{J_A}{J_B}.
\]

## Reproduction

From this directory:

```bash
export PROJECT_PYTHON=/Users/sbilmis/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
"$PROJECT_PYTHON" scripts/run_analysis.py
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noinit -noprompt \
  -script scripts/symbolic_derivation.wl
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noinit -noprompt \
  -script scripts/mathematica_regression.wl
"$PROJECT_PYTHON" scripts/compare_regression.py
"$PROJECT_PYTHON" scripts/build_documents.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd tex/Dsstar_Bsstar_gamma2_paper.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd tex/Dsstar_Bsstar_gamma2_Mathematica_derivation.tex
```

The Python random seed is `20260729`.  The declared scan grids, tolerances,
and retained OPE terms are written to `outputs/csv/run_metadata.csv`.

## Final artifacts and numerical result

The publication-style paper and the detailed Wolfram Language derivation are
in `output/pdf/`.  The latter includes the complete 36-plot atlas and the full
symbolic and numerical Wolfram Language source listings.  The native notebook
is `notebooks/Dsstar_Bsstar_gamma2_derivation.nb`.

The quoted Monte Carlo medians and 16th--84th percentiles are:

| channel | \(g\) | \(\Gamma_\gamma\) [keV] |
|---|---:|---:|
| \(D_{s1}(2460)\) | \(-0.061_{-0.160}^{+0.174}\) | \(0.482_{-0.434}^{+1.609}\) |
| \(D_{s1}(2536)\) | \(-0.017_{-0.069}^{+0.058}\) | \(0.094_{-0.089}^{+0.520}\) |
| \(B_{s1}^{L}\) | \(-0.548_{-0.339}^{+0.320}\) | \(1.589_{-1.298}^{+2.703}\) |
| \(B_{s1}(5830)\) | \(-0.484_{-0.425}^{+0.286}\) | \(2.332_{-1.925}^{+5.735}\) |

These are exploratory predictions.  The \(s_0\) scans are mild, but the
transition Borel-mass variation is about 40% for the charm couplings and more
than 100% for the bottom couplings.  The higher-state A/B components also
partly cancel, producing broad, non-Gaussian width distributions.

Independent Python and native Wolfram Language values pass all 42 regression
entries.  The maximum material relative difference is \(9.13\times10^{-7}\).

## Method boundary

This is a heavy-light external-photon LCSR.  The final transition result uses
hard photon emission plus nonlocal two- and three-particle photon DAs through
twist four.  Ordinary local vacuum-condensate transition terms are excluded.
The often-used heavy-line-emission/local-\(\langle\bar ss\rangle\) term is
retained only in a diagnostic column and is never included in the quoted
couplings.

For the tensor basis current the retained transition OPE is its explicitly
derived leading-power projection,
\(\widehat T_B=m_Q\widehat T_A/(m_Q+m_s)\).  The twist-two trace establishes
this after double Borel transformation; the same leading-power reduction is
applied term by term through twist four.  Subleading tensor-current
\(1/m_Q\) kernels and radiative \(\alpha_s\) corrections are reported as
limitations and covered by a dedicated truncation uncertainty.
