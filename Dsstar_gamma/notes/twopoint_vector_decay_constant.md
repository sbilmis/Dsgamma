# Two-Point Calculation Of \(f_{D_s^\ast}\)

The previous axial \(D_s^\ast\gamma\) benchmark used a temporary
\(f_{D_s^\ast}\) input.  We now calculate it from a local two-point SVZ sum rule
so that the radiative normalization is self-contained.

## Definition

For the vector current

\[
j_\mu=\bar s\gamma_\mu c,
\]

we define

\[
\langle 0|j_\mu|D_s^\ast(q,\epsilon)\rangle
=m_{D_s^\ast}f_{D_s^\ast}\epsilon_\mu .
\]

The two-point correlator is

\[
\Pi_{\mu\nu}(q)=i\int d^4x\,e^{iq\cdot x}
\langle0|T\{j_\mu(x)j_\nu^\dagger(0)\}|0\rangle .
\]

We use the transverse coefficient multiplying \(-g_{\mu\nu}\), following the
standard vector-current convention used in heavy-light SVZ analyses.

## LO Spectral Density

At leading order, neglecting \(m_s\) inside the loop but keeping the strange
channel in the threshold choice, the transverse spectral density is

\[
\rho_T^{\rm LO}(s)
=\frac{1}{8\pi^2}s(1-z)^2(2+z),
\qquad
z=\frac{m_c^2}{s}.
\]

This is the formula quoted in the vector-current setup of Gelhausen,
Khodjamirian, Pivovarov and Rosenthal, arXiv:1305.5432.

The Borel sum rule is

\[
f_{D_s^\ast}^2m_{D_s^\ast}^2e^{-m_{D_s^\ast}^2/M^2}
=\int_{m_c^2}^{s_0}ds\,e^{-s/M^2}\rho_T^{\rm LO}(s).
\]

Equivalently,

\[
f_{D_s^\ast}(M^2,s_0)=
\frac{e^{m_{D_s^\ast}^2/(2M^2)}}{m_{D_s^\ast}}
\left[
\int_{m_c^2}^{s_0}ds\,e^{-s/M^2}\rho_T^{\rm LO}(s)
\right]^{1/2}.
\]

## Window And Result

The script

```bash
python3 Dsstar_gamma/scripts/twopoint_vector_decay_constant.py
```

scans

\[
2.5\le M^2\le4.0~{\rm GeV}^2,
\qquad
s_0=6.5,\ 7.0,\ 7.5~{\rm GeV}^2 .
\]

The output is

```text
Dsstar_gamma/outputs/twopoint_f_Dsstar_vector_LO.csv
Dsstar_gamma/outputs/twopoint_f_Dsstar_vector_LO_summary.txt
```

The scanned envelope is

\[
f_{D_s^\ast}=0.212\ldots0.245~{\rm GeV}.
\]

At the representative point

\[
M^2=3.25~{\rm GeV}^2,\qquad s_0=7.0~{\rm GeV}^2,
\]

we obtain

\[
f_{D_s^\ast}=0.227~{\rm GeV}.
\]

For the radiative calculation we assign

\[
f_{D_s^\ast}=0.227\pm0.028~{\rm GeV},
\]

where the uncertainty combines the window envelope with a conservative
normalization allowance for the missing radiative and higher-dimensional OPE
terms.

## Caveat

This is a leading-order local two-point baseline, not a precision NNLO
determination.  It is sufficient for removing the previous assumed
normalization from the first \(D_s^\ast\gamma\) LCSR pass.  Before a final paper
number, the uncertainty scan should include this \(f_{D_s^\ast}\) uncertainty,
and we should either keep this LO baseline as an internally consistent input or
quote a higher-order literature/lattice value as an external comparison.
