---
title: "Recalculated electromagnetic multipole moments"
subtitle: "Independent derivations, numerical validation and revised analysis"
date: "3 October 2026"
geometry: margin=0.78in
fontsize: 10pt
colorlinks: true
header-includes:
  - \usepackage{booktabs,longtable,amsmath,amssymb}
  - \usepackage{float,needspace,etoolbox}
  - \floatplacement{figure}{H}
  - \pretocmd{\section}{\Needspace{14\baselineskip}}{}{}
  - \pretocmd{\subsubsection}{\Needspace{12\baselineskip}}{}{}
  - \pretocmd{\subsection}{\Needspace{8\baselineskip}}{}{}
---

<!-- BEGIN LATTICE RESULTS -->
# Lattice-normalized results, both twist-2 signs

Updated 3 October 2026 from the blind output frozen in Git commit `54e6958`.
These tables give the lattice-normalized run separately for both signs; the
remaining report records the earlier chi=3.15 baseline and its diagnostics.
The lattice tables supersede that baseline for the requested normalization.
No counterpart results have been used, and an entry-by-entry comparison is pending.

The supplied strange-photon normalization at 1 GeV is
$f_{\perp\gamma,s}=-55.1$ MeV. With
$\langle\bar ss\rangle=0.8(-0.240\,\mathrm{GeV})^3$,
$\chi=4.982277199\,\mathrm{GeV}^{-2}$.
The $+\chi$ label denotes this normalization; $-\chi$ denotes
$f_{\perp\gamma,s}=+55.1$ MeV, obtained by reversing only twist 2.
Chi enters only twist 2 and does not enter the two-point sum rule.
Both implementations were exactly equivalent on every window; all other
sectors and every two-point piece remained identical.

Other inputs are unchanged: zeta1=0.4, zeta2=31/110, kappa=0.2,
kappaPlus=0, f3Gamma=-0.004 GeV$^2$, omegaV=3.8, omegaA=-2.1,
m0$^2$=0.8 GeV$^2$, with the central masses, mixing angles and DA scale of 1 GeV.
The nominal counts are 37, 25, 57, 17; the refitted 90-degree counts are
55, 79, 64, 30. At 90 degrees the lighter state uses $J_A$ and the heavier $J_B$.

$F_M$ comes from $q_\mu\epsilon_\nu$;
$F_D=Q_H-F_M+G_3$ uses the clean projection.
$\mathcal D/e=F_D(\hbar c)^2/m_H^2$.
$F_C$ from $(\epsilon\cdot p)g_{\mu\nu}$ is retained as a check in the window CSV.
Its largest deviation from $Q_H$ is $1.07\times10^{-13}$.

All totals below are medians. Brackets give the 16th and 84th percentiles
of the accepted-window distribution, not statistical errors.
$F_D$ includes $Q_H$; $\mathcal D/e$ is in fm$^2$.

## Nominal mixing angles

| State | Sign | Windows | $\mu/\mu_N$ | $F_D$ | $\mathcal D/e$ [fm$^2$] |
| --- | --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | 37 | 0.032610 | -2.546351 | -0.016391 |
| $D_{s1}(2460)^+$ | $-\chi$ | 37 | 1.169591 | 2.740673 | 0.017642 |
| $D_{s1}(2536)^+$ | $+\chi$ | 25 | 0.298056 | 4.822703 | 0.029219 |
| $D_{s1}(2536)^+$ | $-\chi$ | 25 | 0.193698 | -5.810113 | -0.035202 |
| $B_{s1}(5750)$ | $+\chi$ | 57 | -0.180615 | -1.695005 | -0.001996 |
| $B_{s1}(5750)$ | $-\chi$ | 57 | 0.448094 | 1.169889 | 0.001378 |
| $B_{s1}(5830)$ | $+\chi$ | 17 | 0.379925 | 12.144814 | 0.013917 |
| $B_{s1}(5830)$ | $-\chi$ | 17 | 0.332682 | -24.499228 | -0.028075 |

### $\mu/\mu_N$ window intervals

| State | Sign | Median | 16-84% |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | 0.032610 | [-0.109353, 0.109313] |
| $D_{s1}(2460)^+$ | $-\chi$ | 1.169591 | [1.085413, 1.339248] |
| $D_{s1}(2536)^+$ | $+\chi$ | 0.298056 | [0.272922, 0.325029] |
| $D_{s1}(2536)^+$ | $-\chi$ | 0.193698 | [0.191943, 0.195514] |
| $B_{s1}(5750)$ | $+\chi$ | -0.180615 | [-0.238791, -0.151826] |
| $B_{s1}(5750)$ | $-\chi$ | 0.448094 | [0.397841, 0.537432] |
| $B_{s1}(5830)$ | $+\chi$ | 0.379925 | [0.377483, 0.382594] |
| $B_{s1}(5830)$ | $-\chi$ | 0.332682 | [0.316484, 0.354840] |

### $F_D$ window intervals

| State | Sign | Median | 16-84% |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | -2.546351 | [-3.462775, -2.094731] |
| $D_{s1}(2460)^+$ | $-\chi$ | 2.740673 | [2.136436, 3.892965] |
| $D_{s1}(2536)^+$ | $+\chi$ | 4.822703 | [4.190767, 5.302773] |
| $D_{s1}(2536)^+$ | $-\chi$ | -5.810113 | [-6.354322, -5.171887] |
| $B_{s1}(5750)$ | $+\chi$ | -1.695005 | [-2.182237, -1.457712] |
| $B_{s1}(5750)$ | $-\chi$ | 1.169889 | [0.971854, 1.526393] |
| $B_{s1}(5830)$ | $+\chi$ | 12.144814 | [10.834414, 13.434681] |
| $B_{s1}(5830)$ | $-\chi$ | -24.499228 | [-27.540965, -21.933842] |

### $\mathcal D/e$ [fm$^2$] window intervals

| State | Sign | Median | 16-84% |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | -0.016391 | [-0.022290, -0.013484] |
| $D_{s1}(2460)^+$ | $-\chi$ | 0.017642 | [0.013752, 0.025059] |
| $D_{s1}(2536)^+$ | $+\chi$ | 0.029219 | [0.025391, 0.032128] |
| $D_{s1}(2536)^+$ | $-\chi$ | -0.035202 | [-0.038499, -0.031335] |
| $B_{s1}(5750)$ | $+\chi$ | -0.001996 | [-0.002570, -0.001717] |
| $B_{s1}(5750)$ | $-\chi$ | 0.001378 | [0.001144, 0.001798] |
| $B_{s1}(5830)$ | $+\chi$ | 0.013917 | [0.012416, 0.015396] |
| $B_{s1}(5830)$ | $-\chi$ | -0.028075 | [-0.031561, -0.025135] |

### $\mu/\mu_N$ contribution medians

| State | Sign | Hard | Twist 2 | Twist 3 | Twist 4 | Contact |
| --- | --- | --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | 0.780010 | -0.568491 | -0.094003 | -0.069647 | -0.015260 |
| $D_{s1}(2460)^+$ | $-\chi$ | 0.780010 | 0.568491 | -0.094003 | -0.069647 | -0.015260 |
| $D_{s1}(2536)^+$ | $+\chi$ | 0.638959 | 0.051564 | -0.337387 | -0.156162 | 0.097233 |
| $D_{s1}(2536)^+$ | $-\chi$ | 0.638959 | -0.051564 | -0.337387 | -0.156162 | 0.097233 |
| $B_{s1}(5750)$ | $+\chi$ | 0.133431 | -0.314355 | -0.017663 | 0.017337 | 0.000000 |
| $B_{s1}(5750)$ | $-\chi$ | 0.133431 | 0.314355 | -0.017663 | 0.017337 | 0.000000 |
| $B_{s1}(5830)$ | $+\chi$ | 0.480480 | 0.023798 | -0.329725 | 0.206843 | 0.000000 |
| $B_{s1}(5830)$ | $-\chi$ | 0.480480 | -0.023798 | -0.329725 | 0.206843 | 0.000000 |

### $F_D-Q_H$ contribution medians

| State | Sign | Hard | Twist 2 | Twist 3 | Twist 4 | Contact |
| --- | --- | --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | -2.053195 | -2.643512 | 0.590327 | -0.154038 | 0.701903 |
| $D_{s1}(2460)^+$ | $-\chi$ | -2.053195 | 2.643512 | 0.590327 | -0.154038 | 0.701903 |
| $D_{s1}(2536)^+$ | $+\chi$ | -2.216542 | 5.355061 | 1.791186 | -0.072485 | -1.020116 |
| $D_{s1}(2536)^+$ | $-\chi$ | -2.216542 | -5.355061 | 1.791186 | -0.072485 | -1.020116 |
| $B_{s1}(5750)$ | $+\chi$ | -0.147135 | -1.432447 | -0.387724 | 0.266620 | 0.000000 |
| $B_{s1}(5750)$ | $-\chi$ | -0.147135 | 1.432447 | -0.387724 | 0.266620 | 0.000000 |
| $B_{s1}(5830)$ | $+\chi$ | -5.728128 | 18.334977 | 8.001665 | -8.241479 | 0.000000 |
| $B_{s1}(5830)$ | $-\chi$ | -5.728128 | -18.334977 | 8.001665 | -8.241479 | 0.000000 |

## 90 degrees, refitted pure-current windows

| State | Sign | Windows | $\mu/\mu_N$ | $F_D$ | $\mathcal D/e$ [fm$^2$] |
| --- | --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | 55 | 0.332018 | 0.660839 | 0.004254 |
| $D_{s1}(2460)^+$ | $-\chi$ | 55 | 1.284865 | -1.822083 | -0.011729 |
| $D_{s1}(2536)^+$ | $+\chi$ | 79 | -0.043564 | -2.969843 | -0.017993 |
| $D_{s1}(2536)^+$ | $-\chi$ | 79 | 0.901138 | 3.934440 | 0.023838 |
| $B_{s1}(5750)$ | $+\chi$ | 64 | -0.132233 | 1.817537 | 0.002140 |
| $B_{s1}(5750)$ | $-\chi$ | 64 | 0.457865 | -1.814992 | -0.002137 |
| $B_{s1}(5830)$ | $+\chi$ | 30 | -0.155358 | -3.527316 | -0.004042 |
| $B_{s1}(5830)$ | $-\chi$ | 30 | 0.366107 | 2.366562 | 0.002712 |

### $\mu/\mu_N$ window intervals

| State | Sign | Median | 16-84% |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | 0.332018 | [0.231408, 0.380994] |
| $D_{s1}(2460)^+$ | $-\chi$ | 1.284865 | [1.189955, 1.465082] |
| $D_{s1}(2536)^+$ | $+\chi$ | -0.043564 | [-0.332015, 0.068829] |
| $D_{s1}(2536)^+$ | $-\chi$ | 0.901138 | [0.830043, 1.107031] |
| $B_{s1}(5750)$ | $+\chi$ | -0.132233 | [-0.187257, -0.104411] |
| $B_{s1}(5750)$ | $-\chi$ | 0.457865 | [0.410788, 0.549467] |
| $B_{s1}(5830)$ | $+\chi$ | -0.155358 | [-0.179274, -0.140342] |
| $B_{s1}(5830)$ | $-\chi$ | 0.366107 | [0.341507, 0.402776] |

### $F_D$ window intervals

| State | Sign | Median | 16-84% |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | 0.660839 | [0.492624, 0.725970] |
| $D_{s1}(2460)^+$ | $-\chi$ | -1.822083 | [-2.719671, -1.380208] |
| $D_{s1}(2536)^+$ | $+\chi$ | -2.969843 | [-4.729712, -2.340610] |
| $D_{s1}(2536)^+$ | $-\chi$ | 3.934440 | [2.614303, 7.670882] |
| $B_{s1}(5750)$ | $+\chi$ | 1.817537 | [1.687457, 1.947150] |
| $B_{s1}(5750)$ | $-\chi$ | -1.814992 | [-2.570113, -1.470932] |
| $B_{s1}(5830)$ | $+\chi$ | -3.527316 | [-3.921854, -3.261090] |
| $B_{s1}(5830)$ | $-\chi$ | 2.366562 | [2.065392, 2.859157] |

### $\mathcal D/e$ [fm$^2$] window intervals

| State | Sign | Median | 16-84% |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | 0.004254 | [0.003171, 0.004673] |
| $D_{s1}(2460)^+$ | $-\chi$ | -0.011729 | [-0.017506, -0.008884] |
| $D_{s1}(2536)^+$ | $+\chi$ | -0.017993 | [-0.028656, -0.014181] |
| $D_{s1}(2536)^+$ | $-\chi$ | 0.023838 | [0.015839, 0.046476] |
| $B_{s1}(5750)$ | $+\chi$ | 0.002140 | [0.001987, 0.002293] |
| $B_{s1}(5750)$ | $-\chi$ | -0.002137 | [-0.003027, -0.001732] |
| $B_{s1}(5830)$ | $+\chi$ | -0.004042 | [-0.004494, -0.003737] |
| $B_{s1}(5830)$ | $-\chi$ | 0.002712 | [0.002367, 0.003276] |

### $\mu/\mu_N$ contribution medians

| State | Sign | Hard | Twist 2 | Twist 3 | Twist 4 | Contact |
| --- | --- | --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | 0.839351 | -0.476423 | 0.026448 | -0.057200 | -0.002975 |
| $D_{s1}(2460)^+$ | $-\chi$ | 0.839351 | 0.476423 | 0.026448 | -0.057200 | -0.002975 |
| $D_{s1}(2536)^+$ | $+\chi$ | 0.707351 | -0.470276 | -0.188986 | -0.089902 | 0.000371 |
| $D_{s1}(2536)^+$ | $-\chi$ | 0.707351 | 0.470276 | -0.188986 | -0.089902 | 0.000371 |
| $B_{s1}(5750)$ | $+\chi$ | 0.160582 | -0.294854 | -0.019698 | 0.022629 | 0.000000 |
| $B_{s1}(5750)$ | $-\chi$ | 0.160582 | 0.294854 | -0.019698 | 0.022629 | 0.000000 |
| $B_{s1}(5830)$ | $+\chi$ | 0.114196 | -0.260607 | -0.020594 | 0.012360 | 0.000000 |
| $B_{s1}(5830)$ | $-\chi$ | 0.114196 | 0.260607 | -0.020594 | 0.012360 | 0.000000 |

### $F_D-Q_H$ contribution medians

| State | Sign | Hard | Twist 2 | Twist 3 | Twist 4 | Contact |
| --- | --- | --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | $+\chi$ | -0.283352 | 1.248852 | -1.329886 | 0.005127 | 0.026241 |
| $D_{s1}(2460)^+$ | $-\chi$ | -0.283352 | -1.248852 | -1.329886 | 0.005127 | 0.026241 |
| $D_{s1}(2536)^+$ | $+\chi$ | -2.984153 | -3.452141 | 1.946491 | -0.256859 | 0.786012 |
| $D_{s1}(2536)^+$ | $-\chi$ | -2.984153 | 3.452141 | 1.946491 | -0.256859 | 0.786012 |
| $B_{s1}(5750)$ | $+\chi$ | 1.032532 | 1.807016 | -0.828522 | -0.207140 | 0.000000 |
| $B_{s1}(5750)$ | $-\chi$ | 1.032532 | -1.807016 | -0.828522 | -0.207140 | 0.000000 |
| $B_{s1}(5830)$ | $+\chi$ | -1.204771 | -2.942746 | 0.277191 | 0.343996 | 0.000000 |
| $B_{s1}(5830)$ | $-\chi$ | -1.204771 | 2.942746 | 0.277191 | 0.343996 | 0.000000 |

## Grouping and validation

Hard combines the perturbative heavy and strange lines. Twist 3 combines
the two-particle twist-3 terms with three-particle V and A. Twist 4 combines
two-particle twist 4 with S, S-tilde, T1-T4 and electromagnetic S-gamma and
T4-gamma. Contact is the Wilson-line term.

On all 728 signed rows the five parts add to the complete amplitude:
maximum residual $3.50\times10^{-15}$ for $\mu/\mu_N$ and
$3.91\times10^{-14}$ for $F_D-Q_H$.
The contribution medians need not add to the median total because median
is not a linear operation. The per-window CSV, rather than the rounded
median tables, is the entry-by-entry comparison source.

The positive-power continuum prescription and the nonpositive-power
unsubtracted approximation are unchanged. Doubling three-point quadrature
from 256 to 512 changes the checked observables by at most
$3.59\times10^{-11}$ in absolute units.

The original signed CSVs and their note remain unchanged and are covered by
[the committed SHA-256 manifest](blind_lattice_twist_split/sha256_manifest.json).
The full [window CSV](blind_lattice_twist_split/twist_split_windows.csv),
[summary CSV](blind_lattice_twist_split/twist_split_summary.csv) and
[input and validation note](blind_lattice_twist_split/twist_split_note.md)
contain the reproducible numerical record.

<!-- END LATTICE RESULTS -->

\clearpage

# Earlier baseline: chi=3.15

The calculation has been rerun with repaired OPE coefficients and two-point
residues, a uniform positive-power continuum prescription, the Wilson-line
response, electromagnetic three-particle operators and pole-free projections.
The central setup also enforces the exact twist-4 parameter constraint. All
tables and figures below come from our own symbolic exports and numerical
engine. Counterpart output files are not calculation inputs or acceptance tests.

These are predictions **within the specified LO and twist/OPE truncation**.
The numerical calculation is complete at that truncation. A complete theoretical
error budget, including missing higher-order contributions, is not claimed.

The follow-up independent scan audit confirms that widening the Borel search
adds no accepted points at the adopted thresholds. Widening the threshold
interval does change results: the interval $\sqrt{s_0}-m_H=0.40$--$0.60$ GeV
is a model choice, not a unique consequence of the acceptance cuts. Independent
evaluation of this setup should not be confused with a historically blind
choice of setup. Details are in `blind_definition_scan_audit/AUDIT.md`.

| State | $F_M(0)$ | $\mu/\mu_N$ | $F_D(0)$ | $\mathcal D/e$ [fm$^2$] |
| --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | 0.632386 | 0.241248 | -1.574176 | -0.010133 |
| $D_{s1}(2536)^+$ | 0.752791 | 0.278616 | 2.827267 | 0.017130 |
| $B_{s1}(5750)$ | -0.401926 | -0.065583 | -1.168263 | -0.001376 |
| $B_{s1}(5830)$ | 2.309223 | 0.371700 | 5.405695 | 0.006195 |

Each entry is the median of pointwise results, with the residue evaluated at
that same Borel point and threshold. The physical conversion is

$$
\frac{\mu_H}{\mu_N}=\frac{m_p}{m_H}F_M(0),\qquad
\frac{\mathcal D_H}{e}=
\frac{(\hbar c)^2}{m_H^2}F_D(0).
$$

Six digits are retained for reproducibility, not as an estimate of physical
precision. The following intervals describe the accepted-window distribution;
they are **not statistical errors or a full theoretical uncertainty**.

| State | $\mu/\mu_N$: 16–84% | $F_D$: 16–84% |
| --- | --- | --- |
| $D_{s1}(2460)^+$ | [0.1570, 0.2890] | [-2.1102, -1.3157] |
| $D_{s1}(2536)^+$ | [0.2580, 0.3011] | [2.4636, 3.1460] |
| $B_{s1}(5750)$ | [-0.0961, -0.0504] | [-1.5003, -1.0106] |
| $B_{s1}(5830)$ | [0.3672, 0.3762] | [4.8090, 5.8842] |

# Independently established corrections

## Correlator phase and photon twists

The saved coefficients use $\Pi=-C$, where $C$ is the physical correlator.
The two-particle block had the opposite storage sign to the hard block. Its
entire phase is repaired. The multiline Mathematica definitions of the twist-3
functions also ended before their $\omega_V,\omega_A$ contributions; the complete
right-hand sides are now parenthesized. An explicit four-dimensional trace fixes
the tensor phase. No ABS2020 comparison is used as a correctness criterion.

For three-particle terms the DA weight must be inside the Borel derivative.
Writing $a_i=1/M_i^2$, $u=a_2/(a_1+a_2)$ and
$W_j(u)=\int D\alpha\,dv\,v^j\Phi(\alpha)
\delta(\alpha_{\bar q}+v\alpha_g-u)$, our transform applies

$$
\left(\frac{\partial_{a_1}-\partial_{a_2}}2\right)^k
\left[e^{-(a_1+a_2)m_Q^2}
\frac{(a_1+a_2)^{n-2}}{\Gamma(n)}c(u)W_j(u)\right].
$$

The previous implementation differentiated only the hard factor. That omitted
the product-rule terms and incorrectly made the magnetic $A,V$ contributions
zero. Inverse light-cone factors are treated with the corresponding integrated
DA weights, including their derivatives. The original raw traces are retained.

## Two-point normalization and link response

The light-propagator odd term must be evaluated at $-x$. It reverses the old
$m_s\langle\bar ss\rangle$ contribution. A heavy-line dimension-five term was
also missing. With $m=m_Q$, $d=m_Q+m_s$, $t=M^2$, $S=\langle\bar ss\rangle$,
$X=m_0^2S$ and $E=e^{-m^2/t}$, the repaired terms are

$$
B_{m_sS}=E m_sS
\begin{pmatrix}
 m^2/(2t)&m(m^2/t-1)/(2d)\\
 m(m^2/t-1)/(2d)&m^2(m^2/(2t)-1)/d^2
\end{pmatrix},
$$
$$
B_{5,Q}=-\frac{EX}{12}(m^2/t-1)
\begin{pmatrix}0&1/d\\1/d&2m/d^2\end{pmatrix}.
$$

Both follow from explicit Dirac matrices, independently of the exchanged
two-point values. The pseudoscalar zero-momentum and mixed-condensate identities
provide additional checks of this derivation. Exact bottom-state masses are
used throughout.

The electromagnetic links on the two lines give their net current charge
$Q_H=e_Q-e_s$. Expanding the straight link once yields a uniform integral along
the line. For the **physical** two-point condensate tensor $C_0^{\rm cond}$, our
stored contact block is

$$
\Pi^W_{\mu\nu}=-(e_Q-e_s)\epsilon^\rho\int_0^1dt\,
\frac{\partial C^{\rm cond}_{0,\mu\nu}(P)}{\partial P^\rho}
\bigg|_{P=p+tq}.
$$

The minus sign here is the conversion to $\Pi=-C$. The current-vertex momenta
remain $p$ and $p+q$; only the internal heavy momentum is differentiated. All
five local pieces of the repaired condensate OPE enter. We projected the full
nine-tensor response, including $p_\mu\epsilon_\nu$, before retaining the eight
reported structures. All 16 tensor components reconstruct exactly in every
current channel. This prevents an omitted tensor from contaminating a magnetic
coefficient.

Replacing $\epsilon$ by $q$ turns the line integral into
$C_0^{\rm cond}(p+q)-C_0^{\rm cond}(p)$. This independently establishes the Ward
response. The dependence on $p+tq$ survives the double Borel transform.
Single-leg endpoint terms are a different class of terms.

The electromagnetic $S_\gamma,T_4^\gamma$ terms are obtained by inserting the
electromagnetic field on the heavy line in the same operator traces used for
$S,T_4$. Their charge is $e_Q$. Their weights are integrated and differentiated
independently; no numerical coefficient is copied from the other calculation.

## DA equations of motion and an additional input inconsistency

Our symbolic integrations give the $u^2(1-u)^2$ leading factor in $\mathcal A$
and the $1+2\kappa^+$ factor in $h_\gamma$. The algebraic first line of the
printed $T_4^\gamma$ expression simplifies to the sign multiplying
$(\alpha_{\bar q}-\alpha_q)$; the opposite second-line sign is inconsistent
with that first line. The DA definitions and operator identities are in
[Ball, Braun and Kivel](https://arxiv.org/pdf/hep-ph/0207307), Eqs. (4.32),
(4.39), (4.41), (4.43) and (4.45).

There is a further problem with using all rounded table values as exact inputs:

$$
\zeta_1+11\zeta_2-2\zeta_2^+=\frac72,
\qquad 0.4+11(0.3)-0=3.7\ne3.5.
$$

Our independently integrated first equation of motion leaves

$$
h_\gamma^{\rm derived}(u)-h_\gamma(u)
=-40u(u-1)[1+5u(u-1)]
\left(\zeta_1+11\zeta_2-2\zeta_2^+-\frac72\right).
$$

It equals $-0.5$ at $u=1/2$ for the rounded inputs. Enforcing the constraint
makes all three checked DA identities vanish exactly. Our central convention
holds $\zeta_1=0.4$ and $\zeta_2^+=0$, and derives
$\zeta_2=31/110$. This choice of dependent parameter is explicit, not unique;
the alternative $\zeta_1=0.2,\zeta_2=0.3$ is also calculated below.

Using exact local EM distributions while retaining the truncated gluonic DA
set is not a complete conformal-order variation. Such a replacement is not
assigned a statistical uncertainty here. A consistent higher-order extension
must respect the coupled DA identities.

## Continuum prescription and its limits

For every completed Borel coefficient $c_j t^j e^{-S/t}$ with $j\ge1$ we use

$$
t^j e^{-S/t}\longmapsto\frac1{\Gamma(j)}
\int_S^{s_0}(s-S)^{j-1}e^{-s/t}\,ds
=t^j e^{-S/t}P\left(j,\frac{s_0-S}{t}\right).
$$

The map is applied after all Borel derivatives. For hard terms,
$S=m_Q^2/r+m_s^2/(1-r)$ and the factor remains inside the $r$ integral;
for DA terms $S=m_Q^2$. Nonpositive DA powers are left unsubtracted.
Two-point perturbative terms use their ordinary dispersion integral with the
same threshold; local condensates are retained at their pole support.

The positive-power rule follows directly from its Laplace spectral integral.
Its extension to a complete double-dispersion continuum model remains a
model assumption. In particular, [Belyaev et al., Appendix B](https://arxiv.org/pdf/hep-ph/9410280)
explicitly treat the omission of higher-twist subtraction as an approximation;
zero and negative Borel powers are not proved to be purely diagonal spectral
densities there. We do not claim that the paper establishes an exact universal
rule for every transverse amplitude.

Charge closure is a necessary consistency test, not a uniqueness proof for
the magnetic and quadrupole continuum treatment. A variation that retains the
old DA subtraction but uses the repaired hard map and link term still closes
the charge, while changing the transverse form factors. Its numerical effect
is tabulated below.

# Susceptibility sign and electric charge

## Which susceptibility sign belongs to our conventions?

Our adopted input is **$\chi(1\,\mathrm{GeV})=+3.15\,\mathrm{GeV}^{-2}$**.
This sign belongs to the complete definition

$$
\sigma_{\alpha\beta}=\frac{i}{2}[\gamma_\alpha,\gamma_\beta],\qquad
D_\mu=\partial_\mu-iee_f A_\mu,\qquad
\langle\bar s\sigma_{\alpha\beta}s\rangle_F
=e_s\chi S F_{\alpha\beta},
$$

where the positive elementary charge $e$ is factored out of the last equation
and $S=\langle\bar ss\rangle=-0.0110592\,\mathrm{GeV}^3$.
Thus $\chi S=-0.03483648\,\mathrm{GeV}$ and
$e_s\chi S=+0.01161216\,\mathrm{GeV}$. The sign of the condensate is retained.
These are the conventions of [Ball, Braun and Kivel, Eqs. (2.1), (3.13) and
footnote 7](https://arxiv.org/pdf/hep-ph/0207307).

Negative susceptibilities in another convention need not be contradictory.
For example, defining the same response as $-e_s\chi_{\rm alt}SF$ requires
$\chi_{\rm alt}=-\chi$. [Bali et al., discussion after Eq. (10)](https://arxiv.org/pdf/1209.6015)
explicitly identify the opposite electromagnetic coupling sign relative to BBK.
One must convert the defining operators and coupling together; changing the
number $\chi$ alone in our fixed formulas changes the prediction.

For an outgoing wave $A_\mu=\epsilon_\mu e^{iqx}$,
$F_{\alpha\beta}=-i(\epsilon_\alpha q_\beta-\epsilon_\beta q_\alpha)e^{iqx}$.
This is also the phase used in [Rohrwild, Eqs. (8), (37)](https://arxiv.org/pdf/0708.1405).
The correlator's outer $i$ and heavy propagator give
$iS_Q(k)=(\not k+m_Q)/(m_Q^2-k^2)$. With the tensor Fierz factor $1/8$,
the tensor-DA prefactors are

$$
C:\ -\frac{i}{8}e_s\chi S,\qquad
\Pi=-C:\ +\frac{i}{8}e_s\chi S.
$$

Our independent pure-axial trace is

$$
\operatorname{Tr}[\gamma_\mu\gamma_5(\not k+m_Q)
\gamma_\nu\gamma_5\sigma_{\alpha\beta}]
(\epsilon^\alpha q^\beta-\epsilon^\beta q^\alpha)
=-8im_Q(q_\mu\epsilon_\nu-\epsilon_\mu q_\nu).
$$

Consequently the stored, continuum-subtracted twist-2 coefficient is

$$
\widehat\Pi^{\mathrm{tw2},AA}_{q\epsilon}
=e_s\chi S m_Q\varphi_\gamma(1/2)M^2
\left(e^{-m_Q^2/M^2}-e^{-s_0/M^2}\right)>0.
$$

At $m_Q=1.27$, $M^2=2.05$, $s_0=7.825$ in GeV units and
$\varphi_\gamma(1/2)=1.5$, it is $+0.0196501926792$.
Its contribution to $F_M=-\widehat\Pi_{q\epsilon}/B_H$ is therefore negative
when $B_H>0$. This example is AA; mixing can change a sector's rotated sign.
The trace was checked in Dirac and chiral representations at 32 kinematics
each, with a maximum residual of $1.78\times10^{-14}$.

The earlier error was using the $C$ prefactor in a block stored as $\Pi=-C$.
It was **not a reason to adopt negative $\chi$**. Flipping $\chi$ alone would
only reverse twist 2, leaving the erroneous common phase in twists 3 and 4.
The magnitude of $\chi$ remains a nonperturbative input, not something derived
from the trace or charge normalization.

## Why does $F_C(0)$ equal 1 or 0?

The conserved charge operator gives
$F_C(0)=Q_H=e_Q-e_s$. Our form factor includes the dimensionless hadron charge;
only the elementary $e$ is factored out. Therefore the four values are
$(1,1,0,0)$, at any mixing angle. A neutral meson can still have nonzero magnetic
and quadrupole moments. The integer condition is specific to zero momentum
transfer, not the complete function $F_C(q^2)$.

The charge extraction in our stored convention is

$$
F_C(0)=\frac{\widehat\Pi_g}{2B_H},\qquad
B_H=f_H^2m_H^2e^{-m_H^2/M^2},\qquad
\frac1{M^2}=\frac1{M_1^2}+\frac1{M_2^2}.
$$

The factor 2 follows from $(p+p')\cdot\epsilon=2p\cdot\epsilon$.
The two-point denominator and three-point numerator must use the same inputs,
threshold, effective Borel parameter and rotation. At the retained OPE order,
the separately evaluated current-basis matrices satisfy

$$
\widehat\Pi_g^{Q,\mathrm{hard}}=2e_Q B_{\rm pert},\qquad
\widehat\Pi_g^{s,\mathrm{hard}}=-2e_s B_{\rm pert},\qquad
\widehat\Pi_g^{W}=2Q_H B_{\rm cond}.
$$

Every retained photon-DA sector has zero $g$ coefficient. Hence
$\widehat\Pi_g=2Q_H(B_{\rm pert}+B_{\rm cond})=2Q_HB_H$ before rotation,
and $F_C(0)=Q_H$ after any common rotation. The numerator is computed from
our kernels; it is not replaced by an integer or rescaled to fit one.

For illustration, the historical benchmark points give:

| State | $M^2$ | $s_0$ | Hard / $2B_H$ | Link / $2B_H$ | $F_C(0)$ |
|---|---:|---:|---:|---:|---:|
| $D_{s1}(2460)^+$ | 2.05 | 7.825 | 1.229828 | -0.229828 | 1 |
| $D_{s1}(2536)^+$ | 2.05 | 10.500 | 1.080422 | -0.080422 | 1 |
| $B_{s1}(5750)$ | 6.50 | 44.250 | 0 | 0 | 0 |
| $B_{s1}(5830)$ | 6.50 | 44.275 | 0 | 0 | 0 |

Here $M^2,s_0$ are in GeV$^2$; rows are individual checks, not accepted-window
medians. The maximum unrounded charge residual is $3.80\times10^{-14}$.
The neutral hard entry is a cancellation between nonzero heavy- and
strange-line terms. Omitting the link response, changing only the numerator's
continuum prescription, or using an unrelated residue can spoil the identity.

An externally determined fixed decay constant need not give this exact
cancellation at every Borel point. The exact ratio tests consistency between
our two sum rules; it does not independently establish the physical accuracy
of the decay constants. Likewise, the Ward relation does not determine the
continuum model in every transverse structure.

Explicit runs with positive, zero and negative $\chi$ leave $F_C(0)$ unchanged
while changing magnetic results. **Integer charge does not test the sign of
$\chi$, the correctness of all transverse amplitudes, or uniqueness of the
chosen continuum/window model.** The detailed convention derivation and
machine-readable checks are in `susceptibility_charge_audit/EXPLANATION.md`
and `susceptibility_charge_audit/checks.json`.


# Pole separation and observables

The pseudoscalar contribution at the right current carries $p'_\nu=p_\nu+q_\nu$.
Its most general relevant tensor is

$$
(p_\nu+q_\nu)
[a\epsilon_\mu+b(\epsilon\!\cdot p)p_\mu
 +c(\epsilon\!\cdot p)q_\mu].
$$

It has no $q_\mu\epsilon_\nu$ coefficient, and contributes equally to the
$qq,qp$ coefficients and to $\epsilon q,\epsilon p$. A pseudoscalar at the left
current carries $p_\mu$ and cannot contribute to these clean combinations.
This is an algebraic pole argument; no numerical route agreement is needed
to infer it.

Let $B_H=f_H^2m_H^2e^{-m_H^2/M^2}$ be the repaired rotated two-point sum.
The correlator-oriented spin-one vertex then gives

$$
F_C=\frac{\Pi_g}{2B_H},\qquad
F_M=-\frac{\Pi_{q\epsilon}}{B_H}
=\frac{\Pi_{\epsilon q}-\Pi_{\epsilon p}}{B_H},
$$
$$
G_3=-\frac{m_H^2(\Pi_{qq}-\Pi_{qp})}{B_H},\qquad
F_D=Q_H-F_M+G_3.
$$

The equality of the two magnetic expressions is also tested on the completed
OPE. The older $\epsilon q$ and charge-subtracted $qq$ routes are retained as
diagnostics, but they are not averaged with the pole-free results.

# Inputs, windows and controlled changes

The central inputs are $m_c=1.27$, $m_b=4.18$, $m_s=0.093$ GeV,
$\langle\bar ss\rangle=0.8(-0.240\,{\rm GeV})^3$,
$m_0^2=0.8\,{\rm GeV}^2$, $\chi=3.15\,{\rm GeV}^{-2}$,
$f_{3\gamma}=-0.004\,{\rm GeV}^2$, $\omega_V=3.8$, $\omega_A=-2.1$,
$\kappa=0.2$, $\kappa^+=0$, $\zeta_1=0.4$, $\zeta_2=31/110$,
$\zeta_1^+=\zeta_2^+=\varphi_2=0$. The DA scale is 1 GeV.
The hadron masses are 2.4595, 2.53511, 5.7502 and 5.82911 GeV.
Angles are $26.6^\circ$ for charm and $38.5^\circ$ for bottom.
Conversion constants are $m_p=0.93827208816$ GeV and
$\hbar c=0.1973269804$ GeV fm. Charges are $e_c=2/3$, $e_b=e_s=-1/3$.

The scan uses $\sqrt{s_0}=m_H+\Delta$ with
$\Delta=0.40,0.45,0.50,0.55,0.60$ GeV. The $M^2$ grid is
1.2–5.0 GeV$^2$ in 0.1 steps for charm and 3.0–12.0 GeV$^2$ in 0.25
steps for bottom. Acceptance is decided before computing multipole moments:
$B_H>0$, pole fraction $B_H(s_0)/B_H(\infty)\ge0.5$,
$|B_5|/B_H\le0.1$, and $|m_{\rm SR}/m_H-1|\le0.03$.
Here $m_{\rm SR}^2=-\partial\ln B_H/\partial(1/M^2)$ is evaluated by an
analytic Borel derivative, checked against finite differences.

| State | Points | $M^2$ [GeV$^2$] | $s_0$ [GeV$^2$] | Median $f_H$ [GeV] |
| --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | 37 | 1.60–2.90 | 8.1767–8.7586 | 0.421056 |
| $D_{s1}(2536)^+$ | 25 | 2.30–3.00 | 8.6149–9.8289 | 0.160475 |
| $B_{s1}(5750)$ | 57 | 4.75–9.50 | 38.4425–40.3250 | 0.494055 |
| $B_{s1}(5830)$ | 17 | 5.25–6.75 | 39.4272–41.3335 | 0.083781 |

The row-by-row acceptance decisions, including all rejected candidates, are
recorded in `two_point_window_scan.csv`. These criteria are transparent model
choices, not experimental facts; the computed window spread is therefore kept
distinct from probabilistic input propagation.

## Three-point stability is a separate test

The two-point acceptance conditions do not establish three-point OPE convergence.
The clean magnetic results for $D_{s1}(2460)$ and $B_{s1}(5750)$ retain substantial
Borel dependence over the selected points. The table below gives the median
magnitude of the full retained twist-4 contribution divided by the final form
factor. The twist-4 group contains two-particle twist 4, the gluonic scalar and
tensor operators, and the electromagnetic three-particle terms.

| State | $|F_M^{\rm tw4}/F_M|$ | $|F_D^{\rm tw4}/F_D|$ |
| --- | --- | --- |
| $D_{s1}(2460)^+$ | 0.288 | 0.096 |
| $D_{s1}(2536)^+$ | 0.564 | 0.025 |
| $B_{s1}(5750)$ | 0.264 | 0.226 |
| $B_{s1}(5830)$ | 0.552 | 1.539 |

For $B_{s1}(5830)$ the twist-4 quadrupole contribution exceeds the final result
because of cancellations. Its numerical precision must not be mistaken for a
robust physical uncertainty. The individual groups, ratios to the sum of
absolute group magnitudes, and full ranges are exported in
`three_point_sector_hierarchy.csv` and `three_point_convergence_summary.csv`.
No additional acceptance cut was tuned to force agreement with another result
or to hide these points. More restrictive three-point stability criteria would
define a different analysis and may leave fewer or no accepted points.

## Controlled updates on the original 62 points

The following tables hold the original points fixed and use the repaired
residue point by point. They show the old extraction routes to isolate changes
in the OPE and input setup. The first row already includes the independently
established phase, omega, product-rule and two-point corrections.


### Old $\epsilon q$ magnetic route, $\mu/\mu_N$

| Step | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| Repaired blocks and residue | 0.90427 | 1.14935 | -0.00079 | 0.42800 |
| Uniform positive-power map | 0.48211 | 0.85997 | -0.01430 | 0.36421 |
| Add link response | 0.48146 | 0.86517 | -0.01430 | 0.36421 |
| Add EM three-particle | 0.42421 | 0.78615 | -0.00669 | 0.37463 |
| Correct DA shapes | 0.44026 | 0.59584 | -0.03293 | 0.27125 |
| Change kappa inputs | 0.44395 | 0.60074 | -0.03138 | 0.27295 |
| Enforce zeta constraint | 0.44409 | 0.60035 | -0.03143 | 0.27268 |

### Old $qq$ route, $F_D$

| Step | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| Repaired blocks and residue | -1.07478 | -2.18324 | -1.74856 | 5.86562 |
| Uniform positive-power map | -1.16783 | -1.82697 | -1.68072 | 5.60574 |
| Add link response | -1.28768 | -1.10317 | -1.68072 | 5.60574 |
| Add EM three-particle | -1.59478 | -1.33149 | -1.64566 | 6.12770 |
| Correct DA shapes | -1.95229 | -0.18160 | -1.90932 | 8.61292 |
| Change kappa inputs | -1.92334 | -0.28417 | -1.89367 | 8.41106 |
| Enforce zeta constraint | -1.92304 | -0.28549 | -1.89316 | 8.40708 |

The full file `old62_controlled_updates.csv` also gives the clean routes at
every step. Switching to clean projections and then refitting the windows are
separate operations; a median shift must not be applied to another median.

## Rounded-input comparison and consistent-input sensitivity

For reference, our own calculation with the rounded $\zeta_2=0.3$ gives:

| State | $\mu/\mu_N$ | $F_D$ |
| --- | --- | --- |
| $D_{s1}(2460)^+$ | 0.241807 | -1.577261 |
| $D_{s1}(2536)^+$ | 0.280680 | 2.827215 |
| $B_{s1}(5750)$ | -0.065504 | -1.169136 |
| $B_{s1}(5830)$ | 0.373855 | 5.407271 |

These reproduce the initially relayed rounded central results. That agreement
is not used to validate the calculation. Our central table instead satisfies
the exact DA constraint.

The next variations keep the accepted nominal points fixed and recompute all
affected OPE coefficients and residues. They are model/input sensitivities,
not independent one-sigma errors.


### Clean $\mu/\mu_N$

| Variation | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| Central | 0.24125 | 0.27862 | -0.06558 | 0.37170 |
| Constrain zeta1 instead | 0.24396 | 0.27802 | -0.06537 | 0.37312 |
| Old kappa pair | 0.24330 | 0.25895 | -0.06671 | 0.36080 |
| Old DA subtraction | 0.17568 | 0.24496 | -0.07401 | 0.35805 |

### Clean $F_D$

| Variation | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| Central | -1.57418 | 2.82727 | -1.16826 | 5.40570 |
| Constrain zeta1 instead | -1.58343 | 2.84480 | -1.17246 | 5.46478 |
| Old kappa pair | -1.65272 | 3.12831 | -1.18804 | 6.24863 |
| Old DA subtraction | -1.37072 | 3.01096 | -1.11593 | 5.64881 |

The machine-readable sensitivity table also evaluates $\pm10\%$ changes of
$m_s$, $\langle\bar ss\rangle$, $m_0^2$, $\kappa$ and $\zeta_1$ (with the
constraint maintained), and $\pm1\%$ of $m_Q$. These are derivatives/sensitivity
probes, not assumed measurement uncertainties. Higher-order continuum and OPE
uncertainties remain beyond this truncation.

# Mixing and pure-current analysis

At each angle the OPE and the two-point matrix are rotated together. The
nominal-window curves hold the accepted points fixed so that their angle
dependence is unambiguous. Fixed-residue curves are shown only as a diagnostic.
The pure-current $90^\circ$ results use newly selected windows, with the same
acceptance criteria. At this angle state 1 is axial and state 2 is the
normalized tensor current up to an irrelevant sign.

| State | Refitted points | $\mu/\mu_N$ at $90^\circ$ | $F_D$ at $90^\circ$ |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | 55 | 0.506548 | 0.202758 |
| $D_{s1}(2536)^+$ | 79 | 0.129385 | -1.700287 |
| $B_{s1}(5750)$ | 64 | -0.023798 | 1.118169 |
| $B_{s1}(5830)$ | 30 | -0.059648 | -2.445490 |

The nominal-angle endpoints at $\theta\pm0.6^\circ$ (charm) and
$\theta\pm0.1^\circ$ (bottom), plus the fixed-window $0^\circ$ and $90^\circ$
diagnostics, are in `angle_endpoints.csv`. Pure-current results with their
refitted windows should not be confused with the endpoints of a frozen-window
curve.

The appendix now tabulates all documented extraction routes at nominal mixing,
at $90^\circ$ on frozen nominal windows, and at $90^\circ$ with refitted windows.
It includes the two clean projections and the alternative diagnostic routes.
The latter still exhibit structure dependence. In particular, state 2 is BB
at $90^\circ$, and $P^\mu J_{B\mu}=0$ excludes a pseudoscalar overlap there.
Its remaining route differences therefore cannot all be explained by
pseudoscalar contamination. They remain an unresolved limitation of the
truncated OPE, continuum treatment and pole isolation. Keeping a physical mass
label does not prove that a pure current isolates only that named state.

![Updated residue curves.](figures/decay_constants_vs_angle.png){width=95%}

![Updated magnetic mixing-angle analysis.](figures/magnetic_moments_vs_angle.png){width=95%}

![Updated quadrupole mixing-angle analysis.](figures/quadrupole_moments_vs_angle.png){width=95%}

# Partial uncertainty propagation

We ran 100,000 seeded samples using independent Gaussian photon priors
$\chi=3.15\pm0.30$, $f_{3\gamma}=-0.004\pm0.002$,
$\omega_V=3.8\pm1.8$ and $\omega_A=-2.1\pm1.0$, from
[BBK, Eqs. (3.13) and (4.10)](https://arxiv.org/pdf/hep-ph/0207307), and the
project's angle priors $26.6\pm0.6^\circ$, $38.5\pm0.1^\circ$.
Each sample draws uniformly from the accepted discrete points. The residue is
rotated at the sampled angle. No independent multiplicative residue error is
added, which would double count correlated normalization effects.

The resulting intervals combine only these sampled inputs with the chosen
window distribution. They do not include mass/condensate probability priors,
twist-4 parameter errors, perturbative corrections or omitted higher-twist
terms. The former partial uncertainty results used different photon inputs and
normalization assumptions and are superseded.

| State | Sampled $\mu/\mu_N$ [16–84%] | Sampled $F_D$ [16–84%] |
| --- | --- | --- |
| $D_{s1}(2460)^+$ | 0.2385 [0.1281, 0.3228] | -1.7014 [-2.4360, -0.9426] |
| $D_{s1}(2536)^+$ | 0.3121 [0.0603, 0.4981] | 2.5837 [1.3851, 4.2439] |
| $B_{s1}(5750)$ | -0.0665 [-0.1033, -0.0412] | -1.1651 [-1.5735, -0.9155] |
| $B_{s1}(5830)$ | 0.4045 [0.1484, 0.5943] | 4.4235 [-0.4352, 11.1821] |

The sampling medians need not equal the deterministic central values because
the input distribution and ratios are nonlinear. Results for window-only,
photon-plus-window and angle-plus-window sampling are recorded separately.

# Independent numerical validation

| Check | Maximum residual |
| --- | --- |
| Direct spectral integral versus positive-power map | 9.237e-14 |
| Link charge coefficient versus repaired condensate matrix | 6.939e-18 |
| Charge at randomized masses, thresholds and angles | 4.312e-13 |
| Two clean magnetic projections | 1.021e-14 |
| Two-point Borel derivative, relative | 5.438e-10 |
| All observables, 256 versus 512 quadrature nodes | 3.575e-11 |

The physical identity tests use 32 randomized points; the quadrature test covers
all 136 nominal accepted points. In each current channel every photon-DA sector
has zero charge-structure coefficient. The hard terms plus the link response
give the unit charge separately for each quark line. The charge checks are also
recorded for every nominal point. Exact symbolic checks cover the DA equations
of motion, EM normalization, the full contact tensor reconstruction and the
previous two-point/pseudoscalar identities.

These checks establish the implementation's internal consistency. They do not
prove that the finite OPE and continuum model capture all hadronic effects.

# Artifacts and reproducibility

`calculation_setup.json` contains the input contract. `freeze_sha256.json`
records hashes of the central outputs and exported coefficient files.
`final_pointwise_results.csv` contains every central and pure-current result;
`final_sector_channel_structure.csv` contains the eight-structure, four-channel
sector decomposition on the nominal windows. Additional CSVs contain all
window decisions, controlled changes, angle curves, sensitivity probes and
partial uncertainty outputs.

The executable calculation and verification scripts are in the neighboring
`scripts` directory. The numerical engine reads our symbolic JSON exports and
never imports counterpart code or numerical outputs. The private independent
recalculation source files were not read. Earlier notebooks and step-numbered
tables remain a historical record; this note is the current numerical result.

![Magnetic Borel stability.](figures/magnetic_borel_stability.png){width=95%}

![Quadrupole Borel stability.](figures/quadrupole_borel_stability.png){width=95%}

\clearpage

# Appendix: structure and angle comparisons

These tables cover every extraction route used in the earlier vertex audit,
plus the new pole-free combinations, within the eight available OPE coefficients.
Arbitrary linear recombinations would give infinitely many redundant routes.
All values below were recalculated from our own kernels with the current
EOM-constrained input set. No counterpart values enter.

Let $x_i=\widehat\Pi_i/B_H$, $A=m_H^2$, and $Q=Q_H$.
The labels $pp,qp,pq,qq,g$ include the common factor $(\epsilon\cdot p)$;
$q\epsilon,\epsilon q,\epsilon p$ denote the corresponding unscaled tensors.
Define $G_3=-A(x_{qq}-x_{qp})$.

## Extraction key

| ID | Quantity / formula | Interpretation |
|---|---|---|
| C1 | $F_C=x_g/2$ | Pole-free charge check |
| C2 | $F_C=-Ax_{pp}/2$ | Raw $pp$ diagnostic |
| C3 | $F_C=Q$ | Exact input; not an extraction |
| C4 | $F_C=(-Ax_{qp}+M1)/2$ | $qp,q\epsilon$ closure |
| C5 | $F_C=(-Ax_{qp}+M3)/2$ | $qp,\epsilon q$ closure |
| C6 | $F_C=(-Ax_{qp}+M4)/2$ | $qp,pq$ closure |
| M1 | $F_M=-x_{q\epsilon}$ | Pole-free; central |
| M2 | $F_M=x_{\epsilon q}-x_{\epsilon p}$ | Pole-free cross-check |
| M3 | $F_M=x_{\epsilon q}$ | Raw $\epsilon q$ diagnostic |
| M4 | $F_M=-Ax_{pq}$ | Raw $pq$ diagnostic |
| M5 | $F_M=A(x_{qp}-x_{pp})$ | $pp,qp$ diagnostic |
| M6 | $F_M=x_g+Ax_{qp}$ | $g,qp$ diagnostic |
| M7 | $F_M=2Q+Ax_{qp}$ | $qp$ with exact charge |
| D1 | $F_D=Q-M1+G_3$ | Pole-free; central |
| D2 | $F_D=C1-M2+G_3$ | Pole-free cross-check |
| D3 | $F_D=-Ax_{qq}-Q$ | Raw $qq$ with exact charge |
| D4 | $F_D=-Ax_{qq}-C1$ | $qq,g$ diagnostic |
| D5 | $F_D=-Ax_{qq}-C2$ | $qq,pp$ diagnostic |
| D6 | $F_D=-Ax_{qq}-C4$ | $qq,qp,q\epsilon$ closure |
| D7 | $F_D=-Ax_{qq}-C5$ | $qq,qp,\epsilon q$ closure |
| D8 | $F_D=-Ax_{qq}-C6$ | $qq,qp,pq$ closure |

In formulas referencing another row, M1 means the M1 form factor, not a
Borel mass, and similarly for C1 etc. Diagnostic routes can admit pseudoscalar
contamination and are not averaged into the central prediction. C3 is explicitly
an input, so its integer value is not a successful numerical charge check.

Each number is the median of pointwise extractions, with its own rotated
two-point residue. The CSV summary also gives the 16th and 84th percentiles;
these are window spreads, not complete theoretical errors. Six digits identify
the computation and do not represent physical precision.

At $90^\circ$, $J_1=J_A$ and $J_2=-J_B$: state 1 uses AA and state 2 BB.
The two endpoint signs cancel in the elastic correlator. Mass labels and
hadronic pole assumptions are retained, so these are pure-current diagnostics,
not proof that a pure current isolates only the named physical meson.

The BB current is transverse: $P^\mu J_{B\mu}=0$. It cannot create a
pseudoscalar state. Differences between its diagnostic routes at $90^\circ$
therefore cannot all be attributed to pseudoscalar contamination. The tables
expose residual structure dependence of the truncated OPE, continuum treatment
and pole isolation; its origin has not been fully resolved. Agreement between
the two clean routes is an internal check, not a proof that all extractions
are mutually consistent or that this systematic effect is negligible.

## Clean-angle comparison

Changing the angle on frozen nominal windows and selecting new windows at
$90^\circ$ are distinct operations. Residues are rotated/recomputed in both.
The threshold interval and cuts retain the model limitations documented in
the independent scan audit.

The frozen-window case keeps all nominal points to isolate the angle change;
it is not a newly accepted window. The column "Pass at this angle" reports how
many of those points still satisfy the cuts. The refitted case applies all cuts
again. The independent checks cover 500 scenario-points with 256 and 512
quadrature nodes: charge and clean-route residuals are below $1.1\times10^{-13}$,
and the largest integration change across all routes is $3.6\times10^{-11}$.
These are numerical tolerances, not physical error estimates.

### $\mu/\mu_N$

| State | Nominal | $90^\circ$, frozen | $90^\circ$, refit |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | 0.241248 | 0.476918 | 0.506548 |
| $D_{s1}(2536)^+$ | 0.278616 | 0.210418 | 0.129385 |
| $B_{s1}(5750)$ | -0.065583 | -0.032857 | -0.023798 |
| $B_{s1}(5830)$ | 0.371700 | -0.095018 | -0.059648 |

### $F_D$

| State | Nominal | $90^\circ$, frozen | $90^\circ$, refit |
| --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | -1.574176 | 0.041234 | 0.202758 |
| $D_{s1}(2536)^+$ | 2.827267 | -1.423121 | -1.700287 |
| $B_{s1}(5750)$ | -1.168263 | 1.120063 | 1.118169 |
| $B_{s1}(5830)$ | 5.405695 | -3.127375 | -2.445490 |

## Nominal mixing

| State | Points | Pass at this angle | $M^2$ range | Median $f_H$ [GeV] |
| --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | 37 | 37 | 1.600–2.900 | 0.421056 |
| $D_{s1}(2536)^+$ | 25 | 25 | 2.300–3.000 | 0.160475 |
| $B_{s1}(5750)$ | 57 | 57 | 4.750–9.500 | 0.494055 |
| $B_{s1}(5830)$ | 17 | 17 | 5.250–6.750 | 0.083781 |

### Nominal mixing: $F_C(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| C1 | 1.000000 | 1.000000 | 0.000000 | 0.000000 |
| C2 | 1.037842 | 4.739841 | 0.000000 | 0.000000 |
| C3 | 1.000000 | 1.000000 | 0.000000 | 0.000000 |
| C4 | 0.812738 | -0.554130 | -0.406550 | 2.782338 |
| C5 | 1.080068 | -0.101710 | -0.304999 | 2.569166 |
| C6 | 1.037842 | 4.739841 | 0.000000 | 0.000000 |

### Nominal mixing: $F_M(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| M1 | 0.632386 | 0.752791 | -0.401926 | 2.309223 |
| M2 | 0.632386 | 0.752791 | -0.401926 | 2.309223 |
| M3 | 1.168825 | 1.662996 | -0.200405 | 1.880342 |
| M4 | 1.091496 | 11.331037 | 0.410716 | -3.231498 |
| M5 | 1.091496 | 11.331037 | 0.410716 | -3.231498 |
| M6 | 1.008689 | 3.849282 | 0.410716 | -3.231498 |
| M7 | 1.008689 | 3.849282 | 0.410716 | -3.231498 |

### Nominal mixing: $\mu/\mu_N$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| M1 | 0.241248 | 0.278616 | -0.065583 | 0.371700 |
| M2 | 0.241248 | 0.278616 | -0.065583 | 0.371700 |
| M3 | 0.445894 | 0.615493 | -0.032701 | 0.302666 |
| M4 | 0.416394 | 4.193741 | 0.067017 | -0.520152 |
| M5 | 0.416394 | 4.193741 | 0.067017 | -0.520152 |
| M6 | 0.384804 | 1.424662 | 0.067017 | -0.520152 |
| M7 | 0.384804 | 1.424662 | 0.067017 | -0.520152 |

### Nominal mixing: $F_D(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| D1 | -1.574176 | 2.827267 | -1.168263 | 5.405695 |
| D2 | -1.574176 | 2.827267 | -1.168263 | 5.405695 |
| D3 | -1.949355 | -0.280993 | -1.993503 | 10.969174 |
| D4 | -1.949355 | -0.280993 | -1.993503 | 10.969174 |
| D5 | -1.990758 | -3.982200 | -1.993503 | 10.969174 |
| D6 | -1.761765 | 1.273137 | -1.580856 | 8.185568 |
| D7 | -2.029423 | 0.820718 | -1.682559 | 8.400008 |
| D8 | -1.990758 | -3.982200 | -1.993503 | 10.969174 |

### Nominal mixing: $\mathcal D/e$ [fm$^2$]

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| D1 | -0.010133 | 0.017130 | -0.001376 | 0.006195 |
| D2 | -0.010133 | 0.017130 | -0.001376 | 0.006195 |
| D3 | -0.012548 | -0.001702 | -0.002348 | 0.012570 |
| D4 | -0.012548 | -0.001702 | -0.002348 | 0.012570 |
| D5 | -0.012814 | -0.024127 | -0.002348 | 0.012570 |
| D6 | -0.011340 | 0.007714 | -0.001862 | 0.009380 |
| D7 | -0.013063 | 0.004972 | -0.001981 | 0.009626 |
| D8 | -0.012814 | -0.024127 | -0.002348 | 0.012570 |

## $90^\circ$ on nominal windows

| State | Points | Pass at this angle | $M^2$ range | Median $f_H$ [GeV] |
| --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | 37 | 32 | 1.600–2.900 | 0.249258 |
| $D_{s1}(2536)^+$ | 25 | 21 | 2.300–3.000 | 0.419820 |
| $B_{s1}(5750)$ | 57 | 51 | 4.750–9.500 | 0.321631 |
| $B_{s1}(5830)$ | 17 | 3 | 5.250–6.750 | 0.422464 |

### $90^\circ$ on nominal windows: $F_C(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| C1 | 1.000000 | 1.000000 | 0.000000 | 0.000000 |
| C2 | 2.733221 | 0.967799 | 0.000000 | 0.000000 |
| C3 | 1.000000 | 1.000000 | 0.000000 | 0.000000 |
| C4 | 0.663601 | 0.721525 | -0.881607 | 0.082620 |
| C5 | 0.656635 | 1.059179 | -0.746402 | 0.156981 |
| C6 | 2.733221 | 0.967799 | 0.000000 | 0.000000 |

### $90^\circ$ on nominal windows: $F_M(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| M1 | 1.250150 | 0.568527 | -0.201365 | -0.590306 |
| M2 | 1.250150 | 0.568527 | -0.201365 | -0.590306 |
| M3 | 1.235971 | 1.246236 | 0.069630 | -0.441183 |
| M4 | 5.392601 | 1.079529 | 1.548473 | -0.754520 |
| M5 | 5.392601 | 1.079529 | 1.548473 | -0.754520 |
| M6 | 1.922948 | 1.127009 | 1.548473 | -0.754520 |
| M7 | 1.922948 | 1.127009 | 1.548473 | -0.754520 |

### $90^\circ$ on nominal windows: $\mu/\mu_N$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| M1 | 0.476918 | 0.210418 | -0.032857 | -0.095018 |
| M2 | 0.476918 | 0.210418 | -0.032857 | -0.095018 |
| M3 | 0.471509 | 0.461246 | 0.011362 | -0.071014 |
| M4 | 2.057218 | 0.399545 | 0.252668 | -0.121450 |
| M5 | 2.057218 | 0.399545 | 0.252668 | -0.121450 |
| M6 | 0.733583 | 0.417118 | 0.252668 | -0.121450 |
| M7 | 0.733583 | 0.417118 | 0.252668 | -0.121450 |

### $90^\circ$ on nominal windows: $F_D(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| D1 | 0.041234 | -1.423121 | 1.120063 | -3.127375 |
| D2 | 0.041234 | -1.423121 | 1.120063 | -3.127375 |
| D3 | -0.631564 | -1.982503 | -0.606328 | -2.963237 |
| D4 | -0.631564 | -1.982503 | -0.606328 | -2.963237 |
| D5 | -2.364785 | -1.960934 | -0.606328 | -2.963237 |
| D6 | -0.295165 | -1.702812 | 0.273496 | -3.046713 |
| D7 | -0.288076 | -2.046456 | 0.137789 | -3.123033 |
| D8 | -2.364785 | -1.960934 | -0.606328 | -2.963237 |

### $90^\circ$ on nominal windows: $\mathcal D/e$ [fm$^2$]

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| D1 | 0.000265 | -0.008622 | 0.001319 | -0.003584 |
| D2 | 0.000265 | -0.008622 | 0.001319 | -0.003584 |
| D3 | -0.004065 | -0.012011 | -0.000714 | -0.003396 |
| D4 | -0.004065 | -0.012011 | -0.000714 | -0.003396 |
| D5 | -0.015222 | -0.011881 | -0.000714 | -0.003396 |
| D6 | -0.001900 | -0.010317 | 0.000322 | -0.003491 |
| D7 | -0.001854 | -0.012399 | 0.000162 | -0.003579 |
| D8 | -0.015222 | -0.011881 | -0.000714 | -0.003396 |

## $90^\circ$ with refitted windows

| State | Points | Pass at this angle | $M^2$ range | Median $f_H$ [GeV] |
| --- | --- | --- | --- | --- |
| $D_{s1}(2460)^+$ | 55 | 55 | 1.700–3.500 | 0.251346 |
| $D_{s1}(2536)^+$ | 79 | 79 | 1.300–3.000 | 0.410735 |
| $B_{s1}(5750)$ | 64 | 64 | 5.000–10.500 | 0.320011 |
| $B_{s1}(5830)$ | 30 | 30 | 6.250–9.250 | 0.409471 |

### $90^\circ$ with refitted windows: $F_C(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| C1 | 1.000000 | 1.000000 | 0.000000 | 0.000000 |
| C2 | 2.470842 | 0.962759 | 0.000000 | 0.000000 |
| C3 | 1.000000 | 1.000000 | 0.000000 | 0.000000 |
| C4 | 0.641605 | 0.673179 | -0.838908 | 0.044506 |
| C5 | 0.646684 | 1.099618 | -0.702599 | 0.102508 |
| C6 | 2.470842 | 0.962759 | 0.000000 | 0.000000 |

### $90^\circ$ with refitted windows: $F_M(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| M1 | 1.327817 | 0.349583 | -0.145845 | -0.370571 |
| M2 | 1.327817 | 0.349583 | -0.145845 | -0.370571 |
| M3 | 1.335535 | 1.201434 | 0.127172 | -0.255235 |
| M4 | 5.017385 | 0.933450 | 1.532520 | -0.460251 |
| M5 | 5.017385 | 0.933450 | 1.532520 | -0.460251 |
| M6 | 2.043889 | 1.003225 | 1.532520 | -0.460251 |
| M7 | 2.043889 | 1.003225 | 1.532520 | -0.460251 |

### $90^\circ$ with refitted windows: $\mu/\mu_N$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| M1 | 0.506548 | 0.129385 | -0.023798 | -0.059648 |
| M2 | 0.506548 | 0.129385 | -0.023798 | -0.059648 |
| M3 | 0.509492 | 0.444664 | 0.020751 | -0.041083 |
| M4 | 1.914077 | 0.345480 | 0.250064 | -0.074083 |
| M5 | 1.914077 | 0.345480 | 0.250064 | -0.074083 |
| M6 | 0.779721 | 0.371305 | 0.250064 | -0.074083 |
| M7 | 0.779721 | 0.371305 | 0.250064 | -0.074083 |

### $90^\circ$ with refitted windows: $F_D(0)$

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| D1 | 0.202758 | -1.700287 | 1.118169 | -2.445490 |
| D2 | 0.202758 | -1.700287 | 1.118169 | -2.445490 |
| D3 | -0.517267 | -2.352973 | -0.527701 | -2.361967 |
| D4 | -0.517267 | -2.352973 | -0.527701 | -2.361967 |
| D5 | -2.000620 | -2.321647 | -0.527701 | -2.361967 |
| D6 | -0.157221 | -2.026630 | 0.309740 | -2.403490 |
| D7 | -0.159609 | -2.447370 | 0.173004 | -2.462395 |
| D8 | -2.000620 | -2.321647 | -0.527701 | -2.361967 |

### $90^\circ$ with refitted windows: $\mathcal D/e$ [fm$^2$]

| Route | $D_{s1}(2460)^+$ | $D_{s1}(2536)^+$ | $B_{s1}(5750)$ | $B_{s1}(5830)$ |
| --- | --- | --- | --- | --- |
| D1 | 0.001305 | -0.010302 | 0.001317 | -0.002802 |
| D2 | 0.001305 | -0.010302 | 0.001317 | -0.002802 |
| D3 | -0.003330 | -0.014256 | -0.000621 | -0.002707 |
| D4 | -0.003330 | -0.014256 | -0.000621 | -0.002707 |
| D5 | -0.012878 | -0.014066 | -0.000621 | -0.002707 |
| D6 | -0.001012 | -0.012279 | 0.000365 | -0.002754 |
| D7 | -0.001027 | -0.014828 | 0.000204 | -0.002822 |
| D8 | -0.012878 | -0.014066 | -0.000621 | -0.002707 |
