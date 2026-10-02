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
