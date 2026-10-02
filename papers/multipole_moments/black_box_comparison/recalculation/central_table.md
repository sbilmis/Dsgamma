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

Full intervals and contribution medians: lattice_results.md.
