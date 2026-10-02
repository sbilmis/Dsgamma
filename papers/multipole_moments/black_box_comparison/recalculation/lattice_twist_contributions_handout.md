---
title: "Electromagnetic multipole moments"
subtitle: "Separate perturbative and twist contributions with lattice normalization"
date: "3 October 2026"
geometry: margin=0.7in
fontsize: 10pt
colorlinks: true
header-includes:
  - \usepackage{booktabs,longtable,amsmath,amssymb}
  - \usepackage{needspace}
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhf{}
  - \fancyfoot[C]{\thepage}
  - \renewcommand{\headrulewidth}{0pt}
---

# What does the contact contribution mean?

In these tables, **contact** is the retained electromagnetic Wilson-line
response of the local condensate part of the two-point OPE. A Wilson line
connects the separated charged fields and supplies their gauge phase.
Expanding that line to first order in the external photon field gives an
additional contribution to the three-point correlator.

With the straight-link prescription used in this calculation, the stored
response is

$$
\Pi^{W}_{\mu\nu}=
-Q_H\epsilon^\rho\int_0^1 d\beta\,
\frac{\partial C^{\mathrm{cond}}_{0,\mu\nu}(P;p,p')}{\partial P^\rho}
\bigg|_{P=p+\beta q},\qquad Q_H=e_Q-e_s.
$$

Here $C_0^{\mathrm{cond}}$ is the physical two-point condensate tensor;
the minus sign follows our stored convention $\Pi=-C$.
Only the internal momentum $P$ is differentiated. External current momenta
$p$ and $p'=p+q$ remain fixed. Replacing $\epsilon$ by $q$ turns the integral
into the difference of condensate tensors at the two endpoints, giving the
required electromagnetic Ward response of this block.

This contribution is kept separate from the photon distribution-amplitude
sectors. The perturbative photon-on-quark response is already included in
the perturbative sector and is not added again here. The link response
survives the double Borel transform, so the word "contact" does not mean
that it is discarded from the moments. Ward-response checks support this
prescription; they do not uniquely determine all transverse contributions.

**Why is it zero for Bs?** Every retained link kernel is proportional to
$Q_H$. For these charged Ds currents, $Q_H=1$; for the neutral Bs currents,
$Q_H=0$. The Bs contact entries therefore vanish exactly in this prescription.
The contact block is independent of $\chi$ and is unchanged by reversing twist 2.

# Inputs and reading the tables

The supplied lattice normalization is $f_{\perp\gamma,s}(1\,\mathrm{GeV})=-55.1$
MeV. The $+\chi$ label uses this value; $-\chi$ uses $+55.1$ MeV and flips
only twist 2. The condensate remains
$\langle\bar ss\rangle=0.8(-0.240\,\mathrm{GeV})^3$ in all other terms.
Chi enters only twist 2, not the two-point function or contact block.

Other inputs: $\zeta_1=0.4$, $\zeta_2=31/110$, $\kappa=0.2$,
$\kappa^+=0$, $f_{3\gamma}=-0.004$ GeV$^2$,
$\omega_V=3.8$, $\omega_A=-2.1$, $m_0^2=0.8$ GeV$^2$.
The DA scale is 1 GeV. Nominal angles are $26.6^\circ$ (Ds) and $38.5^\circ$ (Bs).

\clearpage

# Nominal mixing angles

Accepted-window counts in state order: 37, 25, 57, 17. 

Each entry is the median of the indicated contribution. Twist-2 columns give both signs explicitly; the other columns apply to both signs.

## Magnetic moment: $\mu/\mu_N$

| State | Perturbative | Twist 2 (+$\chi$) | Twist 2 (-$\chi$) | Twist 3 | Twist 4 | Contact |
|---|---:|---:|---:|---:|---:|---:|
| $D_{s1}(2460)^+$ | 0.780010 | -0.568491 | 0.568491 | -0.094003 | -0.069647 | -0.015260 |
| $D_{s1}(2536)^+$ | 0.638959 | 0.051564 | -0.051564 | -0.337387 | -0.156162 | 0.097233 |
| $B_{s1}(5750)$ | 0.133431 | -0.314355 | 0.314355 | -0.017663 | 0.017337 | 0.000000 |
| $B_{s1}(5830)$ | 0.480480 | 0.023798 | -0.023798 | -0.329725 | 0.206843 | 0.000000 |

## Quadrupole form factor: $F_D-Q_H$

| State | Perturbative | Twist 2 (+$\chi$) | Twist 2 (-$\chi$) | Twist 3 | Twist 4 | Contact |
|---|---:|---:|---:|---:|---:|---:|
| $D_{s1}(2460)^+$ | -2.053195 | -2.643512 | 2.643512 | 0.590327 | -0.154038 | 0.701903 |
| $D_{s1}(2536)^+$ | -2.216542 | 5.355061 | -5.355061 | 1.791186 | -0.072485 | -1.020116 |
| $B_{s1}(5750)$ | -0.147135 | -1.432447 | 1.432447 | -0.387724 | 0.266620 | 0.000000 |
| $B_{s1}(5830)$ | -5.728128 | 18.334977 | -18.334977 | 8.001665 | -8.241479 | 0.000000 |

## Grouping and normalization

**Perturbative:** heavy- and strange-line photon emission.
**Twist 2:** leading photon DA with lattice normalization.
**Twist 3:** two-particle twist 3 plus three-particle V and A.
**Twist 4:** two-particle twist 4 plus S, S-tilde, T1-T4 and electromagnetic
S-gamma and T4-gamma. **Contact:** Wilson-line condensate response described
on page 1.

The magnetic contributions use $F_M$ from $q_\mu\epsilon_\nu$ and
$\mu/\mu_N=(m_p/m_H)F_M$.
The quadrupole columns give the individual terms in
$F_D-Q_H=-F_M+G_3$ using the clean $G_3$ projection.
To obtain $F_D$, add $Q_H=1$ for Ds and $Q_H=0$ for Bs at each window.
The dimensional quadrupole moment is
$\mathcal D/e=F_D(\hbar c)^2/m_H^2$ in fm$^2$.

**Medians of the parts need not sum to the median total.** Additivity holds
at each window: the largest residual is $3.91\times10^{-14}$.
The six printed decimals are retained for reproducibility, not an estimate
of physical precision. Accepted-window spreads are not a full theoretical
uncertainty; the results use the central LO/twist and continuum prescription.

\clearpage

# 90 degrees: refitted pure-current windows

Accepted-window counts in state order: 55, 79, 64, 30. The lighter state uses $J_A$ and the heavier state uses $J_B$.

Each entry is the median of the indicated contribution. Twist-2 columns give both signs explicitly; the other columns apply to both signs.

## Magnetic moment: $\mu/\mu_N$

| State | Perturbative | Twist 2 (+$\chi$) | Twist 2 (-$\chi$) | Twist 3 | Twist 4 | Contact |
|---|---:|---:|---:|---:|---:|---:|
| $D_{s1}(2460)^+$ | 0.839351 | -0.476423 | 0.476423 | 0.026448 | -0.057200 | -0.002975 |
| $D_{s1}(2536)^+$ | 0.707351 | -0.470276 | 0.470276 | -0.188986 | -0.089902 | 0.000371 |
| $B_{s1}(5750)$ | 0.160582 | -0.294854 | 0.294854 | -0.019698 | 0.022629 | 0.000000 |
| $B_{s1}(5830)$ | 0.114196 | -0.260607 | 0.260607 | -0.020594 | 0.012360 | 0.000000 |

## Quadrupole form factor: $F_D-Q_H$

| State | Perturbative | Twist 2 (+$\chi$) | Twist 2 (-$\chi$) | Twist 3 | Twist 4 | Contact |
|---|---:|---:|---:|---:|---:|---:|
| $D_{s1}(2460)^+$ | -0.283352 | 1.248852 | -1.248852 | -1.329886 | 0.005127 | 0.026241 |
| $D_{s1}(2536)^+$ | -2.984153 | -3.452141 | 3.452141 | 1.946491 | -0.256859 | 0.786012 |
| $B_{s1}(5750)$ | 1.032532 | 1.807016 | -1.807016 | -0.828522 | -0.207140 | 0.000000 |
| $B_{s1}(5830)$ | -1.204771 | -2.942746 | 2.942746 | 0.277191 | 0.343996 | 0.000000 |

## Grouping and normalization

**Perturbative:** heavy- and strange-line photon emission.
**Twist 2:** leading photon DA with lattice normalization.
**Twist 3:** two-particle twist 3 plus three-particle V and A.
**Twist 4:** two-particle twist 4 plus S, S-tilde, T1-T4 and electromagnetic
S-gamma and T4-gamma. **Contact:** Wilson-line condensate response described
on page 1.

The magnetic contributions use $F_M$ from $q_\mu\epsilon_\nu$ and
$\mu/\mu_N=(m_p/m_H)F_M$.
The quadrupole columns give the individual terms in
$F_D-Q_H=-F_M+G_3$ using the clean $G_3$ projection.
To obtain $F_D$, add $Q_H=1$ for Ds and $Q_H=0$ for Bs at each window.
The dimensional quadrupole moment is
$\mathcal D/e=F_D(\hbar c)^2/m_H^2$ in fm$^2$.

**Medians of the parts need not sum to the median total.** Additivity holds
at each window: the largest residual is $3.91\times10^{-14}$.
The six printed decimals are retained for reproducibility, not an estimate
of physical precision. Accepted-window spreads are not a full theoretical
uncertainty; the results use the central LO/twist and continuum prescription.

## Provenance

All entries are read from our signed summary CSV frozen in Git commit
`54e6958`; no counterpart numerical outputs were used. The source definitions
and straight-link audit are in our recalculation report and
`blind_definition_scan_audit/AUDIT.md`. The original CSV SHA-256 hashes are
unchanged. A separate manifest covers this handout and its PDF.
