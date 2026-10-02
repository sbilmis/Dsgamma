# Blind lattice-normalized twist split

Computed with our numerical engine and own symbolic coefficient exports. No counterpart code, numbers, or result files were read. The accepted windows were regenerated from the two-point OPE and matched our archived coordinates: nominal 37, 25, 57, 17 (136 total); refitted theta90 55, 79, 64, 30 (228 total). There are 728 signed window rows and 16 summary rows.

Inputs at DA scale 1 GeV: zeta1=0.4, zeta2=31/110, kappa=0.2, kappaPlus=0, f3Gamma=-0.004 GeV^2, omegaV=3.8, omegaA=-2.1, ssCondensate=0.8*(-0.240)^3=-0.011059199999999998 GeV^3, m0^2=0.8 GeV^2. zeta1Plus=zeta2Plus=phi2=0; ms=0.093 GeV. The supplied lattice normalization -55.1 MeV is used directly (no extra running factor): chi=4.9822771990740753 GeV^-2. +1 means f_perp=-0.0551 GeV; -1 means f_perp=+0.0551 GeV. The labels denote the relative twist-2 sign, not the numerical sign of chi.

Masses/angles/heavy-quark masses: Ds_1(2460): 2.4595 GeV, 26.6 deg, mc=1.27 GeV; Ds_1(2536): 2.53511 GeV, 26.6 deg, mc=1.27 GeV; Bs_1(5750): 5.7502 GeV, 38.5 deg, mb=4.18 GeV; Bs_1(5830): 5.82911 GeV, 38.5 deg, mb=4.18 GeV. theta90 uses J_A for the lighter and J_B for the heavier state, with refitted windows and consistently rotated residues. Q_H=1 for charm and 0 for bottom.

Grouping: hard=pertQ+pertS; tw2=tw2; tw3=two-particle tw3V+tw3A and three-particle V+A; tw4=two-particle tw4A+tw4B and three-particle S, S-tilde, T1–T4, S_gamma, T4_gamma; contact=Wilson-line W. These groups exhaust every evaluated sector.

chi enters ONLY twist 2, and does NOT enter the two-point function. An AST scan of all coefficient exports found chi only in sector tw2. On every window, evaluating chi=0 leaves all other sectors identical and removes tw2; evaluating -chi is exactly identical to multiplying tw2 alone by -1. All two-point matrices, derivatives, full integrals and individual pieces are identical for +chi, -chi and zero chi. Largest sign-equivalence residual: 0; two-point chi-dependence residual: 0. ssCondensate is unchanged in every other term.

Extraction: F_M=-Pi_qeps/N; G3=-m_H^2*(Pi_qq-Pi_qp)/N; F_D=Q_H-F_M+G3. mu/mu_N=0.93827208816 GeV*F_M/m_H; D/(e fm^2)=F_D*(0.1973269804 GeV fm)^2/m_H^2; F_C_g=Pi_g/(2N). N is the rotated two-point Borel residue, and f=sqrt(N*exp(m_H^2/M^2)/m_H^2).

CSV convention: FD_hard, FD_tw2, FD_tw3, FD_tw4, FD_contact sum to FD_without_QH_total. FD_total includes Q_H. mu_over_muN_total is evaluated from the complete amplitude independently of the parts. Summary FD quantiles include Q_H; part medians exclude Q_H. Quantiles use NumPy linear interpolation; these are deterministic accepted-window spreads. Medians of parts need not sum to the median total.

Largest per-window additivity residuals: mu/mu_N 3.4972025275692431e-15; F_D without Q_H 3.907985046680551e-14. Largest overall: 3.907985046680551e-14. Largest FC_g-Q_H residual: 1.0678345879335023e-13. Largest regenerated-vs-archived two-point residue residual: 4.9960036108132044e-16 GeV^4. Doubling three-point Gauss-Legendre quadrature from 256 to 512 changes mu/mu_N, FD and FC_g by at most 3.5822012023345451e-11 in absolute units (common two-point normalization).

Continuum treatment is the central engine's bbkr prescription: positive Borel powers use the diagonal spectral integral, nonpositive DA powers are unsubtracted. This is the same model approximation as the central calculation; no window refit based on three-point results is performed.

The SHA-256 manifest was written after all outputs and checks and before delivery. Source hashes are included for reproducibility. The manifest itself is tracked in the Git commit; self-hashing it within the manifest is not possible.
