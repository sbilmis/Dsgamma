# Finalization review of `main3_savci_222.pdf`

Review date: 24 August 2026  
Scope: attached 22-page PDF, checked against the Mathematica notebooks, batch regressions, and authoritative numerical outputs in `papers/Ds_Bs_gamma`. Existing manuscript sources were not used as the draft and were not edited.

**Editing status.** The current manuscript is online and is not accessible from this workspace. No online manuscript file has been edited. The instructions below are therefore a manual change list to be applied to the online draft. Table II should be retained exactly as written; the regenerated results and figures use those state-specific intervals.

## Executive verdict

The updated symbolic and numerical calculation passes the available independent Mathematica checks, but the PDF is **not ready for submission in its present form**. The transition OPE is now internally consistent, with finite strange-quark-mass corrections retained. The state-specific continuum-threshold intervals printed in Table II are theoretically defensible and have therefore been retained. The abstract, Table III, figures, and central-amplitude discussion must be synchronized with the regenerated outputs that use those intervals.

Recommended referee decision: **major revision**, followed by a short consistency re-review.

## Manual editing sequence for the online draft

Apply the following changes in order so that the formulas, numerical tables, figures, and prose all describe the same calculation.

1. **Keep Table II unchanged.** Retain
   - `Ds1(2460)`: `M^2 = 3.0--4.5 GeV^2`, `s0 = 8.5--9.5 GeV^2`;
   - `Ds1(2536)`: `M^2 = 3.0--4.5 GeV^2`, `s0 = 9.0--10.0 GeV^2`;
   - `Bs1(5750)`: `M^2 = 10.0--14.0 GeV^2`, `s0 = 39.0--41.0 GeV^2`;
   - `Bs1(5830)`: `M^2 = 10.0--14.0 GeV^2`, `s0 = 40.0--42.0 GeV^2`.

2. **Replace the four widths in the abstract** by
   - `45.92^{+10.09}_{-8.51} keV`;
   - `0.55^{+1.73}_{-0.48} keV`;
   - `12.03^{+2.87}_{-2.31} keV`;
   - `8.88^{+2.38}_{-2.01} keV`,
   in the channel order used in the manuscript. Replace the sentence claiming a sub-keV bottom-state suppression: the two bottom widths are comparable within uncertainties, although the higher state still exhibits destructive interference.

3. **Insert the continuum-prescription paragraph immediately after Eq. (26):**

   > In the continuum subtraction, the twist-2, twist-3, two-particle twist-4, gluonic three-particle, and electromagnetic three-particle photon-DA contributions are all multiplied by \(e^{-m_Q^2/M^2}-e^{-s_0/M^2}\). The strange-quark mass is retained in the perturbative threshold \((m_Q+m_s)^2\), the spectral densities, the current normalization, and the mass-dependent coefficient functions; it does not modify this soft-sector continuum factor.

4. **Leave Eq. (27) unchanged.** Its twist-3 and twist-4 exponential differences are correct.

5. **Replace Eq. (31)** by

   \[
   \widehat\Pi_B^{\mathrm{twist3}}
   =\frac{e_s f_{3\gamma}}{m_Q+m_s}
   \left(e^{-m_Q^2/M^2}-e^{-s_0/M^2}\right)
   \left[m_Q^2\psi^{(v)}(u_0)
   -\frac{M^2}{2}\left.\frac{d}{da}
   \{(1-a)\psi^{(v)}(a)\}\right|_{a=u_0}\right].
   \]

   Retain `a=u0` in the derivative evaluation: this is the notation already used in Eq. (31) of the draft.

6. **Update Eq. (B6).** In both the axial and tensor gluonic three-particle lines, replace
   `exp(-m_Q^2/M^2)` by
   `exp(-m_Q^2/M^2)-exp(-s0/M^2)`.

7. **Replace Table III** with the final magnitude percentiles and widths:

   | Channel | `|g| (GeV^-1)` | `Gamma (keV)` |
   |---|---:|---:|
   | `Ds1(2460) -> Ds gamma` | `0.47_{-0.05}^{+0.05}` | `45.92_{-8.51}^{+10.09}` |
   | `Ds1(2536) -> Ds gamma` | `0.04_{-0.03}^{+0.04}` | `0.55_{-0.48}^{+1.73}` |
   | `Bs1(5750) -> Bs gamma` | `0.32_{-0.02}^{+0.02}` | `12.03_{-2.31}^{+2.87}` |
   | `Bs1(5830) -> Bs gamma` | `0.20_{-0.02}^{+0.03}` | `8.88_{-2.01}^{+2.38}` |

8. **Replace the central-amplitude paragraph** by the values evaluated at `M^2=3.75 GeV^2`, `s0=9.0 GeV^2` for `Ds1(2460)`, `s0=9.5 GeV^2` for `Ds1(2536)`, and `theta_Ds=26.6 degrees`:
   - `Ds1(2460)`: `g_A=-0.09604`, `g_B=-0.37832`, `g=-0.47436 GeV^-1`;
   - `Ds1(2536)`: `g_A=-0.51062`, `g_B=+0.48635`, `g=-0.02426 GeV^-1`.

9. **Replace all three figure files in the online project.** Use the regenerated files named `fig1.pdf`, `fig2.pdf`, and `fig3.pdf` from `papers/Ds_Bs_gamma/arxiv/figures/`. Do not keep the old online copies: the filenames are unchanged, so they must be overwritten or re-uploaded. Figure 1 uses the Table II threshold intervals; Figure 2 uses a central accepted sample; Figure 3 uses the same preferred-prior Monte Carlo ensemble as Table III.

10. **Synchronize the remaining numerical statements.** In the literature-comparison table, replace the `This work` row with the four final widths above. Replace the `Ds1(2536)` branching-fraction estimate by `0.60 x 10^-3`, with the approximate 16th--84th percentile interval `0.08 x 10^-3` to `2.46 x 10^-3`. Use the same four widths in the conclusion.

11. **Correct the photon-input attribution.** In Table V, attribute `f_{gamma,s}^perp=-55.1(1.9) MeV` to the lattice reference presently numbered `[28]`; state that the DA shape parameters are from the photon-DA reference presently numbered `[25]`.

12. **Perform a final global search** for the stale values `32.7`, `0.425`, `19.3`, `0.504`, `34.96`, `0.99`, `11.34`, `8.70`, `1.07 x 10^-3`, `-0.07`, `-0.21`, `-0.09`, and `+0.13`. None should remain in the abstract, main text, tables, captions, or conclusion except inside an explicitly labeled historical comparison.

## 1. Calculation and Mathematica audit

The following checks were rerun successfully with Mathematica 15.0 and FeynCalc 10.2.1:

- Current/Wick/trace audit: 19/19 exact checks passed.
- Correlator momentum routing: standard routing and both fermion-line Ward identities passed.
- Hard-loop Feynman parameterization: 4/4 exact denominator checks passed.
- State rotation, residue projection, radiative width, and pure-axial limit: passed.
- Charm AA/AB/BB two-point matrix: Mathematica agrees with the Python central point at the stated regression precision.
- Bottom AA/AB/BB two-point matrix: passed at the stated `10^-6` tolerance.
- Corrected transition numerics after applying the adopted continuum prescription to the complete soft photon-DA sector: all 187 central quantities agree term by term.

The authoritative preferred-prior, lattice-normalized outputs are:

| Channel | Signed `g` median and 16–84% interval (GeV^-1) | `|g|` median and 16–84% interval (GeV^-1) | Width median and 16–84% interval (keV) |
|---|---:|---:|---:|
| Ds1(2460) -> Ds gamma | -0.4676 [-0.5169, -0.4210] | 0.4676 [0.4210, 0.5169] | 45.9157 [37.4092, 56.0037] |
| Ds1(2536) -> Ds gamma | -0.02545 [-0.08186, +0.03342] | 0.04212 [0.01511, 0.08572] | 0.5503 [0.0709, 2.2806] |
| Bs1(5750) -> Bs gamma | +0.3174 [+0.2981, +0.3379] | 0.3174 [0.2981, 0.3379] | 12.0336 [9.7282, 14.9082] |
| Bs1(5830) -> Bs gamma | -0.2043 [-0.2305, -0.1798] | 0.2043 [0.1798, 0.2305] | 8.8775 [6.8673, 11.2598] |

Neither the attached PDF's abstract nor its Table III contains this recalculated result set. The retained continuum windows are centered close to the conventional estimates `s0 = (m_state + 0.5 GeV)^2`: `8.76`, `9.21`, `39.06`, and `40.05 GeV^2`, respectively. The modest extension around these estimates is acceptable because `s0` is an effective duality threshold rather than a measured excitation mass. On the adopted grids, the minimum `(axial spectral proxy, tensor spectral proxy, twist-2 retention)` values are `(0.58, 0.67, 0.78)`, `(0.61, 0.71, 0.81)`, `(0.40, 0.67, 0.79)`, and `(0.41, 0.69, 0.80)`. The bottom axial proxy falls below the common 50% rule of thumb at the least favorable endpoints, so it should be reported transparently as a limitation, but this proxy is not the exact pole fraction of the full transition correlator. In addition, the signed `Ds1(2536)` amplitude crosses zero inside part of the window; its lack of a conventional relative plateau is a physical cancellation effect and must be represented by the broad non-Gaussian uncertainty rather than hidden by choosing a narrower interval.

## 2. Mandatory scientific corrections

### 2.1 Abstract: all four headline widths are stale

Page 1 reports

- 32.7 keV for Ds1(2460),
- 0.425 keV for Ds1(2536),
- 19.3 keV for Bs1(5750), and
- 0.504 keV for Bs1(5830).

Those values belong to an older calculation branch. Replace them with the regenerated preferred-prior, lattice-normalized results:

\[
\Gamma[D_{s1}(2460)\to D_s\gamma]=45.92^{+10.09}_{-8.51}\ \mathrm{keV},
\]
\[
\Gamma[D_{s1}(2536)\to D_s\gamma]=0.55^{+1.73}_{-0.48}\ \mathrm{keV},
\]
\[
\Gamma[B_{s1}(5750)\to B_s\gamma]=12.03^{+2.87}_{-2.31}\ \mathrm{keV},
\]
\[
\Gamma[B_{s1}(5830)\to B_s\gamma]=8.88^{+2.38}_{-2.01}\ \mathrm{keV}.
\]

The qualitative bottom-sector sentence must also be softened. With the current results, the two bottom widths are comparable; the higher state still has destructive interference, but it is not suppressed to the sub-keV value stated in the old abstract.

### 2.2 Continuum prescription: Eq. (27) is correct, but Eqs. (31) and (B6) must be updated

The earlier objection to the twist-3 and twist-4 factors in Eq. (27) is withdrawn. The consistency check confirms that the twist-2, twist-3, two-particle twist-4, gluonic three-particle, and electromagnetic three-particle contributions all contain the continuum-subtracted factor

\[
e^{-m_Q^2/M^2}-e^{-s_0/M^2}.
\]

Consequently, the exponential differences already printed for the axial twist-3 and twist-4 terms in Eq. (27) are correct and should be retained. Finite `m_s` does not alter this soft-sector continuum factor; it is retained in the perturbative threshold `(m_Q+m_s)^2`, spectral densities, current normalization, and mass-dependent coefficient functions.

Two other displayed formulas must be synchronized with the adopted prescription. The tensor twist-3 expression in Eq. (31) should read

\[
\widehat\Pi_B^{\mathrm{twist3}}
=\frac{e_s f_{3\gamma}}{m_Q+m_s}
\left(e^{-m_Q^2/M^2}-e^{-s_0/M^2}\right)
\left[m_Q^2\psi^{(v)}(u_0)
-\frac{M^2}{2}\left.\frac{d}{da}
\{(1-a)\psi^{(v)}(a)\}\right|_{a=u_0}\right].
\]

Likewise, in both lines of Eq. (B6), the factor `exp(-m_Q^2/M^2)` must be replaced by `exp(-m_Q^2/M^2)-exp(-s_0/M^2)`. Equation (B2) already has the required exponential difference.

Add an explicit paragraph after Eq. (26) stating the term-by-term continuum-subtraction rule and the finite-`m_s` extension. This is also the natural place to state the sign convention for `f_{gamma,s}^perp=chi_s<bar s s>`.

Suggested text:

> In the continuum subtraction, the twist-2, twist-3, two-particle twist-4, gluonic three-particle, and electromagnetic three-particle photon-DA contributions are all multiplied by `exp(-m_Q^2/M^2)-exp(-s_0/M^2)`. The strange-quark mass is retained in the perturbative spectral densities, the physical threshold, the interpolating-current normalization, and the relevant coefficient functions.

### 2.3 Page 13 amplitude decomposition is stale

The paragraph quoting approximately `-0.07` and `-0.21` for Ds1(2460), and `-0.09` and `+0.13` for Ds1(2536), does not match the final transition normalization.

At the corrected Mathematica central checkpoints (`M^2 = 3.75 GeV^2`, `s0 = 9.0 GeV^2` for the lower state and `s0 = 9.5 GeV^2` for the higher state, external angle 26.6 degrees), the decomposition is

- Ds1(2460): `g_A = -0.09604`, `g_B = -0.37832`, total `g = -0.47436 GeV^-1`;
- Ds1(2536): `g_A = -0.51062`, `g_B = +0.48635`, total `g = -0.02426 GeV^-1`.

Replace the paragraph with values from one clearly identified central prescription. Do not mix the angle-scan figure's reprojected residues with the fixed central checkpoint. The caption should state whether the two-point residues are reprojected at every plotted angle; the current plotting data do this.

### 2.4 Table III form-factor rounding

Use the magnitude percentiles, not the signed-coupling percentiles, in Table III. The four rounded entries are

\[
|g|=0.47^{+0.05}_{-0.05},\quad
0.04^{+0.04}_{-0.03},\quad
0.32^{+0.02}_{-0.02},\quad
0.20^{+0.03}_{-0.02}\ \mathrm{GeV}^{-1},
\]

in the channel order used in Table III.

### 2.5 Table V source attribution is wrong

The Table V caption says all values are from Ref. [25], while the `f_{gamma,s}^perp` row cites Ref. [29]. The final input `-55.1(1.9) MeV` is the lattice-normalized value from Ref. [28]. Revise the caption to say that the DA shape parameters are from Ref. [25] and the leading-twist normalization is from Ref. [28]. Change the row citation from `[29]` to `[28]`.

## 3. Symbol and convention audit

These points should be corrected globally.

1. **`q` denotes both the photon momentum and the light-quark field.** This is the most serious notation collision. Keep `q` for photon momentum and denote the light field by `q_l`, or simply use `s` throughout because every numerical channel is strange.

2. **Electric charges are ambiguous.** State explicitly that `e_q` and `e_Q` are fractional charges in units of the positron charge `e`. Then use this convention consistently in the electromagnetic current, the background-field propagator, and Eq. (10). As written, “electric charges” plus the extra factor `e` in Eq. (10) can be read as double counting.

3. **Map the basis currents to spectroscopic states.** State explicitly whether `J_A` interpolates `^3P_1` and `J_B` interpolates `^1P_1`, including the phase/sign convention.

4. **Map state labels.** Define `A_1` as the lower physical state and `A_2` as the higher physical state in each flavor sector immediately after Eq. (1) or Eq. (6).

5. **Mass notation varies between `m_i` and uppercase hadron masses in Table I.** Choose one convention. A clean choice is `m_{A_i}` and `m_P` everywhere.

6. **“Decay constants” and “residues” are used interchangeably for `f_i`.** Equation (8) defines `f_i`; call these axial-current decay constants throughout, or define “pole residue `m_i f_i`” separately.

7. **Use the directly normalized photon input in the formula.** Since the numerical analysis samples `f_{gamma,s}^perp = chi_s <sbar s>`, writing the twist-2 term directly in terms of `f_{gamma,s}^perp` would remove a sign/scale ambiguity. If `chi_s` is retained, define the identity and sign convention explicitly.

8. **The tensor current contains `p'^nu`.** Explain whether this is shorthand for a derivative current and give its local operator definition. Otherwise the position-space current in Eqs. (2), (7), and (33) is not fully specified.

9. **Distinguish transition and two-point auxiliary parameters.** The two-point extraction uses its own Borel windows and fitted thresholds. These are not listed in Table II. Add a compact table or a precise paragraph giving the two-point windows, thresholds, pole fractions, and OPE-convergence criteria used to obtain the four `f_i` values.

## 4. Language corrections

The manuscript is generally readable, but the following edits are needed.

- Page 4: “Summing over polarizations of the axial vector meson, for the correlation function from the hadronic side we get,” -> “After summing over the axial-vector polarization, the hadronic representation becomes”.
- Page 7: “Equating ... and performing double Borel transformations ... suppress the contributions” -> “Equating the two representations and applying double Borel transformations suppresses higher-state and continuum contributions.”
- Page 8: “the perturbative spectral density is given as,” -> “the perturbative spectral density is given by”.
- Page 9: “Equation (23) constitute the main theoretical result” -> “Equation (23) constitutes the main theoretical result”.
- Page 10: “If the values of the parameters chi_s and <sbar s> are separately used the uncertainty...” -> “If chi_s and <sbar s> are used separately, the uncertainty in the magnetic susceptibility makes a sizable contribution...”.
- Page 10: revise the mass-source sentence because the unobserved Bs1(5750) mass is from lattice QCD, not the PDG.
- Page 16: “non-perturnative” -> “nonperturbative”.
- Page 16: “we use the photon distribution amplitudes DAs” -> “we use the photon DAs”.
- Page 17, Table V: “parametersare” -> “parameters are”.
- Page 18: “These DAs are given as,” -> “These DAs are given by”.
- Page 18: “parameters appearing in the photon DAs analysis” -> “parameters entering the photon DAs”.
- Use “axial-vector” adjectivally throughout (for example, “axial-vector meson” and “axial-vector current”).
- Use one form consistently: “light-cone QCD sum rules” for the method and “light-cone QCD sum-rule calculation” adjectivally.

Minor numerical consistency edit: use the same regenerated central values, percentile intervals, and significant-digit policy in the abstract, Table III, Eqs. (40)–(41), and Table IV.

## 5. Suggested corrected abstract

We calculate the radiative decay widths of Ds1(2460) -> Ds gamma, Ds1(2536) -> Ds gamma, Bs1(5750) -> Bs gamma, and Bs1(5830) -> Bs gamma using light-cone QCD sum rules, including mixing between the `^1P_1` and `^3P_1` configurations. Using a lattice-normalized leading-twist photon distribution amplitude, we obtain

\[
\Gamma[D_{s1}(2460)\to D_s\gamma]=45.92^{+10.09}_{-8.51}\ \mathrm{keV},\qquad
\Gamma[D_{s1}(2536)\to D_s\gamma]=0.55^{+1.73}_{-0.48}\ \mathrm{keV},
\]
\[
\Gamma[B_{s1}(5750)\to B_s\gamma]=12.03^{+2.87}_{-2.31}\ \mathrm{keV},\qquad
\Gamma[B_{s1}(5830)\to B_s\gamma]=8.88^{+2.38}_{-2.01}\ \mathrm{keV}.
\]

The large hierarchy between the two charm-strange widths arises from constructive interference between the basis-current contributions for Ds1(2460) and destructive interference for Ds1(2536). Destructive interference also occurs for Bs1(5830), although at the preferred bottom-sector mixing angle its width remains comparable to that of the predicted lower Bs1 state. These results demonstrate the sensitivity of heavy-strange axial-vector radiative decays to `^1P_1-^3P_1` mixing.

## 6. Referee-style final remarks

The subject is timely, the use of independent axial and tensor interpolating currents is interesting, and the cancellation pattern is physically valuable. The computational record is stronger than the present PDF suggests: the symbolic routing, Ward identities, two-point matrices, state projection, and final central transition quantities all have independent Mathematica regressions and pass them.

However, the manuscript currently reports incompatible result generations. This is not a cosmetic issue: the abstract's four principal predictions, Table III, the central amplitude decomposition, and the raw continuum factors in Eqs. (31) and (B6) do not describe the now-consistent calculation. Table II's state-specific threshold intervals should be retained, and the numerical results and all three figures should be generated from those same intervals. The paper should therefore undergo a full synchronization pass from the regenerated numerical output before submission.

After the mandatory corrections above, the revised paper should also document (i) the term-by-term continuum prescription for the complete soft photon-DA sector while finite `m_s` corrections are retained, (ii) the two-point windows and quantitative pole/OPE criteria, (iii) the correlation treatment between mixing-angle and residue samples, and (iv) the renormalization-scale treatment of quark masses, condensates, and photon DAs. Subject to those revisions, the work could become suitable for publication.
