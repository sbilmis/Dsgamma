"""Write the current calculation note; exchange reply requires --write-reply."""
import csv
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"recalculation"
STATES=("Ds_1(2460)","Ds_1(2536)","Bs_1(5750)","Bs_1(5830)")
LABELS=(r"$D_{s1}(2460)^+$",r"$D_{s1}(2536)^+$",r"$B_{s1}(5750)$",r"$B_{s1}(5830)$")
def read(name):return list(csv.DictReader((OUT/name).open()))
def table(headers,rows):return "| "+" | ".join(headers)+" |\n| "+" | ".join(["---"]*len(headers))+" |\n"+"\n".join("| "+" | ".join(map(str,r))+" |" for r in rows)+"\n"
def f(r,k,n=6):return f"{float(r[k]):.{n}f}"
def ci(r,k,n=4):return f"[{f(r,k+'_p16',n)}, {f(r,k+'_p84',n)}]"

summary=read("final_summary.csv")
central={r["state"]:r for r in summary if r["scenario"]=="eom_inputs"}
rounded={r["state"]:r for r in summary if r["scenario"]=="rounded_inputs"}
pure={r["state"]:r for r in summary if r["scenario"]=="theta90_refit"}
accepted=read("accepted_points.csv")
updates=read("old62_controlled_summary.csv")
sens=read("model_and_input_sensitivity.csv")
hierarchy=read("three_point_convergence_summary.csv")
angles=read("angle_endpoints.csv")
unc=read("partial_uncertainty_summary.csv")
valid=json.loads((OUT/"independent_validation.json").read_text())
main_table=table(["State","$F_M(0)$",r"$\mu/\mu_N$","$F_D(0)$",r"$\mathcal D/e$ [fm$^2$]"],
 [[label,f(central[s],"FM_qeps_median"),f(central[s],"mu_qeps_median"),f(central[s],"FD_clean_median"),f(central[s],"D_fm2_median")] for s,label in zip(STATES,LABELS)])
window_table=table(["State","Points",r"$M^2$ [GeV$^2$]",r"$s_0$ [GeV$^2$]","Median $f_H$ [GeV]"],
 [[label,central[s]["n_points"],f"{min(float(r['MB2_GeV2']) for r in accepted if r['state']==s):.2f}–{max(float(r['MB2_GeV2']) for r in accepted if r['state']==s):.2f}",
   f"{min(float(r['s0_GeV2']) for r in accepted if r['state']==s):.4f}–{max(float(r['s0_GeV2']) for r in accepted if r['state']==s):.4f}",f(central[s],"f_GeV_median")] for s,label in zip(STATES,LABELS)])
spread_table=table(["State",r"$\mu/\mu_N$: 16–84%",r"$F_D$: 16–84%"],[[label,ci(central[s],"mu_qeps"),ci(central[s],"FD_clean")] for s,label in zip(STATES,LABELS)])

text=r"""---
title: "Recalculated electromagnetic multipole moments"
subtitle: "Independent derivations, numerical validation and revised analysis"
date: "25 September 2026"
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

# Current results

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

"""+main_table+r"""
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

"""+spread_table+r"""
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

"""+window_table+r"""
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

"""+table(["State",r"$|F_M^{\rm tw4}/F_M|$",r"$|F_D^{\rm tw4}/F_D|$"],
 [[label,f(next(r for r in hierarchy if r["state"]==s),"FM_tw4_over_total_abs_median",3),
   f(next(r for r in hierarchy if r["state"]==s),"FD_tw4_over_total_abs_median",3)] for s,label in zip(STATES,LABELS)])+r"""
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

"""
names=[("items1_6_and_repaired_residue","Repaired blocks and residue"),("continuum","Uniform positive-power map"),("contact","Add link response"),("em","Add EM three-particle"),("da_shapes","Correct DA shapes"),("rounded_inputs","Change kappa inputs"),("eom_inputs","Enforce zeta constraint")]
for metric,title in (("mu_epsq_median",r"Old $\epsilon q$ magnetic route, $\mu/\mu_N$"),("FD_qq_QH_median",r"Old $qq$ route, $F_D$")):
    text+="\n### "+title+"\n\n"+table(["Step"]+list(LABELS),[[label]+[f(next(r for r in updates if r["state"]==s and r["scenario"]==sc),metric,5) for s in STATES] for sc,label in names])
text+=r"""
The full file `old62_controlled_updates.csv` also gives the clean routes at
every step. Switching to clean projections and then refitting the windows are
separate operations; a median shift must not be applied to another median.

## Rounded-input comparison and consistent-input sensitivity

For reference, our own calculation with the rounded $\zeta_2=0.3$ gives:

"""+table(["State",r"$\mu/\mu_N$",r"$F_D$"],[[label,f(rounded[s],"mu_qeps_median"),f(rounded[s],"FD_clean_median")] for s,label in zip(STATES,LABELS)])+r"""
These reproduce the initially relayed rounded central results. That agreement
is not used to validate the calculation. Our central table instead satisfies
the exact DA constraint.

The next variations keep the accepted nominal points fixed and recompute all
affected OPE coefficients and residues. They are model/input sensitivities,
not independent one-sigma errors.

"""
variants=[("central","Central"),("constraint_via_zeta1","Constrain zeta1 instead"),("legacy_kappa","Old kappa pair"),("hard_subtraction_only","Old DA subtraction")]
for metric,title in (("mu_qeps_median",r"Clean $\mu/\mu_N$"),("FD_clean_median",r"Clean $F_D$")):
    text+="\n### "+title+"\n\n"+table(["Variation"]+list(LABELS),[[label]+[f(next(r for r in sens if r["state"]==s and r["scenario"]==sc),metric,5) for s in STATES] for sc,label in variants])
text+=r"""
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

"""+table(["State","Refitted points",r"$\mu/\mu_N$ at $90^\circ$",r"$F_D$ at $90^\circ$"],[[label,pure[s]["n_points"],f(pure[s],"mu_qeps_median"),f(pure[s],"FD_clean_median")] for s,label in zip(STATES,LABELS)])+r"""
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

"""
uncrows=[]
for s,label in zip(STATES,LABELS):
    a=next(r for r in unc if r["state"]==s and r["model"]=="photon_angle_and_window" and r["metric"]=="mu_muN")
    b=next(r for r in unc if r["state"]==s and r["model"]=="photon_angle_and_window" and r["metric"]=="FD")
    uncrows.append([label,f"{f(a,'median',4)} [{f(a,'p16',4)}, {f(a,'p84',4)}]",f"{f(b,'median',4)} [{f(b,'p16',4)}, {f(b,'p84',4)}]"])
text+=table(["State",r"Sampled $\mu/\mu_N$ [16–84%]",r"Sampled $F_D$ [16–84%]"],uncrows)+r"""
The sampling medians need not equal the deterministic central values because
the input distribution and ratios are nonlinear. Results for window-only,
photon-plus-window and angle-plus-window sampling are recorded separately.

# Independent numerical validation

"""+table(["Check","Maximum residual"],[
 ["Direct spectral integral versus positive-power map",f"{valid['spectral_map_max_absolute']:.3e}"],
 ["Link charge coefficient versus repaired condensate matrix",f"{valid['contact_two_point_max_absolute']:.3e}"],
 ["Charge at randomized masses, thresholds and angles",f"{valid['charge_max_absolute']:.3e}"],
 ["Two clean magnetic projections",f"{valid['clean_magnetic_routes_max_absolute']:.3e}"],
 ["Two-point Borel derivative, relative",f"{valid['two_point_derivative_max_relative']:.3e}"],
 ["All observables, 256 versus 512 quadrature nodes",f"{valid['quadrature_observable_max_absolute']:.3e}"]])+r"""
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
"""
text=text.replace("# Pole separation and observables",(OUT/"susceptibility_charge_audit/report_section.md").read_text()+"\n# Pole separation and observables",1)
text+="\n\\clearpage\n\n"+(OUT/"structure_tables/report_appendix.md").read_text()
(OUT/"RECALCULATION_REPORT.md").write_text(text)
(OUT/"central_table.md").write_text(main_table+"\nIntervals and assumptions: RECALCULATION_REPORT.md.\n")
latex="\\begin{tabular}{lrrrr}\n\\toprule\nState & $F_M$ & $\\mu/\\mu_N$ & $F_D$ & $\\mathcal D/e$ [fm$^2$] \\\\\n\\midrule\n"
for s,label in zip(STATES,LABELS):latex+=" & ".join([label,f(central[s],"FM_qeps_median"),f(central[s],"mu_qeps_median"),f(central[s],"FD_clean_median"),f(central[s],"D_fm2_median")])+" \\\\\n"
latex+="\\bottomrule\n\\end{tabular}\n"
(OUT/"central_table.tex").write_text(latex)

# Preserve the lattice update whenever the baseline report is regenerated.
if (OUT/'blind_lattice_twist_split/sha256_manifest.json').exists():
    from update_lattice_documents import update
    update()

reply=r"""# Independent recalculation and two methodological clarifications

We have completed our own recalculation. We have not read your private source
directory or used your numerical values as calculation inputs. The raw trace,
Borel, DA-weight and two-point calculations are ours. Our central outputs are
frozen in `../recalculation/freeze_sha256.json`.

The repaired calculation passes independent charge identities (also separated
by quark line), the equality of two clean magnetic projections, direct spectral
integration checks and quadrature refinement. The full nine-tensor link response
was derived with explicit Dirac matrices and reconstructed component by component.

Two points merit discussion before declaring the methodology settled:

1. **The rounded twist-4 inputs are not exactly EOM consistent.** Our direct
   integration of the first BBK (4.41) relation gives

   $$h_\gamma^{\rm derived}-h_\gamma=-40u(u-1)[1+5u(u-1)]
   (\zeta_1+11\zeta_2-2\zeta_2^+-7/2).$$

   With $(0.4,0.3,0)$ the residual is $-0.5$ at $u=1/2$.
   We hold $\zeta_1=0.4$ and use $\zeta_2=31/110$ in our central result.
   The alternative of holding $\zeta_2=0.3$ and using $\zeta_1=0.2$ is a
   separate consistent sensitivity run. Please check this residual independently.

2. **The continuum rule should not be described as an exact universal BBKR
   result.** The positive-power Laplace integral is directly verified. BBKR
   Appendix B treats unsubtracted zero/negative higher-twist powers as an
   approximation; it does not prove all their double spectral densities are
   diagonal. Charge closure is necessary but cannot settle the transverse
   continuum model: our alternative DA-subtraction diagnostic retains exact
   charge while shifting the moments. We use the uniform positive-power map as
   the central model and keep that limitation explicit.

Our finite-conformal-order EM DAs use the same order as the gluonic set. An
exact-local EM replacement alone is not a consistent higher-order uncertainty
estimate without the matching changes to the remaining DA relations.

The original nominal-window criteria give 37,25,57,17 accepted points in our
two-point calculation. Independent pure-current window selection at 90 degrees
gives 55,79,64,30 points. We provide both frozen-window angle curves and refitted
pure-current results.

Our EOM-constrained central values are:

"""+main_table+r"""
The separate rounded-input run reproduces the initially relayed central values.
That comparison was not used as a validation criterion. Please use our frozen
pointwise results and sector table for any further black-box comparison, and
specify whether the rounded or EOM-constrained input set is being tested.

All derivations, assumptions and revised tables are in
`../recalculation/RECALCULATION_REPORT.md`. In particular, window percentiles
are kept separate from the partial photon/angle Monte Carlo intervals. We have
not assigned unsupported uncertainties to missing OPE or conformal terms.
"""
if "--write-reply" in sys.argv:
    (ROOT/"from_claude/codex_reply_4.md").write_text(reply)
print("Wrote independent calculation report and table exports")
