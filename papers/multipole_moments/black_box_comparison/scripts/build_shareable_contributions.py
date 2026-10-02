"""Build a three-page colleague handout from frozen contribution medians."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'recalculation'
FROZEN=OUT/'blind_lattice_twist_split'
PDF=ROOT.parent/'output/pdf/lattice_twist_contributions.pdf'
DOC=OUT/'lattice_twist_contributions_handout.md'

manifest=json.loads((FROZEN/'sha256_manifest.json').read_text())
for name,digest in manifest['outputs'].items():
    assert hashlib.sha256((FROZEN/name).read_bytes()).hexdigest()==digest
rows=list(csv.DictReader((FROZEN/'twist_split_summary.csv').open()))
labels={'Ds_1(2460)':r'$D_{s1}(2460)^+$','Ds_1(2536)':r'$D_{s1}(2536)^+$','Bs_1(5750)':r'$B_{s1}(5750)$','Bs_1(5830)':r'$B_{s1}(5830)$'}
text=r'''---
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

'''
for angle,title,count in [('nominal','Nominal mixing angles','37, 25, 57, 17'),('theta90','90 degrees: refitted pure-current windows','55, 79, 64, 30')]:
    text+='# '+title+'\n\nAccepted-window counts in state order: '+count+'. '
    if angle=='theta90':
        text+=r'The lighter state uses $J_A$ and the heavier state uses $J_B$.'
    text+='\n\nEach entry is the median of the indicated contribution. Twist-2 columns give both signs explicitly; the other columns apply to both signs.\n\n'
    for prefix,title2 in [('mu_over_muN',r'Magnetic moment: $\mu/\mu_N$'),('FD',r'Quadrupole form factor: $F_D-Q_H$')]:
        text+='## '+title2+'\n\n'
        text+='| State | Perturbative | Twist 2 (+$\\chi$) | Twist 2 (-$\\chi$) | Twist 3 | Twist 4 | Contact |\n|---|---:|---:|---:|---:|---:|---:|\n'
        for r in rows:
            if r['angle']!=angle or r['twist2_sign']!='+1':continue
            rev=next(s for s in rows if s['angle']==angle and s['state']==r['state'] and s['twist2_sign']=='-1')
            values=[float(r[prefix+'_hard_median']),float(r[prefix+'_tw2_median']),float(rev[prefix+'_tw2_median'])]+[float(r[prefix+'_'+g+'_median']) for g in ['tw3','tw4','contact']]
            text+='| '+labels[r['state']]+' | '+' | '.join(f'{v:.6f}' for v in values)+' |\n'
        text+='\n'
    text+=r'''## Grouping and normalization

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

'''
    if angle=='nominal':text+='\\clearpage\n\n'
text+=r'''## Provenance

All entries are read from our signed summary CSV frozen in Git commit
`54e6958`; no counterpart numerical outputs were used. The source definitions
and straight-link audit are in our recalculation report and
`blind_definition_scan_audit/AUDIT.md`. The original CSV SHA-256 hashes are
unchanged. A separate manifest covers this handout and its PDF.
'''
DOC.write_text(text)
PDF.parent.mkdir(parents=True,exist_ok=True)
subprocess.run(['pandoc',str(DOC),'--from=markdown+tex_math_dollars','--pdf-engine=xelatex','--output='+str(PDF)],check=True)
files=[DOC,PDF,Path(__file__).resolve()]
freeze={str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(OUT/'lattice_twist_contributions_handout.sha256.json').write_text(json.dumps(freeze,indent=2)+'\n')
print(PDF)
