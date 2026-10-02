"""Render the committed lattice summary into current documentation."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'recalculation'
FROZEN=OUT/'blind_lattice_twist_split'
LABELS={'Ds_1(2460)':r'$D_{s1}(2460)^+$','Ds_1(2536)':r'$D_{s1}(2536)^+$','Bs_1(5750)':r'$B_{s1}(5750)$','Bs_1(5830)':r'$B_{s1}(5830)$'}
BEGIN='<!-- BEGIN LATTICE RESULTS -->'
END='<!-- END LATTICE RESULTS -->'

def table(headers,rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)+'\n'

def fmt(r,k,n=6):return f'{float(r[k]):.{n}f}'

def update():
    manifest=json.loads((FROZEN/'sha256_manifest.json').read_text())
    for name,digest in manifest['outputs'].items():
        assert hashlib.sha256((FROZEN/name).read_bytes()).hexdigest()==digest
    rows=list(csv.DictReader((FROZEN/'twist_split_summary.csv').open()))
    intro=r'''# Lattice-normalized results, both twist-2 signs

Updated 3 October 2026 from the blind output frozen in Git commit `54e6958`.
These tables give the lattice-normalized run separately for both signs; the
remaining report records the earlier chi=3.15 baseline and its diagnostics.
The lattice tables supersede that baseline for the requested normalization.
No counterpart results have been used, and an entry-by-entry comparison is pending.

The supplied strange-photon normalization at 1 GeV is
$f_{\perp\gamma,s}=-55.1$ MeV. With
$\langle\bar ss\rangle=0.8(-0.240\,\mathrm{GeV})^3$,
$\chi=4.982289497\,\mathrm{GeV}^{-2}$.
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

'''
    # Use the exact value from the frozen input rather than a copied decimal.
    chi=-.0551/(.8*(-.240)**3)
    intro=intro.replace('4.982289497',f'{chi:.9f}')
    content=intro
    latex=[]
    central=[]
    for angle,title in [('nominal','Nominal mixing angles'),('theta90','90 degrees, refitted pure-current windows')]:
        data=[r for r in rows if r['angle']==angle]
        totals=table(['State','Sign','Windows',r'$\mu/\mu_N$',r'$F_D$',r'$\mathcal D/e$ [fm$^2$]'],[
            [LABELS[r['state']],r'$+\chi$' if r['twist2_sign']=='+1' else r'$-\chi$',r['n_windows'],fmt(r,'mu_over_muN_median'),fmt(r,'FD_median'),fmt(r,'D_efm2_median')] for r in data])
        content+='## '+title+'\n\n'+totals+'\n'
        central.append('## '+title+'\n\n'+totals)
        for metric,label in [('mu_over_muN',r'$\mu/\mu_N$'),('FD',r'$F_D$'),('D_efm2',r'$\mathcal D/e$ [fm$^2$]')]:
            content+='### '+label+' window intervals\n\n'+table(['State','Sign','Median','16-84%'],[
                [LABELS[r['state']],r'$+\chi$' if r['twist2_sign']=='+1' else r'$-\chi$',fmt(r,metric+'_median'),f"[{fmt(r,metric+'_p16')}, {fmt(r,metric+'_p84')}]" ] for r in data])+'\n'
        for metric,label in [('mu_over_muN',r'$\mu/\mu_N$'),('FD',r'$F_D-Q_H$')]:
            content+='### '+label+' contribution medians\n\n'+table(['State','Sign','Hard','Twist 2','Twist 3','Twist 4','Contact'],[
                [LABELS[r['state']],r'$+\chi$' if r['twist2_sign']=='+1' else r'$-\chi$']+[fmt(r,metric+'_'+g+'_median') for g in ['hard','tw2','tw3','tw4','contact']] for r in data])+'\n'
        latex.append('% '+title+'\n'+r'\begin{tabular}{llrrrr}'+'\n'+r'\toprule'+'\n'+r'State & Sign & Windows & $\mu/\mu_N$ & $F_D$ & $\mathcal D/e$ [fm$^2$] \\'+'\n'+r'\midrule'+'\n'+'\n'.join(' & '.join([LABELS[r['state']],r'$+\chi$' if r['twist2_sign']=='+1' else r'$-\chi$',r['n_windows'],fmt(r,'mu_over_muN_median'),fmt(r,'FD_median'),fmt(r,'D_efm2_median')])+r' \\' for r in data)+'\n'+r'\bottomrule'+'\n'+r'\end{tabular}'+'\n')
    content+=r'''## Grouping and validation

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
'''
    (OUT/'lattice_results.md').write_text(content)
    (OUT/'central_table.md').write_text(intro+'\n'.join(central)+'\nFull intervals and contribution medians: lattice_results.md.\n')
    (OUT/'central_table.tex').write_text('\n\n'.join(latex))
    report=OUT/'RECALCULATION_REPORT.md'
    text=report.read_text()
    if BEGIN in text:
        a=text.index(BEGIN);b=text.index(END,a)+len(END)
        tail=text[b:]
        if tail.startswith('\n\n\\clearpage\n\n'):
            tail=tail[len('\n\n\\clearpage\n\n'):]
        text=text[:a]+tail
    text=text.replace('date: "25 September 2026"','date: "3 October 2026"')
    text=text.replace('# Current results','# Earlier baseline: chi=3.15',1)
    pos=text.index('# Earlier baseline: chi=3.15')
    text=text[:pos]+BEGIN+'\n'+content+'\n'+END+'\n\n\\clearpage\n\n'+text[pos:]
    report.write_text(text)
    for path in [ROOT/'recalculation/README.md',ROOT.parent/'README.md']:
        text=path.read_text()
        marker='<!-- LATTICE UPDATE -->'
        if marker not in text:
            link='lattice_results.md' if path.parent==OUT else 'black_box_comparison/recalculation/lattice_results.md'
            text=marker+'\n**3 October 2026:** [Lattice-normalized results with both twist-2 signs]('+link+'). Includes nominal and refitted 90-degree windows, window intervals and separate contribution medians. Frozen outputs and their SHA-256 manifest remain unchanged.\n\n'+text
            path.write_text(text)
    print('Updated lattice tables, central tables, report and both READMEs.')

if __name__=='__main__':update()
