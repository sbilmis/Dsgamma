"""Make author-review tables and plots from the isolated scan only."""
from __future__ import annotations
import csv
import os
from pathlib import Path

HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.cache/matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME',str(HERE/'.cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import numpy as np


def read(name):
    with (HERE/name).open() as f:
        return list(csv.DictReader(f))


def main():
    rows=read('benchmark_midpoints.csv')
    by={(r['state'],r['mode']):r for r in rows}
    labels={'Ds2460':r'D_{s1}(2460)', 'Ds2536':r'D_{s1}(2536)',
            'Bs5750':r'B_{s1}(5750)', 'Bs5830':r'B_{s1}(5830)'}
    lines=[r'% AUTHOR REVIEW: conditional model-preserving diagnostic, not final referee wording.',
           r'% The off-symmetric manuscript/derivation convention difference remains to be reconciled.',
           r'\begin{table}[ht]',r'\centering',r'\begin{tabular}{lrrrrr}',r'\hline',
           r'Channel & $r_{\rm mass}$ & $g_i(1)$ & $g_i(r_{\rm mass})$ & $\Gamma_i(1)$ & $\Gamma_i(r_{\rm mass})$ \\',
           r' & & \multicolumn{2}{c}{$[\mathrm{GeV}^{-1}]$} & \multicolumn{2}{c}{$[\mathrm{keV}]$} \\',r'\hline']
    for state,label in labels.items():
        a=by[state,'symmetric']; b=by[state,'mass_ratio']
        vals=[float(b['r']),float(a['g_signed']),float(b['g_signed']),
              float(a['Gamma_keV']),float(b['Gamma_keV'])]
        lines.append('$'+label+r'\to '+('D_s' if state.startswith('Ds') else 'B_s')+r'\gamma$ & '+
            f'{vals[0]:.4f} & {vals[1]:+.5f} & {vals[2]:+.5f} & {vals[3]:.3f} & {vals[4]:.3f} '+r'\\')
    lines += [r'\hline',r'\end{tabular}',
        r'\caption{Fixed-input unequal-Borel diagnostic at the working-window midpoints.',
        r'The ratio is $r=M_1^2/M_2^2$ and $r_{\rm mass}=m_{A_i}^2/m_P^2$,',
        r'with the effective $M^2=M_1^2M_2^2/(M_1^2+M_2^2)$ held fixed.',
        r'The existing diagonal hard spectral model and soft-sector continuum prescription are retained.',
        r'The soft fractions follow the stored off-shell derivation, as documented in the calculation notes.',
        r'These point values are not the manuscript Monte Carlo medians or new uncertainty intervals.}',
        r'\label{tab:borel-diagnostic-author-review}',r'\end{table}']
    (HERE/'borel_comparison_author_review.tex').write_text('\n'.join(lines)+'\n')

    plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans',
        'font.size':11,'axes.labelsize':12,'axes.titlesize':12,'axes.linewidth':1.1,
        'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
        'xtick.minor.visible':True,'ytick.minor.visible':True,'svg.fonttype':'none',
        'savefig.bbox':'tight'})
    curves=read('ratio_scan.csv')
    fig,axs=plt.subplots(1,2,figsize=(9.2,3.5))
    for ax,state,title in zip(axs,['Ds2536','Bs5830'],
        [r'(a) $D_{s1}(2536)\to D_s\gamma$',r'(b) $B_{s1}(5830)\to B_s\gamma$']):
        rr=[r for r in curves if r['state']==state]
        x=np.array([float(r['r']) for r in rr])
        y=np.array([float(r['g_midpoint']) for r in rr])
        lo=np.array([float(r['g_min']) for r in rr]); hi=np.array([float(r['g_max']) for r in rr])
        rm=float(rr[0]['r_mass'])
        ax.fill_between(x,lo,hi,color='#0047ab',alpha=.15,lw=0)
        ax.plot(x,y,color='#0047ab',lw=2)
        ax.axvline(1.,color='#777777',ls=':',lw=1.5)
        ax.axvline(rm,color='#b22222',ls='-.',lw=1.5)
        ax.axhline(0,color='black',lw=.7)
        pad=.06*(hi.max()-lo.min())
        ax.set_ylim(lo.min()-pad,hi.max()+pad)
        ax.set_xlim(x.min(),x.max()); ax.set_title(title,pad=10)
        ax.set_xlabel(r'$r=M_1^2/M_2^2$'); ax.set_ylabel(r'$g_2\ [\mathrm{GeV}^{-1}]$')
        ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(which='major',length=5); ax.tick_params(which='minor',length=3)
        ax.text(.97,.04,rf'$r_{{\rm mass}}={rm:.3f}$',ha='right',va='bottom',
                transform=ax.transAxes,color='#b22222',fontsize=10)
    handles=[Line2D([0],[0],color='#0047ab',lw=2,label='Window midpoint'),
        Line2D([0],[0],color='#777777',lw=1.5,ls=':',label=r'$r=1$'),
        Line2D([0],[0],color='#b22222',lw=1.5,ls='-.',label=r'$r=m_A^2/m_P^2$')]
    fig.legend(handles=handles,loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.5,1.04))
    fig.subplots_adjust(left=.085,right=.99,bottom=.18,top=.79,wspace=.29)
    for ext in ['png','svg']:
        fig.savefig(HERE/f'borel_ratio_signed_g2.{ext}',dpi=220)
    plt.close(fig)
    (HERE/'borel_ratio_plot_caption.tex').write_text(
        '% AUTHOR REVIEW ONLY; read README.md before using.\n'
        r'\caption{Signed higher-state couplings as functions of $r=M_1^2/M_2^2$ in the '
        r'fixed-input unequal-Borel diagnostic. Solid curves use the Table~II window midpoints; '
        r'shading spans the sampled effective-$M^2$ and $s_0$ regions and is not a confidence interval. '
        r'The dotted and dash-dotted lines mark the symmetric and mass-based prescriptions. '
        r'The retained diagonal hard spectral model, soft continuum prescription, and off-shell '
        r'fraction convention are specified in the accompanying calculation notes.}'+'\n')
    print('Wrote author-review LaTeX table, PNG/SVG plot, and caption.')


if __name__=='__main__': main()
