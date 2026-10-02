"""Plot frozen signed contributions versus M^2 without mixing s0 trajectories."""
import os
from pathlib import Path
import csv
import hashlib
import json
from collections import defaultdict
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/dsgamma_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'recalculation/blind_lattice_twist_split'
DEST=ROOT.parent/'figures/lattice_twist_borel'
PDF=ROOT.parent/'output/pdf/lattice_twist_contributions_vs_M2.pdf'
STATES=['Ds_1(2460)','Ds_1(2536)','Bs_1(5750)','Bs_1(5830)']
LABELS=[r'$D_{s1}(2460)^+$',r'$D_{s1}(2536)^+$',r'$B_{s1}(5750)$',r'$B_{s1}(5830)$']
MASSES=[2.4595,2.53511,5.7502,5.82911]
PARTS=[('hard','Perturbative','#0072B2'),('tw2','Twist 2','#D55E00'),('tw3','Twist 3','#009E73'),('tw4','Twist 4','#CC79A7'),('contact','Contact','#E69F00'),('total','Total','#111111')]
manifest=json.loads((DATA/'sha256_manifest.json').read_text())
for name,digest in manifest['outputs'].items():
    assert hashlib.sha256((DATA/name).read_bytes()).hexdigest()==digest
rows=list(csv.DictReader((DATA/'twist_split_windows.csv').open()))
for r in rows:
    for k in r:
        if k not in ('angle','state','twist2_sign'):r[k]=float(r[k])
DEST.mkdir(parents=True,exist_ok=True)
PDF.parent.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.labelsize':11,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
files=[]
page_records=[]
with PdfPages(PDF,metadata={'Title':'Lattice-normalized twist contributions versus Borel M squared','Subject':'Frozen accepted windows; both twist-2 signs; nominal and pure-current angles'}) as pdf:
    for angle in ['nominal','theta90']:
        for metric in ['mu_over_muN','FD']:
            for sign in ['+1','-1']:
                fig,axes=plt.subplots(2,2,figsize=(11.7,8.3))
                fig.subplots_adjust(left=.08,right=.98,bottom=.18,top=.82,hspace=.42,wspace=.25)
                angle_text='Nominal mixing angles' if angle=='nominal' else r'$\theta=90^\circ$: refitted pure-current windows'
                metric_text=r'Magnetic moment $\mu/\mu_N$' if metric=='mu_over_muN' else r'Quadrupole form factor $F_D-Q_H$'
                sign_text=r'$+\chi$: $f_\perp=-55.1$ MeV' if sign=='+1' else r'$-\chi$: $f_\perp=+55.1$ MeV'
                fig.suptitle(metric_text+' vs '+r'$M^2$',fontsize=17,y=.975)
                fig.text(.5,.93,angle_text+'  |  '+sign_text,ha='center',fontsize=12)
                highlight_records=[]
                for ax,state,label,mass in zip(axes.flat,STATES,LABELS,MASSES):
                    data=[r for r in rows if r['angle']==angle and r['state']==state and r['twist2_sign']==sign]
                    both=[r for r in rows if r['angle']==angle and r['state']==state]
                    trajectories=defaultdict(list)
                    for r in data:trajectories[r['s0_GeV2']].append(r)
                    highlight=min(trajectories,key=lambda s0:abs(np.sqrt(s0)-mass-.5))
                    delta=np.sqrt(highlight)-mass
                    highlight_records.append(dict(state=state,s0_GeV2=highlight,delta_GeV=float(delta),n_windows=len(data)))
                    def column(part):
                        return (metric+'_'+part) if part!='total' else ('mu_over_muN_total' if metric=='mu_over_muN' else 'FD_without_QH_total')
                    for s0,trajectory in sorted(trajectories.items(),key=lambda item:item[0]==highlight):
                        trajectory.sort(key=lambda r:r['M2_GeV2'])
                        x=[r['M2_GeV2'] for r in trajectory]
                        bold=s0==highlight
                        for part,name,color in PARTS:
                            y=[r[column(part)] for r in trajectory]
                            # Retain isolated accepted windows as points.
                            ax.plot(x,y,color=color,lw=(2.0 if bold else .8),alpha=(1. if bold else .23),
                                    linestyle=('--' if part=='total' else '-'),marker='o' if bold or len(x)==1 else None,
                                    markersize=2.5 if bold else 2,zorder=4 if bold else 2)
                    vals=np.array([r[column(part)] for r in both for part,_,_ in PARTS])
                    lo=min(0.,vals.min());hi=max(0.,vals.max());pad=.07*(hi-lo if hi>lo else 1.)
                    ax.set_ylim(lo-pad,hi+pad)
                    all_x=[r['M2_GeV2'] for r in both]
                    ax.set_xlim(min(all_x)-.04*(max(all_x)-min(all_x)),max(all_x)+.04*(max(all_x)-min(all_x)))
                    ax.axhline(0,color='#999999',lw=.6,zorder=1)
                    ax.grid(True,color='#dddddd',lw=.5,alpha=.6)
                    ax.set_title(label+f'  |  {len(data)} windows\n'+r'Bold: $s_0$'+f' = {highlight:.3f} GeV'+r'$^2$',fontsize=11)
                    ax.set_xlabel(r'$M^2$ [GeV$^2$]')
                    ax.set_ylabel(r'$\mu/\mu_N$' if metric=='mu_over_muN' else r'$F_D-Q_H$')
                handles=[Line2D([0],[0],color=c,lw=2,linestyle='--' if p=='total' else '-',label=n) for p,n,c in PARTS]
                fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.09),ncol=6,frameon=False,fontsize=10)
                fig.text(.08,.065,r'Thick curves: $\sqrt{s_0}-m_H=0.50$ GeV (or nearest accepted threshold). Faint curves: other accepted $s_0$ values.',fontsize=9)
                fig.text(.08,.04,'Only accepted windows are shown; points at different thresholds are never joined. Each point uses its own two-point residue.',fontsize=9)
                fig.text(.08,.015,'Frozen CSV: commit 54e6958. '+('Quadrupole total excludes charge: add Q_H=1 for Ds, 0 for Bs.' if metric=='FD' else 'Both sign pages use identical axis limits for each state.'),fontsize=9)
                pdf.savefig(fig)
                path=DEST/f'{angle}_{metric}_{"plus" if sign=="+1" else "minus"}.png'
                fig.savefig(path,dpi=160)
                files.append(path)
                plt.close(fig)
                page_records.append(dict(angle=angle,metric=metric,twist2_sign=sign,highlighted_thresholds=highlight_records))
files.extend([PDF,Path(__file__).resolve()])
(DEST/'plot_setup.json').write_text(json.dumps(dict(source_sha256=manifest['outputs']['twist_split_windows.csv'],pages=page_records,
    zero_charge_added_to_FD_parts=False,quadrupole_total='FD_without_QH_total',axis_limits='shared between signs for each state and angle; include every accepted threshold'),indent=2)+'\n')
files.append(DEST/'plot_setup.json')
(DEST/'sha256_manifest.json').write_text(json.dumps({str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
print('Created 8-page PDF and 8 PNG figures:',PDF)
