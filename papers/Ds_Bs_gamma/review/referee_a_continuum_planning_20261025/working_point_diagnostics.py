"""Read the preserved scan and distinguish OPE composition from hard retention.

No widths, ensembles, input choices or continuum prescriptions are changed.
The tail integral is only a diagnostic of the current diagonal hard density.
"""
from __future__ import annotations
import csv
import json
import math
import sys
from functools import lru_cache
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'review/referee_a_borel_asymmetry_20261025'
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'scripts'))
from stage1_axial_g1_baseline import rho_p_colangelo_g1
from stage1_tensor_gb_hard_candidate import rho_tensor_hard_candidate


def write_csv(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def infinite_hard(K,inp,n=240,cutoff=80.):
    # s=sth+K*t^2 removes the threshold square root; the tail is exponentially
    # negligible at t=sqrt(cutoff). Recheck with n=400 and cutoff=100.
    x,w=np.polynomial.legendre.leggauss(n)
    t=math.sqrt(cutoff)*(x+1)/2
    w=w*math.sqrt(cutoff)/2
    threshold=(inp['mc']+inp['ms'])**2
    ss=threshold+K[:,None]*t[None,:]**2
    weight=2*K[:,None]*t[None,:]*w[None,:]*np.exp(-ss/K[:,None])
    return np.array([np.sum(weight*f(ss,inp['mc'],inp['ms'],inp['ec'],inp['es']),axis=1)
        for f in (rho_p_colangelo_g1,rho_tensor_hard_candidate)])


@lru_cache(None)
def nodes(n):
    return np.polynomial.legendre.leggauss(n)


def density_roots(func,inp,upper):
    threshold=(inp['mc']+inp['ms'])**2
    grid=np.linspace(threshold+1e-9,upper,10001)
    def f(s): return func(s,inp['mc'],inp['ms'],inp['ec'],inp['es'])
    vals=f(grid);roots=[]
    for j in np.flatnonzero(vals[1:]*vals[:-1]<0):
        a,b=grid[j],grid[j+1];fa=f(a)
        for _ in range(60):
            m=(a+b)/2;fm=f(m)
            if fa*fm<=0:b=m
            else:a=m;fa=fm
        roots.append((a+b)/2)
    return roots


def absolute_hard(K,upper,inp,func,roots,n=96):
    threshold=(inp['mc']+inp['ms'])**2
    bounds=[threshold]+[x for x in roots if threshold<x<upper]+[upper]
    x,w=nodes(n);answer=0.
    for a,b in zip(bounds[:-1],bounds[1:]):
        lo=math.sqrt(max(0.,(a-threshold)/K));hi=math.sqrt((b-threshold)/K)
        t=(hi-lo)*(x+1)/2+lo
        s=threshold+K*t*t
        rho=func(s,inp['mc'],inp['ms'],inp['ec'],inp['es'])
        answer+=(hi-lo)/2*np.sum(w*2*K*t*np.exp(-s/K)*abs(rho))
    return float(answer)


def main():
    data=json.loads((SOURCE/'inputs.json').read_text())
    points=list(csv.DictReader((SOURCE/'benchmark_midpoints.csv').open()))
    ratios=[]; retention=[]; ranges=[]; errors=[]; absolute=[]; absolute_errors=[]
    for row in points:
        output={k:row[k] for k in ['state','mode','r','M2','s0','M1_sq','M2_sq']}
        for base,prefix in [('A','PiA_'),('B','PiB_'),('physical','Pi_')]:
            pieces={k[len(prefix):]:float(v) for k,v in row.items()
                    if k.startswith(prefix) and k!= 'Pi_physical'}
            total=sum(pieces.values())
            output[base+'_pert']=pieces['pert'];output[base+'_total']=total
            output[base+'_pert_over_total']=pieces['pert']/total
            output[base+'_abs_pert_over_sum_abs']=abs(pieces['pert'])/sum(map(abs,pieces.values()))
            tw4=[v for k,v in pieces.items() if k.startswith('tw4_')]
            output[base+'_tw4_signed_over_total']=sum(tw4)/total
            output[base+'_tw4_abs_sum_over_abs_total']=sum(map(abs,tw4))/abs(total)
        ratios.append(output)
    for spec in data['specifications']:
        state=spec['id'];inp=data['QCD_inputs'][state]
        arr=np.load(SOURCE/(state+'_scan.npz'))
        K=arr['M2'];s0=arr['s0'];comp=arr['components']
        j=int(np.argmin(abs(arr['r']-1)))
        hard=comp[0,:,j,:,:]
        infinity=infinite_hard(K,inp)
        control=infinite_hard(K,inp,n=400,cutoff=100.)
        errors.append(float(np.max(abs(infinity-control))))
        idx=len(K)//2;idy=len(s0)//2
        theta=math.radians(spec['theta'])
        wa,wb=(math.sin(theta),math.cos(theta)) if spec['i']==1 else (math.cos(theta),-math.sin(theta))
        rA=hard[0]/infinity[0,:,None];rB=hard[1]/infinity[1,:,None]
        rp=(wa*hard[0]+wb*hard[1])/(wa*infinity[0,:,None]+wb*infinity[1,:,None])
        retention.append(dict(state=state,M2=float(K[idx]),s0=float(s0[idy]),
            A_hard_sub=float(hard[0,idx,idy]),A_hard_infinite=float(infinity[0,idx]),
            R_A=float(rA[idx,idy]),B_hard_sub=float(hard[1,idx,idy]),
            B_hard_infinite=float(infinity[1,idx]),R_B=float(rB[idx,idy]),
            R_physical=float(rp[idx,idy]),
            soft_factor_retention=1-math.exp(-(float(s0[idy])-inp['mc']**2)/float(K[idx]))))
        ranges.append(dict(state=state,M2_min=float(K.min()),M2_max=float(K.max()),
            s0_min=float(s0.min()),s0_max=float(s0.max()),
            R_A_min=float(rA.min()),R_A_max=float(rA.max()),
            R_B_min=float(rB.min()),R_B_max=float(rB.max()),
            R_physical_min=float(rp.min()),R_physical_max=float(rp.max())))
        limits=(inp['mc']+inp['ms'])**2+100*K.max()
        funcs=(rho_p_colangelo_g1,rho_tensor_hard_candidate)
        roots=[density_roots(f,inp,limits) for f in funcs]
        for tag,kk,ss in [('midpoint',float(K[idx]),float(s0[idy])),
                          ('minimum',float(K.max()),float(s0.min())),
                          ('maximum',float(K.min()),float(s0.max()))]:
            record=dict(state=state,point=tag,M2=kk,s0=ss)
            upper=(inp['mc']+inp['ms'])**2+80*kk
            for label,func,zeros in zip(['A','B'],funcs,roots):
                num=absolute_hard(kk,ss,inp,func,zeros)
                den=absolute_hard(kk,upper,inp,func,zeros)
                ctrl=absolute_hard(kk,ss,inp,func,zeros,n=160)/absolute_hard(kk,upper,inp,func,zeros,n=160)
                record['abs_retention_'+label]=num/den
                absolute_errors.append(abs(num/den-ctrl))
            absolute.append(record)
    write_csv('perturbative_to_total.csv',ratios)
    write_csv('hard_retention_midpoints.csv',retention)
    write_csv('hard_retention_windows.csv',ranges)
    write_csv('absolute_hard_retention.csv',absolute)
    (HERE/'numerical_check.json').write_text(json.dumps(dict(
        max_infinite_integral_control_difference=max(errors),
        max_absolute_retention_quadrature_change=max(absolute_errors),
        numerator='previously saved n=1600 diagonal hard integral',
        denominator='same density integrated to infinity, threshold square-root substitution',
        meaning='signed hard-sector retention only; not a positive full-transition pole probability',
        new_widths_or_ensembles_computed=False),indent=2)+'\n')
    for r in retention: print(r['state'],r['R_A'],r['R_B'],r['R_physical'])
    print('Max integral check difference:',max(errors))
    for row in absolute: print('Absolute hard retention:',row)


if __name__=='__main__': main()
