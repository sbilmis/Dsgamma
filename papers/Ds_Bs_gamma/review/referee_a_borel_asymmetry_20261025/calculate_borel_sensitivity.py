"""Isolated unequal-Borel diagnostic within the preserved production model.

The primary branch follows the stored off-shell kernels: the coefficient a
in D=mQ^2-a*pPrime^2-(1-a)*p^2 is 1/(1+r), in BOTH soft sectors.
The hard diagonal spectral model and soft continuum prescription are held
fixed. This is not a derivation of their general two-variable completion.
The manuscript's split argument is evaluated separately as a diagnostic;
it is NOT mixed into the primary results or treated as an uncertainty.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "scripts"))
from corrected_transition_central_python import (
    DS_REFERENCE, ds_inputs, bs_inputs, bs_direct_decay_constants,
)
import rohrwild_transition_exact as prod
from final_stage2_uncertainty_scan import set_photon_shape_params
from stage1_axial_g1_baseline import rho_p_colangelo_g1
from stage1_tensor_gb_hard_candidate import rho_tensor_hard_candidate

COMPONENTS = ("pert", "tw2", "tw3", "tw4_2p", "tw4_3pg", "tw4_3pem")


def write_csv(name, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (HERE / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def write_json(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def lprime(z):
    """Exact symbolic integral of the production-oriented T4gamma density."""
    return -(z-1)*(227+z*(-923+z*(2977+z*(-4523+2352*z))))/10 + 6*math.log(z)


@lru_cache(None)
def convolutions(z, n=64):
    funcs = (prod.f_g_axial, prod.f_g_sigma,
             prod.f_em_axial, prod.f_em_tensor_base)
    return np.array([prod.line_projection(f, z, n=n) for f in funcs] + [lprime(z)])


def rectangle_projection(func, z, n=160):
    """Independent integration of delta(z-b-v*g) by doing v first."""
    x,w = np.polynomial.legendre.leggauss(n)
    a = (1-z)*(x[:,None]+1)/2
    b = z*(x[None,:]+1)/2
    g = 1-a-b
    v = (z-b)/g
    return float(z*(1-z)/4*np.sum(w[:,None]*w[None,:]*func(a,b,g,v)/g))


@lru_cache(None)
def quad_nodes(n):
    return np.polynomial.legendre.leggauss(n)


def hard_grid(K, thresholds, inp, n=1600):
    x,w = quad_nodes(n)
    lo = (inp['mc']+inp['ms'])**2
    out = np.zeros((2,len(K),len(thresholds)))
    for j,s0 in enumerate(thresholds):
        s = (s0-lo)*(x+1)/2+lo
        weight = (s0-lo)*w/2
        for a,func in enumerate((rho_p_colangelo_g1,rho_tensor_hard_candidate)):
            rho = func(s,inp['mc'],inp['ms'],inp['ec'],inp['es'])
            out[a,:,j] = np.sum(np.exp(-s[None,:]/K[:,None])*rho*weight,axis=1)
    return out


def components_grid(K, thresholds, ratios, inp, hard, split=False):
    """Returns [component,basis,r,K,s0] without physical normalization."""
    a = 1/(1+np.asarray(ratios))
    conv = np.array([convolutions(float(1-v if split else v)) for v in a])
    mq,ms,ss,eq,es = (inp[k] for k in ('mc','ms','ss','ec','es'))
    d=mq+ms
    E=np.exp(-mq*mq/K[:,None])-np.exp(-thresholds[None,:]/K[:,None])
    shp=(len(a),len(K),len(thresholds))
    out=np.zeros((6,2,*shp))
    out[0,0]=hard[0][None,:,:]; out[0,1]=hard[1][None,:,:]
    tw2=es*ss*inp['chi']*prod.pda.phi_gamma(a)
    out[1,0]=tw2[:,None,None]*K[None,:,None]*E[None,:,:]
    out[1,1]=mq/d*out[1,0]
    psi=np.array([prod.bar_psi_v(float(v)) for v in a])
    deriv=np.array([prod.bar_psi_v_prime(float(v)) for v in a])
    out[2,0]=(es*inp['f3g']*mq*psi)[:,None,None]*E[None,:,:]
    out[2,1]=es*inp['f3g']/d*(mq*mq*psi[:,None,None]
        -K[None,:,None]/2*(-psi+(1-a)*deriv)[:,None,None])*E[None,:,:]
    t4=es*ss/4*(prod.pda.A_t4(a)+2*prod.pda.B_t4(a))
    out[3,0]=t4[:,None,None]*(1+mq*mq/K)[None,:,None]*E[None,:,:]
    out[3,1]=mq/d*out[3,0]
    out[4,0]=(es*ss*conv[:,0])[:,None,None]*E[None,:,:]
    out[4,1]=(es*ss*mq/d*conv[:,1])[:,None,None]*E[None,:,:]
    out[5,0]=(eq*ss*(conv[:,2]+2*conv[:,4]))[:,None,None]*E[None,:,:]
    out[5,1]=(eq*ss*mq/d*(conv[:,3]+2*conv[:,4]))[:,None,None]*E[None,:,:]
    return out


def project(components, K, ratios, spec, inp):
    t=math.radians(spec['theta'])
    weights=(math.sin(t),math.cos(t)) if spec['i']==1 else (math.cos(t),-math.sin(t))
    pi=weights[0]*components[:,0]+weights[1]*components[:,1]
    r=np.asarray(ratios)
    ma,mp=spec['mass'],inp['m_ds']
    exponent=(ma*ma+mp*mp*r)[:,None]/((1+r)[:,None]*K[None,:])
    norm=np.exp(exponent)*(inp['mc']+inp['ms'])/(ma*spec['f']*mp*mp*inp['f_ds'])
    gc=pi*norm[None,:,:,None]
    g=gc.sum(axis=0)
    photon=(ma*ma-mp*mp)/(2*ma)
    width=(1/137.036)/3*photon**3*1e6*g*g
    return pi,gc,g,width


def specifications():
    fb=bs_direct_decay_constants()
    return [
      dict(id='Ds2460',label='D_s1(2460)',sector='Ds',i=1,mass=2.4595,
           f=DS_REFERENCE['f1'],theta=26.6,Kmin=3.,Kmax=4.5,smin=8.5,smax=9.5),
      dict(id='Ds2536',label='D_s1(2536)',sector='Ds',i=2,mass=2.53512,
           f=DS_REFERENCE['f2'],theta=26.6,Kmin=3.,Kmax=4.5,smin=9.,smax=10.),
      dict(id='Bs5750',label='B_s1(5750)',sector='Bs',i=1,mass=5.750,
           f=fb[0],theta=38.5,Kmin=10.,Kmax=14.,smin=39.,smax=41.),
      dict(id='Bs5830',label='B_s1(5830)',sector='Bs',i=2,mass=5.82873,
           f=fb[1],theta=38.5,Kmin=10.,Kmax=14.,smin=40.,smax=42.),
    ]


def point_row(spec,inp,r,K,s0,comp,pi,gc,g,width,mode):
    row=dict(state=spec['id'],mode=mode,r=r,u0=1/(1+r),
             z3p=(r if mode=='split_argument_diagnostic' else 1)/(1+r),
             M2=K,s0=s0,M1_sq=(1+r)*K,M2_sq=(1+r)*K/r,
             g_signed=float(g),Gamma_keV=float(width),Pi_physical=float(pi.sum()))
    for j,name in enumerate(COMPONENTS):
        row['PiA_'+name]=float(comp[j,0]); row['PiB_'+name]=float(comp[j,1])
        row['Pi_'+name]=float(pi[j]); row['g_'+name]=float(gc[j])
    return row


def check_baseline(specs):
    """Reproduce all saved Borel points, with their original reference masses."""
    lookup={s['id']:s for s in specs}; errors=[]
    mapping={'Ds1_2460':'Ds2460','Ds1_2536':'Ds2536','Bs1_5750':'Bs5750','Bs1_5830':'Bs5830'}
    with (ROOT/'figures/figure1_borel_form_factors/figure1_borel_form_factors.csv').open() as f:
        for row in csv.DictReader(f):
            s=dict(lookup[mapping[row['state']]])
            if s['id']=='Ds2536': s['mass']=DS_REFERENCE['m2']
            if s['id']=='Bs5830': s['mass']=5.82870
            inp=ds_inputs() if s['sector']=='Ds' else bs_inputs(s['mass'])
            K=np.array([float(row['M2'])]); s0=np.array([float(row['s0'])])
            hard=hard_grid(K,s0,inp)
            c=components_grid(K,s0,[1.],inp,hard)
            _,_,g,width=project(c,K,[1.],s,inp)
            errors.append((abs(g.item()-float(row['g_signed'])),
                           abs(width.item()-float(row['Gamma_keV']))))
    return dict(points=len(errors),max_abs_g_error=max(x[0] for x in errors),
                max_abs_width_error_keV=max(x[1] for x in errors))


def numerical_checks(specs):
    projection=[]; lp=[]
    for z in [.34,.376110404998,.390425287006,.458821008661,.465581684937,.5,.54,.625,.66]:
        cv=convolutions(z)
        for j,func in enumerate((prod.f_g_axial,prod.f_g_sigma,prod.f_em_axial,prod.f_em_tensor_base)):
            independent=rectangle_projection(func,z)
            higher=prod.line_projection(func,z,n=112)
            projection.append(dict(z=z,function=func.__name__,line64=float(cv[j]),
                rectangle160=independent,line112=higher,
                delta_rectangle=float(cv[j]-independent),delta_order=float(cv[j]-higher)))
        num=prod.p_line_prime(prod.t4gamma,z,n=44,step=2e-4)
        lp.append(dict(z=z,analytic=lprime(z),numerical=num,difference=lprime(z)-num))
    write_csv('independent_convolution_checks.csv',projection)
    write_csv('analytic_derivative_checks.csv',lp)
    return dict(baseline=check_baseline(specs),
        max_line_vs_rectangle_error=max(abs(x['delta_rectangle']) for x in projection),
        max_line_quadrature_order_error=max(abs(x['delta_order']) for x in projection),
        max_Lprime_analytic_vs_finite_difference_error=max(abs(x['difference']) for x in lp))


def main():
    set_photon_shape_params(-2.1,3.8)
    specs=specifications()
    inputs={s['id']:(ds_inputs() if s['sector']=='Ds' else bs_inputs(s['mass'])) for s in specs}
    write_json('inputs.json',dict(specifications=specs,QCD_inputs=inputs,
        photon_parameters=prod.pda.PARAMS,
        scope='fixed-input, retained diagonal spectral and soft continuum model; no MC refit',
        primary_fraction='u0=1/(1+r), shared by both soft sectors in stored off-shell kernels',
        extra_comparison='only the three-particle argument is reflected; production DA signs retained'))
    print('Checking all saved symmetric points and independent quadratures.',flush=True)
    checks=numerical_checks(specs)
    write_json('numerical_checks.json',checks)
    assert checks['baseline']['max_abs_g_error']<2e-8
    assert checks['baseline']['max_abs_width_error_keV']<1e-5
    assert checks['max_line_vs_rectangle_error']<5e-6
    assert checks['max_line_quadrature_order_error']<1e-10
    assert checks['max_Lprime_analytic_vs_finite_difference_error']<5e-8

    points=[]; windows=[]; curves=[]; envelopes=[]; cross=[]; hardchecks=[]
    for s in specs:
        inp=inputs[s['id']]; rm=s['mass']**2/inp['m_ds']**2
        K=np.linspace(s['Kmin'],s['Kmax'],41)
        s0=np.linspace(s['smin'],s['smax'],21)
        ratios=np.unique(np.r_[np.linspace(.85,1.15*rm,161),
                               np.linspace(1.,rm,101),[1.,rm,.85*rm,1.15]])
        hard=hard_grid(K,s0,inp)
        comp=components_grid(K,s0,ratios,inp,hard)
        pi,gc,g,w=project(comp,K,ratios,s,inp)
        np.savez_compressed(HERE/(s['id']+'_scan.npz'),
            M2=K,s0=s0,r=ratios,components=comp,Pi_components=pi,
            g_components=gc,g_signed=g,Gamma_keV=w)
        mid=(len(K)//2,len(s0)//2)
        for j,r in enumerate(ratios):
            curves.append(dict(state=s['id'],r=float(r),r_mass=rm,u0=1/(1+r),
                g_midpoint=float(g[j,mid[0],mid[1]]),Gamma_midpoint_keV=float(w[j,mid[0],mid[1]]),
                g_min=float(g[j].min()),g_max=float(g[j].max()),
                Gamma_min_keV=float(w[j].min()),Gamma_max_keV=float(w[j].max())))
        for label,r in [('symmetric',1.),('mass_ratio',rm)]:
            j=int(np.argmin(abs(ratios-r)))
            row=point_row(s,inp,r,float(K[mid[0]]),float(s0[mid[1]]),
                comp[:,:,j,mid[0],mid[1]],pi[:,j,mid[0],mid[1]],gc[:,j,mid[0],mid[1]],g[j,mid[0],mid[1]],w[j,mid[0],mid[1]],label)
            points.append(row)
            windows.append(dict(state=s['id'],prescription=label,r=r,
                g_min=float(g[j].min()),g_max=float(g[j].max()),
                Gamma_min_keV=float(w[j].min()),Gamma_max_keV=float(w[j].max())))
            cs=components_grid(K[mid[0]:mid[0]+1],s0[mid[1]:mid[1]+1],[r],inp,
                hard[:,mid[0]:mid[0]+1,mid[1]:mid[1]+1],split=True)
            ps,gcs,gs,ws=project(cs,K[mid[0]:mid[0]+1],[r],s,inp)
            cross.append(point_row(s,inp,r,float(K[mid[0]]),float(s0[mid[1]]),
                cs[:,:,0,0,0],ps[:,0,0,0],gcs[:,0,0,0],gs.item(),ws.item(),
                'split_argument_diagnostic'))
        for name,select in [('between_prescriptions',(ratios>=1)&(ratios<=rm)),
                            ('extended_scan',np.ones(len(ratios),dtype=bool))]:
            idx=np.flatnonzero(select); gg=g[select]; ww=w[select]
            loc=np.unravel_index(np.argmax(ww),ww.shape)
            envelopes.append(dict(state=s['id'],scan=name,rmin=float(ratios[idx[0]]),
                rmax=float(ratios[idx[-1]]),points=int(gg.size),
                g_min=float(gg.min()),g_max=float(gg.max()),
                Gamma_min_sampled_keV=float(ww.min()),Gamma_max_keV=float(ww.max()),
                zero_in_range=bool(gg.min()<=0<=gg.max()),
                max_at_r=float(ratios[idx[loc[0]]]),max_at_M2=float(K[loc[1]]),
                max_at_s0=float(s0[loc[2]])))
        hhigh=hard_grid(K[[0,20,40]],s0[[0,10,20]],inp,n=2400)
        hbase=hard[:,[0,20,40],:][:,:,[0,10,20]]
        hardchecks.append(dict(state=s['id'],max_abs_n1600_vs_2400=float(abs(hhigh-hbase).max())))
        print(s['id'], 'symmetric/mass widths',points[-2]['Gamma_keV'],points[-1]['Gamma_keV'],flush=True)
    write_csv('benchmark_midpoints.csv',points)
    write_csv('benchmark_window_ranges.csv',windows)
    write_csv('ratio_scan.csv',curves)
    write_csv('scan_envelopes.csv',envelopes)
    write_csv('fraction_convention_diagnostic.csv',cross)
    write_csv('hard_quadrature_checks.csv',hardchecks)
    checks['max_hard_quadrature_abs_change']=max(row['max_abs_n1600_vs_2400'] for row in hardchecks)
    checks['all_primary_results_finite']=all(np.isfinite(row['Gamma_midpoint_keV']) for row in curves)
    before=json.loads((HERE/'baseline_sha256.json').read_text())
    changed=[p for p,h in before.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h]
    checks['protected_or_reference_files_changed']=changed
    write_json('numerical_checks.json',checks)
    assert not changed
    print('Saved all scans and checks in',HERE,flush=True)


if __name__=='__main__':
    main()
