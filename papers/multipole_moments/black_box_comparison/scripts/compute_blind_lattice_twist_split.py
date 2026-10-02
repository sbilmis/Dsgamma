"""Blind lattice-normalized twist split from our kernels; no counterpart reads."""
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime, timezone
import ast
import csv
import hashlib
import json
import math
import numpy as np
from recalculation_engine import DEFAULT, STATES, OUT, Engine, two_point, rotate, observables

DEST = OUT / 'blind_lattice_twist_split'
GROUPS = {
    'hard': ('pertQ', 'pertS'),
    'tw2': ('tw2',),
    'tw3': ('tw3V', 'tw3A', '3p_V', '3p_A'),
    'tw4': ('tw4A', 'tw4B', '3p_S', '3p_St', '3p_T1', '3p_T2', '3p_T3', '3p_T4', '3p_Sg', '3p_T4g'),
    'contact': ('W',),
}
CHI = -.0551 / DEFAULT['ssCondensate']
PARAMS = DEFAULT | {'chi': CHI, 'zeta2': 31/110}


def write_csv(name, rows):
    with (DEST/name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def windows(state, theta):
    accepted = []
    for delta in (.4, .45, .5, .55, .6):
        s0 = (state['mass']+delta)**2
        for t in (np.arange(12, 51)/10 if state['mQ']<2 else np.arange(12,49)/4):
            mat, der, full, pieces = two_point(t, s0, state | PARAMS)
            N = rotate(mat, theta, state['state_index'])
            nf = rotate(full, theta, state['state_index'])
            d5 = rotate(pieces['d5_light']+pieces['d5_heavy'], theta, state['state_index'])
            m2 = rotate(der, theta, state['state_index'])/N if N else float('nan')
            mass = math.sqrt(m2) if m2>0 else float('nan')
            if N>0 and nf>0 and N/nf>=.5 and abs(d5)/N<=.1 and abs(mass/state['mass']-1)<=.03:
                accepted.append((float(t), s0, N))
    return accepted


def main():
    DEST.mkdir(exist_ok=True)
    engine = Engine()
    chi_sectors = set()
    sources = [Path(__file__).resolve(), Path(__file__).with_name('recalculation_engine.py')]
    for name in ('hard_two_particle_coefficients.json', 'three_particle_coefficients.json', 'contact_coefficients.json'):
        path = OUT/name
        sources.append(path)
        for block, shape, ch, sector, structure, j, expression in json.loads(path.read_text()):
            if any(isinstance(n, ast.Name) and n.id=='chi' for n in ast.walk(ast.parse(expression, mode='eval'))):
                chi_sectors.add(sector)
    assert chi_sectors == {'tw2'}, chi_sectors
    rows = []
    counts = {}
    max_add = {'mu_over_muN': 0., 'FD_without_QH': 0.}
    max_sign = 0.
    max_two = 0.
    max_old_window = 0.
    max_charge = 0.
    max_refine = 0.
    refine_engine = Engine(quadrature=512)
    for angle, expected, filename in (
        ('nominal', (37,25,57,17), 'accepted_points.csv'),
        ('theta90', (55,79,64,30), 'accepted_points_theta90.csv'),
    ):
        # Own archived windows used only to verify regenerated coordinates.
        archived = list(csv.DictReader((OUT/filename).open()))
        for state, expected_count in zip(STATES, expected):
            theta = state['theta'] if angle=='nominal' else 90.
            state_run = state | {'theta': theta}
            points = windows(state, theta)
            assert len(points)==expected_count, (angle,state,len(points))
            own = {(float(r['MB2_GeV2']),float(r['s0_GeV2'])):r for r in archived if r['state']==state['state']}
            assert set(own)=={(t,s0) for t,s0,N in points}
            counts[angle, state['state']] = len(points)
            for t,s0,N in points:
                max_old_window = max(max_old_window,abs(N-float(own[t,s0]['B_GeV4'])))
                sec,_ = engine.sectors(t,s0,state | PARAMS)
                assert set(sec)=={k for names in GROUPS.values() for k in names}, set(sec)
                reversed_sec,_ = engine.sectors(t,s0,state | PARAMS | {'chi':-CHI})
                zero_sec,_ = engine.sectors(t,s0,state | PARAMS | {'chi':0.})
                for key in sec:
                    expect = -sec[key] if key=='tw2' else sec[key]
                    max_sign = max(max_sign,float(np.max(np.abs(reversed_sec[key]-expect))))
                    assert np.array_equal(zero_sec[key], np.zeros_like(sec[key]) if key=='tw2' else sec[key])
                normal_tp = two_point(t,s0,state | PARAMS)
                reversed_tp = two_point(t,s0,state | PARAMS | {'chi':-CHI})
                zero_tp = two_point(t,s0,state | PARAMS | {'chi':0.})
                for tp in (reversed_tp,zero_tp):
                    for a,b in zip(normal_tp[:3],tp[:3]):
                        max_two = max(max_two,float(np.max(np.abs(a-b))))
                    for k in normal_tp[3]:
                        max_two = max(max_two,float(np.max(np.abs(normal_tp[3][k]-tp[3][k]))))
                refined,_ = refine_engine.sectors(t,s0,state | PARAMS)
                for sign in (1,-1):
                    signed = {k:(sign*v if k=='tw2' else v) for k,v in sec.items()}
                    obs = observables(signed,N,state_run)
                    fine = observables({k:(sign*v if k=='tw2' else v) for k,v in refined.items()},N,state_run)
                    for metric in ('mu_qeps','FD_clean','FC_g'):
                        max_refine = max(max_refine,abs(obs[metric]-fine[metric]))
                    row = dict(angle=angle,twist2_sign=f'{sign:+d}',state=state['state'],M2_GeV2=t,s0_GeV2=s0,
                               f_GeV=math.sqrt(N*math.exp(state['mass']**2/t)/state['mass']**2))
                    mu_parts = {}; fd_parts = {}
                    for group,names in GROUPS.items():
                        mat = sum(signed[k] for k in names)
                        fm = -rotate(mat[:,3],theta,state['state_index'])/N
                        g3 = -state['mass']**2*rotate(mat[:,2]-mat[:,1],theta,state['state_index'])/N
                        mu_parts[group] = .93827208816/state['mass']*fm
                        fd_parts[group] = g3-fm
                    for prefix,parts,total in (
                        ('mu_over_muN',mu_parts,obs['mu_qeps']),
                        ('FD',fd_parts,obs['FD_clean']-(state['eQ']+1/3)),
                    ):
                        row.update({prefix+'_'+k:v for k,v in parts.items()})
                        row[prefix+'_total'] = total
                        residual = abs(sum(parts.values())-total)
                        max_add['mu_over_muN' if prefix=='mu_over_muN' else 'FD_without_QH'] = max(max_add['mu_over_muN' if prefix=='mu_over_muN' else 'FD_without_QH'], residual)
                        assert residual<1e-10
                    # FD_total is reserved for the requested result WITH Q_H.
                    row['FD_without_QH_total'] = row.pop('FD_total')
                    row.update(FC_g=obs['FC_g'],FD_total=obs['FD_clean'],D_efm2=obs['D_fm2'])
                    max_charge = max(max_charge,abs(obs['FC_g']-(state['eQ']+1/3)))
                    rows.append(row)
            print(angle,state['state'],len(points),'windows computed',flush=True)
    assert max_two==0. and max_sign==0.
    summaries = []
    grouped = defaultdict(list)
    for r in rows:
        grouped[r['angle'],r['twist2_sign'],r['state']].append(r)
    for key,data in grouped.items():
        summary = dict(zip(('angle','twist2_sign','state'),key))
        summary['n_windows'] = len(data)
        for label,column in (('mu_over_muN','mu_over_muN_total'),('FD','FD_total'),('D_efm2','D_efm2')):
            for suffix,val in zip(('p16','median','p84'),np.quantile([r[column] for r in data],[.16,.5,.84],method='linear')):
                summary[label+'_'+suffix] = float(val)
        for prefix in ('mu_over_muN','FD'):
            for group in GROUPS:
                summary[prefix+'_'+group+'_median'] = float(np.median([r[prefix+'_'+group] for r in data]))
        summaries.append(summary)
    write_csv('twist_split_windows.csv',rows)
    write_csv('twist_split_summary.csv',summaries)
    note = f'''# Blind lattice-normalized twist split

Computed with our numerical engine and own symbolic coefficient exports. No counterpart code, numbers, or result files were read. The accepted windows were regenerated from the two-point OPE and matched our archived coordinates: nominal 37, 25, 57, 17 (136 total); refitted theta90 55, 79, 64, 30 (228 total). There are {len(rows)} signed window rows and {len(summaries)} summary rows.

Inputs at DA scale 1 GeV: zeta1=0.4, zeta2=31/110, kappa=0.2, kappaPlus=0, f3Gamma=-0.004 GeV^2, omegaV=3.8, omegaA=-2.1, ssCondensate=0.8*(-0.240)^3={DEFAULT['ssCondensate']:.17g} GeV^3, m0^2=0.8 GeV^2. zeta1Plus=zeta2Plus=phi2=0; ms=0.093 GeV. The supplied lattice normalization -55.1 MeV is used directly (no extra running factor): chi={CHI:.17g} GeV^-2. +1 means f_perp=-0.0551 GeV; -1 means f_perp=+0.0551 GeV. The labels denote the relative twist-2 sign, not the numerical sign of chi.

Masses/angles/heavy-quark masses: Ds_1(2460): 2.4595 GeV, 26.6 deg, mc=1.27 GeV; Ds_1(2536): 2.53511 GeV, 26.6 deg, mc=1.27 GeV; Bs_1(5750): 5.7502 GeV, 38.5 deg, mb=4.18 GeV; Bs_1(5830): 5.82911 GeV, 38.5 deg, mb=4.18 GeV. theta90 uses J_A for the lighter and J_B for the heavier state, with refitted windows and consistently rotated residues. Q_H=1 for charm and 0 for bottom.

Grouping: hard=pertQ+pertS; tw2=tw2; tw3=two-particle tw3V+tw3A and three-particle V+A; tw4=two-particle tw4A+tw4B and three-particle S, S-tilde, T1–T4, S_gamma, T4_gamma; contact=Wilson-line W. These groups exhaust every evaluated sector.

chi enters ONLY twist 2, and does NOT enter the two-point function. An AST scan of all coefficient exports found chi only in sector tw2. On every window, evaluating chi=0 leaves all other sectors identical and removes tw2; evaluating -chi is exactly identical to multiplying tw2 alone by -1. All two-point matrices, derivatives, full integrals and individual pieces are identical for +chi, -chi and zero chi. Largest sign-equivalence residual: {max_sign:.17g}; two-point chi-dependence residual: {max_two:.17g}. ssCondensate is unchanged in every other term.

Extraction: F_M=-Pi_qeps/N; G3=-m_H^2*(Pi_qq-Pi_qp)/N; F_D=Q_H-F_M+G3. mu/mu_N=0.93827208816 GeV*F_M/m_H; D/(e fm^2)=F_D*(0.1973269804 GeV fm)^2/m_H^2; F_C_g=Pi_g/(2N). N is the rotated two-point Borel residue, and f=sqrt(N*exp(m_H^2/M^2)/m_H^2).

CSV convention: FD_hard, FD_tw2, FD_tw3, FD_tw4, FD_contact sum to FD_without_QH_total. FD_total includes Q_H. mu_over_muN_total is evaluated from the complete amplitude independently of the parts. Summary FD quantiles include Q_H; part medians exclude Q_H. Quantiles use NumPy linear interpolation; these are deterministic accepted-window spreads. Medians of parts need not sum to the median total.

Largest per-window additivity residuals: mu/mu_N {max_add['mu_over_muN']:.17g}; F_D without Q_H {max_add['FD_without_QH']:.17g}. Largest overall: {max(max_add.values()):.17g}. Largest FC_g-Q_H residual: {max_charge:.17g}. Largest regenerated-vs-archived two-point residue residual: {max_old_window:.17g} GeV^4. Doubling three-point Gauss-Legendre quadrature from 256 to 512 changes mu/mu_N, FD and FC_g by at most {max_refine:.17g} in absolute units (common two-point normalization).

Continuum treatment is the central engine's bbkr prescription: positive Borel powers use the diagonal spectral integral, nonpositive DA powers are unsubtracted. This is the same model approximation as the central calculation; no window refit based on three-point results is performed.

The SHA-256 manifest was written after all outputs and checks and before delivery. Source hashes are included for reproducibility. The manifest itself is tracked in the Git commit; self-hashing it within the manifest is not possible.
'''
    (DEST/'twist_split_note.md').write_text(note)
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(),counterpart_data_used=False,
                    outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DEST.iterdir()) if p.name in ('twist_split_windows.csv','twist_split_summary.csv','twist_split_note.md')},
                    sources={str(p.relative_to(OUT.parents[3])):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (DEST/'sha256_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Outputs frozen; maximum additivity residual',max(max_add.values()),flush=True)

if __name__=='__main__':
    main()
