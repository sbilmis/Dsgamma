"""Fixed-prescription numerical impact of the independently derived d=5 term.

This is a review calculation, not a replacement of accepted ensembles or
published angle inputs. It does not refit continuum thresholds or rerun MC.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(HERE))
from candidate_dimension5 import Inputs, completed_ope_matrix, heavy_line_mixed_component
from twopoint_ds1_matrix_sumrule import ope_matrix, rotation_from_matrix, rotation_matrix


def csv_write(name, rows):
    with (HERE/name).open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def projected(K,s0,theta,state,inp,completion):
    mat=completed_ope_matrix(K,s0,inp)[0] if completion else ope_matrix(K,s0,inp)[0]
    v=rotation_matrix(theta)[state]
    return float(v@mat@v)


def effective_mass(K,s0,theta,state,inp,completion):
    t,h=1/K,2e-4
    low=projected(1/(t-h),s0,theta,state,inp,completion)
    high=projected(1/(t+h),s0,theta,state,inp,completion)
    return math.sqrt(-(math.log(high)-math.log(low))/(2*h))


def main():
    protected = [p for p in (ROOT/"final_manuscript").rglob("*") if p.is_file()]
    protected += [ROOT/"scripts/twopoint_ds1_matrix_sumrule.py",
                  ROOT/"scripts/twopoint_bs1_matrix_sumrule.py",
                  ROOT/"outputs/twopoint_ds1_matrix_mc.csv",
                  ROOT/"outputs/twopoint_bs1_matrix_mc.csv",
                  ROOT/"outputs/final_window_mc_scan.csv"]
    before={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    grids, summaries = [], []
    for name,m,ka,kb,sa,sb,sigma in [
        ("Ds",1.27,3,6,9,11,.6), ("Bs",4.18,9,13,42,44,.1)]:
        rows=[]
        for K in np.linspace(ka,kb,21):
            for s0 in np.linspace(sa,sb,21):
                old_values,new_values=[],[]
                for qq in [-.0165,-.240**3]:
                    inp=Inputs(mc=m,ms=.126,qq=qq,kappa_s=.8,m0_sq=.8)
                    mat,_=ope_matrix(K,s0,inp)
                    old_values.append(rotation_from_matrix(mat)[0])
                    new_values.append(rotation_from_matrix(mat+heavy_line_mixed_component(K,inp))[0])
                rows.append(dict(sector=name,M2_GeV2=float(K),s0_GeV2=float(s0),
                    old_implementation_theta_old_qq_deg=old_values[0],
                    old_implementation_theta_new_qq_deg=old_values[1],
                    completed_d5_theta_old_qq_deg=new_values[0],
                    completed_d5_theta_new_qq_deg=new_values[1],
                    old_implementation_condensate_shift_deg=old_values[1]-old_values[0],
                    completed_d5_condensate_shift_deg=new_values[1]-new_values[0],
                    change_in_sensitivity_deg=(new_values[1]-new_values[0])-(old_values[1]-old_values[0]),
                    quoted_angle_error_deg=sigma))
        grids.extend(rows)
        shifts=[r["completed_d5_condensate_shift_deg"] for r in rows]
        summaries.append(dict(sector=name,midpoint=rows[220],
            completed_d5_shift_min_deg=min(shifts),completed_d5_shift_max_deg=max(shifts),
            maximum_change_in_sensitivity_deg=max(abs(r["change_in_sensitivity_deg"]) for r in rows),
            scope="current exact-mass LO plus local operators through dimension five; not original angle extraction"))
    provenance=json.loads((ROOT/"review/referee1_revision/provenance.json").read_text())
    residues=[]
    for name,row in provenance["selected_two_point_samples"].items():
        inp=Inputs(mc=float(row["mc_GeV" if name=="Ds" else "mb_GeV"]),
            ms=float(row["ms_GeV"]),qq=float(row["ss_GeV3"]),kappa_s=1,
            m0_sq=float(row["m0_sq_GeV2"]),
            mass_low=2.4595 if name=="Ds" else float(row["mass_low_GeV"]),
            mass_high=2.53511 if name=="Ds" else float(row["mass_high_GeV"]))
        K=float(row["M2_GeV2"]); theta=26.6 if name=="Ds" else 38.5
        for state,s0,mass in [(0,float(row["s01_GeV2"]),inp.mass_low),
                              (1,float(row["s02_GeV2"]),inp.mass_high)]:
            a=projected(K,s0,theta,state,inp,False)
            b=projected(K,s0,theta,state,inp,True)
            fold=math.sqrt(a*math.exp(mass**2/K)/mass**2)
            fnew=math.sqrt(b*math.exp(mass**2/K)/mass**2)
            residues.append(dict(sector=name,state=state+1,theta_deg=theta,M2_GeV2=K,s0_GeV2=s0,
                f_old_MeV=1000*fold,f_completed_d5_MeV=1000*fnew,
                f_change_percent=100*(fnew/fold-1),
                fixed_numerator_width_ratio=a/b,
                fixed_numerator_width_change_percent=100*(a/b-1),
                old_effective_mass_GeV=effective_mass(K,s0,theta,state,inp,False),
                completed_d5_effective_mass_GeV=effective_mass(K,s0,theta,state,inp,True)))
    after={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    assert before==after
    csv_write("completed_d5_condensate_grid.csv",grids)
    csv_write("representative_residue_impact.csv",residues)
    report=dict(status="isolated_verified_d5_completion; no production replacement",
        condensate_comparisons=summaries,representative_residue_impact=residues,
        protected_hashes=before,protected_sources_and_outputs_unchanged=True,
        limitations=["No threshold refit or accepted-ensemble regeneration.",
                     "No revised full uncertainties or publication widths are inferred.",
                     "Does not reproduce the original larger-OPE calculation."])
    (HERE/"numerical_impact.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report|{"protected_hashes":f"{len(before)} files verified unchanged"},indent=2))


if __name__=="__main__":
    main()
