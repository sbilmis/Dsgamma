"""Literal transcription of published Eqs. (A1),(A2), NOT original code.

K=M^2; q=<sbar s>; z=m0^2; G2=<g^2G^2>; G3=<g^3G^3>.
The known current-definition typo is not used to rescale these expressions.
The second angle below is explicitly a label-assignment hypothesis motivated
by the independently checked ml=0 leading terms; it is not an adopted repair.
"""
from pathlib import Path
import csv
import json
import math
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/"scripts"))
from twopoint_ds1_matrix_sumrule import Inputs,local_components


def printed(K, s0, m, l=.126, q=-.0165*.8, z=.8,
            G2=4*math.pi**2*.012, G3=.57, nodes=256):
    G4=G2**2
    p=math.pi**2
    a1 = [
        m**5*(4*G4*m**2+24*z*l*m*G3)*q/(2**9*3**4*K**7),
        -m**3*((52*G4-72*G3*z)*m**2+48*z*l*m*(4*G3+G2*m**2))*q/(2**9*3**4*K**6),
        m*(2*m*(m*(14*G4-51*G3*z-12*G2*z*m**2)+6*l*G3*(z+2*m**2)-48*G2*z*l*m**2))*q/(2**8*3**3*K**5),
        -(5*G2-6*m*(7*z*m-2*m**3-2*l*(6*z-m**2)))*m*G2*q/(2**6*3**3*K**4),
        6*m*G3*(23*z*m-8*m**3+l*(30*z+8*m**2))*q/(2**8*3**3*m*K**4),
        -(G2*(4*z*m-6*m**3-l*(4*z-7*m**2))-12*z*l*m**4)*q/(144*K**3),
        G3*(2*z+12*l*m+7*m**2)*q/(288*m*K**3),
        -(2*m*(G2-6*z*m**2)-3*l*(G2-4*z*m**2))*q/(48*K**2),
        G3*(l*m**2+2*p*q)/(192*m*K**2*p),
        -m*(G2*l+24*p*(z+2*l*m)*q)/(96*K*p),
        -G3*(l+2*m)/(384*m*K*p),
        (G2*(l+m)-24*p*q*(z-2*l*m+4*m**2))/(96*m*p),
        -G3*(3*l+8*m)/(384*m**3*p),
    ]
    a2 = [
        -m**3*(14*G4*m**2+96*G3*z*l*m)*q/(2**8*3**4*K**6),
        -m*(-(22*G4-36*G3*z)*m**2-7*z*l*m*(17*G3+4*G2*m**2))*q/(2**7*3**3*K**5),
        -G2*m*(25*G2-6*m*(5*z*m-l*(14*z-4*m**2)))*q/(2**6*3**3*K**4),
        G3*(5*(G2+z*m**2)-l*m*(23*z-10*m**2))*q/(2**6*3**2*m*K**4),
        (3*G2*(G2-6*z*m**2)+l*m*G2*(20*z-39*m**2)+72*z*l*m**5)*q/(432*m*K**3),
        G3*(108*z*m-l*(20*z+81*m**2))*q/(2**6*3**3*m**2*K**3),
        -l*(3*G3-4*m**2*(3*G2-14*z*m**2))*q/(96*m**2*K**2),
        (2*G2-3*m*(4*z*m-l*(z+4*m**2)))*q/(12*m*K),
        -G3*(2*l-m)/(96*m*K*p),
        -(35*G2*m**2-2*l*(G2*m+8*p*(z+9*m**2)*q))/(96*m**2*p),
        -G3*(5*l-2*m)/(96*m**3*p),
    ]
    x,w = np.polynomial.legendre.leggauss(nodes)
    s=(s0-m**2)*x/2+(s0+m**2)/2
    weights=w*(s0-m**2)/2*np.exp(-s/K)
    a1_integrand = -(G2/(32*s*p)+3*(s-m**2)*(s*m-m**3-l*(s+m**2))/(8*m*s*p))
    a2_integrand = G2*(4*s-m**2)/(48*s**2*p)-G3*(2*s+m**2)/(384*m**2*s**2*p) \
                  -(s-m**2)**2*(s**2-6*s*l*m-m**4)/(8*p*m**2*s**2)
    f=math.exp(-m**2/K)
    return dict(A1=f*sum(a1)+float(weights@a1_integrand),
                A2=f*sum(a2)+float(weights@a2_integrand),
                A1_local=f*sum(a1), A2_local=f*sum(a2),
                A1_integral=float(weights@a1_integrand), A2_integral=float(weights@a2_integrand))


def angle(ab, difference):
    return math.degrees(.5*math.atan2(-2*ab, difference)) % 90


def main():
    rows=[]
    for sector,m,K,s0,published in [("Ds",1.27,4.5,10,26.6),("Bs",4.18,11,43,38.5)]:
        for mode, extras in [("full_printed",{}), ("without_gluon_condensates",dict(G2=0,G3=0)),
                             ("perturbative_only",dict(q=0,G2=0,G3=0))]:
            result=printed(K,s0,m,**extras)
            finer=printed(K,s0,m,nodes=512,**extras)
            assert max(abs(result[k]-finer[k]) for k in result)<1e-11
            rows.append(dict(sector=sector,mode=mode,M2=K,s0=s0,**result,
                        literal_labels_angle=angle(result["A2"],result["A1"]),
                        LO_label_hypothesis_angle=angle(-result["A1"],result["A2"]),
                        quoted_angle=published))
    with (HERE/"printed_appendix_midpoint_check.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    local_checks=[]
    for m,K,s0 in [(1.27,4.5,10),(4.18,11,43)]:
        q=-.0165*.8;z=.8
        full_local=printed(K,s0,m,l=0,q=q,z=z,G2=0,G3=0)
        d3_only=printed(K,s0,m,l=0,q=q,z=0,G2=0,G3=0)
        a1d5=full_local["A1_local"]-d3_only["A1_local"]
        a2d5=full_local["A2_local"]-d3_only["A2_local"]
        current=local_components(K,Inputs(mc=m,ms=0,qq=q,kappa_s=1,m0_sq=z))["d5_mixed"]
        current_difference=current[0,0]-current[1,1]
        assert np.isclose(a1d5,-current[0,1],rtol=1e-12)
        assert np.isclose(a2d5,2*current_difference,rtol=1e-12)
        local_checks.append(dict(m=m,M2=K,printed_A1_d5=a1d5,current_minus_AB_d5=-current[0,1],
                                 printed_A2_d5=a2d5,current_D_d5=current_difference,
                                 ratio_A2_to_current_D=a2d5/current_difference))
    (HERE/"printed_appendix_massless_d5_check.json").write_text(json.dumps(local_checks,indent=2)+"\n")
    print(json.dumps(rows,indent=2))


if __name__=="__main__":
    main()
