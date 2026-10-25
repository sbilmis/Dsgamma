"""Independent NumPy Dirac/color traces with boosted external momenta."""
from pathlib import Path
import csv
import json
import math

import numpy as np

HERE = Path(__file__).resolve().parent
eye2, zero2, eye4 = np.eye(2), np.zeros((2, 2)), np.eye(4)
pauli = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
         np.diag([1, -1])]
gamma = [np.block([[eye2, zero2], [zero2, -eye2]])]
gamma += [np.block([[zero2, p], [-p, zero2]]) for p in pauli]
gamma5 = 1j * gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]
eta = np.array([1, -1, -1, -1])
metric = np.diag(eta)
sigma = np.array([[0.5j*(a@b-b@a) for b in gamma] for a in gamma])

generators = []
for i, j in [(0, 1), (0, 2), (1, 2)]:
    symmetric, antisymmetric = np.zeros((3,3), complex), np.zeros((3,3), complex)
    symmetric[i,j] = symmetric[j,i] = 0.5
    antisymmetric[i,j], antisymmetric[j,i] = -0.5j, 0.5j
    generators += [symmetric, antisymmetric]
generators += [np.diag([1, -1, 0])/2, np.diag([1, 1, -2])/(2*math.sqrt(3))]
color_trace = sum(np.trace(t@t) for t in generators).real
sigma_trace = sum(eta[a]*eta[b]*np.trace(sigma[a,b]@sigma[a,b])
                  for a in range(4) for b in range(4)).real
assert abs(color_trace-4) < 1e-13 and abs(sigma_trace-48) < 1e-13


def traced(m, d, s, spatial):
    p = np.r_[math.sqrt(s+np.dot(spatial, spatial)), spatial]
    p_lower = eta*p
    projection = (-metric + np.outer(p_lower, p_lower)/s)/3
    numerator = sum(p_lower[i]*gamma[i] for i in range(4)) + m*eye4
    currents = [
        [gamma[i]@gamma5 for i in range(4)],
        [1j*sum(sigma[i,j]*p_lower[j] for j in range(4))@gamma5/d
         for i in range(4)],
    ]
    out = np.zeros((2,2), complex)
    for i in range(2):
        for j in range(2):
            for mu in range(4):
                for nu in range(4):
                    right = gamma[0]@currents[j][nu].conj().T@gamma[0]
                    for a in range(4):
                        for b in range(4):
                            heavy_G = -(sigma[a,b]@numerator+numerator@sigma[a,b])/(4*(s-m*m)**2)
                            light_open_G = -eta[a]*eta[b]*sigma[a,b]/192
                            out[i,j] += projection[mu,nu]*color_trace*np.trace(
                                currents[i][mu]@heavy_G@right@light_open_G)
    assert np.max(np.abs(out.imag)) < 1e-11
    return out.real


def main():
    rows = []
    for m in [1.27, 4.18]:
        for ms in [0, .126]:
            d = m+ms
            for factor in [.4, .8, 1.4]:
                s = factor*m*m
                expected = -s/(12*d*(m*m-s)**2)*np.array([[0,1],[1,2*m/d]])
                for boost in [np.zeros(3), math.sqrt(s)*np.array([.31,-.17,.43])]:
                    matrix = traced(m,d,s,boost)
                    error = float(np.max(np.abs(matrix-expected)))
                    assert error < 1e-10
                    rows.append(dict(m=m, ms_in_normalization=ms, s=s,
                        boosted=bool(np.any(boost)), maximum_absolute_error=error,
                        AA=matrix[0,0], AB=matrix[0,1], BA=matrix[1,0], BB=matrix[1,1]))
    with (HERE/"covariant_trace_check.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    result=dict(points=len(rows),color_trace=color_trace,sigma_trace=sigma_trace,
        maximum_absolute_error=max(r["maximum_absolute_error"] for r in rows),
        checks="Independent gamma/color contraction, transverse projection, rest/boosted momentum",
        all_pass=True)
    (HERE/"covariant_trace_check.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
