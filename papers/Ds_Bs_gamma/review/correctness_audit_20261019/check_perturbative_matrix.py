"""Independent 4x4 Dirac traces, plus the printed 2025 massless-light limit.

Run from any directory. Outputs stay in this audit folder. This tests the
current LO spectral matrix, not the original 2025 numerical implementation.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from twopoint_ds1_matrix_sumrule import Inputs, spectral_matrix

eye2, zero2 = np.eye(2), np.zeros((2, 2))
pauli = [np.array([[0, 1], [1, 0]]),
         np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
gamma = [np.block([[eye2, zero2], [zero2, -eye2]])]
gamma += [np.block([[zero2, sig], [-sig, zero2]]) for sig in pauli]
gamma5 = 1j * gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]


def from_trace(s, mh, ml):
    """Cut loop: antiquark (qslash-ml), quark (kslash+mh), Dirac adjoint."""
    energy = math.sqrt(s)
    momentum = math.sqrt((s - (mh + ml)**2) * (s - (mh - ml)**2)) / (2 * energy)
    el = (s + ml**2 - mh**2) / (2 * energy)
    eh = (s + mh**2 - ml**2) / (2 * energy)
    anti = el * gamma[0] + momentum * gamma[3] - ml * np.eye(4)
    quark = eh * gamma[0] - momentum * gamma[3] + mh * np.eye(4)
    axial = [g @ gamma5 for g in gamma[1:]]
    tensor = [1j * (0.5j * (g @ gamma[0] - gamma[0] @ g))
              @ gamma5 * energy / (mh + ml) for g in gamma[1:]]
    currents = [axial, tensor]
    trace = np.empty((2, 2))
    for a in range(2):
        for b in range(2):
            raw = sum(np.trace(anti @ currents[a][i] @ quark @ gamma[0]
                               @ currents[b][i].conj().T @ gamma[0])
                      for i in range(3)) / 3
            assert abs(raw.imag) < 1e-10
            trace[a, b] = raw.real
    # Nc * sqrt(lambda)/(16*pi^2*s), including the cut phase space.
    return trace * 3 * (2 * energy * momentum) / (16 * math.pi**2 * s)


def main():
    records = []
    max_relative = 0.0
    for mh, ml in [(1.27, 0), (1.27, .093), (1.27, .126),
                   (4.18, 0), (4.18, .093), (4.18, .126)]:
        for multiplier in [1.01, 1.2, 2, 4]:
            s = multiplier * (mh + ml)**2
            direct = from_trace(s, mh, ml)
            implemented = spectral_matrix(np.array([s]), Inputs(mc=mh, ms=ml))[0]
            rel = float(np.max(np.abs(direct - implemented)) / np.max(np.abs(direct)))
            max_relative = max(max_relative, rel)
            assert rel < 1e-11
            assert np.min(np.linalg.eigvalsh(direct)) >= -1e-12
            records.append(dict(mh=mh, ml=ml, s=s, relative_error=rel,
                                rho_AA=direct[0, 0], rho_AB=direct[0, 1], rho_BB=direct[1, 1]))
    with (HERE / "perturbative_trace_check.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)

    # Transcribed from final published PDF, Eqs. (A1),(A2), pp. 6-7.
    # Set ml=0 and all condensates=0; remove exp(-s/M2) from each integrand.
    for mh, s in [(1.27, 5), (4.18, 30)]:
        rho = from_trace(s, mh, 0)
        printed_A1 = -3 * (s - mh**2)**2 / (8 * math.pi**2 * s)
        printed_A2 = -(s - mh**2)**2 * (s**2 - mh**4) / (8 * math.pi**2 * mh**2 * s**2)
        assert np.isclose(printed_A1, -rho[0, 1], rtol=1e-12)
        assert np.isclose(printed_A2, rho[0, 0] - rho[1, 1], rtol=1e-12)
    sources = [ROOT / "scripts/twopoint_ds1_matrix_sumrule.py",
               HERE / "references/PhysRevD.111.014005.pdf",
               HERE / "references/article_R1.tex"]
    result = dict(points=len(records), max_relative_error=max_relative,
                  status="Current LO matrix passes independent cut-loop trace check",
                  printed_appendix_massless_limit={
                      "A1_label": "AA-BB", "A1_matches": "minus current AB",
                      "A2_label": "AB", "A2_matches": "current AA-BB",
                      "limitation": "Printed-label mismatch is not proof of an error in original numerical code"},
                  sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sources})
    (HERE / "perturbative_trace_check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
