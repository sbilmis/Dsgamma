"""Saddle-level reduction of the x_alpha gamma_beta Stage-2 kernels.

The raw FeynCalc kernels from step7 contain k.x, k.p and k.q.  Before doing the
full Fourier/Borel derivative treatment, it is useful to impose the light-cone
momentum routing

    k = p + a q,      q^2 = 0,

where a = alpha_qbar + v alpha_g is the photon momentum fraction entering the
heavy-quark line.  This script reduces the ratio checks from step7 under that
substitution and, optionally, the heavy-line on-shell relation

    k^2 = m_Q^2 = p^2 + 2 a p.q.

The result is diagnostic only.  It shows why the x_alpha gamma_beta term cannot
be absorbed into the simple m_Q/(m_Q+m_s) trace ratio before the derivative
operation associated with x_alpha/(q.x) is performed.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


def main():
    ratios = {
        "S": "(p^2 + (a+1) p.q) / m_Q^2",
        "T1": "(p^2 + (3a+2) p.q) / (3 m_Q^2)",
        "T2": "(p^2 + (2-a) p.q) / m_Q^2",
        "T3": "(p^2 + 2a p.q) / m_Q^2",
    }
    ratios_onshell = {
        "S": "(p^2 + (a+1) p.q) / (p^2 + 2a p.q)",
        "T1": "(p^2 + (3a+2) p.q) / (3[p^2 + 2a p.q])",
        "T2": "(p^2 + (2-a) p.q) / (p^2 + 2a p.q)",
        "T3": "1",
    }

    lines = [
        "Saddle reduction of x_alpha gamma_beta ratio checks",
        "===================================================",
        "Substitution: k = p + a q, q^2 = 0, so",
        "  k.q = p.q, k.p = p^2 + a p.q, (k.x-p.x)/(q.x)=a.",
        "",
        "Ratios are defined as B_kernel / [(mQ/(mQ+ms)) A_kernel].",
        "They would equal 1 if the simple trace-ratio rule held.",
        "",
    ]
    for name, expr in ratios.items():
        lines.append(f"{name}: {expr}")
        lines.append(f"{name} with k^2=mQ^2=p^2+2 a p.q: {ratios_onshell[name]}")
    lines.extend(
        [
            "",
            "Conclusion: the x_alpha gamma_beta term has nontrivial a and",
            "virtuality dependence at the raw saddle level.  Its contribution",
            "must be obtained by applying the Fourier/Borel derivative rules,",
            "not by multiplying the axial F1 term by mQ/(mQ+ms).",
        ]
    )

    path = OUT / "step8_xgamma_saddle_reduction.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
