"""Double-pole reduction checks for the tensor-current hard term.

The tensor hard numerator contains p^2/(p.q).  In the diagonal double
dispersion form used by the axial benchmark, the perturbative invariant is

    rho(s) / [(s - p^2) (s - P^2)],       P = p + q.

Only the part proportional to this double pole survives the double Borel
projection.  Therefore external virtualities in a numerator are replaced by the
dispersion variable s up to single-pole terms:

    p^2 / [(s-p^2)(s-P^2)]
      = s / [(s-p^2)(s-P^2)] - 1/(s-P^2).

The last term has no pole in p^2 and is annihilated by the double Borel
transform in (p^2, P^2).  The same identity holds for P^2.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


def main():
    lines = [
        "Double-pole reduction for tensor hard numerator",
        "================================================",
        "Let P^2=(p+q)^2 and use the diagonal spectral form",
        "  1/[(s-p^2)(s-P^2)].",
        "",
        "Then",
        "  p^2/[(s-p^2)(s-P^2)]",
        "    = s/[(s-p^2)(s-P^2)] - 1/(s-P^2),",
        "",
        "and",
        "  P^2/[(s-p^2)(s-P^2)]",
        "    = s/[(s-p^2)(s-P^2)] - 1/(s-p^2).",
        "",
        "The second term in each identity is a single pole.  It is removed by",
        "the double Borel projection in p^2 and P^2.  Therefore, for the",
        "double-pole sum rule, p^2 -> s and P^2 -> s inside numerator factors.",
        "",
        "Conclusion: the p^2 -> s replacement used in",
        "scripts/stage1_tensor_gb_hard_candidate.py is justified at the",
        "double-pole level.  Denominator-cancellation and single-pole pieces",
        "should remain excluded from the final double-Borel sum rule.",
    ]
    path = OUT / "step4_double_pole_p2_reduction.txt"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
