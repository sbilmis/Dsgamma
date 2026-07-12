#!/usr/bin/env python3
"""Support audit for pseudoscalar Bc gamma G2 contact rows."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"

p2, P2 = sp.symbols("p2 P2")
pq = (P2 - p2) / 2
t1, t2, t3 = sp.symbols("t1 t2 t3", positive=True)


def dot(a, b):
    ap, aq = a
    bp, bq = b
    return sp.expand(ap * bp * p2 + (ap * bq + aq * bp) * pq)


def support_delta(topology: str, support: tuple[int, ...]) -> sp.Expr:
    if topology == "c_photon":
        shifts = {1: (0, 1), 2: (0, 0), 3: (-1, 0)}
    elif topology == "bbar_photon":
        shifts = {1: (0, 0), 2: (-1, 0), 3: (-1, 1)}
    else:
        raise ValueError(topology)
    params = {1: t1, 2: t2, 3: t3}
    total = sum(params[i] for i in support)
    if total == 0:
        return sp.Integer(0)
    summed = (
        sum(params[i] * shifts[i][0] for i in support),
        sum(params[i] * shifts[i][1] for i in support),
    )
    raw = sum(params[i] * dot(shifts[i], shifts[i]) for i in support)
    return sp.factor(sp.expand(raw - dot(summed, summed) / total))


def activity(delta: sp.Expr) -> str:
    has_p2 = sp.diff(delta, p2) != 0
    has_P2 = sp.diff(delta, P2) != 0
    if has_p2 and has_P2:
        return "double_borel_active"
    if has_p2 or has_P2:
        return "single_virtuality"
    return "no_external_virtuality"


def refined_class(support: str, act: str, delta: str) -> str:
    compact = delta.replace(" ", "")
    if act == "no_external_virtuality":
        return "no_external_virtuality"
    if act == "single_virtuality":
        return "single_physical_virtuality"
    if support == "{1,3}" and "-P2+2*p2" in compact:
        return "crossed_single_combination"
    return "ordinary_double_channel_candidate"


def main() -> None:
    rows = list(csv.DictReader((OUT / "step3_ps_g2_denominator_reduction.csv").open()))
    out_rows = []
    detail = Counter()
    refined = Counter()
    for row in rows:
        eff = tuple(int(row[f"eff_n{i}"]) for i in (1, 2, 3))
        support = tuple(i + 1 for i, value in enumerate(eff) if value > 0)
        sector = "all_positive" if all(value > 0 for value in eff) else "contact_or_derivative"
        delta = support_delta(row["topology"], support)
        act = activity(delta)
        support_text = "{" + ",".join(map(str, support)) + "}"
        ref = refined_class(support_text, act, str(delta)) if sector == "contact_or_derivative" else "ordinary_all_positive"
        if sector == "contact_or_derivative":
            detail[(row["topology"], row["current"], row["class"], support_text, act, ref)] += 1
            refined[ref] += 1
        out_rows.append(
            {
                "topology": row["topology"],
                "class": row["class"],
                "current": row["current"],
                "object": row["object"],
                "eff": str(eff),
                "support": support_text,
                "sector": sector,
                "borel_activity": act,
                "delta_external": str(delta),
                "refined_borel_support": ref,
            }
        )

    with (OUT / "step3_ps_g2_contact_borel_admissibility.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    summary_rows = [
        {
            "topology": k[0],
            "current": k[1],
            "class": k[2],
            "support": k[3],
            "borel_activity": k[4],
            "refined_borel_support": k[5],
            "count": v,
        }
        for k, v in sorted(detail.items())
    ]
    with (OUT / "step3_ps_g2_contact_borel_admissibility_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "Bc1 -> Bc gamma refined contact Borel-admissibility audit",
        "=========================================================",
        "",
        "Contact/derivative rows by refined support:",
    ]
    for key, count in sorted(refined.items()):
        lines.append(f"- {key}: {count}")
    ordinary = refined["ordinary_double_channel_candidate"]
    lines.extend(
        [
            "",
            f"Ordinary two-channel contact rows: {ordinary}",
            "Rows with crossed_single_combination depend on (p-q)^2 and are not",
            "added to the standard physical double-Borel sum.",
        ]
    )
    (OUT / "step3_ps_g2_contact_borel_admissibility.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
