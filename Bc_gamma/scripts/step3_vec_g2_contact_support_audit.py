#!/usr/bin/env python3
"""Classify vector G2 contact rows by double-Borel support.

The denominator-reduction CSV tells us the effective powers of the three
propagator denominators.  A term with at least one nonpositive effective power
is in the contact/derivative sector.  Not every such row is equally dangerous:
after the denominator cancellation, some remaining denominator supports depend
on both external virtualities, while others depend on only one virtuality and
therefore vanish under a double Borel transform.

This script works at the level of denominator support only.  It does not yet
evaluate derivative contact terms.  It answers the reproducibility question:
"Which contact rows can still contribute to the double-Borel sum rule?"
"""

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


def dot(vec_a: tuple[sp.Expr, sp.Expr], vec_b: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    """Dot product for vectors represented as a*p + b*q."""

    ap, aq = vec_a
    bp, bq = vec_b
    return sp.expand(ap * bp * p2 + (ap * bq + aq * bp) * pq)


def support_delta(topology: str, support: tuple[int, ...]) -> sp.Expr:
    """Return the external-invariant part of the shifted denominator.

    The support entries are 1,2,3 for denominators whose effective powers are
    positive.  Mass constants are irrelevant for deciding Borel activity and
    are omitted.
    """

    if topology == "c_photon":
        # D1=(k+q)^2-mc^2, D2=k^2-mc^2, D3=(k-p)^2-mb^2.
        shifts = {1: (sp.Integer(0), sp.Integer(1)), 2: (sp.Integer(0), sp.Integer(0)), 3: (sp.Integer(-1), sp.Integer(0))}
    elif topology == "bbar_photon":
        # D1=k^2-mc^2, D2=(k-p)^2-mb^2, D3=(k-p+q)^2-mb^2.
        shifts = {1: (sp.Integer(0), sp.Integer(0)), 2: (sp.Integer(-1), sp.Integer(0)), 3: (sp.Integer(-1), sp.Integer(1))}
    else:
        raise ValueError(topology)

    params = {1: t1, 2: t2, 3: t3}
    total = sum(params[i] for i in support)
    if total == 0:
        return sp.Integer(0)
    sum_a = (
        sum(params[i] * shifts[i][0] for i in support),
        sum(params[i] * shifts[i][1] for i in support),
    )
    raw = sum(params[i] * dot(shifts[i], shifts[i]) for i in support)
    shifted = sp.expand(raw - dot(sum_a, sum_a) / total)
    return sp.factor(shifted)


def activity(delta: sp.Expr) -> str:
    has_p2 = sp.diff(delta, p2) != 0
    has_P2 = sp.diff(delta, P2) != 0
    if has_p2 and has_P2:
        return "double_borel_active"
    if has_p2 or has_P2:
        return "single_virtuality"
    return "no_external_virtuality"


def main() -> None:
    rows = list(csv.DictReader((OUT / "step3_vec_g2_denominator_reduction.csv").open()))
    out_rows: list[dict[str, object]] = []
    support_counter = Counter()
    contact_counter = Counter()
    for row in rows:
        eff = tuple(int(row[f"eff_n{i}"]) for i in (1, 2, 3))
        support = tuple(i + 1 for i, value in enumerate(eff) if value > 0)
        sector = "all_positive" if all(value > 0 for value in eff) else "contact_or_derivative"
        delta = support_delta(row["topology"], support)
        act = activity(delta)
        key = (row["topology"], row["current"], row["class"], sector, support, act)
        support_counter[key] += 1
        if sector == "contact_or_derivative":
            contact_counter[(support, act)] += 1
        out_rows.append(
            {
                "topology": row["topology"],
                "class": row["class"],
                "current": row["current"],
                "object": row["object"],
                "eff": str(eff),
                "support": "{" + ",".join(map(str, support)) + "}",
                "sector": sector,
                "borel_activity": act,
                "delta_external": str(delta),
            }
        )

    with (OUT / "step3_vec_g2_contact_support_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    active_rows = [
        row for row in out_rows
        if row["sector"] == "contact_or_derivative"
        and row["borel_activity"] == "double_borel_active"
    ]
    with (OUT / "step3_vec_g2_contact_active_candidates.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(active_rows[0].keys()))
        writer.writeheader()
        writer.writerows(active_rows)

    summary_rows = []
    for key, count in sorted(support_counter.items(), key=lambda item: (str(item[0]), item[1])):
        topology, current, cls, sector, support, act = key
        summary_rows.append(
            {
                "topology": topology,
                "current": current,
                "class": cls,
                "sector": sector,
                "support": "{" + ",".join(map(str, support)) + "}",
                "borel_activity": act,
                "count": count,
                "delta_external": str(support_delta(topology, support)),
            }
        )
    with (OUT / "step3_vec_g2_contact_support_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    contact_total = sum(contact_counter.values())
    active_total = sum(count for (support, act), count in contact_counter.items() if act == "double_borel_active")
    single_total = sum(count for (support, act), count in contact_counter.items() if act == "single_virtuality")
    none_total = sum(count for (support, act), count in contact_counter.items() if act == "no_external_virtuality")

    lines = [
        "Bc1 -> Bc* gamma contact-support audit",
        "======================================",
        "",
        f"Contact/derivative rows: {contact_total}",
        f"Double-Borel-active contact rows: {active_total}",
        f"Single-virtuality contact rows: {single_total}",
        f"No-external-virtuality contact rows: {none_total}",
        "",
        "Contact rows by remaining denominator support:",
    ]
    for (support, act), count in sorted(contact_counter.items(), key=lambda item: (item[0][1], item[0][0])):
        lines.append(f"- support {support}, {act}: {count}")
    lines.extend(
        [
            "",
            "Interpretation:",
            "- Rows with no_external_virtuality vanish under the double Borel transform.",
            "- Rows with single_virtuality are subtraction/contact terms for one channel",
            "  and vanish under a double Borel transform in both virtualities.",
            "- Rows marked double_borel_active remain candidates for derivative-contact",
            "  Borel maps and must be treated explicitly.",
            "- The active candidate rows are written separately to",
            "  step3_vec_g2_contact_active_candidates.csv.",
        ]
    )
    (OUT / "step3_vec_g2_contact_support_audit.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
