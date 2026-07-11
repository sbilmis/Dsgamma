#!/usr/bin/env python3
"""Summarize the denominator-reduced radiative G2 term inventory."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def main() -> None:
    rows = list(csv.DictReader((OUT / "step3_ps_g2_denominator_reduction.csv").open()))
    for row in rows:
        row["eff"] = tuple(int(row[f"eff_n{i}"]) for i in (1, 2, 3))
        row["num"] = tuple(int(row[f"num_d{i}"]) for i in (1, 2, 3))
        row["borel_active_all_positive"] = all(x > 0 for x in row["eff"])
        row["has_contact_or_derivative"] = any(x <= 0 for x in row["eff"])

    by_sector = Counter()
    by_eff = Counter()
    by_topology = Counter()
    for row in rows:
        sector = "all_positive" if row["borel_active_all_positive"] else "contact_or_derivative"
        by_sector[sector] += 1
        by_eff[row["eff"]] += 1
        by_topology[(row["topology"], row["class"], row["current"], sector)] += 1

    summary_rows = []
    for key, count in sorted(by_topology.items()):
        topology, cls, current, sector = key
        summary_rows.append(
            {
                "topology": topology,
                "class": cls,
                "current": current,
                "sector": sector,
                "count": count,
            }
        )
    with (OUT / "step3_ps_g2_denominator_reduction_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "Bc1 -> Bc gamma denominator-reduced G2 inventory",
        "================================================",
        "",
        f"Total reduced terms: {len(rows)}",
        f"All-positive effective denominator powers: {by_sector['all_positive']}",
        f"Contact/derivative sector terms: {by_sector['contact_or_derivative']}",
        "",
        "Most common effective denominator powers:",
    ]
    for eff, count in by_eff.most_common(15):
        lines.append(f"- {eff}: {count}")
    lines.extend(
        [
            "",
            "Interpretation:",
            "- The all-positive sector can be represented by ordinary Schwinger",
            "  parameters with positive powers.",
            "- The contact/derivative sector is not optional: it is the majority of",
            "  the reduced expression and requires derivative double-Borel terms.",
            "- Therefore a result based only on the all-positive sector would not be",
            "  the full explicit radiative G2 correction.",
        ]
    )
    (OUT / "step3_ps_g2_denominator_reduction_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
