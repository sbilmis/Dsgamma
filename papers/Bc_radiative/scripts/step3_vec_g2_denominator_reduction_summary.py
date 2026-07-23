#!/usr/bin/env python3
"""Summarize the vector denominator-reduced radiative G2 term inventory."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def main() -> None:
    rows = list(csv.DictReader((OUT / "step3_vec_g2_denominator_reduction.csv").open()))
    for row in rows:
        row["eff"] = tuple(int(row[f"eff_n{i}"]) for i in (1, 2, 3))
        row["num"] = tuple(int(row[f"num_d{i}"]) for i in (1, 2, 3))
        row["borel_active_all_positive"] = all(x > 0 for x in row["eff"])
        row["has_contact_or_derivative"] = any(x <= 0 for x in row["eff"])

    by_sector = Counter()
    by_eff = Counter()
    by_topology = Counter()
    by_current = Counter()
    by_class = Counter()
    for row in rows:
        sector = "all_positive" if row["borel_active_all_positive"] else "contact_or_derivative"
        by_sector[sector] += 1
        by_eff[row["eff"]] += 1
        by_topology[(row["topology"], row["class"], row["current"], sector)] += 1
        by_current[(row["current"], sector)] += 1
        by_class[(row["class"], sector)] += 1

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
    with (OUT / "step3_vec_g2_denominator_reduction_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "Bc1 -> Bc* gamma denominator-reduced G2 inventory",
        "=================================================",
        "",
        f"Total reduced terms: {len(rows)}",
        f"All-positive effective denominator powers: {by_sector['all_positive']}",
        f"Contact/derivative sector terms: {by_sector['contact_or_derivative']}",
        "",
        "By current:",
    ]
    for key, count in sorted(by_current.items()):
        lines.append(f"- current {key[0]}, {key[1]}: {count}")
    lines.append("")
    lines.append("By field-insertion class:")
    for key, count in sorted(by_class.items()):
        lines.append(f"- {key[0]}, {key[1]}: {count}")
    lines.append("")
    lines.append("Most common effective denominator powers:")
    for eff, count in by_eff.most_common(15):
        lines.append(f"- {eff}: {count}")
    lines.extend(
        [
            "",
            "Interpretation:",
            "- The all-positive sector can be represented by ordinary Schwinger",
            "  parameters with positive powers.",
            "- The contact/derivative sector is sizeable, so it must be audited",
            "  before one quotes a full explicit vector radiative G2 correction.",
            "- The refined contact audit is performed separately by",
            "  step3_vec_g2_contact_borel_admissibility.py; in the current",
            "  reduction it finds no ordinary two-channel contact contribution.",
            "- This inventory is the completed algebraic reduction of the vector",
            "  G2 numerator; the remaining step is analytic/numeric evaluation of",
            "  the all-positive radiative G2 sector.",
        ]
    )
    (OUT / "step3_vec_g2_denominator_reduction_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
