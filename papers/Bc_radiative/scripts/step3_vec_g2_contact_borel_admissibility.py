#!/usr/bin/env python3
"""Refine the vector G2 contact audit by physical double-Borel support.

The first contact-support audit labels a contact row as "double-Borel active"
whenever the reduced denominator contains both p^2 and P^2.  That condition is
necessary but not sufficient for an ordinary two-channel sum rule: a reduced
two-denominator term can depend only on the crossed momentum (p-q)^2 =
2 p^2 - P^2.  Such a term is not a simultaneous pole in the final and initial
hadron channels.

This script classifies the contact/derivative rows into ordinary physical
support, single-virtuality support, crossed single-combination support, and
no-external-virtuality support.  It is still a support-level audit; it does not
evaluate the full numerator/contact derivative maps.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def classify(row: dict[str, str]) -> str:
    """Return the refined Borel-support class for one reduced row."""

    support = row["support"]
    activity = row["borel_activity"]
    delta = row["delta_external"].replace(" ", "")

    if activity == "no_external_virtuality":
        return "no_external_virtuality"
    if activity == "single_virtuality":
        return "single_physical_virtuality"

    # The only first-pass active contact rows found so far are the anti-bottom
    # photon {d1,d3} rows.  Their denominator is controlled by (p-q)^2, i.e.
    # 2 p^2 - P^2, not by independent p^2 and P^2 pole variables.
    if support == "{1,3}" and "-P2+2*p2" in delta:
        return "crossed_single_combination"

    # Keep this branch for future extensions.  A genuine ordinary double-pole
    # contact term would need independent physical support in both channels.
    return "ordinary_double_channel_candidate"


def main() -> None:
    rows = list(csv.DictReader((OUT / "step3_vec_g2_contact_support_audit.csv").open()))
    contact_rows = [row for row in rows if row["sector"] == "contact_or_derivative"]

    refined_rows: list[dict[str, str]] = []
    counter: Counter[str] = Counter()
    detailed: Counter[tuple[str, str, str, str, str]] = Counter()

    for row in contact_rows:
        refined = classify(row)
        counter[refined] += 1
        detailed[(row["topology"], row["current"], row["class"], row["support"], refined)] += 1
        out = dict(row)
        out["refined_borel_support"] = refined
        refined_rows.append(out)

    with (OUT / "step3_vec_g2_contact_borel_admissibility.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(refined_rows[0].keys()))
        writer.writeheader()
        writer.writerows(refined_rows)

    summary_rows = [
        {
            "topology": topology,
            "current": current,
            "class": cls,
            "support": support,
            "refined_borel_support": refined,
            "count": count,
        }
        for (topology, current, cls, support, refined), count in sorted(detailed.items())
    ]
    with (OUT / "step3_vec_g2_contact_borel_admissibility_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "Bc1 -> Bc* gamma refined contact Borel-admissibility audit",
        "==========================================================",
        "",
        "This is a refinement of step3_vec_g2_contact_support_audit.py.",
        "The first audit counted rows as active if the reduced denominator",
        "contained both p^2 and P^2.  Here we additionally require ordinary",
        "two-channel support, not only a crossed single combination.",
        "",
        "Contact/derivative rows by refined support:",
    ]
    for key, count in sorted(counter.items()):
        lines.append(f"- {key}: {count}")

    crossed = counter["crossed_single_combination"]
    ordinary = counter["ordinary_double_channel_candidate"]
    lines.extend(
        [
            "",
            "Interpretation:",
            f"- {crossed} rows depend on the crossed virtuality (p-q)^2 = 2 p^2 - P^2.",
            "  They are not ordinary simultaneous initial/final hadron poles.",
            f"- {ordinary} contact rows have ordinary two-channel support at this",
            "  denominator-support level.",
            "- Therefore the previously isolated 131 first-pass active candidates",
            "  are crossed-contact candidates rather than direct physical double-pole",
            "  terms.  They should be kept in the audit trail, but they are not added",
            "  to the physical double-Borel numerical sum unless a separate crossed",
            "  contact prescription is introduced.",
        ]
    )

    text = "\n".join(lines) + "\n"
    (OUT / "step3_vec_g2_contact_borel_admissibility.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
