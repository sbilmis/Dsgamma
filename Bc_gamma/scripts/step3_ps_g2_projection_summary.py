#!/usr/bin/env python3
"""Summarize the generated Bc1 -> Bc gamma G2 projection chunks."""

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


SOURCES = [
    ("single_line_G2", "c_photon", OUT / "step3_ps_g2_single_line_cphoton.txt"),
    ("single_line_G2", "bbar_photon", OUT / "step3_ps_g2_single_line_bbarphoton.txt"),
    ("open_pair_GG", "c_photon", OUT / "step3_ps_g2_open_pair_cphoton.txt"),
    ("open_pair_GG", "bbar_photon", OUT / "step3_ps_g2_open_pair_bbarphoton.txt"),
]


def parse_source(kind: str, photon_topology: str, path: Path):
    rows = []
    pattern_single = re.compile(r"current=(?P<cur>[AB]), segment=(?P<a>[^,]+).*term_count=(?P<n>\d+)")
    pattern_pair = re.compile(r"current=(?P<cur>[AB]), pair=\((?P<a>[^,]+),(?P<b>[^)]+)\).*term_count=(?P<n>\d+)")
    for line in path.read_text().splitlines():
        match = pattern_single.search(line) or pattern_pair.search(line)
        if not match:
            continue
        rows.append(
            {
                "class": kind,
                "photon_topology": photon_topology,
                "current": match.group("cur"),
                "segment_a": match.group("a"),
                "segment_b": match.groupdict().get("b") or "",
                "term_count": int(match.group("n")),
                "source_file": str(path.relative_to(ROOT)),
            }
        )
    return rows


def main():
    rows = []
    for kind, photon_topology, path in SOURCES:
        rows.extend(parse_source(kind, photon_topology, path))
    out_csv = OUT / "step3_ps_g2_projection_summary.csv"
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    total = sum(r["term_count"] for r in rows)
    by_class = {}
    for row in rows:
        by_class[row["class"]] = by_class.get(row["class"], 0) + row["term_count"]
    lines = [
        "Bc1 -> Bc gamma G2 projection summary",
        "======================================",
        "",
        f"Projection rows: {len(rows)}",
        f"Total projected term count: {total}",
    ]
    for key, value in sorted(by_class.items()):
        lines.append(f"{key}: {value} projected terms")
    lines.extend(
        [
            "",
            "Status:",
            "- Numerator projection is complete for the dimension-4 topology inventory.",
            "- Double-Borel reduction/integration is not yet complete.",
            "- The next step is reducing each row to denominator powers in (p^2, P^2).",
        ]
    )
    out_txt = OUT / "step3_ps_g2_projection_summary.txt"
    out_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_txt}")


if __name__ == "__main__":
    main()
