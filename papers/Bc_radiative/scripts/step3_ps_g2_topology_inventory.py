#!/usr/bin/env python3
"""Topology inventory for Bc1 -> Bc gamma dimension-4 gluon condensates.

This is not a numerical OPE contribution.  It enumerates the background-gluon
insertions that must be derived for the radiative three-point correlator and
writes a machine-readable checklist.  The point is to avoid running one huge
FeynCalc job before we know the topology count.
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


SEGMENTS = {
    "c_photon": [
        ("c_after_photon", "k+q", "mc", "D_c1=(k+q)^2-m_c^2"),
        ("c_before_photon", "k", "mc", "D_c2=k^2-m_c^2"),
        ("b_spectator", "k-p", "mb", "D_c3=(k-p)^2-m_b^2"),
    ],
    "bbar_photon": [
        ("c_spectator", "k", "mc", "D_b1=k^2-m_c^2"),
        ("b_before_photon", "k-p", "mb", "D_b2=(k-p)^2-m_b^2"),
        ("b_after_photon", "k-p+q", "mb", "D_b3=(k-p+q)^2-m_b^2"),
    ],
}


def make_rows():
    rows = []
    for photon_topology, segments in SEGMENTS.items():
        for seg_name, momentum, mass, denominator in segments:
            rows.append(
                {
                    "class": "single_line_G2",
                    "photon_topology": photon_topology,
                    "segment_a": seg_name,
                    "segment_b": "",
                    "momentum_a": momentum,
                    "momentum_b": "",
                    "mass_a": mass,
                    "mass_b": "",
                    "denominator_a": denominator,
                    "denominator_b": "",
                    "description": "two background gluons contracted on one heavy propagator segment",
                }
            )
        for i, seg_a in enumerate(segments):
            for seg_b in segments[i + 1 :]:
                rows.append(
                    {
                        "class": "cross_line_GG" if seg_a[2] != seg_b[2] else "same_flavor_two_segment_GG",
                        "photon_topology": photon_topology,
                        "segment_a": seg_a[0],
                        "segment_b": seg_b[0],
                        "momentum_a": seg_a[1],
                        "momentum_b": seg_b[1],
                        "mass_a": seg_a[2],
                        "mass_b": seg_b[2],
                        "denominator_a": seg_a[3],
                        "denominator_b": seg_b[3],
                        "description": "one open background gluon on each listed propagator segment",
                    }
                )
    return rows


def main():
    rows = make_rows()
    path = OUT / "step3_ps_g2_topology_inventory.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = [
        "Bc1 -> Bc gamma G2 topology inventory",
        "======================================",
        "",
        f"Total topology rows: {len(rows)}",
        f"Single-line G2 rows: {sum(r['class'] == 'single_line_G2' for r in rows)}",
        f"Open-gluon pair rows: {sum(r['class'] != 'single_line_G2' for r in rows)}",
        "",
        "Interpretation:",
        "- This is a derivation checklist, not a numerical contribution.",
        "- Each row must be projected onto the E1 tensor and reduced to a double-Borel form.",
        "- The radiative G2 correction should not be borrowed from the two-point Bc-mixing OPE.",
    ]
    txt = OUT / "step3_ps_g2_topology_inventory.txt"
    txt.write_text("\n".join(summary) + "\n")
    print("\n".join(summary))
    print(f"Wrote {path}")
    print(f"Wrote {txt}")


if __name__ == "__main__":
    main()
