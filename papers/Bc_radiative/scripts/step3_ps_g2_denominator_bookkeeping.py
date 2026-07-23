#!/usr/bin/env python3
"""Denominator bookkeeping for Bc1 -> Bc gamma radiative G2 terms.

The FeynCalc step3 scripts project the Dirac numerators.  This script records
the denominator powers and the Schwinger-parameter quadratic form for each
projected row.  It intentionally does not produce a numerical condensate
correction; its purpose is to make the next double-Borel reduction auditable.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


@dataclass(frozen=True)
class Topology:
    name: str
    segments: tuple[str, str, str]
    masses: tuple[str, str, str]
    shifts: tuple[str, str, str]
    schwinger_total: str
    phi_p2_pq: str
    phi_p2_P2: str
    mass_term: str


TOPOLOGIES = {
    "c_photon": Topology(
        name="c_photon",
        segments=("c_after_photon", "c_before_photon", "b_spectator"),
        masses=("m_c", "m_c", "m_b"),
        shifts=("k+q", "k", "k-p"),
        schwinger_total="A = a + b + r",
        phi_p2_pq="Phi = r(a+b)/A p^2 + 2 a r/A (p.q) - (a+b)m_c^2 - r m_b^2",
        phi_p2_P2="Phi = r b/A p^2 + a r/A P^2 - (a+b)m_c^2 - r m_b^2",
        mass_term="(a+b)m_c^2 + r m_b^2",
    ),
    "bbar_photon": Topology(
        name="bbar_photon",
        segments=("c_spectator", "b_before_photon", "b_after_photon"),
        masses=("m_c", "m_b", "m_b"),
        shifts=("k", "k-p", "k-p+q"),
        schwinger_total="A = a + b + r",
        phi_p2_pq=(
            "Phi = a(b+r)/A p^2 + 2 r(b+r)/A (p.q) "
            "- a m_c^2 - (b+r)m_b^2"
        ),
        phi_p2_P2=(
            "Phi = (b+r)(a-r)/A p^2 + r(b+r)/A P^2 "
            "- a m_c^2 - (b+r)m_b^2"
        ),
        mass_term="a m_c^2 + (b+r)m_b^2",
    ),
}


def powers_for_single_line(segments: tuple[str, str, str], selected: str) -> dict[str, int]:
    return {segment: (4 if segment == selected else 1) for segment in segments}


def powers_for_open_pair(
    segments: tuple[str, str, str], first: str, second: str
) -> dict[str, int]:
    return {segment: (2 if segment in {first, second} else 1) for segment in segments}


def rows():
    for topo in TOPOLOGIES.values():
        for current in ("A", "B"):
            for selected in topo.segments:
                yield {
                    "topology": topo.name,
                    "current": current,
                    "class": "single_line_G2",
                    "object": selected,
                    **{f"power_{segment}": powers_for_single_line(topo.segments, selected)[segment] for segment in topo.segments},
                    "segments": "; ".join(topo.segments),
                    "masses": "; ".join(topo.masses),
                    "shifts": "; ".join(topo.shifts),
                    "schwinger_total": topo.schwinger_total,
                    "phi_p2_pq": topo.phi_p2_pq,
                    "phi_p2_P2": topo.phi_p2_P2,
                    "mass_term": topo.mass_term,
                }
            for i, first in enumerate(topo.segments):
                for second in topo.segments[i + 1 :]:
                    yield {
                        "topology": topo.name,
                        "current": current,
                        "class": "open_pair_GG",
                        "object": f"{first},{second}",
                        **{f"power_{segment}": powers_for_open_pair(topo.segments, first, second)[segment] for segment in topo.segments},
                        "segments": "; ".join(topo.segments),
                        "masses": "; ".join(topo.masses),
                        "shifts": "; ".join(topo.shifts),
                        "schwinger_total": topo.schwinger_total,
                        "phi_p2_pq": topo.phi_p2_pq,
                        "phi_p2_P2": topo.phi_p2_P2,
                        "mass_term": topo.mass_term,
                    }


def main() -> None:
    table = list(rows())
    fieldnames = [
        "topology",
        "current",
        "class",
        "object",
        "power_c_after_photon",
        "power_c_before_photon",
        "power_b_spectator",
        "power_c_spectator",
        "power_b_before_photon",
        "power_b_after_photon",
        "segments",
        "masses",
        "shifts",
        "schwinger_total",
        "phi_p2_pq",
        "phi_p2_P2",
        "mass_term",
    ]
    csv_path = OUT / "step3_ps_g2_denominator_bookkeeping.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in table:
            writer.writerow(row)

    txt_path = OUT / "step3_ps_g2_denominator_bookkeeping.txt"
    lines = [
        "Bc1 -> Bc gamma G2 denominator bookkeeping",
        "===========================================",
        "",
        f"Rows: {len(table)} current-split denominator rows",
        "",
        "Schwinger parameters:",
        "- c_photon: a=c_after_photon, b=c_before_photon, r=b_spectator",
        "- bbar_photon: a=c_spectator, b=b_before_photon, r=b_after_photon",
        "",
        "Quadratic forms after loop-momentum shift:",
        f"- c_photon, p.q form: {TOPOLOGIES['c_photon'].phi_p2_pq}",
        f"- c_photon, P^2 form: {TOPOLOGIES['c_photon'].phi_p2_P2}",
        f"- bbar_photon, p.q form: {TOPOLOGIES['bbar_photon'].phi_p2_pq}",
        f"- bbar_photon, P^2 form: {TOPOLOGIES['bbar_photon'].phi_p2_P2}",
        "",
        "Power convention:",
        "- single_line_G2: selected propagator power is 4; the other two powers are 1.",
        "- open_pair_GG: selected propagator powers are 2 and 2; the remaining power is 1.",
        "",
        "Status:",
        "- This file attaches denominator powers to the projected numerator rows.",
        "- The next step is applying the double-Borel delta constraints to these",
        "  Schwinger-parameter kernels and integrating the remaining parameter.",
    ]
    txt_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
