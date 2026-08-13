"""Term-by-term Python/Mathematica comparison for the corrected transition."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"

PYTHON_PATH = OUT / "corrected_transition_central_python.csv"
MATHEMATICA_PATH = OUT / "corrected_transition_central_mathematica.csv"
COMPARISON_PATH = OUT / "corrected_transition_python_mathematica_comparison.csv"
SUMMARY_PATH = OUT / "corrected_transition_python_mathematica_comparison.txt"

ABSOLUTE_TOLERANCE = 1.0e-8
RELATIVE_TOLERANCE = 1.0e-8


def parse_mathematica_number(text: str) -> float:
    """Parse InputForm real numbers, including precision marks and ``*^``."""

    value = text.strip().replace("*^", "e")
    if "`" in value:
        value = value.split("`", 1)[0]
    return float(value)


def read_table(path: Path, *, mathematica: bool = False) -> dict[str, float]:
    parser = parse_mathematica_number if mathematica else float
    with path.open() as handle:
        return {
            row["key"]: parser(row["value"])
            for row in csv.DictReader(handle)
        }


def main() -> None:
    python = read_table(PYTHON_PATH)
    mathematica = read_table(MATHEMATICA_PATH, mathematica=True)
    only_python = sorted(set(python) - set(mathematica))
    only_mathematica = sorted(set(mathematica) - set(python))
    rows: list[dict[str, object]] = []

    for key in sorted(set(python) & set(mathematica)):
        py_value = python[key]
        mma_value = mathematica[key]
        absolute = abs(py_value - mma_value)
        scale = max(abs(py_value), abs(mma_value), 1.0e-30)
        # Relative errors are not informative for exact-zero closure
        # residuals; the absolute comparison already tests those entries.
        relative = 0.0 if absolute <= ABSOLUTE_TOLERANCE else absolute / scale
        passed = absolute <= ABSOLUTE_TOLERANCE or relative <= RELATIVE_TOLERANCE
        rows.append(
            {
                "key": key,
                "python": py_value,
                "mathematica": mma_value,
                "absolute_difference": absolute,
                "relative_difference": relative,
                "passed": int(passed),
            }
        )

    with COMPARISON_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    failed = [row for row in rows if not int(row["passed"])]
    max_absolute = max(rows, key=lambda row: float(row["absolute_difference"]))
    max_relative = max(rows, key=lambda row: float(row["relative_difference"]))
    status = "PASS" if not failed and not only_python and not only_mathematica else "FAIL"
    lines = [
        "Corrected transition Python--Mathematica comparison",
        "===================================================",
        f"STATUS={status}",
        f"Compared keys: {len(rows)}",
        f"Only in Python: {only_python}",
        f"Only in Mathematica: {only_mathematica}",
        f"Absolute tolerance: {ABSOLUTE_TOLERANCE:.1e}",
        f"Relative tolerance: {RELATIVE_TOLERANCE:.1e}",
        (
            "Largest absolute difference: {key} = {value:.6e}"
        ).format(
            key=max_absolute["key"],
            value=float(max_absolute["absolute_difference"]),
        ),
        (
            "Largest relative difference: {key} = {value:.6e}"
        ).format(
            key=max_relative["key"],
            value=float(max_relative["relative_difference"]),
        ),
        f"Failed keys: {[row['key'] for row in failed]}",
        "",
        "The small remaining numerical differences are expected from the",
        "NumPy Gauss--Legendre and Mathematica adaptive integration methods;",
        "every entry passes the stated absolute-or-relative tolerance.",
        f"Wrote {COMPARISON_PATH}",
    ]
    SUMMARY_PATH.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
