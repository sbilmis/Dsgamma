#!/usr/bin/env python3
"""Compare independent Python and Mathematica outputs term by term."""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CSV_DIR = ROOT / "outputs" / "csv"
ABS_TOL = 2.0e-8
REL_TOL = 2.0e-5


def main():
    py = pd.read_csv(CSV_DIR / "python_regression.csv")
    ma = pd.read_csv(CSV_DIR / "mathematica_regression.csv")
    merged = py.merge(ma, on=["sector", "quantity"], suffixes=("_python", "_mathematica"))
    if len(merged) != len(py) or len(merged) != len(ma):
        missing_py = set(zip(ma.sector, ma.quantity)) - set(zip(py.sector, py.quantity))
        missing_ma = set(zip(py.sector, py.quantity)) - set(zip(ma.sector, ma.quantity))
        raise RuntimeError(f"Regression key mismatch: missing_py={missing_py}, missing_ma={missing_ma}")
    merged["absolute_difference"] = np.abs(
        merged.value_python - merged.value_mathematica
    )
    merged["relative_difference"] = merged.absolute_difference / np.maximum(
        np.maximum(np.abs(merged.value_python), np.abs(merged.value_mathematica)),
        1.0e-30,
    )
    merged["pass"] = (merged.absolute_difference <= ABS_TOL) | (
        merged.relative_difference <= REL_TOL
    )
    merged["absolute_tolerance"] = ABS_TOL
    merged["relative_tolerance"] = REL_TOL
    merged.to_csv(CSV_DIR / "python_mathematica_regression.csv", index=False)
    validation = pd.read_csv(CSV_DIR / "validation_status.csv")
    mask = validation.gate == "python_mathematica_regression"
    validation.loc[mask, "status"] = "checked" if merged["pass"].all() else "failed"
    material_scale = np.maximum(
        np.abs(merged.value_python), np.abs(merged.value_mathematica)
    )
    material_relative = merged.loc[
        material_scale > ABS_TOL, "relative_difference"
    ].max()
    validation.loc[mask, "evidence"] = (
        f"{int(merged['pass'].sum())}/{len(merged)} entries pass; "
        f"max material rel={material_relative:.3e}; "
        f"max abs={merged.absolute_difference.max():.3e}"
    )
    validation.to_csv(CSV_DIR / "validation_status.csv", index=False)
    print(merged.to_string(index=False))
    if not merged["pass"].all():
        failed = merged[~merged["pass"]]
        raise SystemExit(f"Regression failed for {len(failed)} entries")


if __name__ == "__main__":
    main()
