"""Compare independent Wolfram transcription with the Python reconstruction."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
grid = list(csv.DictReader((HERE / "condensate_shift_grid.csv").open()))
wolfram = json.loads((HERE / "wolfram_high_precision.json").read_text())
keys = ("theta_old_deg", "theta_new_deg", "delta_theta_deg",
        "A1_old", "A1_new", "A2_old", "A2_new")
errors = {key: 0.0 for key in keys}
for row in wolfram:
    matches = [r for r in grid if r["sector"] == row["sector"]
               and abs(float(r["M2_GeV2"])-row["M2_GeV2"]) < 1e-12
               and abs(float(r["s0_GeV2"])-row["s0_GeV2"]) < 1e-12]
    assert len(matches) == 1
    for key in keys:
        error = abs(float(matches[0][key])-row[key])
        errors[key] = max(errors[key], error)
        assert error < (1e-9 if key.endswith("deg") else 1e-11), (row, key, error)
assert len(wolfram) == 18
summary = json.loads((HERE / "comparison_summary.json").read_text())
for path, expected in summary["hashes"].items():
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected, path
result = dict(
    points=len(wolfram), matched_scalar_values=len(wolfram)*len(keys),
    maximum_absolute_errors=errors, all_pass=True,
    method="Independent Wolfram algebra grouped by condensate power; 50-digit adaptive quadrature",
    limitation="Numerical/transcription agreement does not establish the original program or resolve the OPE audit",
    manuscript_and_input_sources_unchanged=True,
)
(HERE / "verification.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result, indent=2))
