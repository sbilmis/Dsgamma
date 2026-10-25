"""Paired condensate-only sensitivity check; not a published-baseline rerun.

All generated files stay beside this script. The primary calculation uses
the explicitly marked A1/A2 assignment reconstruction in the earlier audit.
The original numerical program and the mixed-condensate discrepancy remain
unresolved. A second, current-LO+d3+d5 calculation is a robustness diagnostic.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = ROOT / "review/correctness_audit_20261019"
sys.path.insert(0, str(AUDIT))
from check_printed_appendix import printed, angle
from twopoint_ds1_matrix_sumrule import Inputs, ope_matrix, rotation_from_matrix

QQ_OLD = -0.0165
QQ_NEW = -(0.240**3)
KAPPA = 0.8
SS_OLD = KAPPA * QQ_OLD
SS_NEW = KAPPA * QQ_NEW
FIXED = dict(l=0.126, z=0.8, G2=4 * math.pi**2 * 0.012, G3=0.57)
SECTORS = [
    dict(sector="Ds", m=1.27, M2_min=3.0, M2_max=6.0,
         s0_min=9.0, s0_max=11.0, published_angle=26.6, quoted_sigma=0.6),
    dict(sector="Bs", m=4.18, M2_min=9.0, M2_max=13.0,
         s0_min=42.0, s0_max=44.0, published_angle=38.5, quoted_sigma=0.1),
]


def write_csv(name, rows):
    with (HERE / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(K, s0, m, ss, nodes=256, **override):
    return printed(K, s0, m, q=ss, nodes=nodes, **(FIXED | override))


def theta(values):
    # A1 = -AB, A2 = AA-BB: assignment inferred from the checked LO terms.
    return angle(-values["A1"], values["A2"])


def paired_row(sector, K, s0, nodes=256):
    old = evaluate(K, s0, sector["m"], SS_OLD, nodes)
    new = evaluate(K, s0, sector["m"], SS_NEW, nodes)
    a, b = theta(old), theta(new)
    Nold, Nnew = 2 * old["A1"], 2 * new["A1"]
    Dold, Dnew = old["A2"], new["A2"]
    independent_difference = math.degrees(0.5 * math.atan2(
        Nnew * Dold - Nold * Dnew, Dnew * Dold + Nnew * Nold))
    assert abs((b - a) - independent_difference) < 5e-12
    return dict(
        sector=sector["sector"], M2_GeV2=float(K), s0_GeV2=float(s0),
        theta_old_deg=a, theta_new_deg=b, delta_theta_deg=b-a,
        quoted_sigma_deg=sector["quoted_sigma"],
        delta_over_quoted_sigma=(b-a)/sector["quoted_sigma"],
        delta_percent_of_reconstructed_old_angle=100*(b-a)/a,
        A1_old=old["A1"], A1_new=new["A1"],
        A2_old=old["A2"], A2_new=new["A2"],
    )


def main():
    protected = sorted(p for p in (ROOT / "final_manuscript").rglob("*") if p.is_file())
    sources = [
        Path(__file__).resolve(), AUDIT / "check_printed_appendix.py",
        AUDIT / "references/article_R1.tex",
        AUDIT / "references/PhysRevD.111.014005.pdf",
        ROOT / "scripts/twopoint_ds1_matrix_sumrule.py",
        ROOT / "scripts/twopoint_bs1_matrix_sumrule.py",
        ROOT / "scripts/mixing_angle_inputs.py",
    ]
    before = {str(p.relative_to(ROOT)): digest(p) for p in protected + sources}
    rows, midpoints, summaries = [], [], []
    quadrature_checks, current_checks, decompositions, envelopes = [], [], [], []
    affine_errors = []
    for sec in SECTORS:
        masses = np.linspace(sec["M2_min"], sec["M2_max"], 21)
        thresholds = np.linspace(sec["s0_min"], sec["s0_max"], 21)
        this_grid = [paired_row(sec, K, s0) for K in masses for s0 in thresholds]
        rows.extend(this_grid)
        midK, mids = float(masses[10]), float(thresholds[10])
        midpoint = paired_row(sec, midK, mids)
        midpoints.append(midpoint)
        lower = min(this_grid, key=lambda r: r["delta_theta_deg"])
        upper = max(this_grid, key=lambda r: r["delta_theta_deg"])
        array = np.array([r["delta_theta_deg"] for r in this_grid]).reshape(21, 21)
        summaries.append(dict(
            sector=sec["sector"], grid_shape=[21, 21],
            midpoint=midpoint,
            delta_min_deg=lower["delta_theta_deg"],
            delta_max_deg=upper["delta_theta_deg"],
            minimum_at_M2_s0=[lower["M2_GeV2"], lower["s0_GeV2"]],
            maximum_at_M2_s0=[upper["M2_GeV2"], upper["s0_GeV2"]],
            delta_decreases_with_M2_on_grid=bool(np.all(np.diff(array, axis=0) < 0)),
            delta_decreases_with_s0_on_grid=bool(np.all(np.diff(array, axis=1) < 0)),
            theta_old_grid_min_deg=min(r["theta_old_deg"] for r in this_grid),
            theta_old_grid_max_deg=max(r["theta_old_deg"] for r in this_grid),
            theta_new_grid_min_deg=min(r["theta_new_deg"] for r in this_grid),
            theta_new_grid_max_deg=max(r["theta_new_deg"] for r in this_grid),
        ))
        # Nine matched boundary/midpoint locations in each sector.
        for K in (masses[0], masses[10], masses[-1]):
            for s0 in (thresholds[0], thresholds[10], thresholds[-1]):
                coarse = paired_row(sec, K, s0, nodes=256)
                fine = paired_row(sec, K, s0, nodes=512)
                error = max(abs(coarse[key]-fine[key]) for key in
                            ("theta_old_deg", "theta_new_deg", "delta_theta_deg"))
                assert error < 1e-9
                quadrature_checks.append(dict(sector=sec["sector"], M2_GeV2=float(K),
                    s0_GeV2=float(s0), maximum_angle_error_deg=error))
                values = []
                for qq in (QQ_OLD, QQ_NEW):
                    inp = Inputs(mc=sec["m"], ms=FIXED["l"], qq=qq,
                                 kappa_s=KAPPA, m0_sq=FIXED["z"])
                    matrix, _ = ope_matrix(K, s0, inp)
                    values.append(rotation_from_matrix(matrix)[0])
                current_checks.append(dict(sector=sec["sector"], M2_GeV2=float(K),
                    s0_GeV2=float(s0), theta_old_deg=values[0], theta_new_deg=values[1],
                    delta_theta_deg=values[1]-values[0],
                    status="current_LO_d3_d5_diagnostic_not_original_implementation"))

        # Separate perturbative, q-independent gluon and all ss-dependent terms.
        perturbative = evaluate(midK, mids, sec["m"], 0, G2=0, G3=0)
        zero = evaluate(midK, mids, sec["m"], 0)
        old = evaluate(midK, mids, sec["m"], SS_OLD)
        new = evaluate(midK, mids, sec["m"], SS_NEW)
        half = evaluate(midK, mids, sec["m"], (SS_OLD+SS_NEW)/2)
        affine_errors.append(max(abs(half[k]-(old[k]+new[k])/2) for k in ("A1", "A2")))
        for invariant in ("A1", "A2"):
            decompositions.append(dict(sector=sec["sector"], invariant=invariant,
                M2_GeV2=midK, s0_GeV2=mids,
                perturbative=perturbative[invariant],
                gluon_only=zero[invariant]-perturbative[invariant],
                ss_dependent_old=old[invariant]-zero[invariant],
                ss_dependent_new=new[invariant]-zero[invariant],
                total_old=old[invariant], total_new=new[invariant]))
        for label, center, endpoints in (
            ("earlier_qq", QQ_OLD, (-0.018, -0.015)),
            ("updated_qq", QQ_NEW, (-(0.250**3), -(0.230**3))),
        ):
            endpoint_angles = [theta(evaluate(midK, mids, sec["m"], KAPPA*q)) for q in endpoints]
            envelopes.append(dict(sector=sec["sector"], input_set=label,
                M2_GeV2=midK, s0_GeV2=mids, qq_center_GeV3=center,
                qq_lower_GeV3=endpoints[0], qq_upper_GeV3=endpoints[1],
                theta_central_deg=theta(evaluate(midK, mids, sec["m"], KAPPA*center)),
                theta_endpoint_min_deg=min(endpoint_angles),
                theta_endpoint_max_deg=max(endpoint_angles),
                interpretation="qq_only_endpoint_envelope_not_total_uncertainty"))

    assert max(affine_errors) < 1e-12
    after = {str(p.relative_to(ROOT)): digest(p) for p in protected + sources}
    assert before == after
    write_csv("condensate_shift_grid.csv", rows)
    write_csv("midpoint_comparison.csv", midpoints)
    write_csv("quadrature_checks.csv", quadrature_checks)
    write_csv("current_truncation_crosscheck.csv", current_checks)
    write_csv("signed_midpoint_contributions.csv", decompositions)
    write_csv("qq_only_endpoint_envelopes.csv", envelopes)
    result = dict(
        status="conditional_appendix_reconstruction_not_validated_original_rerun",
        primary_formula="tan(2 theta)=2 F_A1/F_A2, branch 0<=theta<90 degrees",
        assignment="F_A1=-AB, F_A2=AA-BB, motivated by independently checked leading terms",
        unresolved="Original numerical implementation and mixed-condensate coefficient discrepancy",
        scope="central light-condensate change only; all other inputs and windows fixed",
        qq_scale_GeV=1.0, qq_old_GeV3=QQ_OLD, qq_new_GeV3=QQ_NEW,
        ss_old_GeV3=SS_OLD, ss_new_GeV3=SS_NEW, kappa_s=KAPPA,
        mixed_old_GeV5=FIXED["z"]*SS_OLD, mixed_new_GeV5=FIXED["z"]*SS_NEW,
        decrease_in_condensate_magnitude_percent=100*(1-QQ_NEW/QQ_OLD),
        fixed_inputs=FIXED, sectors=SECTORS, results=summaries,
        maximum_256_vs_512_node_angle_difference_deg=max(r["maximum_angle_error_deg"] for r in quadrature_checks),
        maximum_affine_condensate_residual=max(affine_errors),
        protected_manuscript_file_count=len(protected), protected_and_source_hashes_unchanged=True,
        hashes=before,
        limitations=[
            "Grid extrema are sampled extrema, not a continuous optimization proof.",
            "Quoted-error comparisons do not propagate or validate the original uncertainties.",
            "No scale running, change of strange mass, or change of kappa uncertainty is included.",
            "No numerical replacement of the published external angle inputs is made.",
        ],
    )
    (HERE / "comparison_summary.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"midpoints": midpoints, "grid_summary": summaries,
                      "quadrature_max_deg": result["maximum_256_vs_512_node_angle_difference_deg"],
                      "protected_files_unchanged": True}, indent=2))


if __name__ == "__main__":
    main()
