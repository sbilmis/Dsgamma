"""Prepare Referee 1 figures using the existing Figure 2 two-point inputs.

Re-evaluate the residues at all 161 original angles and verify agreement with
the saved Figure 2 CSV before writing the new figures. No existing scientific
outputs or manuscript files are changed. Run from any working directory with
Python, NumPy and Matplotlib; SciPy is not required.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIGURES = HERE / "figures"
FIGURES.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(HERE / ".cache"))
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

from twopoint_ds1_matrix_sumrule import Inputs, projected_ope

SOURCE_CSV = ROOT / "figures/figure2_mixing_angle_widths/figure2_mixing_angle_widths.csv"
SOURCE_BUILDER = ROOT / "figures/figure2_mixing_angle_widths/build_figure2_data.py"
PREFERRED_ANGLE = {"Ds": 26.6, "Bs": 38.5}
HQ_ANGLE = 35.3


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def readonly_manifest():
    return {
        str(path.relative_to(ROOT)): sha256(path)
        for path in sorted((ROOT / "final_manuscript").rglob("*"))
        if path.is_file()
    }


def representative_sample(path):
    """Same median/MAD selection as the existing Figure 2 builder."""
    accepted = [r for r in read_csv(path) if int(r.get("accepted", "1")) == 1]
    candidates = (
        "theta_deg", "mc_GeV", "mb_GeV", "ms_GeV", "ss_GeV3",
        "m0_sq_GeV2", "M2_GeV2", "s01_GeV2", "s02_GeV2",
        "mass_low_GeV", "mass_high_GeV",
    )
    keys = [key for key in candidates if key in accepted[0] and accepted[0][key] != ""]
    values = np.asarray([[float(row[key]) for key in keys] for row in accepted])
    center = np.median(values, axis=0)
    scale = np.median(np.abs(values - center), axis=0)
    scale = np.where(scale > 1e-12, scale, 1.0)
    index = int(np.argmin(np.sum(((values - center) / scale) ** 2, axis=1)))
    return accepted[index]


def inputs_for_sample(sample, sector):
    return Inputs(
        mc=float(sample["mc_GeV"] if sector == "Ds" else sample["mb_GeV"]),
        ms=float(sample["ms_GeV"]),
        qq=float(sample["ss_GeV3"]),
        kappa_s=1.0,
        m0_sq=float(sample["m0_sq_GeV2"]),
        mass_low=2.4595 if sector == "Ds" else float(sample["mass_low_GeV"]),
        mass_high=2.53511 if sector == "Ds" else float(sample["mass_high_GeV"]),
    )


def calculate(sample, sector, theta):
    inp = inputs_for_sample(sample, sector)
    m2 = float(sample["M2_GeV2"])
    constants = []
    for channel, mass, threshold in (
        (0, inp.mass_low, float(sample["s01_GeV2"])),
        (1, inp.mass_high, float(sample["s02_GeV2"])),
    ):
        projection, _ = projected_ope(m2, threshold, theta, channel, inp)
        if projection <= 0 or not math.isfinite(projection):
            raise ValueError(f"Nonpositive or nonfinite projection: {sector}, {theta}, {channel}")
        constants.append(math.sqrt(projection * math.exp(mass**2 / m2) / mass**2))
    f1, f2 = constants
    return dict(sector=sector, theta_deg=theta, f1_MeV=1000*f1,
                f2_MeV=1000*f2, f2_over_f1=f2/f1)


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def format_panel(ax, sector, ratio=False):
    ax.set_title(r"(a) Charm-strange, $D_s$" if sector == "Ds" else r"(b) Bottom-strange, $B_s$", fontsize=11)
    ax.set_xlim(25, 45)
    ax.set_xlabel(r"$\theta$ [degrees]", fontsize=11)
    ax.set_ylabel(r"$f_2(\theta)/f_1(\theta)$" if ratio else r"$f_i(\theta)$ [MeV]", fontsize=11)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=9)
    ax.grid(False)
    ax.axvline(PREFERRED_ANGLE[sector], color="#b22222", ls="-.", lw=1.2, zorder=0)
    ax.axvline(HQ_ANGLE, color="0.48", ls=":", lw=1.2, zorder=0)
    ax.text(PREFERRED_ANGLE[sector]+0.3, 0.04,
            rf"${PREFERRED_ANGLE[sector]:.1f}^\circ$", color="#9d2020",
            transform=ax.get_xaxis_transform(), fontsize=8.5, va="bottom")
    ax.text(HQ_ANGLE-0.3, 0.04, r"$35.3^\circ$", color="0.35",
            transform=ax.get_xaxis_transform(), fontsize=8.5, va="bottom", ha="right")


def make_figures(rows):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "mathtext.fontset": "dejavusans", "pdf.fonttype": 42,
        "ps.fonttype": 42, "axes.linewidth": 0.8,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.65))
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.17, top=0.79, wspace=0.30)
    for ax, sector in zip(axes, ("Ds", "Bs")):
        selected = [r for r in rows if r["sector"] == sector]
        angles = [r["theta_deg"] for r in selected]
        ax.plot(angles, [r["f1_MeV"] for r in selected], color="#0047ab", lw=2)
        ax.plot(angles, [r["f2_MeV"] for r in selected], color="#111111", ls="--", lw=2)
        ax.set_ylim(0, 600)
        ax.yaxis.set_major_locator(MultipleLocator(100))
        format_panel(ax, sector)
    handles = [Line2D([0], [0], color="#0047ab", lw=2, label=r"$f_1$: lower state"),
               Line2D([0], [0], color="#111111", lw=2, ls="--", label=r"$f_2$: higher state")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.54, 1.0), fontsize=10, handlelength=3, columnspacing=3)
    fig.savefig(FIGURES / "decay_constants_vs_theta.pdf")
    fig.savefig(FIGURES / "decay_constants_vs_theta.png", dpi=240)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.25))
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.18, top=0.88, wspace=0.30)
    for ax, sector in zip(axes, ("Ds", "Bs")):
        selected = [r for r in rows if r["sector"] == sector]
        ax.plot([r["theta_deg"] for r in selected],
                [r["f2_over_f1"] for r in selected], color="#0047ab", lw=2)
        # Common scale facilitates direct comparison of the sectors.
        ax.set_ylim(0.10, 0.55)
        ax.yaxis.set_major_locator(MultipleLocator(0.10))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        format_panel(ax, sector, ratio=True)
    fig.savefig(FIGURES / "decay_constant_ratio_vs_theta.pdf")
    fig.savefig(FIGURES / "decay_constant_ratio_vs_theta.png", dpi=240)
    plt.close(fig)


def main():
    protected_before = readonly_manifest()
    source = read_csv(SOURCE_CSV)
    samples, sources, curves, benchmarks = {}, {}, [], []
    errors = []
    for sector in ("Ds", "Bs"):
        sample_path = ROOT / "outputs" / f"twopoint_{sector.lower()}1_matrix_mc.csv"
        sample = representative_sample(sample_path)
        samples[sector] = sample
        sources[str(sample_path.relative_to(ROOT))] = sha256(sample_path)
        by_angle = {}
        for row in source:
            if row["sector"] == sector:
                by_angle.setdefault(float(row["theta_deg"]), []).append(row)
        assert len(by_angle) == 161, (sector, len(by_angle))
        for theta, original in sorted(by_angle.items()):
            result = calculate(sample, sector, theta)
            curves.append(result)
            assert len(original) == 2
            for row in original:
                assert math.isclose(float(row["two_point_source_theta_deg"]), float(sample["theta_deg"]), abs_tol=1e-12)
                for key in ("f1", "f2"):
                    error = abs(result[f"{key}_MeV"] / 1000 - float(row[f"{key}_GeV"]))
                    errors.append(error)
                    assert error < 1e-12, (sector, theta, key, error)
        for theta in sorted({25.0, PREFERRED_ANGLE[sector], HQ_ANGLE, 45.0}):
            benchmarks.append(calculate(sample, sector, theta))
    write_csv(HERE / "decay_constants_theta.csv", curves)
    write_csv(HERE / "benchmark_values.csv", benchmarks)
    make_figures(curves)
    protected_after = readonly_manifest()
    assert protected_before == protected_after, "Protected manuscript files changed"
    provenance = {
        "purpose": "Referee 1 item 3: fixed-input angle scans consistent with existing Figure 2",
        "interpretation": "Rotated-current residues, not unconditionally weak axial-current decay constants",
        "uncertainty": "No uncertainty bands; not marginalized posterior predictions",
        "source_sha256": {str(SOURCE_CSV.relative_to(ROOT)): sha256(SOURCE_CSV),
                          str(SOURCE_BUILDER.relative_to(ROOT)): sha256(SOURCE_BUILDER),
                          "scripts/twopoint_ds1_matrix_sumrule.py": sha256(ROOT / "scripts/twopoint_ds1_matrix_sumrule.py"),
                          **sources},
        "selected_two_point_samples": samples,
        "angle_grid": {"min_degrees": 25, "max_degrees": 45, "points_per_sector": 161},
        "maximum_absolute_agreement_error_GeV": max(errors),
        "number_of_residue_comparisons_with_existing_csv": len(errors),
        "final_manuscript_unchanged": protected_before == protected_after,
        "final_manuscript_sha256": protected_before,
    }
    (HERE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"maximum_absolute_agreement_error_GeV": max(errors),
                      "benchmarks": benchmarks, "final_manuscript_unchanged": True}, indent=2))


if __name__ == "__main__":
    main()
