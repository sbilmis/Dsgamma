<!-- LATTICE UPDATE -->
**3 October 2026:** [Lattice-normalized results with both twist-2 signs](lattice_results.md). Includes nominal and refitted 90-degree windows, window intervals and separate contribution medians. Frozen outputs and their SHA-256 manifest remain unchanged.

# Current calculation

[Report](RECALCULATION_REPORT.md) · [central table](central_table.md) · [LaTeX table](central_table.tex) · [setup](calculation_setup.json)

**Follow-up verification:** [definition and scan audit](blind_definition_scan_audit/AUDIT.md). The link prescription has an independent coordinate-space derivation. Expanded scans support the original Borel bounds, but do not uniquely select the threshold interval. The tables remain conditional on that interval; independent evaluation of a shared setup is not a historically blind choice of setup.

[Susceptibility sign and charge normalization](susceptibility_charge_audit/EXPLANATION.md) gives the explicit convention chain, fresh Dirac-trace checks and charge-sector decomposition. In our convention chi is positive; integer F_C does not test its sign.

[Complete structure and angle tables](structure_tables/STRUCTURE_TABLES.md) now cover all previously documented extraction routes plus the clean combinations: nominal mixing, 90 degrees on nominal windows, and 90 degrees with refitted windows. The clean refitted 90-degree results were already calculated; the full side-by-side structure tables are new. Both this comparison and the sign/charge explanation are incorporated in the main report and PDF.

The alternative routes are diagnostics and are not averaged into the central prediction. Differences persist even in the transverse BB channel at 90 degrees, where a pseudoscalar overlap is absent. Their full explanation remains unresolved. The accepted-window spreads do not quantify this structure dependence.

The primary results use zeta2=31/110, obtained by imposing the exact twist-4 parameter relation at zeta1=0.4 and zeta2Plus=0. Rounded inputs are a separately identified run. The report explains the continuum model, the OPE truncation and the limitations of the two-point window criteria for three-point stability.

## Data

- `final_pointwise_results.csv`: all nominal and refitted pure-current moments, form factors and residues.
- `final_summary.csv`: medians and deterministic 16–84% window intervals.
- `structure_tables/all_routes_pointwise.csv`: 500 recomputed scenario-points, every route, residue, eight stored amplitudes and acceptance diagnostics.
- `structure_tables/all_routes_summary.csv`: all route medians, 16th/84th percentiles and extrema; `windows.csv` records counts and ranges.
- `structure_tables/structure_tables.tex`: 15 LaTeX tables, with the route key in the companion Markdown report.
- `structure_tables/validation.json`: pure-pole algebra, pseudoscalar cancellation, charge closure, clean-route agreement and quadrature refinement.
- `final_sector_channel_structure.csv`: every retained sector, channel and structure on nominal windows.
- `two_point_window_scan.csv`: accepted and rejected candidates, with every criterion.
- `old62_controlled_updates.csv`: isolated updates on the original points.
- `four_point_sector_benchmark.csv`: our sector-by-sector exports at the four historical benchmark points, with each input and continuum choice labeled.
- `old62_angle_pointwise.csv`: every original point evaluated over the mixing-angle scan.
- `angle_curves.csv`, `angle_endpoints.csv`: rotated and fixed-residue angle diagnostics.
- `model_and_input_sensitivity.csv`: deterministic parameter and continuum variations.
- `three_point_sector_hierarchy.csv`: twist hierarchy and cancellation diagnostics.
- `partial_uncertainty_summary.csv`: conditional photon/angle/window sampling; assumptions in `partial_uncertainty_model.json`.
- `independent_validation.json`, `da_eom_checks.wl`, `contact_tensor_reconstruction_checks.wl`: independent tests, without counterpart data.

The eight figures are available in `figures` as PDF and PNG.

## Reproduction

Numerical scripts require Python 3.10+ and NumPy. Plotting additionally requires Matplotlib. Report generation uses Pandoc, XeLaTeX and Poppler. Symbolic regeneration requires Wolfram Language and FeynCalc.

From the repository root, run these scripts in order with a Python environment containing the required packages:

```sh
python papers/multipole_moments/black_box_comparison/scripts/run_recalculation.py
python papers/multipole_moments/black_box_comparison/scripts/audit_three_point_convergence.py
python papers/multipole_moments/black_box_comparison/scripts/validate_recalculation.py
python papers/multipole_moments/black_box_comparison/scripts/recalculate_partial_uncertainty.py
python papers/multipole_moments/black_box_comparison/scripts/export_recalculation_exchange.py
python papers/multipole_moments/black_box_comparison/scripts/plot_recalculation.py
python papers/multipole_moments/black_box_comparison/scripts/tabulate_all_structures.py
python papers/multipole_moments/black_box_comparison/scripts/write_recalculation_report.py
bash papers/multipole_moments/black_box_comparison/scripts/build_recalculation_pdf.sh
```

Report generation preserves the historical exchange reply by default; passing `--write-reply` explicitly also regenerates that reply. The sign/charge section is maintained in `susceptibility_charge_audit/report_section.md`.

The numerical engine reads the three `*_coefficients.json` files generated from our own symbolic kernels. To regenerate them, execute these neighboring Wolfram scripts in order:

1. `derive_epsp_hard.wl` and `derive_epsp_two.wl`.
2. `export_recalculation_kernels.wl`.
3. `build_parametric_da_weights.wl`, then `export_parametric_three_particle.wl`.
4. `derive_contact_recalculation.wl`, then `export_contact_coefficients.wl`.
5. `check_da_equations_of_motion.wl`.

All source paths are explicit and restricted to our own derivations. The old 62 points come from our frozen `stage1` input table. Acceptance of new points is calculated from our two-point sum rule within a prescribed search domain. That does not independently establish the domain or cutoff choices. Neither counterpart CSVs nor counterpart source files feed the numerical engine.

`freeze_sha256.json` records the central calculation. `delivery_manifest.json` additionally records the final tables, figures, report and derivation scripts.
