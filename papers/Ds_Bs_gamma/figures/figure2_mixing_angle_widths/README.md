# Figure 2: Mixing-Angle Widths

This folder contains the reproducible draft for the mixing-angle dependence of
the radiative widths.

## Files

- `build_figure2_data.py` builds `figure2_mixing_angle_widths.csv`.
- `plot_figure2.wl` creates the Mathematica/MaTeX PDF and PNG.
- `plot_figure2_preview.py` creates the Python comparison PDF and PNG.

## Central Choices

- Angle scan: `25 deg <= theta <= 45 deg`.
- Reference line: `theta = 35.3 deg`.
- `D_s` sector: `M^2 = 3.75 GeV^2`, with `s0 = 9.0, 9.5 GeV^2`.
- `B_s` sector: `M^2 = 12.0 GeV^2`, with `s0 = 40.0, 41.0 GeV^2`.
- Physical-current residues are reprojected from the accepted AA/AB/BB
  two-point matrices at each plotted angle.
- The plots use no background grid and use line styles chosen to remain
  separable in grayscale.
