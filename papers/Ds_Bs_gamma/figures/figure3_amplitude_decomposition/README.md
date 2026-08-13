# Figure 3: Amplitude Decomposition

This folder contains the signed physical-current amplitude decomposition for
the four radiative channels.

Current recommendation: keep the plot outputs as diagnostics and use
`amplitude_decomposition_table.tex` in the manuscript instead of a Figure 3
graphic.  The table communicates the sign pattern more directly.

## Files

- `build_figure3_data.py` builds `figure3_amplitude_decomposition.csv`.
- `plot_figure3.wl` creates the Mathematica/MaTeX PDF and PNG.
- `plot_figure3_preview.py` creates the Python comparison PDF and PNG.
- `amplitude_decomposition_table.tex` is the recommended manuscript replacement.

## Central Choices

- Mixing angle: `theta = 35.3 deg`.
- Transition points match Figure 2:
  - `D_s` sector: `M^2 = 3.75 GeV^2`, with `s0 = 9.0, 9.5 GeV^2`.
  - `B_s` sector: `M^2 = 12.0 GeV^2`, with `s0 = 40.0, 41.0 GeV^2`.
- The plotted quantities are the two signed components entering the physical
  amplitude and their sum.
