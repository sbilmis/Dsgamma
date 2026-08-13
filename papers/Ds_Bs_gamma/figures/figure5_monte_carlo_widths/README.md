# Figure 5: Monte Carlo Width Distributions

This folder contains the preferred lattice-normalized Monte Carlo width
distributions with the mixing angle varied over `25 deg <= theta <= 45 deg`.

## Files

- `build_figure5_data.py` builds the Monte Carlo data and summary CSV files.
- `plot_figure5_preview.py` creates the Python comparison PDF and PNG.

## Varied Inputs

- Transition Borel parameter `M^2`.
- Continuum threshold `s0`.
- Quark masses, hadron masses, decay constants, condensates, photon-DA shape
  parameters, and lattice-normalized `f_gamma,s^perp`.
- Accepted two-point residue samples.
- Mixing angle over `25 deg <= theta <= 45 deg`.

The quoted central values should be medians with 16--84 percentile intervals.
Gaussian fits are diagnostic only and should not define the paper values for
the cancellation-sensitive channels.
