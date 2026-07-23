# \(D_s^\ast\gamma\) and \(B_s^\ast\gamma\) paper

Light-cone QCD sum-rule analysis of

- `Ds1(2460,2536) -> Ds* gamma`, and
- `Bs1(5750,5830) -> Bs* gamma`.

The vector final state requires its own hadronic parametrization, Lorentz
structures, decay-width formula, and vector decay constants. It reuses the
mixed-current and photon-DA conventions of the pseudoscalar-final-state paper.

## Layout

- `Dsstar_gamma/` - charm-strange scripts, notes, and outputs.
- `Bsstar_gamma/` - bottom-strange scripts, notes, and outputs.
- `scripts/` - analyses comparing or combining both vector channels.
- `outputs/` - combined tables and cross-channel diagnostics.
- `notes/` - combined analysis notes.
- `manuscript/` - standalone PRD-style paper draft.

Run scripts from `papers/Dsstar_Bsstar_gamma/`. The shared photon-DA module is
in `../../shared/photon_da.py`, and the manuscript uses the bibliography in
`../../shared/bibliography/`.
