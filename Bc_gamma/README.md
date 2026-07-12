# Bc1 radiative transitions

Workspace for a separate paper on

- `Bc1(6743) -> Bc gamma`
- `Bc1(6743) -> Bc* gamma`
- `Bc1(6750) -> Bc gamma`
- `Bc1(6750) -> Bc* gamma`

This channel is not a direct copy of the `Ds`/`Bs` light-cone calculation.
Both valence quarks are heavy, so the leading baseline should be a hard-photon
external-field / three-point QCD sum-rule calculation with photon emission from
the charm and bottom lines. The photon-DA and light-quark-condensate terms used
for strange spectators are not present as leading nonperturbative input here.

The notes in `notes/` track the physics decisions and the derivation. Scripts in
`scripts/` are for reproducible numerical checks and tables.

Current controlled pseudoscalar-channel driver:

- `scripts/bc_ps_complete_analysis.py`

This script recomputes the \(B_{c1}\to B_c\gamma\) hard-photon baseline with
physical two-point residues \(f_1\) and \(f_2\), exports stability/Monte Carlo
tables, and writes publication-style PDF plots into `outputs/`.  The
\(B_c^\ast\gamma\) channel now uses the standard vector-current normalization;
its contact-support audit finds no ordinary two-channel contact contribution.
The \(B_c\gamma\) widths are the controlled predictions.  The \(B_c^\ast\gamma\)
widths include the explicitly integrated all-positive dimension-four \(G^2\)
sector and have completed contact-support and tensor-current sign audits.  The
large \(6750/6743\) vector hierarchy is therefore a convention-fixed
mixed-current result in the adopted two-point state assignment, driven by
destructive interference in the lower state and constructive interference in
the higher state.
