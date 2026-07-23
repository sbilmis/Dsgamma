# Radiative heavy-meson projects

This repository contains three paper-level projects:

| Project | Channels | Method |
|---|---|---|
| [`papers/Ds_Bs_gamma`](papers/Ds_Bs_gamma/) | \(D_{s1},B_{s1}\to D_s,B_s\,\gamma\) | light-cone QCD sum rules, pseudoscalar final state |
| [`papers/Dsstar_Bsstar_gamma`](papers/Dsstar_Bsstar_gamma/) | \(D_{s1},B_{s1}\to D_s^\ast,B_s^\ast\,\gamma\) | light-cone QCD sum rules, vector final state |
| [`papers/Bc_radiative`](papers/Bc_radiative/) | \(B_{c1}\to B_c^{(\ast)}\gamma\) | double-heavy external-field / three-point QCD sum rules |

Reusable conventions, photon distribution amplitudes, and the common PRD
bibliography are kept in [`shared`](shared/). Each paper owns its manuscripts,
notes, scripts, and generated outputs.

Run calculation scripts from the relevant paper directory so the Mathematica
scripts place outputs in that paper's `outputs/` tree. Python scripts resolve
their project and shared paths from their file locations.

The repository-level `.venv/` remains the common Python environment.
