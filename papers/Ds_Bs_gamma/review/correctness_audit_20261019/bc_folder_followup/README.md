# Follow-up on the author-supplied mixing-code folder

19 October 2026. The author supplied `/Users/sbilmis/developer/physics/Bc_mixing`
and recalled Mathematica and Python, but did not remember the original
heavy-strange calculation's filename.

## What the supplied folder contains

`BcMixingMomentum.wl` and its notebook implement the **heavy-heavy Bc**
correlator for bbar c. The retained OPE uses the perturbative loop and gluon
condensates. The parameter lists contain mb, mc, G2 and G3; they do not contain
the strange-quark and mixed quark-gluon condensate terms needed for the
original Ds/Bs analysis. The coordinate-space implementations have the same
Bc scope. The saved publication MC summaries give means around 43.30 and
43.34 degrees in the two retained versions. These are Bc outputs, not the
26.6/38.5-degree heavy-strange outputs.

The Bc repository history starts on 13 May 2026. Its README cites the earlier
heavy-light study as methodology. No positively identified copy of the
original 2025 heavy-strange numerical calculation was found in the supplied
folder, including its retained archives. This is a statement about the files
inspected, not proof that the original code is unavailable elsewhere.

## Additional algebraic check completed

`check_bc_perturbative_mapping.wl` loads the supplied package read-only and
compares its actual perturbative density functions against the current
heavy-strange polynomials. The flavor-role mapping is mb → ms for the
antiquark and mc → mQ for the quark.

For **each** of AA, AB and BB, the exact ratio is 3. Thus the two perturbative
matrices agree up to one common factor, which cancels in the perturbative
mixing-angle ratio. This concerns the overall spin-projection normalization;
it is not a change to the tensor-current factor `1/(mQ+ms)`.

See `bc_perturbative_mapping.json` and `wolfram_run.log`. The result supports
the current leading-order matrix, but does not resolve the earlier appendix's
mixed-condensate difference. The Bc gluon-condensate expansion must not be
used as the full heavy-light OPE merely by substituting a strange-quark mass.

Reproduce from the present project root:

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noinit -noprompt -script review/correctness_audit_20261019/bc_folder_followup/check_bc_perturbative_mapping.wl
```

## Earlier D/B paper project found nearby

The related paper folder is:

`/Users/sbilmis/developer/physics/LateX/papers/mixing_of_D_and_B_mesons`

It contains the earlier manuscript, a `Submission/` directory, figures and
20 precomputed `.dat` curves. `Submission/article_R1.tex` quotes the published
28.2, 26.6, 38.6 and 38.5-degree results. The available sources
`Figures/data_points/fig1_D10_D10prime_data/fig1.py` and `fig1_me.nb` are plotting
routines importing precomputed data; no sum-rule integration or angle
evaluation was found in them. The repository's history for this paper also
lists only these two computational source paths.

The available raw curves belong to an earlier stage: for example the Ds data
at s0=10 begin near 17.11 degrees and the Bs data at s0=43 near 21.37 degrees.
The top-level older manuscript quotes corresponding older-angle results.
These raw curves are not a validated baseline for the final paper. No cause
for the historical change is inferred, and the author-confirmed current
normalization is not reopened. See `historical_plot_data_inventory.csv`.

## Search scope and preservation

- Inspected the supplied folder's source and documentation, including ignored
  archive paths; a content pass covered 167 candidate text/source/data files.
- Searched 1,536 Mathematica/Python-style source paths in the available
  developer tree using paper identifiers, angle-output names and relevant
  physics terms. This included generic filenames, not just names containing
  “mixing”. No original heavy-strange implementation was identified.
- Checked 475 additional Jupyter/other-source paths for the same physics
  terms in cell sources; four unrelated course notebooks were not parseable.
- Recorded hashes of the 31 directly relevant source/data files in
  `source_manifest.json`.
- Made no edits or moves in either supplied research folder. The Bc tracked
  working tree remained clean. All 11 protected present-manuscript hashes
  remain unchanged. New files are confined to the present audit directory,
  with status notes updated in the present project.

## Remaining limitation

This lead strengthens the perturbative cross-check and locates the earlier
paper's document project. It does not yet recover the calculation that
produced the final Ds/Bs central values and uncertainties. The original
implementation, or a complete independent reconstruction of its heavy-light
OPE and numerical prescription, is still needed before claiming that the
published angle baseline and its updated-input sensitivity are verified.
