# Manuscript Folder

This folder is now an archived local snapshot.  It is not the active paper
editing target.  The active calculation and redo-analysis record is in
`../notes/Ds1_to_Ds_gamma_LCSR_notes.tex`; the paper itself will be edited
elsewhere.

## Archived Drafts

- `dsbs_radiative_lcsr_polished.tex` - archived polished draft snapshot.
- `dsbs_radiative_lcsr_polished.pdf` - compiled polished draft.
- `dsbs_radiative_lcsr.tex` - shorter compact draft kept for comparison.
- `references.bib` - bibliography used by both manuscript drafts.

## Build If Needed

From this folder:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error dsbs_radiative_lcsr_polished.tex
```

The manuscript figures and tables are loaded from `../outputs/`.

## Current Physics Normalization

The bottom-strange rows use Pullin-Zwicky basis-current inputs followed by the
mixed-current diagonalization closure.  The current preferred observed-state
result is

```text
Gamma[Bs1(5830) -> Bs gamma] = 0.342 [0.160, 1.026] keV
```

with the lattice photon normalization.
