# Manuscript Folder

Use this folder when editing the paper.

## Main Drafts

- `dsbs_radiative_lcsr_polished.tex` - current polished paper draft.
- `dsbs_radiative_lcsr_polished.pdf` - compiled polished draft.
- `dsbs_radiative_lcsr.tex` - shorter compact draft kept for comparison.
- `references.bib` - bibliography used by both manuscript drafts.

## Build

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
