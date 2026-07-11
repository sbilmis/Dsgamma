# Draft: \(D_s^\ast\gamma\) and \(B_s^\ast\gamma\)

This folder contains a standalone PRD-style draft for the vector-final-state
radiative transitions

\[
D_{s1}(2460,2536)\to D_s^\ast\gamma,\qquad
B_{s1}(5750,5830)\to B_s^\ast\gamma .
\]

It is intentionally separate from the current \(D_s\gamma\) manuscript.  The
draft uses the same bibliography file as the first paper:

```text
../draft_prd/all.bib
```

Compile from this folder with

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

The figures and comparison tables are read from the existing calculation
outputs under `../outputs`, `../Dsstar_gamma/outputs`, and
`../Bsstar_gamma/outputs`.
