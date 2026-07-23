# Mathematica notebooks

Open `DsBs_gamma_symbolic_derivation.nb` in Mathematica and evaluate the input
cells from top to bottom (`Evaluation` -> `Evaluate Notebook`), or run them one
at a time with `Shift+Enter`.

The notebook covers:

1. conventions and axial/tensor current vertices;
2. the standard correlator and Wick signs;
3. free, vacuum-condensate, electromagnetic, and one-gluon propagator terms;
4. honest implementation status for the condensate contributions;
5. heavy- and strange-line propagator routing;
6. fermion-line Ward identities;
7. hard Dirac traces and E1 triangle cores;
8. shifted Feynman-parameter denominators; and
9. two-particle soft Fierz traces;
10. physical-state projection and residue-weighted couplings;
11. mixed-current decay constants and the exact pure-axial Colangelo limit;
12. Colangelo-style axial OPE components and the final mixed coupling assembly.

It is generated from `../scripts/build_dsbs_symbolic_notebook.wl`.  Edit the
generator when changing notebook source cells, then rebuild with

```bash
/Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
  -script scripts/build_dsbs_symbolic_notebook.wl
```

from the `papers/Ds_Bs_gamma/` directory.  The four `mathematica_*.wl` files
in `../scripts/` are the corresponding noninteractive regression checks.
