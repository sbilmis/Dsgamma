# Ds1 -> Ds* gamma LCSR extension

This folder is a separate workspace for the vector-final-state channels

\[
  D_{s1}(2460)\to D_s^\ast\gamma,\qquad
  D_{s1}(2536)\to D_s^\ast\gamma,
\]

and, later if desired,

\[
  B_{s1}(5750)\to B_s^\ast\gamma,\qquad
  B_{s1}(5830)\to B_s^\ast\gamma.
\]

The existing `D_s gamma` analysis should remain stable while this extension is
developed.  The scripts and notes here may import common inputs and photon DA
utilities from the parent project, but final outputs for this channel should be
written inside this folder.

## Recommendation

Treat this as a second paper or a companion paper unless the calculation turns
out to be almost identical numerically.  The vector final state changes the
hadronic parametrization, Lorentz structures, decay-width formula, and input
decay constants.  It is related enough to reuse infrastructure, but different
enough that adding it to the current paper would make that paper less focused.

## Immediate Tasks

1. Define the vector-final-state correlator with a \(D_s^\ast\) interpolating
   current.
2. Choose the gauge-invariant Lorentz structure used to extract the form
   factor.
3. Derive the hadronic representation and width formula for
   \(1^+\to1^-\gamma\).
4. Reuse the photon-DA/OPE machinery where possible, but redo the Dirac traces
   because the final current changes.
5. Build a controlled Stage-1 baseline, then add the same Stage-2 pieces as in
   the \(D_s\gamma\) study.
6. Only after the baseline is stable, decide whether the results belong in the
   current paper as an added section/table or deserve a separate article.

