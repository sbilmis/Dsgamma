# Angle-Only Benchmark Diagnostic

This note records a diagnostic check of whether the differences between our
widths and representative literature values can be explained only by changing
the physical mixing angle.

The script

```bash
python3 scripts/angle_benchmark_scan.py
```

reads the final Monte Carlo rows, reconstructs the two reduced unmixed
amplitudes from

\[
  f_1 g_1=\sin\theta\,A+\cos\theta\,B,\qquad
  f_2 g_2=\cos\theta\,A-\sin\theta\,B,
\]

and then rotates \(A\) and \(B\) over \(0^\circ\le\theta\le90^\circ\), keeping
all other sampled QCD inputs fixed.  The output is stored in

```text
outputs/angle_benchmark_scan.csv
```

## Main Results

For the pseudoscalar final state, \(D_s\gamma\), angle variation alone cannot
reproduce the large Bondar--Milstein value for
\(D_{s1}(2460)\to D_s\gamma\).  With the preferred lattice photon input, the
best two-channel compromise occurs near \(\theta\simeq5^\circ\), but it gives
only

\[
  \Gamma[D_{s1}(2460)\to D_s\gamma]\simeq 27.4~{\rm keV},\qquad
  \Gamma[D_{s1}(2536)\to D_s\gamma]\simeq 7.18~{\rm keV}.
\]

Thus the discrepancy in this channel is not a mixing-angle effect; it reflects
a different underlying normalization/dynamics of the transition amplitude.

For the vector final state, \(D_s^\ast\gamma\), the situation is different.
With the preferred lattice photon input, the best simultaneous comparison to
the Bondar--Milstein pair
\[
  (104~{\rm keV},\,29~{\rm keV})
\]
is obtained near \(\theta\simeq37.8^\circ\), giving

\[
  \Gamma[D_{s1}(2460)\to D_s^\ast\gamma]\simeq 89.1~{\rm keV},\qquad
  \Gamma[D_{s1}(2536)\to D_s^\ast\gamma]\simeq 28.3~{\rm keV}.
\]

This is inside the usual physical angle region and supports the interpretation
that the vector-channel hierarchy is largely an interference/mixing effect.

For the bottom-strange vector final state, the angle-only scan shows mixed
behavior.  The \(B_{s1}(5830)\to B_s^\ast\gamma\) Bondar--Milstein value
\(\sim41\) keV can be reproduced near \(\theta\simeq59.4^\circ\), outside the
usual charm-strange angle window.  The \(B_{s1}(5750)\to B_s^\ast\gamma\)
target \(\sim20\) keV is not reached by changing the angle in the preferred
lattice-photon setup; the closest median is about \(11.0\) keV.  This suggests
that the bottom-sector discrepancy is not mainly a small angle-choice issue.

## Interpretation

The discrepancies are therefore channel dependent:

- \(D_s\gamma\): angle variation does not explain the large
  \(D_{s1}(2460)\) width in Bondar--Milstein.
- \(D_s^\ast\gamma\): the preferred lattice normalization plus a physical
  angle around \(38^\circ\) gives a pattern close to Bondar--Milstein.
- \(B_s^\ast\gamma\): angle variation helps one state but not the full pair,
  so remaining differences likely come from model-dependent radial wave
  functions, heavy-quark scaling, and normalization choices.

This diagnostic should be used only as a consistency check.  It keeps our QCD
input ensemble fixed and rotates the already extracted reduced amplitudes; it
does not refit the two-point decay constants or the continuum thresholds at
each angle.
